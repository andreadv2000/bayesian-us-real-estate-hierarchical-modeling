# Lightweight DIC-only refit for M0, M1, and selected heteroskedastic M2.
#
# rjags::dic.samples() requires a live JAGS model object, so DIC cannot be
# recomputed from the saved output/mcmc/*.rds posterior draws alone.  This
# script refits each model only far enough to create a live post-warmup JAGS
# object and sample DIC.  It does not monitor or save posterior parameter
# chains, and it does not overwrite the production posterior draws from
# scripts/05_fit_models.R.
#
# Environment variables:
#   DIC_N_CHAINS  default 4; must be at least 2 for rjags::dic.samples()
#   DIC_N_ADAPT   default 2000
#   DIC_N_BURN    default 5000
#   DIC_N_ITER    default 5000
#   DIC_SEED      default 20260623
#   DIC_M2_SIGMA_PRIOR_CENTER       default 0.5
#   DIC_M2_GAMMA0_PRIOR_SCALE       default 0.5
#   DIC_M2_GAMMA_PRIOR_SCALE        default 0.2
#   DIC_M2_SIGMA_STATE_PRIOR_SCALE  default 0.2
#   DIC_DRY_RUN   default false; if true, validate inputs and stop before JAGS

required_dic_packages <- c("here", "readr", "rjags", "tibble")
missing_dic_packages <- required_dic_packages[!vapply(
  required_dic_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_dic_packages) > 0L) {
  stop(
    "Install the package(s) required for lightweight DIC: ",
    paste(missing_dic_packages, collapse = ", "),
    call. = FALSE
  )
}

source(here::here("scripts", "12_m2_heteroskedastic_helpers.R"))
source(here::here("scripts", "07_model_comparison.R"))

read_dic_positive_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.integer(default))
  }

  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) {
    stop("Environment variable `", name, "` must be a positive integer.", call. = FALSE)
  }
  parsed
}

read_dic_positive_number <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.numeric(default))
  }

  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed) || parsed <= 0) {
    stop("Environment variable `", name, "` must be a positive number.", call. = FALSE)
  }
  parsed
}

read_dic_flag <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, unset = "")))
  if (!nzchar(value)) {
    return(default)
  }
  if (value %in% c("1", "true", "yes")) {
    return(TRUE)
  }
  if (value %in% c("0", "false", "no")) {
    return(FALSE)
  }
  stop("Environment variable `", name, "` must be true or false.", call. = FALSE)
}

timestamp_message <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

augment_zip_centering_data <- function(data_list) {
  required_data <- c("N", "S", "Z", "state", "zip")
  missing_data <- setdiff(required_data, names(data_list))
  if (length(missing_data) > 0L) {
    stop(
      "Cannot build ZIP centering data because JAGS data are missing: ",
      paste(missing_data, collapse = ", "),
      call. = FALSE
    )
  }

  zip_state <- vapply(seq_len(data_list$Z), function(zip_id) {
    states <- unique(as.integer(data_list$state[data_list$zip == zip_id]))
    states <- states[!is.na(states)]
    if (length(states) != 1L) {
      stop("Each ZIP must map to exactly one state for M2 ZIP centering.", call. = FALSE)
    }
    states
  }, integer(1L))
  n_zips_by_state <- tabulate(zip_state, nbins = data_list$S)
  if (any(n_zips_by_state < 1L)) {
    stop("Every state must have at least one ZIP for M2 ZIP centering.", call. = FALSE)
  }

  k_max <- max(n_zips_by_state)
  zip_in_state <- matrix(1L, nrow = data_list$S, ncol = k_max)
  zip_in_state_weight <- matrix(0, nrow = data_list$S, ncol = k_max)
  for (state_id in seq_len(data_list$S)) {
    state_zips <- which(zip_state == state_id)
    zip_in_state[state_id, seq_along(state_zips)] <- state_zips
    zip_in_state_weight[state_id, seq_along(state_zips)] <- 1
  }

  data_list$K <- as.integer(k_max)
  data_list$zip_state <- as.integer(zip_state)
  data_list$n_zips_by_state <- as.integer(n_zips_by_state)
  data_list$zip_in_state <- zip_in_state
  data_list$zip_in_state_weight <- zip_in_state_weight
  data_list
}

select_jags_data <- function(data_list, required_names) {
  missing_names <- setdiff(required_names, names(data_list))
  if (length(missing_names) > 0L) {
    stop(
      "JAGS data are missing required item(s): ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }
  data_list[required_names]
}

make_dic_inits <- function(data_list, model_file, n_chains, seed) {
  model_text <- paste(readLines(model_file, warn = FALSE), collapse = "\n")
  set.seed(seed)
  lapply(seq_len(n_chains), function(chain_id) {
    values <- list(
      alpha = stats::rnorm(1L, 0, 0.2),
      beta = stats::rnorm(data_list$P, 0, 0.1),
      sigma = stats::runif(1L, 0.5, 1.5),
      gamma0 = stats::rnorm(
        1L,
        if ("log_sigma_prior_mean" %in% names(data_list)) {
          data_list$log_sigma_prior_mean
        } else {
          log(0.5)
        },
        0.1
      ),
      gamma = if ("P_sigma" %in% names(data_list)) {
        stats::rnorm(data_list$P_sigma, 0, 0.03)
      } else {
        numeric(0)
      },
      sd_sigma_state = stats::runif(1L, 0.02, 0.2),
      sd_state = stats::runif(1L, 0.1, 0.8),
      sd_zip = stats::runif(1L, 0.1, 0.8),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = as.integer(seed + 1000L * chain_id)
    )
    if (!grepl("sigma ~", model_text, fixed = TRUE)) {
      values$sigma <- NULL
    }
    if (!grepl("gamma0", model_text, fixed = TRUE)) {
      values$gamma0 <- NULL
    }
    if (!grepl("gamma[", model_text, fixed = TRUE)) {
      values$gamma <- NULL
    }
    if (!grepl("sd_sigma_state", model_text, fixed = TRUE)) {
      values$sd_sigma_state <- NULL
    }
    if (!grepl("sd_state", model_text, fixed = TRUE)) {
      values$sd_state <- NULL
    }
    if (!grepl("sd_zip", model_text, fixed = TRUE)) {
      values$sd_zip <- NULL
    }
    values
  })
}

fit_live_dic_model <- function(spec, n_chains, n_adapt, n_burn, seed) {
  if (!file.exists(spec$model_file)) {
    stop("JAGS model file not found: ", spec$model_file, call. = FALSE)
  }

  inits <- make_dic_inits(
    data_list = spec$data,
    model_file = spec$model_file,
    n_chains = n_chains,
    seed = seed
  )
  model <- rjags::jags.model(
    file = spec$model_file,
    data = spec$data,
    inits = inits,
    n.chains = n_chains,
    n.adapt = n_adapt,
    quiet = TRUE
  )
  stats::update(model, n.iter = n_burn, progress.bar = "none")
  list(model = model, inits = inits)
}

summarise_dic_object <- function(model_name, dic_object) {
  d_bar <- sum(as.numeric(dic_object$deviance))
  p_d <- sum(as.numeric(dic_object$penalty))
  if (!all(is.finite(c(d_bar, p_d)))) {
    stop("DIC calculation returned a non-finite deviance or penalty.", call. = FALSE)
  }
  tibble::tibble(model = model_name, Dbar = d_bar, pD = p_d, DIC = d_bar + p_d)
}

run_lightweight_dic <- function(
    analysis_objects_path = here::here("data", "processed", "analysis_objects.rds"),
    tables_dir = here::here("output", "tables"),
    logs_dir = here::here("output", "logs"),
    n_chains = read_dic_positive_integer("DIC_N_CHAINS", 4L),
    n_adapt = read_dic_positive_integer("DIC_N_ADAPT", 2000L),
    n_burn = read_dic_positive_integer("DIC_N_BURN", 5000L),
    n_iter = read_dic_positive_integer("DIC_N_ITER", 5000L),
    seed = read_dic_positive_integer("DIC_SEED", 20260623L),
    m2_sigma_prior_center = read_dic_positive_number("DIC_M2_SIGMA_PRIOR_CENTER", 0.5),
    m2_gamma0_prior_scale = read_dic_positive_number("DIC_M2_GAMMA0_PRIOR_SCALE", 0.5),
    m2_gamma_prior_scale = read_dic_positive_number("DIC_M2_GAMMA_PRIOR_SCALE", 0.2),
    m2_sigma_state_prior_scale = read_dic_positive_number("DIC_M2_SIGMA_STATE_PRIOR_SCALE", 0.2),
    dry_run = read_dic_flag("DIC_DRY_RUN", default = FALSE)) {
  if (!file.exists(analysis_objects_path)) {
    stop(
      "Prepared analysis objects not found: ", analysis_objects_path,
      ". Run scripts/03_prepare_jags_data.R first.",
      call. = FALSE
    )
  }
  if (any(c(n_chains, n_adapt, n_burn, n_iter, seed) < 1L)) {
    stop("DIC settings must be positive integers.", call. = FALSE)
  }
  if (n_chains < 2L) {
    stop(
      "DIC_N_CHAINS must be at least 2 because rjags::dic.samples() ",
      "requires two or more parallel chains.",
      call. = FALSE
    )
  }

  analysis_objects <- readRDS(analysis_objects_path)
  if (!"jags_train" %in% names(analysis_objects)) {
    stop("Analysis objects do not contain the frozen JAGS training data.", call. = FALSE)
  }

  jags_train <- analysis_objects$jags_train
  jags_train_m2 <- prepare_heteroskedastic_m2_data(
    analysis_objects = analysis_objects,
    sigma_prior_center = m2_sigma_prior_center,
    gamma0_prior_scale = m2_gamma0_prior_scale,
    gamma_prior_scale = m2_gamma_prior_scale,
    sigma_state_prior_scale = m2_sigma_state_prior_scale
  )
  model_specs <- list(
    list(
      key = "M0",
      label = "M0: pooled",
      model_file = here::here("models", "m0_pooled.jags"),
      data = select_jags_data(jags_train, c("N", "P", "y", "X"))
    ),
    list(
      key = "M1",
      label = "M1: state",
      model_file = here::here("models", "m1_state.jags"),
      data = select_jags_data(jags_train, c("N", "P", "y", "X", "S", "state"))
    ),
    list(
      key = "M2",
      label = "M2: heteroskedastic state + ZIP",
      model_file = here::here("models", "m2_heteroskedastic_state_sigma.jags"),
      data = select_jags_data(
        jags_train_m2,
        c(
          "N", "P", "y", "X", "S", "Z", "state", "zip",
          "K", "zip_state", "n_zips_by_state",
          "zip_in_state", "zip_in_state_weight",
          "X_sigma", "P_sigma", "log_sigma_prior_mean",
          "tau_gamma0_prior", "tau_gamma_prior",
          "tau_sigma_state_prior"
        )
      )
    )
  )

  for (spec in model_specs) {
    if (!file.exists(spec$model_file)) {
      stop("JAGS model file not found: ", spec$model_file, call. = FALSE)
    }
  }

  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  table_path <- file.path(tables_dir, "model_comparison_dic.csv")
  partial_path <- file.path(tables_dir, "model_comparison_dic_partial.csv")
  settings_path <- file.path(tables_dir, "model_comparison_dic_settings.csv")
  runtime_path <- file.path(tables_dir, "model_comparison_dic_runtime.csv")
  log_path <- file.path(logs_dir, "model_comparison_dic.log")

  settings <- tibble::tibble(
    setting = c(
      "workflow", "n_chains", "n_adapt", "n_burn", "n_dic_iter",
      "seed", "m2_sigma_prior_center", "m2_gamma0_prior_scale",
      "m2_gamma_prior_scale", "m2_sigma_state_prior_scale",
      "analysis_objects_path"
    ),
    value = c(
      "lightweight_dic_only_refit",
      as.character(n_chains),
      as.character(n_adapt),
      as.character(n_burn),
      as.character(n_iter),
      as.character(seed),
      as.character(m2_sigma_prior_center),
      as.character(m2_gamma0_prior_scale),
      as.character(m2_gamma_prior_scale),
      as.character(m2_sigma_state_prior_scale),
      analysis_objects_path
    )
  )
  readr::write_csv(settings, settings_path)

  timestamp_message(
    "Starting lightweight DIC-only refit. Settings: chains=", n_chains,
    ", adapt=", n_adapt,
    ", burn=", n_burn,
    ", DIC iter=", n_iter,
    ". Production posterior draw files will not be overwritten."
  )

  if (dry_run) {
    timestamp_message(
      "DIC_DRY_RUN=true; validated data/model inputs and wrote settings to ",
      settings_path, ". No JAGS models were fitted."
    )
    return(invisible(list(settings = settings, model_specs = model_specs)))
  }

  dic_rows <- list()
  runtime_rows <- list()
  log_lines <- c(
    "DIC calculation from rjags::dic.samples(type = \"pD\").",
    paste("Workflow:", "lightweight DIC-only refit"),
    paste("rjags version:", as.character(utils::packageVersion("rjags"))),
    paste("DIC iterations per model:", n_iter),
    paste("Adaptation iterations per model:", n_adapt),
    paste("Warm-up iterations per model:", n_burn),
    paste("Chains per model:", n_chains),
    paste("Seed:", seed),
    "Production posterior summaries, diagnostics, predictions, and PPC use output/mcmc/*.rds from scripts/05_fit_models.R.",
    "This DIC run creates live JAGS objects only for rjags::dic.samples(); it does not save posterior parameter chains.",
    "Dbar is posterior mean deviance; pD is the effective parameter count; DIC = Dbar + pD.",
    "All models use the same frozen training observations.",
    ""
  )

  for (index in seq_along(model_specs)) {
    spec <- model_specs[[index]]
    timestamp_message("Fitting live ", spec$label, " object for DIC.")
    started_at <- Sys.time()
    fit <- fit_live_dic_model(
      spec = spec,
      n_chains = n_chains,
      n_adapt = n_adapt,
      n_burn = n_burn,
      seed = seed + 10000L * index
    )
    fit_seconds <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))

    timestamp_message("Sampling DIC for ", spec$label, ".")
    dic_started_at <- Sys.time()
    dic_object <- sample_dic(fit$model, n_iter = n_iter)
    dic_seconds <- as.numeric(difftime(Sys.time(), dic_started_at, units = "secs"))

    dic_rows[[spec$key]] <- summarise_dic_object(spec$label, dic_object)
    runtime_rows[[spec$key]] <- tibble::tibble(
      model = spec$key,
      label = spec$label,
      fit_seconds = fit_seconds,
      dic_seconds = dic_seconds,
      total_seconds = fit_seconds + dic_seconds
    )
    partial_table <- do.call(rbind, dic_rows)
    readr::write_csv(partial_table, partial_path)

    log_lines <- c(
      log_lines,
      paste0("Printed DIC object for ", spec$label, ":"),
      utils::capture.output(print(dic_object)),
      ""
    )

    timestamp_message(
      "Finished ", spec$label, " DIC stage in ",
      round(fit_seconds + dic_seconds, 1),
      " seconds. Partial DIC table: ", partial_path, "."
    )

    rm(fit, dic_object)
    invisible(gc())
  }

  dic_table <- do.call(rbind, dic_rows)
  dic_table <- dic_table[order(dic_table$DIC), , drop = FALSE]
  dic_table$delta_DIC <- dic_table$DIC - min(dic_table$DIC)
  readr::write_csv(dic_table, table_path)

  runtime_table <- do.call(rbind, runtime_rows)
  readr::write_csv(runtime_table, runtime_path)

  writeLines(log_lines, con = log_path, useBytes = TRUE)
  timestamp_message(
    "Lightweight DIC complete. Saved final table to ", table_path,
    " and log to ", log_path, "."
  )

  invisible(list(
    table = dic_table,
    runtime = runtime_table,
    settings = settings,
    table_path = table_path,
    log_path = log_path
  ))
}

if (sys.nframe() == 0L) {
  run_lightweight_dic()
}
