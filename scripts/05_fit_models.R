# Fit the M0, M1, and centered M2 models on the frozen training sample.
#
# This production fitting stage is checkpointed: each model's draws are saved
# immediately after that model finishes, before DIC is attempted.  Interrupted
# M2 or DIC runs therefore do not erase completed M0/M1 work.

required_fit_packages <- c("coda", "here", "rjags")
missing_fit_packages <- required_fit_packages[!vapply(
  required_fit_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_fit_packages) > 0L) {
  stop(
    "Install the package(s) required for model fitting: ",
    paste(missing_fit_packages, collapse = ", "),
    call. = FALSE
  )
}

source(here::here("scripts", "00_setup.R"))

read_positive_integer <- function(name, default) {
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

read_flag <- function(name, default = FALSE) {
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

read_model_selection <- function(name, default = c("M0", "M1", "M2")) {
  value <- toupper(trimws(Sys.getenv(name, unset = "")))
  if (!nzchar(value) || value %in% c("ALL", "*")) {
    return(default)
  }

  requested <- unique(strsplit(gsub("\\s+", "", value), "[,;]+")[[1L]])
  requested <- requested[nzchar(requested)]
  valid <- c("M0", "M1", "M2")
  invalid <- setdiff(requested, valid)
  if (length(invalid) > 0L || length(requested) == 0L) {
    stop(
      name, " must be ALL or a comma-separated subset of M0, M1, M2.",
      call. = FALSE
    )
  }
  valid[valid %in% requested]
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

fit_settings <- list(
  n_chains = read_positive_integer("FIT_N_CHAINS", 4L),
  n_adapt = read_positive_integer("FIT_N_ADAPT", 2000L),
  n_burn = read_positive_integer("FIT_N_BURN", 5000L),
  n_iter = read_positive_integer("FIT_N_ITER", 10000L),
  thin = read_positive_integer("FIT_THIN", 1L),
  seed = read_positive_integer("FIT_SEED", 20260621L)
)
selected_models <- read_model_selection("FIT_MODELS")
force_all <- read_flag("FIT_FORCE", default = FALSE)
run_dic <- read_flag("FIT_RUN_DIC", default = TRUE)
dic_iter <- read_positive_integer("FIT_DIC_ITER", 5000L)

fit_jags_model <- function(model_file, data_list, parameters,
                           n_chains = 4L, n_adapt = 2000L,
                           n_burn = 5000L, n_iter = 10000L,
                           thin = 1L, seed = 20260621L) {
  if (!file.exists(model_file)) {
    stop("JAGS model file not found: ", model_file, call. = FALSE)
  }
  if (!is.character(parameters) || length(parameters) == 0L) {
    stop("Provide at least one parameter to monitor.", call. = FALSE)
  }
  if (any(c(n_chains, n_adapt, n_burn, n_iter, thin) < 1L) ||
      any(c(n_chains, n_adapt, n_burn, n_iter, thin) !=
        as.integer(c(n_chains, n_adapt, n_burn, n_iter, thin)))) {
    stop("MCMC settings must be positive integers.", call. = FALSE)
  }

  required_data <- c("N", "P", "y", "X")
  missing_data <- setdiff(required_data, names(data_list))
  if (length(missing_data) > 0L) {
    stop(
      "JAGS data are missing required item(s): ",
      paste(missing_data, collapse = ", "),
      call. = FALSE
    )
  }

  set.seed(seed)
  make_inits <- function(chain_id) {
    list(
      alpha = stats::rnorm(1L, 0, 0.2),
      beta = stats::rnorm(data_list$P, 0, 0.1),
      sigma = stats::runif(1L, 0.5, 1.5),
      sd_state = stats::runif(1L, 0.1, 0.8),
      sd_zip = stats::runif(1L, 0.1, 0.8),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = as.integer(seed + 1000L * chain_id)
    )
  }

  inits <- lapply(seq_len(n_chains), make_inits)
  model_text <- paste(readLines(model_file, warn = FALSE), collapse = "\n")
  if (!grepl("sd_state", model_text, fixed = TRUE)) {
    inits <- lapply(inits, function(values) values[setdiff(names(values), "sd_state")])
  }
  if (!grepl("sd_zip", model_text, fixed = TRUE)) {
    inits <- lapply(inits, function(values) values[setdiff(names(values), "sd_zip")])
  }

  model <- rjags::jags.model(
    file = model_file,
    data = data_list,
    inits = inits,
    n.chains = n_chains,
    n.adapt = n_adapt,
    quiet = TRUE
  )
  stats::update(model, n.iter = n_burn, progress.bar = "none")
  draws <- rjags::coda.samples(
    model = model,
    variable.names = parameters,
    n.iter = n_iter,
    thin = thin,
    progress.bar = "none"
  )

  list(model = model, draws = draws, inits = inits)
}

analysis_path <- here::here("data", "processed", "analysis_objects.rds")
if (!file.exists(analysis_path)) {
  stop(
    "Prepared analysis objects not found: ", analysis_path,
    ". Run scripts/03_prepare_jags_data.R first.",
    call. = FALSE
  )
}

analysis_objects <- readRDS(analysis_path)
if (!"jags_train" %in% names(analysis_objects)) {
  stop("Analysis objects do not contain the frozen JAGS training data.", call. = FALSE)
}
jags_train_m2 <- augment_zip_centering_data(analysis_objects$jags_train)

mcmc_output_dir <- here::here("output", "mcmc")
tables_dir <- here::here("output", "tables")
dir.create(mcmc_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

model_specs <- list(
  M0 = list(
    label = "M0: pooled",
    model_file = here::here("models", "m0_pooled.jags"),
    data = analysis_objects$jags_train,
    parameters = c("alpha", "beta", "sigma"),
    output_path = file.path(mcmc_output_dir, "m0_draws.rds")
  ),
  M1 = list(
    label = "M1: state sum-to-zero",
    model_file = here::here("models", "m1_state.jags"),
    data = analysis_objects$jags_train,
    parameters = c("alpha", "beta", "sigma", "sd_state", "u"),
    output_path = file.path(mcmc_output_dir, "m1_draws.rds")
  ),
  M2 = list(
    label = "M2: centered state + within-state centered ZIP",
    model_file = here::here("models", "m2_state_zip.jags"),
    data = jags_train_m2,
    parameters = c("alpha", "beta", "sigma", "sd_state", "sd_zip", "u", "v"),
    output_path = file.path(mcmc_output_dir, "m2_draws.rds")
  )
)

timestamp_message(
  "Starting production fit stage for ", paste(selected_models, collapse = ", "),
  ". Settings: chains=", fit_settings$n_chains,
  ", adapt=", fit_settings$n_adapt,
  ", burn=", fit_settings$n_burn,
  ", iter=", fit_settings$n_iter,
  ", thin=", fit_settings$thin,
  "."
)

fit_results <- list()
fit_runtime <- list()
for (model_id in selected_models) {
  spec <- model_specs[[model_id]]
  force_model <- force_all || read_flag(paste0("FIT_FORCE_", model_id), default = FALSE)
  model_mtime <- file.info(spec$model_file)$mtime
  output_mtime <- if (file.exists(spec$output_path)) file.info(spec$output_path)$mtime else as.POSIXct(NA)
  output_current <- file.exists(spec$output_path) &&
    !is.na(output_mtime) &&
    !is.na(model_mtime) &&
    output_mtime >= model_mtime

  if (output_current && !force_model) {
    timestamp_message(
      "Skipping ", spec$label, " because current draws already exist at ",
      spec$output_path, ". Set FIT_FORCE=true or FIT_FORCE_", model_id,
      "=true to refit."
    )
    fit_runtime[[model_id]] <- data.frame(
      model = model_id,
      label = spec$label,
      status = "skipped_existing_draws",
      output_path = spec$output_path,
      elapsed_seconds = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }
  if (file.exists(spec$output_path) && !output_current && !force_model) {
    timestamp_message(
      "Refitting ", spec$label, " because ", spec$output_path,
      " is older than its model file."
    )
  }

  timestamp_message("Starting ", spec$label, ".")
  started_at <- Sys.time()
  fit <- fit_jags_model(
    model_file = spec$model_file,
    data_list = spec$data,
    parameters = spec$parameters,
    n_chains = fit_settings$n_chains,
    n_adapt = fit_settings$n_adapt,
    n_burn = fit_settings$n_burn,
    n_iter = fit_settings$n_iter,
    thin = fit_settings$thin,
    seed = fit_settings$seed
  )
  saveRDS(fit$draws, spec$output_path)
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
  timestamp_message(
    "Finished ", spec$label, " in ", round(elapsed_seconds, 1),
    " seconds. Saved draws to ", spec$output_path, "."
  )

  fit_results[[model_id]] <- fit
  fit_runtime[[model_id]] <- data.frame(
    model = model_id,
    label = spec$label,
    status = "fitted",
    output_path = spec$output_path,
    elapsed_seconds = elapsed_seconds,
    stringsAsFactors = FALSE
  )
}

runtime_table <- do.call(rbind, fit_runtime)
runtime_path <- file.path(tables_dir, "model_fit_runtime.csv")
utils::write.csv(runtime_table, runtime_path, row.names = FALSE)
timestamp_message("Saved fit runtime log to ", runtime_path, ".")

if (run_dic) {
  required_dic_models <- c("M0", "M1", "M2")
  missing_live_fits <- setdiff(required_dic_models, names(fit_results))
  if (length(missing_live_fits) == 0L) {
    timestamp_message(
      "Calculating DIC after all model draws have been saved. DIC iterations per model: ",
      dic_iter, "."
    )
    source(here::here("scripts", "07_model_comparison.R"))
    dic_results <- calculate_and_write_dic(
      fit_m0 = fit_results$M0,
      fit_m1 = fit_results$M1,
      fit_m2 = fit_results$M2,
      n_iter = dic_iter
    )
    timestamp_message(
      "Saved DIC comparison to ",
      here::here("output", "tables", "model_comparison_dic.csv"),
      "."
    )
  } else {
    timestamp_message(
      "Skipping DIC because rjags::dic.samples() requires live JAGS objects ",
      "for M0, M1, and M2 in the current R session. Missing live fit(s): ",
      paste(missing_live_fits, collapse = ", "),
      ". Draw files are saved/available. To force a same-session DIC run, ",
      "set FIT_MODELS=ALL and FIT_FORCE=true."
    )
  }
} else {
  timestamp_message("Skipping DIC because FIT_RUN_DIC=false.")
}

timestamp_message("Production fit stage complete.")
