# Experimental fixed-nu Student-t residual M2 caveat check.
#
# This script keeps the residual degrees of freedom fixed and uses the
# normal-gamma scale-mixture representation of a Student-t likelihood. It tests
# the same heavy-tailed residual idea as `scripts/14_m2_student_t_experiment.R`,
# but avoids the very slow unknown-nu JAGS updates seen in smoke testing.

source(here::here("scripts", "14_m2_student_t_experiment.R"))

student_t_fixed_read_positive_number <- function(name, default) {
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

student_t_fixed_format_number <- function(value) {
  gsub("\\.", "p", format(value, trim = TRUE, scientific = FALSE))
}

student_t_fixed_default_label <- function(nu, geo_sd_prior_scale) {
  paste0(
    "nu", student_t_fixed_format_number(nu),
    "_geo", student_t_fixed_format_number(geo_sd_prior_scale)
  )
}

student_t_fixed_read_label <- function(name, nu, geo_sd_prior_scale) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) {
    return(student_t_fixed_default_label(nu, geo_sd_prior_scale))
  }
  if (!grepl("^[A-Za-z0-9_-]+$", value)) {
    stop(
      "Environment variable `", name,
      "` may contain only letters, numbers, underscores, and hyphens.",
      call. = FALSE
    )
  }
  value
}

make_student_t_fixed_inits <- function(data_list, n_chains, seed) {
  set.seed(seed)
  lapply(seq_len(n_chains), function(chain_id) {
    list(
      alpha = stats::rnorm(1L, 0, 0.2),
      beta = stats::rnorm(data_list$P, 0, 0.1),
      sigma = stats::runif(1L, 0.4, 1.2),
      sd_state = stats::runif(1L, 0.1, 0.8),
      sd_zip = stats::runif(1L, 0.1, 0.8),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = as.integer(seed + 1000L * chain_id)
    )
  })
}

student_t_fixed_required_parameters <- function(n_predictors, n_states, n_zips) {
  c(
    "alpha", "sigma", "sd_state", "sd_zip",
    student_t_indexed_parameter_names("beta", n_predictors),
    student_t_indexed_parameter_names("u", n_states),
    student_t_indexed_parameter_names("v", n_zips)
  )
}

validate_student_t_fixed_draw_matrix <- function(draw_matrix, analysis_objects) {
  X <- as.matrix(analysis_objects$X_train)
  draw_matrix <- as.matrix(draw_matrix)
  if (nrow(draw_matrix) == 0L || ncol(draw_matrix) == 0L ||
      is.null(colnames(draw_matrix)) || any(!is.finite(draw_matrix))) {
    stop("Fixed-nu Student-t M2 posterior draws are empty, unnamed, or non-finite.", call. = FALSE)
  }
  required_parameters <- student_t_fixed_required_parameters(
    n_predictors = ncol(X),
    n_states = length(analysis_objects$state_levels),
    n_zips = length(analysis_objects$zip_levels)
  )
  missing_parameters <- setdiff(required_parameters, colnames(draw_matrix))
  if (length(missing_parameters) > 0L) {
    stop(
      "Fixed-nu Student-t M2 draws are missing parameter(s): ",
      paste(missing_parameters, collapse = ", "),
      call. = FALSE
    )
  }
  if (any(draw_matrix[, "sigma"] <= 0) ||
      any(draw_matrix[, "sd_state"] <= 0) ||
      any(draw_matrix[, "sd_zip"] <= 0)) {
    stop("Fixed-nu Student-t M2 draws contain invalid scale values.", call. = FALSE)
  }
  invisible(TRUE)
}

fit_student_t_fixed_m2 <- function(model_file, data_list, parameters, settings,
                                   partial_draws_path = NULL) {
  if (!file.exists(model_file)) {
    stop("Fixed-nu Student-t M2 model file not found: ", model_file, call. = FALSE)
  }
  data_list <- augment_zip_centering_data(data_list)
  data_list$nu <- settings$nu
  data_list$tau_state_prior <- 1 / (settings$geo_sd_prior_scale^2)
  data_list$tau_zip_prior <- 1 / (settings$geo_sd_prior_scale^2)
  inits <- make_student_t_fixed_inits(data_list, settings$n_chains, settings$seed)

  student_t_timestamp_message("Compiling fixed-nu Student-t M2 with nu=", settings$nu, ".")
  model <- rjags::jags.model(
    file = model_file,
    data = data_list,
    inits = inits,
    n.chains = settings$n_chains,
    n.adapt = 0L,
    quiet = TRUE
  )
  student_t_timestamp_message("Compiled fixed-nu Student-t M2.")
  run_student_t_adapt_chunks(model, settings$n_adapt, settings$chunk_iter)
  run_student_t_update_chunks(model, settings$n_burn, settings$chunk_iter, "burn-in")
  draws <- sample_student_t_coda_chunks(
    model = model,
    parameters = parameters,
    total_iter = settings$n_iter,
    thin = settings$thin,
    chunk_iter = settings$chunk_iter,
    partial_draws_path = partial_draws_path
  )
  list(model = model, draws = draws, inits = inits)
}

add_fixed_nu_column <- function(draw_matrix, nu) {
  draw_matrix <- as.matrix(draw_matrix)
  if ("nu" %in% colnames(draw_matrix)) {
    draw_matrix[, "nu"] <- nu
    return(draw_matrix)
  }
  cbind(draw_matrix, nu = rep(nu, nrow(draw_matrix)))
}

write_student_t_fixed_comparison_tables <- function(heldout_metrics, ppc_summary,
                                                    tables_dir, output_stem) {
  production_metrics_path <- file.path(tables_dir, "heldout_predictive_metrics.csv")
  if (file.exists(production_metrics_path)) {
    production_metrics <- readr::read_csv(production_metrics_path, show_col_types = FALSE)
    production_m2 <- production_metrics[production_metrics$model == "M2", , drop = FALSE]
    if (nrow(production_m2) == 1L) {
      readr::write_csv(
        rbind(production_m2, heldout_metrics),
        file.path(tables_dir, paste0(output_stem, "_vs_m2_heldout_comparison.csv"))
      )
    }
  }

  production_ppc_path <- file.path(tables_dir, "posterior_predictive_check_summary.csv")
  if (file.exists(production_ppc_path)) {
    production_ppc <- readr::read_csv(production_ppc_path, show_col_types = FALSE)
    production_m2 <- production_ppc[production_ppc$model == "M2", , drop = FALSE]
    key_statistics <- c(
      "overall_iqr_standardized_log_price",
      "overall_q90_standardized_log_price",
      "zip_mean_variance_standardized_log_price"
    )
    comparison <- rbind(
      production_m2[production_m2$statistic %in% key_statistics, ],
      ppc_summary[ppc_summary$statistic %in% key_statistics, ]
    )
    readr::write_csv(
      comparison,
      file.path(tables_dir, paste0(output_stem, "_vs_m2_ppc_comparison.csv"))
    )
  }
}

run_student_t_fixed_m2_experiment <- function(
    analysis_objects_path = here::here("data", "processed", "analysis_objects.rds"),
    model_file = here::here("models", "m2_student_t_fixed_nu_mixture.jags"),
    mcmc_output_dir = here::here("output", "mcmc"),
    tables_dir = here::here("output", "tables"),
    figures_dir = here::here("output", "figures")) {
  fixed_nu <- student_t_fixed_read_positive_number("STUDENT_T_FIXED_NU_M2_NU", 4)
  geo_sd_prior_scale <- student_t_fixed_read_positive_number(
    "STUDENT_T_FIXED_NU_M2_GEO_SD_PRIOR_SCALE",
    0.3
  )
  settings <- list(
    n_chains = student_t_read_positive_integer("STUDENT_T_FIXED_NU_M2_N_CHAINS", 4L),
    n_adapt = student_t_read_positive_integer("STUDENT_T_FIXED_NU_M2_N_ADAPT", 2000L),
    n_burn = student_t_read_positive_integer("STUDENT_T_FIXED_NU_M2_N_BURN", 5000L),
    n_iter = student_t_read_positive_integer("STUDENT_T_FIXED_NU_M2_N_ITER", 10000L),
    thin = student_t_read_positive_integer("STUDENT_T_FIXED_NU_M2_THIN", 1L),
    seed = student_t_read_positive_integer("STUDENT_T_FIXED_NU_M2_SEED", 20260702L),
    chunk_iter = student_t_read_positive_integer("STUDENT_T_FIXED_NU_M2_CHUNK_ITER", 1000L),
    max_heldout_draws = student_t_read_positive_integer("STUDENT_T_FIXED_NU_M2_MAX_HELDOUT_DRAWS", 2000L),
    max_ppc_draws = student_t_read_positive_integer("STUDENT_T_FIXED_NU_M2_MAX_PPC_DRAWS", 1000L),
    force = student_t_read_flag("STUDENT_T_FIXED_NU_M2_FORCE", default = FALSE),
    nu = fixed_nu,
    geo_sd_prior_scale = geo_sd_prior_scale,
    label = student_t_fixed_read_label(
      "STUDENT_T_FIXED_NU_M2_LABEL",
      fixed_nu,
      geo_sd_prior_scale
    )
  )
  if (settings$nu <= 2) {
    stop("STUDENT_T_FIXED_NU_M2_NU must be greater than 2 for finite residual variance.", call. = FALSE)
  }

  if (!file.exists(analysis_objects_path)) {
    stop(
      "Prepared analysis objects not found: ", analysis_objects_path,
      ". Run scripts/03_prepare_jags_data.R first.",
      call. = FALSE
    )
  }
  dir.create(mcmc_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  analysis_objects <- readRDS(analysis_objects_path)
  validate_student_t_analysis_objects(analysis_objects)
  parameters <- c("alpha", "beta", "sigma", "sd_state", "sd_zip", "u", "v")
  output_stem <- paste0("m2_student_t_fixednu_", settings$label)
  model_label <- paste0("M2_student_t_fixednu_", settings$label)
  draws_path <- file.path(mcmc_output_dir, paste0(output_stem, "_draws.rds"))
  partial_draws_path <- file.path(mcmc_output_dir, paste0(output_stem, "_draw_chunks_partial.rds"))
  runtime_path <- file.path(tables_dir, paste0(output_stem, "_runtime.csv"))

  readr::write_csv(
    data.frame(setting = names(settings), value = unlist(settings), row.names = NULL),
    file.path(tables_dir, paste0(output_stem, "_settings.csv"))
  )
  readr::write_csv(
    data.frame(
      residual_distribution = "Student-t scale mixture",
      nu = settings$nu,
      geo_sd_prior_scale = settings$geo_sd_prior_scale,
      label = settings$label,
      stringsAsFactors = FALSE
    ),
    file.path(tables_dir, paste0(output_stem, "_terms.csv"))
  )

  freshness_mtime <- max(
    file.info(c(model_file, here::here("scripts", "15_m2_student_t_fixed_nu_experiment.R")))$mtime,
    na.rm = TRUE
  )
  draws_current <- file.exists(draws_path) &&
    file.info(draws_path)$mtime >= freshness_mtime

  fit_started_at <- Sys.time()
  if (draws_current && !settings$force) {
    student_t_timestamp_message(
      "Skipping fixed-nu Student-t M2 refit because current draws already exist at ",
      draws_path, ". Set STUDENT_T_FIXED_NU_M2_FORCE=true to refit."
    )
    draws <- readRDS(draws_path)
    fit_status <- "skipped_existing_draws"
    fit_seconds <- NA_real_
  } else {
    student_t_timestamp_message(
      "Starting fixed-nu Student-t M2 fit. Settings: chains=", settings$n_chains,
      ", adapt=", settings$n_adapt,
      ", burn=", settings$n_burn,
      ", iter=", settings$n_iter,
      ", thin=", settings$thin,
      ", chunk=", settings$chunk_iter,
      ", nu=", settings$nu,
      ", geo_sd_prior_scale=", settings$geo_sd_prior_scale,
      "."
    )
    fit <- fit_student_t_fixed_m2(
      model_file = model_file,
      data_list = analysis_objects$jags_train,
      parameters = parameters,
      settings = settings,
      partial_draws_path = partial_draws_path
    )
    draws <- fit$draws
    saveRDS(draws, draws_path)
    fit_seconds <- as.numeric(difftime(Sys.time(), fit_started_at, units = "secs"))
    fit_status <- "fitted"
    student_t_timestamp_message(
      "Saved fixed-nu Student-t M2 draws to ", draws_path,
      " after ", round(fit_seconds, 1), " seconds."
    )
  }

  draw_matrix <- as.matrix(posterior::as_draws_matrix(draws))
  validate_student_t_fixed_draw_matrix(draw_matrix, analysis_objects)

  convergence_started_at <- Sys.time()
  convergence <- summarise_student_t_convergence(draws)
  convergence_path <- file.path(tables_dir, paste0(output_stem, "_posterior_summary.csv"))
  readr::write_csv(convergence, convergence_path)
  convergence_seconds <- as.numeric(difftime(Sys.time(), convergence_started_at, units = "secs"))
  student_t_timestamp_message(
    "Saved fixed-nu Student-t M2 convergence summary to ", convergence_path,
    ". Converged ", sum(convergence$converged), "/", nrow(convergence), " parameters."
  )

  eval_draw_matrix <- add_fixed_nu_column(draw_matrix, settings$nu)
  set.seed(settings$seed + 101L)
  heldout_draw_ids <- sort(sample(
    seq_len(nrow(eval_draw_matrix)),
    size = min(settings$max_heldout_draws, nrow(eval_draw_matrix)),
    replace = FALSE
  ))
  heldout <- evaluate_student_t_holdout(
    draw_matrix = eval_draw_matrix,
    analysis_objects = analysis_objects,
    draw_ids = heldout_draw_ids
  )
  heldout$metrics$model <- model_label
  heldout$predictions$model <- model_label
  heldout_metrics_path <- file.path(tables_dir, paste0(output_stem, "_heldout_predictive_metrics.csv"))
  heldout_predictions_path <- file.path(tables_dir, paste0(output_stem, "_heldout_predictive_predictions.csv"))
  readr::write_csv(heldout$metrics, heldout_metrics_path)
  readr::write_csv(heldout$predictions, heldout_predictions_path)

  set.seed(settings$seed + 201L)
  ppc_draw_ids <- sort(sample(
    seq_len(nrow(eval_draw_matrix)),
    size = min(settings$max_ppc_draws, nrow(eval_draw_matrix)),
    replace = FALSE
  ))
  selected_state_ids <- select_ppc_state_ids(analysis_objects$state_levels)
  ppc <- run_student_t_ppc(
    draw_matrix = eval_draw_matrix,
    analysis_objects = analysis_objects,
    draw_ids = ppc_draw_ids,
    selected_state_ids = selected_state_ids
  )
  ppc$summary$model <- model_label
  ppc$replicated_statistics$model <- model_label
  ppc$observed_statistics$model <- model_label
  ppc_summary_path <- file.path(tables_dir, paste0(output_stem, "_posterior_predictive_check_summary.csv"))
  ppc_replicates_path <- file.path(tables_dir, paste0(output_stem, "_posterior_predictive_check_replicates.csv"))
  ppc_figure_path <- file.path(figures_dir, paste0(output_stem, "_posterior_predictive_checks.png"))
  readr::write_csv(ppc$summary, ppc_summary_path)
  readr::write_csv(ppc$replicated_statistics, ppc_replicates_path)
  save_ppc_plot(
    replicated_statistics = ppc$replicated_statistics,
    observed_statistics = ppc$observed_statistics,
    model = paste0("M2 Student-t residual fixed nu=", settings$nu),
    figure_path = ppc_figure_path
  )
  write_student_t_fixed_comparison_tables(
    heldout$metrics,
    ppc$summary,
    tables_dir,
    output_stem
  )

  runtime <- data.frame(
    stage = c("fit_or_load_draws", "convergence_summary", "heldout_ppc"),
    status = c(fit_status, "completed", "completed"),
    elapsed_seconds = c(fit_seconds, convergence_seconds, NA_real_),
    stringsAsFactors = FALSE
  )
  readr::write_csv(runtime, runtime_path)
  student_t_timestamp_message(
    "Fixed-nu Student-t M2 experiment complete. Key outputs: ", convergence_path,
    ", ", heldout_metrics_path, ", ", ppc_summary_path, "."
  )
  invisible(list(
    draws_path = draws_path,
    convergence = convergence,
    heldout = heldout,
    ppc = ppc
  ))
}

if (sys.nframe() == 0L) {
  run_student_t_fixed_m2_experiment()
}
