# Planned sensitivity analyses (implementation guide point 18).
#
# Run this only after the production M2 fit has converged and the normal-model
# posterior predictive checks have been reviewed.  The default workflow refits
# M2 with wider, proper geographic-SD priors and compares posterior and
# held-out conclusions with the production fit.  The lot-size and Student-t
# branches are deliberately gated by the conditions specified in the guide.

required_sensitivity_packages <- c("coda", "here", "posterior", "readr", "rjags")
missing_sensitivity_packages <- required_sensitivity_packages[!vapply(
  required_sensitivity_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_sensitivity_packages) > 0L) {
  stop(
    "Install the package(s) required for planned sensitivity analyses: ",
    paste(missing_sensitivity_packages, collapse = ", "),
    call. = FALSE
  )
}

# Reuse the held-out metric calculation so the sensitivity conclusion uses the
# same RMSE, MAE, coverage, and stable log predictive density definitions.
comparison_script <- here::here("scripts", "07_model_comparison.R")
if (!file.exists(comparison_script)) {
  stop("Held-out comparison helper not found: ", comparison_script, call. = FALSE)
}
source(comparison_script, local = FALSE)

validate_sensitivity_integer <- function(value, name) {
  valid <- is.numeric(value) &&
    length(value) == 1L &&
    is.finite(value) &&
    value >= 1L &&
    value %% 1 == 0
  if (!valid) {
    stop(name, " must be one positive integer.", call. = FALSE)
  }
  as.integer(value)
}

validate_mcmc_settings <- function(n_chains, n_adapt, n_burn, n_iter, thin) {
  settings <- c(n_chains, n_adapt, n_burn, n_iter, thin)
  setting_names <- c("n_chains", "n_adapt", "n_burn", "n_iter", "thin")
  if (any(!is.finite(settings)) || any(settings < 1L) ||
      any(settings != as.integer(settings))) {
    stop(
      paste(setting_names, collapse = ", "),
      " must all be positive integers.",
      call. = FALSE
    )
  }
  stats::setNames(as.integer(settings), setting_names)
}

validate_m2_analysis_objects <- function(analysis_objects) {
  required_objects <- c(
    "jags_train", "X_test", "test", "predictor_names", "state_levels", "zip_levels"
  )
  missing_objects <- setdiff(required_objects, names(analysis_objects))
  if (length(missing_objects) > 0L) {
    stop(
      "Analysis objects are missing: ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

make_m2_initial_values <- function(data_list, n_chains, seed) {
  lapply(seq_len(n_chains), function(chain_id) {
    list(
      alpha = stats::rnorm(1L, 0, 0.2),
      beta = stats::rnorm(data_list$P, 0, 0.1),
      sigma = stats::runif(1L, 0.5, 1.5),
      sd_state = stats::runif(1L, 0.1, 0.8),
      sd_zip = stats::runif(1L, 0.1, 0.8),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = as.integer(seed + 1000L * chain_id)
    )
  })
}

fit_m2_sensitivity_model <- function(model_file, data_list, parameters,
                                     n_chains = 4L, n_adapt = 2000L,
                                     n_burn = 5000L, n_iter = 10000L,
                                     thin = 1L, seed = 20260622L) {
  if (!file.exists(model_file)) {
    stop("JAGS model file not found: ", model_file, call. = FALSE)
  }
  required_data <- c("N", "P", "S", "Z", "y", "X", "state", "zip")
  missing_data <- setdiff(required_data, names(data_list))
  if (length(missing_data) > 0L) {
    stop(
      "JAGS data are missing: ", paste(missing_data, collapse = ", "),
      call. = FALSE
    )
  }
  settings <- validate_mcmc_settings(n_chains, n_adapt, n_burn, n_iter, thin)
  seed <- validate_sensitivity_integer(seed, "seed")
  if (!is.character(parameters) || length(parameters) == 0L) {
    stop("Provide at least one monitored parameter.", call. = FALSE)
  }

  set.seed(seed)
  inits <- make_m2_initial_values(data_list, settings[["n_chains"]], seed)
  model <- rjags::jags.model(
    file = model_file,
    data = data_list,
    inits = inits,
    n.chains = settings[["n_chains"]],
    n.adapt = settings[["n_adapt"]],
    quiet = TRUE
  )
  stats::update(model, n.iter = settings[["n_burn"]], progress.bar = "none")
  draws <- rjags::coda.samples(
    model = model,
    variable.names = parameters,
    n.iter = settings[["n_iter"]],
    thin = settings[["thin"]],
    progress.bar = "none"
  )
  list(model = model, draws = draws, inits = inits, settings = settings)
}

summarise_convergence <- function(draws) {
  draws_array <- posterior::as_draws_array(draws)
  if (posterior::ndraws(draws_array) == 0L || posterior::nchains(draws_array) < 2L) {
    stop("Sensitivity diagnostics require non-empty draws from at least two chains.", call. = FALSE)
  }
  summary <- posterior::summarise_draws(
    draws_array,
    mean,
    median,
    sd,
    posterior::rhat,
    posterior::ess_bulk,
    posterior::ess_tail,
    posterior::mcse_mean
  )
  summary <- as.data.frame(summary)
  names(summary) <- sub("^posterior::", "", names(summary))
  summary$mcse_mean_sd_ratio <- summary$mcse_mean / summary$sd
  summary$rhat_pass <- is.finite(summary$rhat) & summary$rhat < 1.01
  summary$bulk_ess_pass <- is.finite(summary$ess_bulk) & summary$ess_bulk > 400
  summary$tail_ess_pass <- is.finite(summary$ess_tail) & summary$ess_tail > 400
  summary$mcse_pass <- is.finite(summary$mcse_mean_sd_ratio) &
    summary$mcse_mean_sd_ratio < 0.05
  summary$converged <- summary$rhat_pass & summary$bulk_ess_pass &
    summary$tail_ess_pass & summary$mcse_pass
  summary
}

summarise_sensitivity_parameters <- function(draw_matrix, analysis_objects,
                                             fit_label) {
  draw_matrix <- as.matrix(draw_matrix)
  beta_parameters <- indexed_parameter_names(
    "beta", length(analysis_objects$predictor_names)
  )
  parameters <- c(beta_parameters, "sd_state", "sd_zip")
  missing_parameters <- setdiff(parameters, colnames(draw_matrix))
  if (length(missing_parameters) > 0L) {
    stop(
      fit_label, " draws are missing: ", paste(missing_parameters, collapse = ", "),
      call. = FALSE
    )
  }

  labels <- c(analysis_objects$predictor_names, "State SD", "ZIP SD")
  data.frame(
    fit = fit_label,
    parameter = parameters,
    label = labels,
    posterior_median = vapply(
      parameters, function(parameter) stats::median(draw_matrix[, parameter]), numeric(1)
    ),
    lower_95 = vapply(
      parameters,
      function(parameter) stats::quantile(draw_matrix[, parameter], 0.025, names = FALSE),
      numeric(1)
    ),
    upper_95 = vapply(
      parameters,
      function(parameter) stats::quantile(draw_matrix[, parameter], 0.975, names = FALSE),
      numeric(1)
    ),
    row.names = NULL
  )
}

compare_sensitivity_parameters <- function(production_summary, sensitivity_summary) {
  production <- production_summary[, c("parameter", "label", "posterior_median")]
  sensitivity <- sensitivity_summary[, c("parameter", "posterior_median")]
  names(production)[3L] <- "production_posterior_median"
  names(sensitivity)[2L] <- "wider_prior_posterior_median"
  comparison <- merge(production, sensitivity, by = "parameter", sort = FALSE)
  comparison$median_difference_wider_minus_production <-
    comparison$wider_prior_posterior_median - comparison$production_posterior_median
  comparison
}

evaluate_m2_prior_sensitivity <- function(production_draws, sensitivity_draws,
                                          analysis_objects,
                                          max_predictive_draws = 2000L,
                                          seed = 20260622L) {
  production_draws <- as.matrix(posterior::as_draws_matrix(production_draws))
  sensitivity_draws <- as.matrix(posterior::as_draws_matrix(sensitivity_draws))
  max_predictive_draws <- validate_sensitivity_integer(
    max_predictive_draws, "max_predictive_draws"
  )
  seed <- validate_sensitivity_integer(seed, "seed")
  common_available_draws <- min(nrow(production_draws), nrow(sensitivity_draws))
  n_draws <- min(max_predictive_draws, common_available_draws)
  if (n_draws < 1L) {
    stop("Both M2 fits need at least one posterior draw.", call. = FALSE)
  }
  set.seed(seed)
  draw_ids <- sort(sample(seq_len(common_available_draws), size = n_draws))
  n_test <- nrow(as.matrix(analysis_objects$X_test))
  standard_normal_draws <- matrix(
    stats::rnorm(n_test * n_draws), nrow = n_test, ncol = n_draws
  )

  production <- evaluate_model_on_holdout(
    model = "M2",
    draw_matrix = production_draws,
    analysis_objects = analysis_objects,
    draw_ids = draw_ids,
    standard_normal_draws = standard_normal_draws
  )
  sensitivity <- evaluate_model_on_holdout(
    model = "M2",
    draw_matrix = sensitivity_draws,
    analysis_objects = analysis_objects,
    draw_ids = draw_ids,
    standard_normal_draws = standard_normal_draws
  )
  production$metrics$fit <- "Production M2"
  sensitivity$metrics$fit <- "M2 with wider geographic-SD priors"
  production$predictions$fit <- production$metrics$fit
  sensitivity$predictions$fit <- sensitivity$metrics$fit

  list(
    metrics = rbind(production$metrics, sensitivity$metrics),
    predictions = rbind(production$predictions, sensitivity$predictions),
    draw_ids = draw_ids
  )
}

run_planned_sensitivity_analyses <- function(
    analysis_objects_path = here::here("data", "processed", "analysis_objects.rds"),
    production_draws_path = here::here("output", "mcmc", "m2_draws.rds"),
    model_file = here::here("models", "m2_wide_sd_priors.jags"),
    mcmc_output_dir = here::here("output", "mcmc"),
    tables_dir = here::here("output", "tables"),
    n_chains = 4L, n_adapt = 2000L, n_burn = 5000L, n_iter = 10000L,
    thin = 1L, max_predictive_draws = 2000L, seed = 20260622L) {
  if (!file.exists(analysis_objects_path)) {
    stop(
      "Prepared analysis objects not found: ", analysis_objects_path,
      ". Run scripts/03_prepare_jags_data.R first.",
      call. = FALSE
    )
  }
  if (!file.exists(production_draws_path)) {
    stop(
      "Converged production M2 draws not found: ", production_draws_path,
      ". Run scripts/05_fit_models.R and scripts/06_diagnostics.R first.",
      call. = FALSE
    )
  }
  analysis_objects <- readRDS(analysis_objects_path)
  validate_m2_analysis_objects(analysis_objects)
  production_draws <- readRDS(production_draws_path)
  parameters <- c("alpha", "beta", "sigma", "sd_state", "sd_zip", "u", "v")
  sensitivity_fit <- fit_m2_sensitivity_model(
    model_file = model_file,
    data_list = analysis_objects$jags_train,
    parameters = parameters,
    n_chains = n_chains,
    n_adapt = n_adapt,
    n_burn = n_burn,
    n_iter = n_iter,
    thin = thin,
    seed = seed
  )

  dir.create(mcmc_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  sensitivity_draws_path <- file.path(mcmc_output_dir, "m2_wide_sd_prior_draws.rds")
  saveRDS(sensitivity_fit$draws, sensitivity_draws_path)

  convergence <- summarise_convergence(sensitivity_fit$draws)
  convergence_path <- file.path(tables_dir, "m2_wide_sd_prior_diagnostics.csv")
  readr::write_csv(convergence, convergence_path)
  if (!all(convergence$converged)) {
    failed_parameters <- convergence$variable[!convergence$converged]
    stop(
      "The wider-prior M2 fit did not meet convergence criteria. Investigate: ",
      paste(failed_parameters, collapse = ", "),
      ". Diagnostics saved to ", convergence_path,
      call. = FALSE
    )
  }

  production_matrix <- as.matrix(posterior::as_draws_matrix(production_draws))
  sensitivity_matrix <- as.matrix(posterior::as_draws_matrix(sensitivity_fit$draws))
  production_summary <- summarise_sensitivity_parameters(
    production_matrix, analysis_objects, "Production M2"
  )
  sensitivity_summary <- summarise_sensitivity_parameters(
    sensitivity_matrix, analysis_objects, "M2 with wider geographic-SD priors"
  )
  parameter_summary <- rbind(production_summary, sensitivity_summary)
  parameter_comparison <- compare_sensitivity_parameters(
    production_summary, sensitivity_summary
  )
  performance <- evaluate_m2_prior_sensitivity(
    production_draws = production_draws,
    sensitivity_draws = sensitivity_fit$draws,
    analysis_objects = analysis_objects,
    max_predictive_draws = max_predictive_draws,
    seed = seed
  )

  parameter_summary_path <- file.path(tables_dir, "m2_prior_sensitivity_parameter_summary.csv")
  parameter_comparison_path <- file.path(
    tables_dir, "m2_prior_sensitivity_parameter_comparison.csv"
  )
  predictive_metrics_path <- file.path(tables_dir, "m2_prior_sensitivity_predictive_metrics.csv")
  predictive_predictions_path <- file.path(
    tables_dir, "m2_prior_sensitivity_predictive_predictions.csv"
  )
  readr::write_csv(parameter_summary, parameter_summary_path)
  readr::write_csv(parameter_comparison, parameter_comparison_path)
  readr::write_csv(performance$metrics, predictive_metrics_path)
  readr::write_csv(performance$predictions, predictive_predictions_path)

  lot_size_status <- if ("log_acre_lot_z" %in% analysis_objects$predictor_names) {
    "Not run: log acre lot is already in the core M0-M2 design matrix."
  } else {
    paste(
      "Not run automatically: the core model omitted log acre lot, so a separately",
      "labelled complete-case M2 is required before comparing its changed target population."
    )
  }
  student_t_status <- paste(
    "Not run: no clear normal-residual posterior-predictive failure was confirmed.",
    "models/m2_student_t_residual.jags is available only if that evidence is",
    "documented; then refit on the frozen sample and rerun diagnostics, held-out",
    "evaluation, and posterior predictive checks."
  )
  decision_log <- data.frame(
    analysis = c("Prior sensitivity", "Lot-size complete-case sensitivity", "Student-t residual extension"),
    status = c(
      "Run: M2 with wider proper geographic-SD priors (half-normal scale 0.6 vs 0.3).",
      lot_size_status,
      student_t_status
    ),
    stringsAsFactors = FALSE
  )
  decision_log_path <- file.path(tables_dir, "planned_sensitivity_analysis_log.csv")
  readr::write_csv(decision_log, decision_log_path)

  message(
    "Saved converged wider-prior M2 draws to ", sensitivity_draws_path,
    " and prior-sensitivity comparisons to ", tables_dir, "."
  )
  invisible(list(
    sensitivity_draws_path = sensitivity_draws_path,
    convergence = convergence,
    parameter_summary = parameter_summary,
    parameter_comparison = parameter_comparison,
    predictive_metrics = performance$metrics,
    decision_log = decision_log
  ))
}

if (sys.nframe() == 0L) {
  run_planned_sensitivity_analyses()
}
