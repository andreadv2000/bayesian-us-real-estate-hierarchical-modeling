# Planned sensitivity analyses (implementation guide point 18).
#
# This script is retained as a legacy normal-residual M2 prior-sensitivity
# check.  On the caveat-fix branch, the official M2 is heteroskedastic; this
# script still compares the original normal-residual M2 with its wider-prior
# variant and should be interpreted as baseline sensitivity evidence.

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

read_sensitivity_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.integer(default))
  }

  parsed <- suppressWarnings(as.integer(value))
  validate_sensitivity_integer(parsed, name)
}

read_sensitivity_flag <- function(name, default = FALSE) {
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

chunk_sizes <- function(total_iter, chunk_iter) {
  total_iter <- validate_sensitivity_integer(total_iter, "total_iter")
  chunk_iter <- validate_sensitivity_integer(chunk_iter, "chunk_iter")
  starts <- seq.int(1L, total_iter, by = chunk_iter)
  as.integer(pmin(chunk_iter, total_iter - starts + 1L))
}

stall_hint <- function(phase, chunk_iter) {
  timestamp_message(
    "Entering ", phase, " chunk of up to ", chunk_iter,
    " JAGS iteration(s). If this message does not advance for a long time ",
    "and Rscript CPU stays at 0, the process is likely stalled inside this ",
    "JAGS call; interrupting will keep any previously saved final draw file, ",
    "but not the currently active chunk."
  )
}

run_update_chunks <- function(model, total_iter, chunk_iter, phase) {
  chunks <- chunk_sizes(total_iter, chunk_iter)
  for (chunk_id in seq_along(chunks)) {
    stall_hint(phase, chunks[[chunk_id]])
    started_at <- Sys.time()
    stats::update(model, n.iter = chunks[[chunk_id]], progress.bar = "none")
    elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
    timestamp_message(
      "Finished ", phase, " chunk ", chunk_id, "/", length(chunks),
      " (", chunks[[chunk_id]], " iteration(s)) in ",
      round(elapsed, 1), " seconds."
    )
  }
  invisible(model)
}

run_adapt_chunks <- function(model, total_iter, chunk_iter) {
  chunks <- chunk_sizes(total_iter, chunk_iter)
  for (chunk_id in seq_along(chunks)) {
    stall_hint("adaptation", chunks[[chunk_id]])
    started_at <- Sys.time()
    adapted <- rjags::adapt(
      object = model,
      n.iter = chunks[[chunk_id]],
      end.adaptation = chunk_id == length(chunks),
      progress.bar = "none"
    )
    elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
    timestamp_message(
      "Finished adaptation chunk ", chunk_id, "/", length(chunks),
      " (", chunks[[chunk_id]], " iteration(s)) in ",
      round(elapsed, 1), " seconds. Adapted=", adapted, "."
    )
  }
  invisible(model)
}

combine_mcmc_chunks <- function(draw_chunks) {
  if (length(draw_chunks) == 0L) {
    stop("No MCMC chunks were sampled.", call. = FALSE)
  }
  n_chains <- length(draw_chunks[[1L]])
  chains <- lapply(seq_len(n_chains), function(chain_id) {
    combined <- do.call(rbind, lapply(draw_chunks, function(chunk) chunk[[chain_id]]))
    if (!inherits(combined, "mcmc")) {
      combined <- coda::mcmc(combined)
    }
    combined
  })
  coda::mcmc.list(chains)
}

sample_coda_chunks <- function(model, parameters, total_iter, thin,
                               chunk_iter, partial_path = NULL) {
  chunks <- chunk_sizes(total_iter, chunk_iter)
  draw_chunks <- vector("list", length(chunks))
  for (chunk_id in seq_along(chunks)) {
    stall_hint("sampling", chunks[[chunk_id]])
    started_at <- Sys.time()
    draw_chunks[[chunk_id]] <- rjags::coda.samples(
      model = model,
      variable.names = parameters,
      n.iter = chunks[[chunk_id]],
      thin = thin,
      progress.bar = "none"
    )
    elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
    timestamp_message(
      "Finished sampling chunk ", chunk_id, "/", length(chunks),
      " (", chunks[[chunk_id]], " iteration(s)) in ",
      round(elapsed, 1), " seconds."
    )
    if (!is.null(partial_path)) {
      saveRDS(draw_chunks[seq_len(chunk_id)], partial_path)
      timestamp_message("Saved partial sensitivity draw chunks to ", partial_path, ".")
    }
  }
  combine_mcmc_chunks(draw_chunks)
}

fit_m2_sensitivity_model <- function(model_file, data_list, parameters,
                                     n_chains = 4L, n_adapt = 2000L,
                                     n_burn = 5000L, n_iter = 10000L,
                                     thin = 1L, seed = 20260622L,
                                     chunk_iter = 1000L,
                                     partial_draws_path = NULL) {
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
  chunk_iter <- validate_sensitivity_integer(chunk_iter, "chunk_iter")
  seed <- validate_sensitivity_integer(seed, "seed")
  if (!is.character(parameters) || length(parameters) == 0L) {
    stop("Provide at least one monitored parameter.", call. = FALSE)
  }
  data_list <- augment_zip_centering_data(data_list)

  set.seed(seed)
  inits <- make_m2_initial_values(data_list, settings[["n_chains"]], seed)
  timestamp_message("Compiling wider-prior M2 sensitivity JAGS model.")
  model <- rjags::jags.model(
    file = model_file,
    data = data_list,
    inits = inits,
    n.chains = settings[["n_chains"]],
    n.adapt = 0L,
    quiet = TRUE
  )
  timestamp_message("Compiled wider-prior M2 sensitivity JAGS model.")
  run_adapt_chunks(
    model = model,
    total_iter = settings[["n_adapt"]],
    chunk_iter = chunk_iter
  )
  run_update_chunks(
    model = model,
    total_iter = settings[["n_burn"]],
    chunk_iter = chunk_iter,
    phase = "burn-in"
  )
  draws <- sample_coda_chunks(
    model = model,
    parameters = parameters,
    total_iter = settings[["n_iter"]],
    thin = settings[["thin"]],
    chunk_iter = chunk_iter,
    partial_path = partial_draws_path
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

legacy_m2_linear_predictor <- function(draw, X_test, state_id, zip_id) {
  beta_names <- indexed_parameter_names("beta", ncol(X_test))
  state_names <- indexed_parameter_names("u", max(state_id))
  zip_names <- indexed_parameter_names("v", max(zip_id))
  as.vector(
    draw["alpha"] +
      X_test %*% unname(draw[beta_names]) +
      unname(draw[state_names])[state_id] +
      unname(draw[zip_names])[zip_id]
  )
}

legacy_update_log_density_accumulator <- function(log_max, scaled_sum, log_density) {
  new_log_max <- pmax(log_max, log_density)
  list(
    log_max = new_log_max,
    scaled_sum =
      exp(log_max - new_log_max) * scaled_sum +
      exp(log_density - new_log_max)
  )
}

evaluate_legacy_m2_on_holdout <- function(draw_matrix, analysis_objects,
                                          draw_ids, standard_normal_draws,
                                          fit_label) {
  draw_matrix <- as.matrix(draw_matrix)
  X_test <- as.matrix(analysis_objects$X_test)
  test_y <- as.numeric(analysis_objects$test$y)
  state_id <- as.integer(analysis_objects$test$state_id)
  zip_id <- as.integer(analysis_objects$test$zip_id_num)
  n_test <- length(test_y)
  required_parameters <- c(
    "alpha", "sigma", "sd_zip",
    indexed_parameter_names("beta", ncol(X_test)),
    indexed_parameter_names("u", length(analysis_objects$state_levels)),
    indexed_parameter_names("v", length(analysis_objects$zip_levels))
  )
  missing_parameters <- setdiff(required_parameters, colnames(draw_matrix))
  if (length(missing_parameters) > 0L) {
    stop(
      fit_label, " legacy M2 draws are missing: ",
      paste(missing_parameters, collapse = ", "),
      call. = FALSE
    )
  }

  posterior_mean_sum <- numeric(n_test)
  predictive_draws <- matrix(NA_real_, nrow = n_test, ncol = length(draw_ids))
  log_density_max <- numeric(n_test)
  log_density_sum <- numeric(n_test)
  for (draw_position in seq_along(draw_ids)) {
    draw <- as.numeric(draw_matrix[draw_ids[[draw_position]], , drop = TRUE])
    names(draw) <- colnames(draw_matrix)
    mu <- legacy_m2_linear_predictor(draw, X_test, state_id, zip_id)
    posterior_mean_sum <- posterior_mean_sum + mu
    predictive_draws[, draw_position] <-
      mu + draw["sigma"] * standard_normal_draws[, draw_position]
    log_density <- stats::dnorm(test_y, mean = mu, sd = draw["sigma"], log = TRUE)
    if (draw_position == 1L) {
      log_density_max <- log_density
      log_density_sum <- rep(1, n_test)
    } else {
      accumulator <- legacy_update_log_density_accumulator(
        log_density_max,
        log_density_sum,
        log_density
      )
      log_density_max <- accumulator$log_max
      log_density_sum <- accumulator$scaled_sum
    }
  }

  posterior_mean <- posterior_mean_sum / length(draw_ids)
  lower_95 <- apply(
    predictive_draws, 1L, stats::quantile,
    probs = 0.025, names = FALSE, type = 8L
  )
  upper_95 <- apply(
    predictive_draws, 1L, stats::quantile,
    probs = 0.975, names = FALSE, type = 8L
  )
  log_predictive_density <- log_density_max + log(log_density_sum) - log(length(draw_ids))
  covered_95 <- test_y >= lower_95 & test_y <= upper_95

  predictions <- data.frame(
    model = "legacy_normal_M2",
    fit = fit_label,
    test_row = seq_len(n_test),
    state_id = state_id,
    zip_id_num = zip_id,
    observed_standardized_log_price = test_y,
    posterior_mean_prediction = posterior_mean,
    lower_95_predictive = lower_95,
    upper_95_predictive = upper_95,
    covered_by_95_predictive_interval = covered_95,
    log_predictive_density = log_predictive_density,
    row.names = NULL
  )
  metrics <- data.frame(
    model = "legacy_normal_M2",
    fit = fit_label,
    n_test = n_test,
    n_posterior_draws = length(draw_ids),
    RMSE = sqrt(mean((posterior_mean - test_y)^2)),
    MAE = mean(abs(posterior_mean - test_y)),
    coverage_95 = mean(covered_95),
    total_lppd = sum(log_predictive_density),
    MLPD = mean(log_predictive_density),
    row.names = NULL
  )
  list(metrics = metrics, predictions = predictions)
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

  production <- evaluate_legacy_m2_on_holdout(
    draw_matrix = production_draws,
    analysis_objects = analysis_objects,
    draw_ids = draw_ids,
    standard_normal_draws = standard_normal_draws,
    fit_label = "Production legacy normal M2"
  )
  sensitivity <- evaluate_legacy_m2_on_holdout(
    draw_matrix = sensitivity_draws,
    analysis_objects = analysis_objects,
    draw_ids = draw_ids,
    standard_normal_draws = standard_normal_draws,
    fit_label = "Legacy normal M2 with wider geographic-SD priors"
  )

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
    n_chains = read_sensitivity_integer("SENS_N_CHAINS", 4L),
    n_adapt = read_sensitivity_integer("SENS_N_ADAPT", 2000L),
    n_burn = read_sensitivity_integer("SENS_N_BURN", 5000L),
    n_iter = read_sensitivity_integer("SENS_N_ITER", 10000L),
    thin = read_sensitivity_integer("SENS_THIN", 1L),
    max_predictive_draws = read_sensitivity_integer("SENS_MAX_PREDICTIVE_DRAWS", 2000L),
    seed = read_sensitivity_integer("SENS_SEED", 20260622L),
    chunk_iter = read_sensitivity_integer("SENS_CHUNK_ITER", 1000L),
    force = read_sensitivity_flag("SENS_FORCE", default = FALSE),
    strict_convergence = read_sensitivity_flag("SENS_STRICT_CONVERGENCE", default = FALSE)) {
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
  dir.create(mcmc_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  sensitivity_draws_path <- file.path(mcmc_output_dir, "m2_wide_sd_prior_draws.rds")
  partial_draws_path <- file.path(mcmc_output_dir, "m2_wide_sd_prior_draw_chunks_partial.rds")
  settings_path <- file.path(tables_dir, "m2_prior_sensitivity_settings.csv")
  runtime_path <- file.path(tables_dir, "m2_prior_sensitivity_runtime.csv")
  convergence_path <- file.path(tables_dir, "m2_wide_sd_prior_diagnostics.csv")
  convergence_status_path <- file.path(tables_dir, "m2_wide_sd_prior_convergence_status.csv")

  settings_table <- data.frame(
    setting = c(
      "n_chains", "n_adapt", "n_burn", "n_iter", "thin",
      "max_predictive_draws", "seed", "chunk_iter", "force",
      "strict_convergence"
    ),
    value = c(
      n_chains, n_adapt, n_burn, n_iter, thin,
      max_predictive_draws, seed, chunk_iter, force, strict_convergence
    ),
    stringsAsFactors = FALSE
  )
  readr::write_csv(settings_table, settings_path)
  timestamp_message("Saved sensitivity settings to ", settings_path, ".")

  analysis_objects <- readRDS(analysis_objects_path)
  validate_m2_analysis_objects(analysis_objects)
  production_draws <- readRDS(production_draws_path)
  parameters <- c("alpha", "beta", "sigma", "sd_state", "sd_zip", "u", "v")

  freshness_mtime <- max(
    file.info(c(model_file, here::here("scripts", "11_sensitivity_analyses.R")))$mtime,
    na.rm = TRUE
  )
  output_mtime <- if (file.exists(sensitivity_draws_path)) {
    file.info(sensitivity_draws_path)$mtime
  } else {
    as.POSIXct(NA)
  }
  output_current <- file.exists(sensitivity_draws_path) &&
    !is.na(output_mtime) &&
    output_mtime >= freshness_mtime

  fit_started_at <- Sys.time()
  if (output_current && !force) {
    timestamp_message(
      "Skipping wider-prior M2 refit because current sensitivity draws exist at ",
      sensitivity_draws_path, ". Set SENS_FORCE=true to refit."
    )
    sensitivity_draws <- readRDS(sensitivity_draws_path)
    fit_status <- "skipped_existing_draws"
    fit_seconds <- NA_real_
  } else {
    if (file.exists(sensitivity_draws_path) && !force) {
      timestamp_message(
        "Refitting wider-prior M2 because existing sensitivity draws are older ",
        "than the model/script."
      )
    }
    timestamp_message(
      "Starting wider-prior M2 sensitivity fit. Settings: chains=", n_chains,
      ", adapt=", n_adapt, ", burn=", n_burn, ", iter=", n_iter,
      ", thin=", thin, ", chunk=", chunk_iter, "."
    )
    timestamp_message(
      "Stall monitor: during long JAGS calls, heartbeat messages appear after ",
      "each chunk. If no heartbeat appears for much longer than a previous ",
      "similar chunk and Rscript CPU stays at 0, treat it as stalled."
    )
    sensitivity_fit <- fit_m2_sensitivity_model(
      model_file = model_file,
      data_list = analysis_objects$jags_train,
      parameters = parameters,
      n_chains = n_chains,
      n_adapt = n_adapt,
      n_burn = n_burn,
      n_iter = n_iter,
      thin = thin,
      seed = seed,
      chunk_iter = chunk_iter,
      partial_draws_path = partial_draws_path
    )
    sensitivity_draws <- sensitivity_fit$draws
    saveRDS(sensitivity_draws, sensitivity_draws_path)
    fit_seconds <- as.numeric(difftime(Sys.time(), fit_started_at, units = "secs"))
    fit_status <- "fitted"
    timestamp_message(
      "Saved wider-prior M2 sensitivity draws immediately to ",
      sensitivity_draws_path, " after ", round(fit_seconds, 1), " seconds."
    )
  }

  convergence_started_at <- Sys.time()
  timestamp_message("Summarising wider-prior M2 convergence diagnostics.")
  convergence <- summarise_convergence(sensitivity_draws)
  beta_parameters <- indexed_parameter_names("beta", length(analysis_objects$predictor_names))
  primary_parameters <- c("alpha", beta_parameters, "sigma", "sd_state", "sd_zip")
  convergence$primary_parameter <- convergence$variable %in% primary_parameters
  readr::write_csv(convergence, convergence_path)
  convergence_seconds <- as.numeric(difftime(Sys.time(), convergence_started_at, units = "secs"))

  failed_parameters <- convergence$variable[!convergence$converged]
  failed_primary_parameters <- convergence$variable[
    convergence$primary_parameter & !convergence$converged
  ]
  convergence_status <- data.frame(
    n_parameters = nrow(convergence),
    n_converged = sum(convergence$converged),
    n_failed = length(failed_parameters),
    n_primary_parameters = sum(convergence$primary_parameter),
    n_primary_failed = length(failed_primary_parameters),
    strict_convergence = strict_convergence,
    stringsAsFactors = FALSE
  )
  readr::write_csv(convergence_status, convergence_status_path)
  timestamp_message(
    "Saved sensitivity diagnostics to ", convergence_path,
    ". Converged ", convergence_status$n_converged, "/",
    convergence_status$n_parameters, " parameters; primary failures: ",
    convergence_status$n_primary_failed, "."
  )
  if (length(failed_parameters) > 0L && strict_convergence) {
    stop(
      "The wider-prior M2 fit did not meet convergence criteria. Investigate: ",
      paste(failed_parameters, collapse = ", "),
      ". Diagnostics saved to ", convergence_path,
      call. = FALSE
    )
  } else if (length(failed_parameters) > 0L) {
    timestamp_message(
      "Continuing despite non-strict sensitivity convergence failures. ",
      "Set SENS_STRICT_CONVERGENCE=true to make these failures stop the script."
    )
  }

  production_matrix <- as.matrix(posterior::as_draws_matrix(production_draws))
  sensitivity_matrix <- as.matrix(posterior::as_draws_matrix(sensitivity_draws))
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
  timestamp_message("Evaluating held-out sensitivity performance.")
  performance <- evaluate_m2_prior_sensitivity(
    production_draws = production_draws,
    sensitivity_draws = sensitivity_draws,
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

  runtime_table <- data.frame(
    stage = c("fit_or_load_draws", "convergence_summary"),
    status = c(fit_status, "completed"),
    elapsed_seconds = c(fit_seconds, convergence_seconds),
    stringsAsFactors = FALSE
  )
  readr::write_csv(runtime_table, runtime_path)

  lot_size_status <- if ("log_acre_lot_z" %in% analysis_objects$predictor_names) {
    "Not run: log acre lot is already in the core M0-M2 design matrix."
  } else {
    paste(
      "Not run automatically: the core model omitted log acre lot, so a separately",
      "labelled complete-case M2 is required before comparing its changed target population."
    )
  }
  student_t_status <- paste(
    "Not run: the caveat-fix branch selected the heteroskedastic M2 instead.",
    "The discarded Student-t experiment scripts/models were removed after they",
    "failed to provide a better official M2 candidate."
  )
  decision_log <- data.frame(
    analysis = c("Prior sensitivity", "Lot-size complete-case sensitivity", "Student-t residual extension"),
    status = c(
      "Run: legacy normal-residual M2 with wider proper geographic-SD priors (half-normal scale 0.6 vs 0.3).",
      lot_size_status,
      student_t_status
    ),
    stringsAsFactors = FALSE
  )
  decision_log_path <- file.path(tables_dir, "planned_sensitivity_analysis_log.csv")
  readr::write_csv(decision_log, decision_log_path)

  timestamp_message(
    "Saved wider-prior M2 draws to ", sensitivity_draws_path,
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
