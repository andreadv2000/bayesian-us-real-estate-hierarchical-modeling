# Selected heteroskedastic residual M2 fit and checks.
#
# This model keeps the M2 mean structure and geographic effects, but replaces
# the single residual SD with observation-specific residual SDs. The log
# residual SD depends on property predictors plus a centered state-level
# residual-scale effect.  On the caveat-fix branch, this is the selected M2
# candidate rather than a discarded experiment.

source(here::here("scripts", "12_m2_heteroskedastic_helpers.R"))

hetero_read_positive_number <- function(name, default) {
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

hetero_read_label <- function(name, default) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) {
    return(default)
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

prepare_hetero_objects <- function(analysis_objects, settings) {
  validate_student_t_analysis_objects(analysis_objects)
  X_sigma_train <- as.matrix(analysis_objects$X_train)
  X_sigma_test <- as.matrix(analysis_objects$X_test)
  if (ncol(X_sigma_train) != length(analysis_objects$predictor_names) ||
      ncol(X_sigma_test) != length(analysis_objects$predictor_names) ||
      any(!is.finite(X_sigma_train)) ||
      any(!is.finite(X_sigma_test))) {
    stop("Saved design matrices are inconsistent with sigma predictors.", call. = FALSE)
  }

  updated <- analysis_objects
  updated$jags_train$X_sigma <- X_sigma_train
  updated$jags_train$P_sigma <- as.integer(ncol(X_sigma_train))
  updated$jags_train$log_sigma_prior_mean <- log(settings$sigma_prior_center)
  updated$jags_train$tau_gamma0_prior <- 1 / (settings$gamma0_prior_scale^2)
  updated$jags_train$tau_gamma_prior <- 1 / (settings$gamma_prior_scale^2)
  updated$jags_train$tau_sigma_state_prior <- 1 / (settings$sigma_state_prior_scale^2)
  updated$X_sigma_train <- X_sigma_train
  updated$X_sigma_test <- X_sigma_test
  updated$sigma_predictor_names <- analysis_objects$predictor_names
  updated
}

make_hetero_inits <- function(data_list, n_chains, seed, settings) {
  set.seed(seed)
  lapply(seq_len(n_chains), function(chain_id) {
    list(
      alpha = stats::rnorm(1L, 0, 0.2),
      beta = stats::rnorm(data_list$P, 0, 0.1),
      gamma0 = stats::rnorm(1L, log(settings$sigma_prior_center), 0.1),
      gamma = stats::rnorm(data_list$P_sigma, 0, 0.03),
      sd_sigma_state = stats::runif(1L, 0.02, 0.2),
      sd_state = stats::runif(1L, 0.1, 0.8),
      sd_zip = stats::runif(1L, 0.1, 0.8),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = as.integer(seed + 1000L * chain_id)
    )
  })
}

hetero_required_parameters <- function(n_predictors, n_states, n_zips) {
  c(
    "alpha", "gamma0", "sd_sigma_state", "sd_state", "sd_zip",
    student_t_indexed_parameter_names("beta", n_predictors),
    student_t_indexed_parameter_names("gamma", n_predictors),
    student_t_indexed_parameter_names("h", n_states),
    student_t_indexed_parameter_names("u", n_states),
    student_t_indexed_parameter_names("v", n_zips)
  )
}

validate_hetero_draw_matrix <- function(draw_matrix, analysis_objects) {
  X <- as.matrix(analysis_objects$X_train)
  draw_matrix <- as.matrix(draw_matrix)
  if (nrow(draw_matrix) == 0L || ncol(draw_matrix) == 0L ||
      is.null(colnames(draw_matrix)) || any(!is.finite(draw_matrix))) {
    stop("Heteroskedastic M2 posterior draws are empty, unnamed, or non-finite.", call. = FALSE)
  }
  required_parameters <- hetero_required_parameters(
    n_predictors = ncol(X),
    n_states = length(analysis_objects$state_levels),
    n_zips = length(analysis_objects$zip_levels)
  )
  missing_parameters <- setdiff(required_parameters, colnames(draw_matrix))
  if (length(missing_parameters) > 0L) {
    stop(
      "Heteroskedastic M2 draws are missing parameter(s): ",
      paste(missing_parameters, collapse = ", "),
      call. = FALSE
    )
  }
  if (any(draw_matrix[, "sd_sigma_state"] <= 0) ||
      any(draw_matrix[, "sd_state"] <= 0) ||
      any(draw_matrix[, "sd_zip"] <= 0)) {
    stop("Heteroskedastic M2 draws contain invalid scale values.", call. = FALSE)
  }
  invisible(TRUE)
}

fit_hetero_m2 <- function(model_file, data_list, parameters, settings,
                          partial_draws_path = NULL) {
  if (!file.exists(model_file)) {
    stop("Heteroskedastic M2 model file not found: ", model_file, call. = FALSE)
  }
  data_list <- augment_zip_centering_data(data_list)
  inits <- make_hetero_inits(data_list, settings$n_chains, settings$seed, settings)

  hetero_timestamp_message("Compiling heteroskedastic M2.")
  model <- rjags::jags.model(
    file = model_file,
    data = data_list,
    inits = inits,
    n.chains = settings$n_chains,
    n.adapt = 0L,
    quiet = TRUE
  )
  hetero_timestamp_message("Compiled heteroskedastic M2.")
  run_hetero_adapt_chunks(model, settings$n_adapt, settings$chunk_iter)
  run_hetero_update_chunks(model, settings$n_burn, settings$chunk_iter, "burn-in")
  draws <- sample_hetero_coda_chunks(
    model = model,
    parameters = parameters,
    total_iter = settings$n_iter,
    thin = settings$thin,
    chunk_iter = settings$chunk_iter,
    partial_draws_path = partial_draws_path
  )
  list(model = model, draws = draws, inits = inits)
}

hetero_timestamp_message <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

run_hetero_adapt_chunks <- function(model, total_iter, chunk_iter) {
  chunks <- student_t_chunk_sizes(total_iter, chunk_iter)
  for (chunk_id in seq_along(chunks)) {
    hetero_timestamp_message("Starting heteroskedastic M2 adaptation chunk ", chunk_id, "/", length(chunks), ".")
    started_at <- Sys.time()
    adapted <- rjags::adapt(
      object = model,
      n.iter = chunks[[chunk_id]],
      end.adaptation = chunk_id == length(chunks),
      progress.bar = "none"
    )
    hetero_timestamp_message(
      "Finished heteroskedastic M2 adaptation chunk ", chunk_id, "/", length(chunks),
      " in ", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 1),
      " seconds. Adapted=", adapted, "."
    )
  }
  invisible(model)
}

run_hetero_update_chunks <- function(model, total_iter, chunk_iter, phase) {
  chunks <- student_t_chunk_sizes(total_iter, chunk_iter)
  for (chunk_id in seq_along(chunks)) {
    hetero_timestamp_message("Starting heteroskedastic M2 ", phase, " chunk ", chunk_id, "/", length(chunks), ".")
    started_at <- Sys.time()
    stats::update(model, n.iter = chunks[[chunk_id]], progress.bar = "none")
    hetero_timestamp_message(
      "Finished heteroskedastic M2 ", phase, " chunk ", chunk_id, "/", length(chunks),
      " in ", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 1),
      " seconds."
    )
  }
  invisible(model)
}

combine_hetero_mcmc_chunks <- function(draw_chunks) {
  if (length(draw_chunks) == 0L) {
    stop("No MCMC chunks were sampled.", call. = FALSE)
  }
  n_chains <- length(draw_chunks[[1L]])
  coda::mcmc.list(lapply(seq_len(n_chains), function(chain_id) {
    combined <- do.call(rbind, lapply(draw_chunks, function(chunk) chunk[[chain_id]]))
    if (inherits(combined, "mcmc")) {
      combined
    } else {
      coda::mcmc(combined)
    }
  }))
}

sample_hetero_coda_chunks <- function(model, parameters, total_iter, thin,
                                      chunk_iter, partial_draws_path = NULL) {
  chunks <- student_t_chunk_sizes(total_iter, chunk_iter)
  draw_chunks <- vector("list", length(chunks))
  for (chunk_id in seq_along(chunks)) {
    hetero_timestamp_message("Starting heteroskedastic M2 sampling chunk ", chunk_id, "/", length(chunks), ".")
    started_at <- Sys.time()
    draw_chunks[[chunk_id]] <- rjags::coda.samples(
      model = model,
      variable.names = parameters,
      n.iter = chunks[[chunk_id]],
      thin = thin,
      progress.bar = "none"
    )
    hetero_timestamp_message(
      "Finished heteroskedastic M2 sampling chunk ", chunk_id, "/", length(chunks),
      " in ", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 1),
      " seconds."
    )
    if (!is.null(partial_draws_path)) {
      saveRDS(draw_chunks[seq_len(chunk_id)], partial_draws_path)
      hetero_timestamp_message("Saved partial heteroskedastic M2 draw chunks to ", partial_draws_path, ".")
    }
  }
  combine_hetero_mcmc_chunks(draw_chunks)
}

linear_predictor_hetero <- function(draw, X, state_id, zip_id) {
  X <- as.matrix(X)
  beta_names <- student_t_indexed_parameter_names("beta", ncol(X))
  state_names <- student_t_indexed_parameter_names("u", max(state_id))
  zip_names <- student_t_indexed_parameter_names("v", max(zip_id))
  as.vector(
    draw["alpha"] +
      X %*% unname(draw[beta_names]) +
      unname(draw[state_names])[state_id] +
      unname(draw[zip_names])[zip_id]
  )
}

residual_sigma_hetero <- function(draw, X_sigma, state_id) {
  X_sigma <- as.matrix(X_sigma)
  gamma_names <- student_t_indexed_parameter_names("gamma", ncol(X_sigma))
  state_sigma_names <- student_t_indexed_parameter_names("h", max(state_id))
  log_sigma <- as.vector(
    draw["gamma0"] +
      X_sigma %*% unname(draw[gamma_names]) +
      unname(draw[state_sigma_names])[state_id]
  )
  exp(log_sigma)
}

simulate_hetero_replicated_outcome <- function(draw, X, X_sigma, state_id, zip_id) {
  mu <- linear_predictor_hetero(draw, X, state_id, zip_id)
  sigma <- residual_sigma_hetero(draw, X_sigma, state_id)
  stats::rnorm(length(mu), mean = mu, sd = sigma)
}

evaluate_hetero_holdout <- function(draw_matrix, analysis_objects, draw_ids) {
  validate_hetero_draw_matrix(draw_matrix, analysis_objects)
  X_test <- as.matrix(analysis_objects$X_test)
  X_sigma_test <- as.matrix(analysis_objects$X_sigma_test)
  test_y <- as.numeric(analysis_objects$test$y)
  state_id <- as.integer(analysis_objects$test$state_id)
  zip_id <- as.integer(analysis_objects$test$zip_id_num)

  n_test <- length(test_y)
  n_draws <- length(draw_ids)
  posterior_mean_sum <- numeric(n_test)
  predictive_draws <- matrix(NA_real_, nrow = n_test, ncol = n_draws)
  log_density_max <- numeric(n_test)
  log_density_sum <- numeric(n_test)

  for (draw_position in seq_along(draw_ids)) {
    draw <- as.numeric(draw_matrix[draw_ids[[draw_position]], , drop = TRUE])
    names(draw) <- colnames(draw_matrix)
    mu <- linear_predictor_hetero(draw, X_test, state_id, zip_id)
    sigma <- residual_sigma_hetero(draw, X_sigma_test, state_id)
    posterior_mean_sum <- posterior_mean_sum + mu
    predictive_draws[, draw_position] <- stats::rnorm(n_test, mean = mu, sd = sigma)
    log_density <- stats::dnorm(test_y, mean = mu, sd = sigma, log = TRUE)
    if (draw_position == 1L) {
      log_density_max <- log_density
      log_density_sum <- rep(1, n_test)
    } else {
      accumulator <- student_t_update_log_density_accumulator(
        log_density_max,
        log_density_sum,
        log_density
      )
      log_density_max <- accumulator$log_max
      log_density_sum <- accumulator$scaled_sum
    }
  }

  posterior_mean <- posterior_mean_sum / n_draws
  lower_95 <- apply(
    predictive_draws, 1L, stats::quantile,
    probs = 0.025, names = FALSE, type = 8L
  )
  upper_95 <- apply(
    predictive_draws, 1L, stats::quantile,
    probs = 0.975, names = FALSE, type = 8L
  )
  log_predictive_density <- log_density_max + log(log_density_sum) - log(n_draws)
  covered_95 <- test_y >= lower_95 & test_y <= upper_95

  predictions <- data.frame(
    model = "M2_heteroskedastic",
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
    model = "M2_heteroskedastic",
    n_test = n_test,
    n_posterior_draws = n_draws,
    RMSE = sqrt(mean((posterior_mean - test_y)^2)),
    MAE = mean(abs(posterior_mean - test_y)),
    coverage_95 = mean(covered_95),
    total_lppd = sum(log_predictive_density),
    MLPD = mean(log_predictive_density),
    row.names = NULL
  )
  list(metrics = metrics, predictions = predictions)
}

run_hetero_ppc <- function(draw_matrix, analysis_objects, draw_ids,
                           selected_state_ids) {
  validate_hetero_draw_matrix(draw_matrix, analysis_objects)
  X_train <- as.matrix(analysis_objects$X_train)
  X_sigma_train <- as.matrix(analysis_objects$X_sigma_train)
  jags_train <- analysis_objects$jags_train
  state_id <- as.integer(jags_train$state)
  zip_id <- as.integer(jags_train$zip)

  observed_statistics <- calculate_ppc_statistics(
    y = jags_train$y,
    state_id = state_id,
    zip_id = zip_id,
    state_levels = analysis_objects$state_levels,
    selected_state_ids = selected_state_ids
  )
  replicated_statistics <- do.call(rbind, lapply(seq_along(draw_ids), function(index) {
    draw <- as.numeric(draw_matrix[draw_ids[[index]], , drop = TRUE])
    names(draw) <- colnames(draw_matrix)
    y_rep <- simulate_hetero_replicated_outcome(
      draw = draw,
      X = X_train,
      X_sigma = X_sigma_train,
      state_id = state_id,
      zip_id = zip_id
    )
    statistics <- calculate_ppc_statistics(
      y = y_rep,
      state_id = state_id,
      zip_id = zip_id,
      state_levels = analysis_objects$state_levels,
      selected_state_ids = selected_state_ids
    )
    statistics$replication <- index
    statistics$model <- "M2_heteroskedastic"
    statistics
  }))
  observed_statistics$model <- "M2_heteroskedastic"
  summary <- summarise_ppc_statistics(observed_statistics, replicated_statistics)
  summary$model <- "M2_heteroskedastic"
  list(
    summary = summary,
    replicated_statistics = replicated_statistics,
    observed_statistics = observed_statistics
  )
}

write_hetero_comparison_tables <- function(heldout_metrics, ppc_summary,
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

run_hetero_m2_experiment <- function(
    analysis_objects_path = here::here("data", "processed", "analysis_objects.rds"),
    model_file = here::here("models", "m2_heteroskedastic_state_sigma.jags"),
    mcmc_output_dir = here::here("output", "mcmc"),
    tables_dir = here::here("output", "tables"),
    figures_dir = here::here("output", "figures")) {
  settings <- list(
    n_chains = student_t_read_positive_integer("HETERO_M2_N_CHAINS", 4L),
    n_adapt = student_t_read_positive_integer("HETERO_M2_N_ADAPT", 2000L),
    n_burn = student_t_read_positive_integer("HETERO_M2_N_BURN", 5000L),
    n_iter = student_t_read_positive_integer("HETERO_M2_N_ITER", 10000L),
    thin = student_t_read_positive_integer("HETERO_M2_THIN", 1L),
    seed = student_t_read_positive_integer("HETERO_M2_SEED", 20260703L),
    chunk_iter = student_t_read_positive_integer("HETERO_M2_CHUNK_ITER", 1000L),
    max_heldout_draws = student_t_read_positive_integer("HETERO_M2_MAX_HELDOUT_DRAWS", 2000L),
    max_ppc_draws = student_t_read_positive_integer("HETERO_M2_MAX_PPC_DRAWS", 1000L),
    force = student_t_read_flag("HETERO_M2_FORCE", default = FALSE),
    sigma_prior_center = hetero_read_positive_number("HETERO_M2_SIGMA_PRIOR_CENTER", 0.5),
    gamma0_prior_scale = hetero_read_positive_number("HETERO_M2_GAMMA0_PRIOR_SCALE", 0.5),
    gamma_prior_scale = hetero_read_positive_number("HETERO_M2_GAMMA_PRIOR_SCALE", 0.2),
    sigma_state_prior_scale = hetero_read_positive_number("HETERO_M2_SIGMA_STATE_PRIOR_SCALE", 0.2),
    label = hetero_read_label("HETERO_M2_LABEL", "xall_state")
  )

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

  analysis_objects <- prepare_hetero_objects(readRDS(analysis_objects_path), settings)
  parameters <- c(
    "alpha", "beta", "gamma0", "gamma", "sd_sigma_state",
    "sd_state", "sd_zip", "h", "u", "v"
  )
  output_stem <- paste0("m2_heteroskedastic_", settings$label)
  model_label <- paste0("M2_heteroskedastic_", settings$label)
  draws_path <- file.path(mcmc_output_dir, paste0(output_stem, "_draws.rds"))
  partial_draws_path <- file.path(mcmc_output_dir, paste0(output_stem, "_draw_chunks_partial.rds"))
  runtime_path <- file.path(tables_dir, paste0(output_stem, "_runtime.csv"))

  readr::write_csv(
    data.frame(setting = names(settings), value = unlist(settings), row.names = NULL),
    file.path(tables_dir, paste0(output_stem, "_settings.csv"))
  )
  readr::write_csv(
    data.frame(
      sigma_predictor = analysis_objects$sigma_predictor_names,
      stringsAsFactors = FALSE
    ),
    file.path(tables_dir, paste0(output_stem, "_sigma_terms.csv"))
  )

  freshness_mtime <- max(
    file.info(c(model_file, here::here("scripts", "17_m2_heteroskedastic_experiment.R")))$mtime,
    na.rm = TRUE
  )
  draws_current <- file.exists(draws_path) &&
    file.info(draws_path)$mtime >= freshness_mtime

  fit_started_at <- Sys.time()
  if (draws_current && !settings$force) {
    hetero_timestamp_message(
      "Skipping heteroskedastic M2 refit because current draws already exist at ",
      draws_path, ". Set HETERO_M2_FORCE=true to refit."
    )
    draws <- readRDS(draws_path)
    fit_status <- "skipped_existing_draws"
    fit_seconds <- NA_real_
  } else {
    hetero_timestamp_message(
      "Starting heteroskedastic M2 fit. Settings: chains=", settings$n_chains,
      ", adapt=", settings$n_adapt,
      ", burn=", settings$n_burn,
      ", iter=", settings$n_iter,
      ", thin=", settings$thin,
      ", chunk=", settings$chunk_iter,
      ", gamma_prior_scale=", settings$gamma_prior_scale,
      ", sigma_state_prior_scale=", settings$sigma_state_prior_scale,
      "."
    )
    fit <- fit_hetero_m2(
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
    hetero_timestamp_message(
      "Saved heteroskedastic M2 draws to ", draws_path,
      " after ", round(fit_seconds, 1), " seconds."
    )
  }

  draw_matrix <- as.matrix(posterior::as_draws_matrix(draws))
  validate_hetero_draw_matrix(draw_matrix, analysis_objects)

  convergence_started_at <- Sys.time()
  convergence <- summarise_student_t_convergence(draws)
  convergence_path <- file.path(tables_dir, paste0(output_stem, "_posterior_summary.csv"))
  readr::write_csv(convergence, convergence_path)
  convergence_seconds <- as.numeric(difftime(Sys.time(), convergence_started_at, units = "secs"))
  hetero_timestamp_message(
    "Saved heteroskedastic M2 convergence summary to ", convergence_path,
    ". Converged ", sum(convergence$converged), "/", nrow(convergence), " parameters."
  )

  set.seed(settings$seed + 101L)
  heldout_draw_ids <- sort(sample(
    seq_len(nrow(draw_matrix)),
    size = min(settings$max_heldout_draws, nrow(draw_matrix)),
    replace = FALSE
  ))
  heldout <- evaluate_hetero_holdout(
    draw_matrix = draw_matrix,
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
    seq_len(nrow(draw_matrix)),
    size = min(settings$max_ppc_draws, nrow(draw_matrix)),
    replace = FALSE
  ))
  selected_state_ids <- select_ppc_state_ids(analysis_objects$state_levels)
  ppc <- run_hetero_ppc(
    draw_matrix = draw_matrix,
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
    model = "M2 heteroskedastic residual",
    figure_path = ppc_figure_path
  )
  write_hetero_comparison_tables(heldout$metrics, ppc$summary, tables_dir, output_stem)

  runtime <- data.frame(
    stage = c("fit_or_load_draws", "convergence_summary", "heldout_ppc"),
    status = c(fit_status, "completed", "completed"),
    elapsed_seconds = c(fit_seconds, convergence_seconds, NA_real_),
    stringsAsFactors = FALSE
  )
  readr::write_csv(runtime, runtime_path)
  hetero_timestamp_message(
    "Heteroskedastic M2 experiment complete. Key outputs: ", convergence_path,
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
  run_hetero_m2_experiment()
}
