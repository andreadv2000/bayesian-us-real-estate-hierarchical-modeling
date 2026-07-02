# Experimental state-varying house-size slope M2 caveat check.
#
# This script tests whether allowing the log-house-size coefficient to vary by
# state improves the posterior-predictive over-dispersion caveat. It is kept
# separate from the production pipeline and writes `m2_state_slope_*` outputs.

required_state_slope_packages <- c("coda", "ggplot2", "here", "posterior", "readr", "rjags")
missing_state_slope_packages <- required_state_slope_packages[!vapply(
  required_state_slope_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_state_slope_packages) > 0L) {
  stop(
    "Install the package(s) required for the state-slope M2 experiment: ",
    paste(missing_state_slope_packages, collapse = ", "),
    call. = FALSE
  )
}

source(here::here("scripts", "00_setup.R"))
source(here::here("scripts", "08_posterior_predictive_checks.R"))

slope_read_positive_integer <- function(name, default) {
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

slope_read_flag <- function(name, default = FALSE) {
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

slope_timestamp_message <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

slope_indexed_parameter_names <- function(prefix, count) {
  paste0(prefix, "[", seq_len(count), "]")
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

prepare_state_slope_objects <- function(analysis_objects,
                                        slope_predictor = "log_house_size_z") {
  required_objects <- c(
    "jags_train", "train", "test", "X_train", "X_test",
    "predictor_names", "state_levels", "zip_levels"
  )
  missing_objects <- setdiff(required_objects, names(analysis_objects))
  if (length(missing_objects) > 0L) {
    stop(
      "Analysis objects are missing: ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  X_train <- as.matrix(analysis_objects$X_train)
  X_test <- as.matrix(analysis_objects$X_test)
  size_col <- match(slope_predictor, colnames(X_train))
  if (is.na(size_col)) {
    stop("Slope predictor not found in X_train: ", slope_predictor, call. = FALSE)
  }
  if (nrow(X_train) != analysis_objects$jags_train$N ||
      ncol(X_train) != length(analysis_objects$predictor_names) ||
      ncol(X_test) != length(analysis_objects$predictor_names) ||
      any(!is.finite(X_train)) ||
      any(!is.finite(X_test))) {
    stop("Saved design matrices are inconsistent with the JAGS objects.", call. = FALSE)
  }

  updated <- analysis_objects
  updated$jags_train$size_col <- as.integer(size_col)
  updated$state_slope_predictor <- slope_predictor
  updated$state_slope_size_col <- as.integer(size_col)
  updated
}

make_state_slope_inits <- function(data_list, n_chains, seed) {
  set.seed(seed)
  lapply(seq_len(n_chains), function(chain_id) {
    list(
      alpha = stats::rnorm(1L, 0, 0.2),
      beta = stats::rnorm(data_list$P, 0, 0.1),
      sigma = stats::runif(1L, 0.5, 1.5),
      sd_state = stats::runif(1L, 0.1, 0.8),
      sd_zip = stats::runif(1L, 0.1, 0.8),
      sd_slope = stats::runif(1L, 0.02, 0.4),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = as.integer(seed + 1000L * chain_id)
    )
  })
}

slope_chunk_sizes <- function(total_iter, chunk_iter) {
  starts <- seq.int(1L, total_iter, by = chunk_iter)
  as.integer(pmin(chunk_iter, total_iter - starts + 1L))
}

run_adapt_chunks <- function(model, total_iter, chunk_iter) {
  chunks <- slope_chunk_sizes(total_iter, chunk_iter)
  for (chunk_id in seq_along(chunks)) {
    slope_timestamp_message("Starting state-slope M2 adaptation chunk ", chunk_id, "/", length(chunks), ".")
    started_at <- Sys.time()
    adapted <- rjags::adapt(
      object = model,
      n.iter = chunks[[chunk_id]],
      end.adaptation = chunk_id == length(chunks),
      progress.bar = "none"
    )
    slope_timestamp_message(
      "Finished state-slope M2 adaptation chunk ", chunk_id, "/", length(chunks),
      " in ", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 1),
      " seconds. Adapted=", adapted, "."
    )
  }
  invisible(model)
}

run_update_chunks <- function(model, total_iter, chunk_iter, phase) {
  chunks <- slope_chunk_sizes(total_iter, chunk_iter)
  for (chunk_id in seq_along(chunks)) {
    slope_timestamp_message("Starting state-slope M2 ", phase, " chunk ", chunk_id, "/", length(chunks), ".")
    started_at <- Sys.time()
    stats::update(model, n.iter = chunks[[chunk_id]], progress.bar = "none")
    slope_timestamp_message(
      "Finished state-slope M2 ", phase, " chunk ", chunk_id, "/", length(chunks),
      " in ", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 1),
      " seconds."
    )
  }
  invisible(model)
}

combine_mcmc_chunks <- function(draw_chunks) {
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

sample_coda_chunks <- function(model, parameters, total_iter, thin, chunk_iter,
                               partial_draws_path = NULL) {
  chunks <- slope_chunk_sizes(total_iter, chunk_iter)
  draw_chunks <- vector("list", length(chunks))
  for (chunk_id in seq_along(chunks)) {
    slope_timestamp_message("Starting state-slope M2 sampling chunk ", chunk_id, "/", length(chunks), ".")
    started_at <- Sys.time()
    draw_chunks[[chunk_id]] <- rjags::coda.samples(
      model = model,
      variable.names = parameters,
      n.iter = chunks[[chunk_id]],
      thin = thin,
      progress.bar = "none"
    )
    slope_timestamp_message(
      "Finished state-slope M2 sampling chunk ", chunk_id, "/", length(chunks),
      " in ", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 1),
      " seconds."
    )
    if (!is.null(partial_draws_path)) {
      saveRDS(draw_chunks[seq_len(chunk_id)], partial_draws_path)
      slope_timestamp_message("Saved partial state-slope M2 draw chunks to ", partial_draws_path, ".")
    }
  }
  combine_mcmc_chunks(draw_chunks)
}

fit_state_slope_m2 <- function(model_file, data_list, parameters, settings,
                               partial_draws_path = NULL) {
  if (!file.exists(model_file)) {
    stop("State-slope M2 model file not found: ", model_file, call. = FALSE)
  }
  data_list <- augment_zip_centering_data(data_list)
  inits <- make_state_slope_inits(data_list, settings$n_chains, settings$seed)

  slope_timestamp_message("Compiling state-slope M2.")
  model <- rjags::jags.model(
    file = model_file,
    data = data_list,
    inits = inits,
    n.chains = settings$n_chains,
    n.adapt = 0L,
    quiet = TRUE
  )
  slope_timestamp_message("Compiled state-slope M2.")
  run_adapt_chunks(model, settings$n_adapt, settings$chunk_iter)
  run_update_chunks(model, settings$n_burn, settings$chunk_iter, "burn-in")
  draws <- sample_coda_chunks(
    model = model,
    parameters = parameters,
    total_iter = settings$n_iter,
    thin = settings$thin,
    chunk_iter = settings$chunk_iter,
    partial_draws_path = partial_draws_path
  )
  list(model = model, draws = draws, inits = inits)
}

state_slope_required_parameters <- function(n_predictors, n_states, n_zips) {
  c(
    "alpha", "sigma", "sd_state", "sd_zip", "sd_slope",
    slope_indexed_parameter_names("beta", n_predictors),
    slope_indexed_parameter_names("u", n_states),
    slope_indexed_parameter_names("v", n_zips),
    slope_indexed_parameter_names("w", n_states)
  )
}

validate_state_slope_draw_matrix <- function(draw_matrix, analysis_objects) {
  X <- as.matrix(analysis_objects$X_train)
  draw_matrix <- as.matrix(draw_matrix)
  if (nrow(draw_matrix) == 0L || ncol(draw_matrix) == 0L ||
      is.null(colnames(draw_matrix)) || any(!is.finite(draw_matrix))) {
    stop("State-slope M2 posterior draws are empty, unnamed, or non-finite.", call. = FALSE)
  }
  required_parameters <- state_slope_required_parameters(
    n_predictors = ncol(X),
    n_states = length(analysis_objects$state_levels),
    n_zips = length(analysis_objects$zip_levels)
  )
  missing_parameters <- setdiff(required_parameters, colnames(draw_matrix))
  if (length(missing_parameters) > 0L) {
    stop(
      "State-slope M2 draws are missing parameter(s): ",
      paste(missing_parameters, collapse = ", "),
      call. = FALSE
    )
  }
  if (any(draw_matrix[, "sigma"] <= 0) ||
      any(draw_matrix[, "sd_state"] <= 0) ||
      any(draw_matrix[, "sd_zip"] <= 0) ||
      any(draw_matrix[, "sd_slope"] <= 0)) {
    stop("State-slope M2 draws contain invalid scale values.", call. = FALSE)
  }
  invisible(TRUE)
}

linear_predictor_state_slope <- function(draw, X, state_id, zip_id, size_col) {
  X <- as.matrix(X)
  beta_names <- slope_indexed_parameter_names("beta", ncol(X))
  state_names <- slope_indexed_parameter_names("u", max(state_id))
  zip_names <- slope_indexed_parameter_names("v", max(zip_id))
  slope_names <- slope_indexed_parameter_names("w", max(state_id))
  as.vector(
    draw["alpha"] +
      X %*% unname(draw[beta_names]) +
      unname(draw[slope_names])[state_id] * X[, size_col] +
      unname(draw[state_names])[state_id] +
      unname(draw[zip_names])[zip_id]
  )
}

slope_update_log_density_accumulator <- function(log_max, scaled_sum, log_density) {
  new_log_max <- pmax(log_max, log_density)
  list(
    log_max = new_log_max,
    scaled_sum =
      scaled_sum * exp(log_max - new_log_max) +
      exp(log_density - new_log_max)
  )
}

evaluate_state_slope_holdout <- function(draw_matrix, analysis_objects, draw_ids,
                                         standard_normal_draws) {
  validate_state_slope_draw_matrix(draw_matrix, analysis_objects)
  X_test <- as.matrix(analysis_objects$X_test)
  test_y <- as.numeric(analysis_objects$test$y)
  state_id <- as.integer(analysis_objects$test$state_id)
  zip_id <- as.integer(analysis_objects$test$zip_id_num)
  size_col <- analysis_objects$state_slope_size_col

  n_test <- length(test_y)
  n_draws <- length(draw_ids)
  posterior_mean_sum <- numeric(n_test)
  predictive_draws <- matrix(NA_real_, nrow = n_test, ncol = n_draws)
  log_density_max <- numeric(n_test)
  log_density_sum <- numeric(n_test)

  for (draw_position in seq_along(draw_ids)) {
    draw <- as.numeric(draw_matrix[draw_ids[[draw_position]], , drop = TRUE])
    names(draw) <- colnames(draw_matrix)
    mu <- linear_predictor_state_slope(draw, X_test, state_id, zip_id, size_col)
    posterior_mean_sum <- posterior_mean_sum + mu
    predictive_draws[, draw_position] <-
      mu + draw["sigma"] * standard_normal_draws[, draw_position]
    log_density <- stats::dnorm(test_y, mean = mu, sd = draw["sigma"], log = TRUE)
    if (draw_position == 1L) {
      log_density_max <- log_density
      log_density_sum <- rep(1, n_test)
    } else {
      accumulator <- slope_update_log_density_accumulator(
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
    model = "M2_state_slope",
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
    model = "M2_state_slope",
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

simulate_state_slope_replicated_outcome <- function(draw, X, state_id, zip_id, size_col) {
  mu <- linear_predictor_state_slope(draw, X, state_id, zip_id, size_col)
  stats::rnorm(length(mu), mean = mu, sd = draw["sigma"])
}

run_state_slope_ppc <- function(draw_matrix, analysis_objects, draw_ids,
                                selected_state_ids) {
  validate_state_slope_draw_matrix(draw_matrix, analysis_objects)
  X_train <- as.matrix(analysis_objects$X_train)
  jags_train <- analysis_objects$jags_train
  state_id <- as.integer(jags_train$state)
  zip_id <- as.integer(jags_train$zip)
  size_col <- analysis_objects$state_slope_size_col

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
    y_rep <- simulate_state_slope_replicated_outcome(
      draw = draw,
      X = X_train,
      state_id = state_id,
      zip_id = zip_id,
      size_col = size_col
    )
    statistics <- calculate_ppc_statistics(
      y = y_rep,
      state_id = state_id,
      zip_id = zip_id,
      state_levels = analysis_objects$state_levels,
      selected_state_ids = selected_state_ids
    )
    statistics$replication <- index
    statistics$model <- "M2_state_slope"
    statistics
  }))
  observed_statistics$model <- "M2_state_slope"
  summary <- summarise_ppc_statistics(observed_statistics, replicated_statistics)
  summary$model <- "M2_state_slope"
  list(
    summary = summary,
    replicated_statistics = replicated_statistics,
    observed_statistics = observed_statistics
  )
}

summarise_state_slope_convergence <- function(draws) {
  draws_array <- posterior::as_draws_array(draws)
  summary <- as.data.frame(posterior::summarise_draws(
    draws_array,
    mean,
    median,
    sd,
    posterior::quantile2,
    posterior::rhat,
    posterior::ess_bulk,
    posterior::ess_tail,
    posterior::mcse_mean
  ))
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

write_state_slope_comparison_tables <- function(heldout_metrics, ppc_summary, tables_dir) {
  production_metrics_path <- file.path(tables_dir, "heldout_predictive_metrics.csv")
  if (file.exists(production_metrics_path)) {
    production_metrics <- readr::read_csv(production_metrics_path, show_col_types = FALSE)
    production_m2 <- production_metrics[production_metrics$model == "M2", , drop = FALSE]
    if (nrow(production_m2) == 1L) {
      readr::write_csv(
        rbind(production_m2, heldout_metrics),
        file.path(tables_dir, "m2_state_slope_vs_m2_heldout_comparison.csv")
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
      file.path(tables_dir, "m2_state_slope_vs_m2_ppc_comparison.csv")
    )
  }
}

run_state_slope_m2_experiment <- function(
    analysis_objects_path = here::here("data", "processed", "analysis_objects.rds"),
    model_file = here::here("models", "m2_state_varying_size_slope.jags"),
    mcmc_output_dir = here::here("output", "mcmc"),
    tables_dir = here::here("output", "tables"),
    figures_dir = here::here("output", "figures")) {
  settings <- list(
    n_chains = slope_read_positive_integer("STATE_SLOPE_M2_N_CHAINS", 4L),
    n_adapt = slope_read_positive_integer("STATE_SLOPE_M2_N_ADAPT", 2000L),
    n_burn = slope_read_positive_integer("STATE_SLOPE_M2_N_BURN", 5000L),
    n_iter = slope_read_positive_integer("STATE_SLOPE_M2_N_ITER", 10000L),
    thin = slope_read_positive_integer("STATE_SLOPE_M2_THIN", 1L),
    seed = slope_read_positive_integer("STATE_SLOPE_M2_SEED", 20260702L),
    chunk_iter = slope_read_positive_integer("STATE_SLOPE_M2_CHUNK_ITER", 1000L),
    max_heldout_draws = slope_read_positive_integer("STATE_SLOPE_M2_MAX_HELDOUT_DRAWS", 2000L),
    max_ppc_draws = slope_read_positive_integer("STATE_SLOPE_M2_MAX_PPC_DRAWS", 1000L),
    force = slope_read_flag("STATE_SLOPE_M2_FORCE", default = FALSE)
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

  analysis_objects <- prepare_state_slope_objects(readRDS(analysis_objects_path))
  parameters <- c("alpha", "beta", "sigma", "sd_state", "sd_zip", "sd_slope", "u", "v", "w")
  draws_path <- file.path(mcmc_output_dir, "m2_state_slope_draws.rds")
  partial_draws_path <- file.path(mcmc_output_dir, "m2_state_slope_draw_chunks_partial.rds")
  runtime_path <- file.path(tables_dir, "m2_state_slope_runtime.csv")

  readr::write_csv(
    data.frame(setting = names(settings), value = unlist(settings), row.names = NULL),
    file.path(tables_dir, "m2_state_slope_settings.csv")
  )
  readr::write_csv(
    data.frame(
      slope_predictor = analysis_objects$state_slope_predictor,
      size_col = analysis_objects$state_slope_size_col,
      stringsAsFactors = FALSE
    ),
    file.path(tables_dir, "m2_state_slope_terms.csv")
  )

  freshness_mtime <- max(
    file.info(c(model_file, here::here("scripts", "13_m2_state_slope_experiment.R")))$mtime,
    na.rm = TRUE
  )
  draws_current <- file.exists(draws_path) &&
    file.info(draws_path)$mtime >= freshness_mtime

  fit_started_at <- Sys.time()
  if (draws_current && !settings$force) {
    slope_timestamp_message(
      "Skipping state-slope M2 refit because current draws already exist at ",
      draws_path, ". Set STATE_SLOPE_M2_FORCE=true to refit."
    )
    draws <- readRDS(draws_path)
    fit_status <- "skipped_existing_draws"
    fit_seconds <- NA_real_
  } else {
    slope_timestamp_message(
      "Starting state-slope M2 fit. Settings: chains=", settings$n_chains,
      ", adapt=", settings$n_adapt,
      ", burn=", settings$n_burn,
      ", iter=", settings$n_iter,
      ", thin=", settings$thin,
      ", chunk=", settings$chunk_iter,
      "."
    )
    fit <- fit_state_slope_m2(
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
    slope_timestamp_message(
      "Saved state-slope M2 draws to ", draws_path,
      " after ", round(fit_seconds, 1), " seconds."
    )
  }

  draw_matrix <- as.matrix(posterior::as_draws_matrix(draws))
  validate_state_slope_draw_matrix(draw_matrix, analysis_objects)

  convergence_started_at <- Sys.time()
  convergence <- summarise_state_slope_convergence(draws)
  convergence_path <- file.path(tables_dir, "m2_state_slope_posterior_summary.csv")
  readr::write_csv(convergence, convergence_path)
  convergence_seconds <- as.numeric(difftime(Sys.time(), convergence_started_at, units = "secs"))
  slope_timestamp_message(
    "Saved state-slope M2 convergence summary to ", convergence_path,
    ". Converged ", sum(convergence$converged), "/", nrow(convergence), " parameters."
  )

  set.seed(settings$seed + 101L)
  heldout_draw_ids <- sort(sample(
    seq_len(nrow(draw_matrix)),
    size = min(settings$max_heldout_draws, nrow(draw_matrix)),
    replace = FALSE
  ))
  standard_normal_draws <- matrix(
    stats::rnorm(nrow(analysis_objects$test) * length(heldout_draw_ids)),
    nrow = nrow(analysis_objects$test),
    ncol = length(heldout_draw_ids)
  )
  heldout <- evaluate_state_slope_holdout(
    draw_matrix = draw_matrix,
    analysis_objects = analysis_objects,
    draw_ids = heldout_draw_ids,
    standard_normal_draws = standard_normal_draws
  )
  heldout_metrics_path <- file.path(tables_dir, "m2_state_slope_heldout_predictive_metrics.csv")
  heldout_predictions_path <- file.path(tables_dir, "m2_state_slope_heldout_predictive_predictions.csv")
  readr::write_csv(heldout$metrics, heldout_metrics_path)
  readr::write_csv(heldout$predictions, heldout_predictions_path)

  set.seed(settings$seed + 201L)
  ppc_draw_ids <- sort(sample(
    seq_len(nrow(draw_matrix)),
    size = min(settings$max_ppc_draws, nrow(draw_matrix)),
    replace = FALSE
  ))
  selected_state_ids <- select_ppc_state_ids(analysis_objects$state_levels)
  ppc <- run_state_slope_ppc(
    draw_matrix = draw_matrix,
    analysis_objects = analysis_objects,
    draw_ids = ppc_draw_ids,
    selected_state_ids = selected_state_ids
  )
  ppc_summary_path <- file.path(tables_dir, "m2_state_slope_posterior_predictive_check_summary.csv")
  ppc_replicates_path <- file.path(tables_dir, "m2_state_slope_posterior_predictive_check_replicates.csv")
  ppc_figure_path <- file.path(figures_dir, "m2_state_slope_posterior_predictive_checks.png")
  readr::write_csv(ppc$summary, ppc_summary_path)
  readr::write_csv(ppc$replicated_statistics, ppc_replicates_path)
  save_ppc_plot(
    replicated_statistics = ppc$replicated_statistics,
    observed_statistics = ppc$observed_statistics,
    model = "M2 state-varying size slope",
    figure_path = ppc_figure_path
  )
  write_state_slope_comparison_tables(heldout$metrics, ppc$summary, tables_dir)

  runtime <- data.frame(
    stage = c("fit_or_load_draws", "convergence_summary", "heldout_ppc"),
    status = c(fit_status, "completed", "completed"),
    elapsed_seconds = c(fit_seconds, convergence_seconds, NA_real_),
    stringsAsFactors = FALSE
  )
  readr::write_csv(runtime, runtime_path)
  slope_timestamp_message(
    "State-slope M2 experiment complete. Key outputs: ", convergence_path,
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
  run_state_slope_m2_experiment()
}
