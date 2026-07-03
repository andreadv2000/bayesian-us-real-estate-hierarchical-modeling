# Experimental contaminated-normal residual M2 caveat check.
#
# This model keeps the production M2 geography but splits residual noise into a
# narrow core component and a wider outlier component. The goal is to reduce the
# overall IQR/q90 PPC over-dispersion without destroying ZIP heterogeneity.

source(here::here("scripts", "15_m2_student_t_fixed_nu_experiment.R"))

contam_read_positive_number <- function(name, default) {
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

contam_read_label <- function(name, default) {
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

make_contam_inits <- function(data_list, n_chains, seed, settings) {
  set.seed(seed)
  lapply(seq_len(n_chains), function(chain_id) {
    list(
      alpha = stats::rnorm(1L, 0, 0.2),
      beta = stats::rnorm(data_list$P, 0, 0.1),
      sigma_core = stats::runif(1L, 0.2, 0.8),
      sigma_extra = stats::runif(1L, 0.2, 1.2),
      p_outlier = stats::rbeta(
        1L,
        shape1 = settings$p_outlier_alpha,
        shape2 = settings$p_outlier_beta
      ),
      sd_state = stats::runif(1L, 0.1, 0.8),
      sd_zip = stats::runif(1L, 0.1, 0.8),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = as.integer(seed + 1000L * chain_id)
    )
  })
}

contam_required_parameters <- function(n_predictors, n_states, n_zips) {
  c(
    "alpha", "sigma_core", "sigma_outlier", "p_outlier", "sd_state", "sd_zip",
    student_t_indexed_parameter_names("beta", n_predictors),
    student_t_indexed_parameter_names("u", n_states),
    student_t_indexed_parameter_names("v", n_zips)
  )
}

validate_contam_draw_matrix <- function(draw_matrix, analysis_objects) {
  X <- as.matrix(analysis_objects$X_train)
  draw_matrix <- as.matrix(draw_matrix)
  if (nrow(draw_matrix) == 0L || ncol(draw_matrix) == 0L ||
      is.null(colnames(draw_matrix)) || any(!is.finite(draw_matrix))) {
    stop("Contaminated-normal M2 posterior draws are empty, unnamed, or non-finite.", call. = FALSE)
  }
  required_parameters <- contam_required_parameters(
    n_predictors = ncol(X),
    n_states = length(analysis_objects$state_levels),
    n_zips = length(analysis_objects$zip_levels)
  )
  missing_parameters <- setdiff(required_parameters, colnames(draw_matrix))
  if (length(missing_parameters) > 0L) {
    stop(
      "Contaminated-normal M2 draws are missing parameter(s): ",
      paste(missing_parameters, collapse = ", "),
      call. = FALSE
    )
  }
  if (any(draw_matrix[, "sigma_core"] <= 0) ||
      any(draw_matrix[, "sigma_outlier"] <= 0) ||
      any(draw_matrix[, "sigma_outlier"] <= draw_matrix[, "sigma_core"]) ||
      any(draw_matrix[, "sd_state"] <= 0) ||
      any(draw_matrix[, "sd_zip"] <= 0) ||
      any(draw_matrix[, "p_outlier"] <= 0) ||
      any(draw_matrix[, "p_outlier"] >= 1)) {
    stop("Contaminated-normal M2 draws contain invalid scale or mixture-probability values.", call. = FALSE)
  }
  invisible(TRUE)
}

fit_contam_m2 <- function(model_file, data_list, parameters, settings,
                          partial_draws_path = NULL) {
  if (!file.exists(model_file)) {
    stop("Contaminated-normal M2 model file not found: ", model_file, call. = FALSE)
  }
  data_list <- augment_zip_centering_data(data_list)
  data_list$p_outlier_alpha <- settings$p_outlier_alpha
  data_list$p_outlier_beta <- settings$p_outlier_beta
  data_list$tau_sigma_extra_prior <- 1 / (settings$sigma_extra_prior_scale^2)
  inits <- make_contam_inits(data_list, settings$n_chains, settings$seed, settings)

  student_t_timestamp_message("Compiling contaminated-normal M2.")
  model <- rjags::jags.model(
    file = model_file,
    data = data_list,
    inits = inits,
    n.chains = settings$n_chains,
    n.adapt = 0L,
    quiet = TRUE
  )
  student_t_timestamp_message("Compiled contaminated-normal M2.")
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

linear_predictor_contam <- function(draw, X, state_id, zip_id) {
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

contam_log_density <- function(y, mu, sigma_core, sigma_outlier, p_outlier) {
  core_log <- log1p(-p_outlier) +
    stats::dnorm(y, mean = mu, sd = sigma_core, log = TRUE)
  outlier_log <- log(p_outlier) +
    stats::dnorm(y, mean = mu, sd = sigma_outlier, log = TRUE)
  max_log <- pmax(core_log, outlier_log)
  max_log + log(exp(core_log - max_log) + exp(outlier_log - max_log))
}

simulate_contam_replicated_outcome <- function(draw, X, state_id, zip_id) {
  mu <- linear_predictor_contam(draw, X, state_id, zip_id)
  outlier_flag <- stats::rbinom(length(mu), size = 1L, prob = draw["p_outlier"])
  sigma <- ifelse(outlier_flag == 1L, draw["sigma_outlier"], draw["sigma_core"])
  stats::rnorm(length(mu), mean = mu, sd = sigma)
}

evaluate_contam_holdout <- function(draw_matrix, analysis_objects, draw_ids) {
  validate_contam_draw_matrix(draw_matrix, analysis_objects)
  X_test <- as.matrix(analysis_objects$X_test)
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
    mu <- linear_predictor_contam(draw, X_test, state_id, zip_id)
    posterior_mean_sum <- posterior_mean_sum + mu
    predictive_draws[, draw_position] <- simulate_contam_replicated_outcome(
      draw = draw,
      X = X_test,
      state_id = state_id,
      zip_id = zip_id
    )
    log_density <- contam_log_density(
      y = test_y,
      mu = mu,
      sigma_core = draw["sigma_core"],
      sigma_outlier = draw["sigma_outlier"],
      p_outlier = draw["p_outlier"]
    )
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
    model = "M2_contaminated_normal",
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
    model = "M2_contaminated_normal",
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

run_contam_ppc <- function(draw_matrix, analysis_objects, draw_ids,
                           selected_state_ids) {
  validate_contam_draw_matrix(draw_matrix, analysis_objects)
  X_train <- as.matrix(analysis_objects$X_train)
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
    y_rep <- simulate_contam_replicated_outcome(
      draw = draw,
      X = X_train,
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
    statistics$model <- "M2_contaminated_normal"
    statistics
  }))
  observed_statistics$model <- "M2_contaminated_normal"
  summary <- summarise_ppc_statistics(observed_statistics, replicated_statistics)
  summary$model <- "M2_contaminated_normal"
  list(
    summary = summary,
    replicated_statistics = replicated_statistics,
    observed_statistics = observed_statistics
  )
}

write_contam_comparison_tables <- function(heldout_metrics, ppc_summary,
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

run_contam_m2_experiment <- function(
    analysis_objects_path = here::here("data", "processed", "analysis_objects.rds"),
    model_file = here::here("models", "m2_contaminated_normal.jags"),
    mcmc_output_dir = here::here("output", "mcmc"),
    tables_dir = here::here("output", "tables"),
    figures_dir = here::here("output", "figures")) {
  settings <- list(
    n_chains = student_t_read_positive_integer("CONTAM_M2_N_CHAINS", 4L),
    n_adapt = student_t_read_positive_integer("CONTAM_M2_N_ADAPT", 2000L),
    n_burn = student_t_read_positive_integer("CONTAM_M2_N_BURN", 5000L),
    n_iter = student_t_read_positive_integer("CONTAM_M2_N_ITER", 10000L),
    thin = student_t_read_positive_integer("CONTAM_M2_THIN", 1L),
    seed = student_t_read_positive_integer("CONTAM_M2_SEED", 20260703L),
    chunk_iter = student_t_read_positive_integer("CONTAM_M2_CHUNK_ITER", 1000L),
    max_heldout_draws = student_t_read_positive_integer("CONTAM_M2_MAX_HELDOUT_DRAWS", 2000L),
    max_ppc_draws = student_t_read_positive_integer("CONTAM_M2_MAX_PPC_DRAWS", 1000L),
    force = student_t_read_flag("CONTAM_M2_FORCE", default = FALSE),
    p_outlier_alpha = contam_read_positive_number("CONTAM_M2_P_OUTLIER_ALPHA", 2),
    p_outlier_beta = contam_read_positive_number("CONTAM_M2_P_OUTLIER_BETA", 20),
    sigma_extra_prior_scale = contam_read_positive_number("CONTAM_M2_SIGMA_EXTRA_PRIOR_SCALE", 1),
    label = contam_read_label("CONTAM_M2_LABEL", "beta2_20_extra1")
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

  analysis_objects <- readRDS(analysis_objects_path)
  validate_student_t_analysis_objects(analysis_objects)
  parameters <- c(
    "alpha", "beta", "sigma_core", "sigma_outlier",
    "p_outlier", "sd_state", "sd_zip", "u", "v"
  )
  output_stem <- paste0("m2_contaminated_normal_", settings$label)
  model_label <- paste0("M2_contaminated_normal_", settings$label)
  draws_path <- file.path(mcmc_output_dir, paste0(output_stem, "_draws.rds"))
  partial_draws_path <- file.path(mcmc_output_dir, paste0(output_stem, "_draw_chunks_partial.rds"))
  runtime_path <- file.path(tables_dir, paste0(output_stem, "_runtime.csv"))

  readr::write_csv(
    data.frame(setting = names(settings), value = unlist(settings), row.names = NULL),
    file.path(tables_dir, paste0(output_stem, "_settings.csv"))
  )
  readr::write_csv(
    data.frame(
      residual_distribution = "Contaminated normal",
      p_outlier_prior = paste0(
        "Beta(", settings$p_outlier_alpha, ", ", settings$p_outlier_beta, ")"
      ),
      sigma_extra_prior_scale = settings$sigma_extra_prior_scale,
      label = settings$label,
      stringsAsFactors = FALSE
    ),
    file.path(tables_dir, paste0(output_stem, "_terms.csv"))
  )

  freshness_mtime <- max(
    file.info(c(model_file, here::here("scripts", "16_m2_contaminated_normal_experiment.R")))$mtime,
    na.rm = TRUE
  )
  draws_current <- file.exists(draws_path) &&
    file.info(draws_path)$mtime >= freshness_mtime

  fit_started_at <- Sys.time()
  if (draws_current && !settings$force) {
    student_t_timestamp_message(
      "Skipping contaminated-normal M2 refit because current draws already exist at ",
      draws_path, ". Set CONTAM_M2_FORCE=true to refit."
    )
    draws <- readRDS(draws_path)
    fit_status <- "skipped_existing_draws"
    fit_seconds <- NA_real_
  } else {
    student_t_timestamp_message(
      "Starting contaminated-normal M2 fit. Settings: chains=", settings$n_chains,
      ", adapt=", settings$n_adapt,
      ", burn=", settings$n_burn,
      ", iter=", settings$n_iter,
      ", thin=", settings$thin,
      ", chunk=", settings$chunk_iter,
      ", p_outlier_prior=Beta(", settings$p_outlier_alpha,
      ", ", settings$p_outlier_beta, "), sigma_extra_prior_scale=",
      settings$sigma_extra_prior_scale,
      "."
    )
    fit <- fit_contam_m2(
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
      "Saved contaminated-normal M2 draws to ", draws_path,
      " after ", round(fit_seconds, 1), " seconds."
    )
  }

  draw_matrix <- as.matrix(posterior::as_draws_matrix(draws))
  validate_contam_draw_matrix(draw_matrix, analysis_objects)

  convergence_started_at <- Sys.time()
  convergence <- summarise_student_t_convergence(draws)
  convergence_path <- file.path(tables_dir, paste0(output_stem, "_posterior_summary.csv"))
  readr::write_csv(convergence, convergence_path)
  convergence_seconds <- as.numeric(difftime(Sys.time(), convergence_started_at, units = "secs"))
  student_t_timestamp_message(
    "Saved contaminated-normal M2 convergence summary to ", convergence_path,
    ". Converged ", sum(convergence$converged), "/", nrow(convergence), " parameters."
  )

  set.seed(settings$seed + 101L)
  heldout_draw_ids <- sort(sample(
    seq_len(nrow(draw_matrix)),
    size = min(settings$max_heldout_draws, nrow(draw_matrix)),
    replace = FALSE
  ))
  heldout <- evaluate_contam_holdout(
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
  ppc <- run_contam_ppc(
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
    model = "M2 contaminated-normal residual",
    figure_path = ppc_figure_path
  )
  write_contam_comparison_tables(heldout$metrics, ppc$summary, tables_dir, output_stem)

  runtime <- data.frame(
    stage = c("fit_or_load_draws", "convergence_summary", "heldout_ppc"),
    status = c(fit_status, "completed", "completed"),
    elapsed_seconds = c(fit_seconds, convergence_seconds, NA_real_),
    stringsAsFactors = FALSE
  )
  readr::write_csv(runtime, runtime_path)
  student_t_timestamp_message(
    "Contaminated-normal M2 experiment complete. Key outputs: ", convergence_path,
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
  run_contam_m2_experiment()
}
