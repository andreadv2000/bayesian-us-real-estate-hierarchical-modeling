# Calculate and record DIC plus held-out performance for the three fitted
# JAGS models.
#
# This file defines functions only. Source it from scripts/05_fit_models.R while
# the live JAGS model objects are available: JAGS model objects cannot be
# restored from an RDS file and then used with rjags::dic.samples().

required_comparison_packages <- c(
  "dplyr", "here", "posterior", "readr", "rjags", "tibble"
)
missing_comparison_packages <- required_comparison_packages[!vapply(
  required_comparison_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_comparison_packages) > 0L) {
  stop(
    "Install the package(s) required for DIC comparison: ",
    paste(missing_comparison_packages, collapse = ", "),
    call. = FALSE
  )
}

validate_dic_iterations <- function(n_iter) {
  valid <- is.numeric(n_iter) &&
    length(n_iter) == 1L &&
    is.finite(n_iter) &&
    n_iter >= 1L &&
    n_iter %% 1 == 0
  if (!valid) {
    stop("n_iter must be one positive integer.", call. = FALSE)
  }
  as.integer(n_iter)
}

sample_dic <- function(jags_model, n_iter = 5000L) {
  if (!inherits(jags_model, "jags")) {
    stop("jags_model must be a live rjags model object.", call. = FALSE)
  }

  rjags::dic.samples(
    model = jags_model,
    n.iter = validate_dic_iterations(n_iter),
    type = "pD",
    progress.bar = "none"
  )
}

summarise_dic <- function(jags_model, n_iter = 5000L) {
  dic_object <- sample_dic(jags_model, n_iter = n_iter)
  d_bar <- sum(as.numeric(dic_object$deviance))
  p_d <- sum(as.numeric(dic_object$penalty))

  if (!all(is.finite(c(d_bar, p_d)))) {
    stop("DIC calculation returned a non-finite deviance or penalty.", call. = FALSE)
  }

  tibble::tibble(Dbar = d_bar, pD = p_d, DIC = d_bar + p_d)
}

calculate_and_write_dic <- function(fit_m0, fit_m1, fit_m2,
                                    n_iter = 5000L,
                                    table_path = here::here(
                                      "output", "tables", "model_comparison_dic.csv"
                                    ),
                                    log_path = here::here(
                                      "output", "logs", "model_comparison_dic.log"
                                    )) {
  fit_objects <- list(
    "M0: pooled" = fit_m0,
    "M1: state" = fit_m1,
    "M2: state + ZIP" = fit_m2
  )
  invalid_fits <- names(fit_objects)[!vapply(
    fit_objects,
    function(fit) is.list(fit) && inherits(fit$model, "jags"),
    FUN.VALUE = logical(1)
  )]
  if (length(invalid_fits) > 0L) {
    stop(
      "DIC requires live JAGS fits for: ",
      paste(invalid_fits, collapse = ", "),
      ". Run scripts/05_fit_models.R in this R session.",
      call. = FALSE
    )
  }

  n_iter <- validate_dic_iterations(n_iter)
  dic_objects <- lapply(fit_objects, function(fit) {
    sample_dic(fit$model, n_iter = n_iter)
  })
  dic_table <- dplyr::bind_rows(Map(
    function(model_name, dic_object) {
      d_bar <- sum(as.numeric(dic_object$deviance))
      p_d <- sum(as.numeric(dic_object$penalty))
      tibble::tibble(model = model_name, Dbar = d_bar, pD = p_d, DIC = d_bar + p_d)
    },
    names(dic_objects),
    dic_objects
  )) |>
    dplyr::arrange(DIC) |>
    dplyr::mutate(delta_DIC = DIC - min(DIC))

  dic_components <- unlist(dic_table[c("Dbar", "pD", "DIC")], use.names = FALSE)
  if (any(!is.finite(dic_components))) {
    stop("DIC comparison returned invalid values.", call. = FALSE)
  }

  dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(dic_table, table_path)

  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  log_lines <- c(
    "DIC calculation from rjags::dic.samples(type = \"pD\").",
    paste("rjags version:", as.character(utils::packageVersion("rjags"))),
    paste("DIC iterations per model:", n_iter),
    "Dbar is posterior mean deviance; pD is the effective parameter count; DIC = Dbar + pD.",
    "All models were fitted on the same frozen training observations in scripts/05_fit_models.R.",
    "",
    "Printed DIC object for M0: pooled:",
    utils::capture.output(print(dic_objects[["M0: pooled"]]))
  )
  writeLines(log_lines, con = log_path, useBytes = TRUE)

  list(table = dic_table, dic_objects = dic_objects)
}

# -------------------------------------------------------------------------
# Held-out predictive evaluation (implementation guide point 17)
# -------------------------------------------------------------------------
#
# DIC needs live JAGS model objects and is therefore calculated by
# scripts/05_fit_models.R.  Held-out evaluation instead uses the saved draws,
# so this section can be run only after production chains have passed the
# convergence checks.  It deliberately uses the same test rows, draw indices,
# and standard-normal innovations for M0, M1, and M2.

validate_positive_integer <- function(value, name) {
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

indexed_parameter_names <- function(prefix, count) {
  paste0(prefix, "[", seq_len(count), "]")
}

required_prediction_parameters <- function(model, n_predictors, n_states,
                                           n_zips) {
  model <- match.arg(model, c("M0", "M1", "M2"))
  parameters <- c("alpha", "sigma", indexed_parameter_names("beta", n_predictors))
  if (model %in% c("M1", "M2")) {
    parameters <- c(parameters, "sd_state", indexed_parameter_names("u", n_states))
  }
  if (identical(model, "M2")) {
    parameters <- c(parameters, "sd_zip", indexed_parameter_names("v", n_zips))
  }
  parameters
}

validate_heldout_inputs <- function(model, draw_matrix, analysis_objects) {
  model <- match.arg(model, c("M0", "M1", "M2"))
  required_objects <- c("test", "X_test", "state_levels", "zip_levels")
  missing_objects <- setdiff(required_objects, names(analysis_objects))
  if (length(missing_objects) > 0L) {
    stop(
      "Analysis objects are missing: ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  test <- analysis_objects$test
  required_test_columns <- c("y", "state_id", "zip_id_num")
  missing_test_columns <- setdiff(required_test_columns, names(test))
  if (length(missing_test_columns) > 0L) {
    stop(
      "Held-out test data are missing: ",
      paste(missing_test_columns, collapse = ", "),
      call. = FALSE
    )
  }

  X_test <- as.matrix(analysis_objects$X_test)
  test_y <- as.numeric(test$y)
  state_id <- as.integer(test$state_id)
  zip_id <- as.integer(test$zip_id_num)
  if (nrow(X_test) == 0L ||
      nrow(X_test) != length(test_y) ||
      nrow(X_test) != length(state_id) ||
      nrow(X_test) != length(zip_id) ||
      any(!is.finite(X_test)) || any(!is.finite(test_y)) ||
      anyNA(state_id) || anyNA(zip_id) ||
      any(state_id < 1L) || any(zip_id < 1L) ||
      any(state_id > length(analysis_objects$state_levels)) ||
      any(zip_id > length(analysis_objects$zip_levels))) {
    stop(
      "Held-out outcomes, covariates, or geography indices are inconsistent.",
      call. = FALSE
    )
  }

  draw_matrix <- as.matrix(draw_matrix)
  if (nrow(draw_matrix) == 0L || ncol(draw_matrix) == 0L ||
      is.null(colnames(draw_matrix)) || any(!is.finite(draw_matrix))) {
    stop("Posterior draws are empty, unnamed, or contain non-finite values.", call. = FALSE)
  }
  required_parameters <- required_prediction_parameters(
    model = model,
    n_predictors = ncol(X_test),
    n_states = length(analysis_objects$state_levels),
    n_zips = length(analysis_objects$zip_levels)
  )
  missing_parameters <- setdiff(required_parameters, colnames(draw_matrix))
  if (length(missing_parameters) > 0L) {
    stop(
      model, " draws are missing required parameter(s): ",
      paste(missing_parameters, collapse = ", "),
      call. = FALSE
    )
  }
  if (any(draw_matrix[, "sigma"] <= 0)) {
    stop(model, " draws contain a non-positive residual SD.", call. = FALSE)
  }
  if (identical(model, "M2") && any(draw_matrix[, "sd_zip"] <= 0)) {
    stop("M2 draws contain a non-positive ZIP SD.", call. = FALSE)
  }

  invisible(TRUE)
}

linear_predictor_for_draw <- function(model, draw, X_test, state_id, zip_id) {
  model <- match.arg(model, c("M0", "M1", "M2"))
  beta_names <- indexed_parameter_names("beta", ncol(X_test))
  mu <- as.vector(draw["alpha"] + X_test %*% draw[beta_names])

  if (model %in% c("M1", "M2")) {
    state_names <- indexed_parameter_names("u", max(state_id))
    mu <- mu + unname(draw[state_names])[state_id]
  }
  if (identical(model, "M2")) {
    zip_names <- indexed_parameter_names("v", max(zip_id))
    mu <- mu + unname(draw[zip_names])[zip_id]
  }
  mu
}

# Calculate a running log-sum-exp without storing an N_test-by-M likelihood
# matrix.  That keeps the held-out calculation stable and limits memory use.
update_log_density_accumulator <- function(log_max, scaled_sum, log_density) {
  new_log_max <- pmax(log_max, log_density)
  list(
    log_max = new_log_max,
    scaled_sum =
      exp(log_max - new_log_max) * scaled_sum +
      exp(log_density - new_log_max)
  )
}

evaluate_model_on_holdout <- function(model, draw_matrix, analysis_objects,
                                      draw_ids, standard_normal_draws) {
  model <- match.arg(model, c("M0", "M1", "M2"))
  validate_heldout_inputs(model, draw_matrix, analysis_objects)
  draw_matrix <- as.matrix(draw_matrix)
  draw_ids <- as.integer(draw_ids)
  n_draws <- validate_positive_integer(length(draw_ids), "draw_ids length")
  if (anyNA(draw_ids) || any(draw_ids < 1L) ||
      any(draw_ids > nrow(draw_matrix)) || anyDuplicated(draw_ids)) {
    stop("draw_ids must be unique, valid rows of draw_matrix.", call. = FALSE)
  }

  test <- analysis_objects$test
  X_test <- as.matrix(analysis_objects$X_test)
  test_y <- as.numeric(test$y)
  state_id <- as.integer(test$state_id)
  zip_id <- as.integer(test$zip_id_num)
  n_test <- length(test_y)

  standard_normal_draws <- as.matrix(standard_normal_draws)
  if (!identical(dim(standard_normal_draws), c(n_test, n_draws)) ||
      any(!is.finite(standard_normal_draws))) {
    stop(
      "standard_normal_draws must be a finite N_test-by-N_draws matrix.",
      call. = FALSE
    )
  }

  posterior_mean_sum <- numeric(n_test)
  predictive_draws <- matrix(NA_real_, nrow = n_test, ncol = n_draws)
  log_density_max <- numeric(n_test)
  log_density_sum <- numeric(n_test)

  for (draw_position in seq_len(n_draws)) {
    draw <- as.numeric(draw_matrix[draw_ids[draw_position], , drop = TRUE])
    names(draw) <- colnames(draw_matrix)
    mu <- linear_predictor_for_draw(model, draw, X_test, state_id, zip_id)
    posterior_mean_sum <- posterior_mean_sum + mu
    predictive_draws[, draw_position] <-
      mu + draw["sigma"] * standard_normal_draws[, draw_position]

    log_density <- stats::dnorm(test_y, mean = mu, sd = draw["sigma"], log = TRUE)
    if (draw_position == 1L) {
      log_density_max <- log_density
      log_density_sum <- rep(1, n_test)
    } else {
      accumulator <- update_log_density_accumulator(
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
  log_predictive_density <-
    log_density_max + log(log_density_sum) - log(n_draws)
  covered_95 <- test_y >= lower_95 & test_y <= upper_95

  predictions <- data.frame(
    model = model,
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
    model = model,
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

run_heldout_predictive_evaluation <- function(
    draws_paths = list(
      M0 = here::here("output", "mcmc", "m0_draws.rds"),
      M1 = here::here("output", "mcmc", "m1_draws.rds"),
      M2 = here::here("output", "mcmc", "m2_draws.rds")
    ),
    analysis_objects_path = here::here("data", "processed", "analysis_objects.rds"),
    tables_dir = here::here("output", "tables"),
    max_predictive_draws = 2000L,
    seed = 20260621L) {
  expected_models <- c("M0", "M1", "M2")
  if (is.null(names(draws_paths)) || !identical(names(draws_paths), expected_models)) {
    stop("draws_paths must be a named list ordered as M0, M1, M2.", call. = FALSE)
  }
  max_predictive_draws <- validate_positive_integer(
    max_predictive_draws,
    "max_predictive_draws"
  )
  seed <- validate_positive_integer(seed, "seed")
  if (!file.exists(analysis_objects_path)) {
    stop(
      "Prepared analysis objects not found: ", analysis_objects_path,
      ". Run scripts/03_prepare_jags_data.R first.",
      call. = FALSE
    )
  }
  missing_draw_files <- unlist(draws_paths, use.names = TRUE)
  missing_draw_files <- missing_draw_files[!file.exists(missing_draw_files)]
  if (length(missing_draw_files) > 0L) {
    stop(
      "Converged posterior draws not found for ",
      paste(names(missing_draw_files), collapse = ", "),
      ". Run scripts/05_fit_models.R and scripts/06_diagnostics.R first.",
      call. = FALSE
    )
  }

  analysis_objects <- readRDS(analysis_objects_path)
  draw_matrices <- lapply(draws_paths, function(path) {
    as.matrix(posterior::as_draws_matrix(readRDS(path)))
  })
  common_available_draws <- min(vapply(draw_matrices, nrow, integer(1)))
  n_draws <- min(max_predictive_draws, common_available_draws)
  if (n_draws < 1L) {
    stop("At least one posterior draw is required for every model.", call. = FALSE)
  }

  set.seed(seed)
  draw_ids <- sort(sample(seq_len(common_available_draws), size = n_draws))
  n_test <- nrow(as.matrix(analysis_objects$X_test))
  standard_normal_draws <- matrix(
    stats::rnorm(n_test * n_draws),
    nrow = n_test,
    ncol = n_draws
  )

  results <- Map(
    function(model, draw_matrix) {
      evaluate_model_on_holdout(
        model = model,
        draw_matrix = draw_matrix,
        analysis_objects = analysis_objects,
        draw_ids = draw_ids,
        standard_normal_draws = standard_normal_draws
      )
    },
    names(draw_matrices),
    draw_matrices
  )
  metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
  predictions <- do.call(rbind, lapply(results, `[[`, "predictions"))
  metadata <- data.frame(
    seed = seed,
    n_test = n_test,
    common_available_draws = common_available_draws,
    selected_posterior_draws = n_draws,
    draw_indices = paste(draw_ids, collapse = ";"),
    stringsAsFactors = FALSE
  )

  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  metrics_path <- file.path(tables_dir, "heldout_predictive_metrics.csv")
  predictions_path <- file.path(tables_dir, "heldout_predictive_predictions.csv")
  metadata_path <- file.path(tables_dir, "heldout_predictive_evaluation_metadata.csv")
  readr::write_csv(metrics, metrics_path)
  readr::write_csv(predictions, predictions_path)
  readr::write_csv(metadata, metadata_path)

  message(
    "Saved held-out predictive metrics for M0, M1, and M2 to ", metrics_path,
    ". The evaluation used the same ", n_draws,
    " posterior-draw indices and predictive innovations for every model."
  )
  invisible(list(
    metrics = metrics,
    predictions = predictions,
    metadata = metadata,
    draw_ids = draw_ids,
    metrics_path = metrics_path,
    predictions_path = predictions_path,
    metadata_path = metadata_path
  ))
}

if (sys.nframe() == 0L) {
  run_heldout_predictive_evaluation()
}
