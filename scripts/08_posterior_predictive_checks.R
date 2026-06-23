# Posterior predictive checks for the M0, M1, and M2 fitted models.
#
# This script implements point 16 of IMPLEMENTATION_GUIDE.md.  It compares
# replicated, standardised log-price outcomes against the observed training
# outcomes.  The checks are diagnostics of specific distributional and
# geographic features; their tail probabilities are not mechanical pass/fail
# tests.

required_ppc_packages <- c("ggplot2", "here", "posterior", "readr")
missing_ppc_packages <- required_ppc_packages[!vapply(
  required_ppc_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_ppc_packages) > 0L) {
  stop(
    "Install the package(s) required for posterior predictive checks: ",
    paste(missing_ppc_packages, collapse = ", "),
    call. = FALSE
  )
}

# Reuse the M2 predictive simulator so posterior predictions and posterior
# predictive checks generate outcomes from exactly the same fitted hierarchy.
prediction_script <- here::here("scripts", "10_posterior_predictions.R")
if (!file.exists(prediction_script)) {
  stop("Posterior prediction helper not found: ", prediction_script, call. = FALSE)
}
source(prediction_script, local = FALSE)

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

required_indexed_parameter_names <- function(prefix, count) {
  paste0(prefix, "[", seq_len(count), "]")
}

model_parameter_names <- function(model, n_predictors, n_states, n_zips) {
  model <- match.arg(model, c("M0", "M1", "M2"))
  parameters <- c(
    "alpha", "sigma",
    required_indexed_parameter_names("beta", n_predictors)
  )
  if (model %in% c("M1", "M2")) {
    parameters <- c(
      parameters,
      required_indexed_parameter_names("u", n_states)
    )
  }
  if (identical(model, "M2")) {
    parameters <- c(
      parameters,
      "sd_zip",
      required_indexed_parameter_names("v", n_zips)
    )
  }
  parameters
}

validate_ppc_inputs <- function(draw_matrix, analysis_objects, model) {
  model <- match.arg(model, c("M0", "M1", "M2"))
  required_objects <- c("jags_train", "X_train", "state_levels", "zip_levels")
  missing_objects <- setdiff(required_objects, names(analysis_objects))
  if (length(missing_objects) > 0L) {
    stop(
      "Analysis objects are missing: ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  jags_train <- analysis_objects$jags_train
  required_jags_fields <- c("y", "X", "state", "zip")
  missing_jags_fields <- setdiff(required_jags_fields, names(jags_train))
  if (length(missing_jags_fields) > 0L) {
    stop(
      "JAGS training data are missing: ",
      paste(missing_jags_fields, collapse = ", "),
      call. = FALSE
    )
  }

  X <- as.matrix(analysis_objects$X_train)
  y <- as.numeric(jags_train$y)
  state <- as.integer(jags_train$state)
  zip <- as.integer(jags_train$zip)
  if (nrow(X) < 2L ||
      length(y) != nrow(X) ||
      length(state) != nrow(X) ||
      length(zip) != nrow(X) ||
      any(!is.finite(X)) ||
      any(!is.finite(y)) ||
      anyNA(state) ||
      anyNA(zip) ||
      any(state < 1L) ||
      any(zip < 1L)) {
    stop("Training outcomes, covariates, and geography indices are inconsistent.", call. = FALSE)
  }
  if (length(unique(zip)) < 2L) {
    stop("Posterior predictive checks require at least two observed ZIP codes.", call. = FALSE)
  }
  if (max(state) > length(analysis_objects$state_levels) ||
      max(zip) > length(analysis_objects$zip_levels)) {
    stop("Training geography indices exceed the saved geography lookups.", call. = FALSE)
  }
  if (!identical(unname(X), unname(as.matrix(jags_train$X)))) {
    stop("X_train does not match the JAGS training design matrix.", call. = FALSE)
  }

  draw_matrix <- as.matrix(draw_matrix)
  if (nrow(draw_matrix) == 0L || ncol(draw_matrix) == 0L ||
      is.null(colnames(draw_matrix)) || any(!is.finite(draw_matrix))) {
    stop("Posterior draws are empty, unnamed, or contain non-finite values.", call. = FALSE)
  }
  required_parameters <- model_parameter_names(
    model = model,
    n_predictors = ncol(X),
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
  if (any(draw_matrix[, "sigma"] <= 0) ||
      (identical(model, "M2") && any(draw_matrix[, "sd_zip"] <= 0))) {
    stop(model, " draws contain a non-positive scale parameter.", call. = FALSE)
  }

  invisible(TRUE)
}

# State summaries are selected from the sorted state lookup, not their observed
# outcomes.  Supplying `selected_states` fixes a project-specific alternative.
select_ppc_state_ids <- function(state_levels, selected_states = NULL,
                                 max_states = 4L) {
  state_levels <- as.character(state_levels)
  max_states <- validate_positive_integer(max_states, "max_states")
  if (length(state_levels) == 0L || anyDuplicated(state_levels)) {
    stop("state_levels must contain unique state labels.", call. = FALSE)
  }

  if (is.null(selected_states)) {
    n_selected <- min(max_states, length(state_levels))
    selected_ids <- unique(as.integer(round(seq(
      from = 1L,
      to = length(state_levels),
      length.out = n_selected
    ))))
  } else {
    selected_ids <- match(as.character(selected_states), state_levels)
    if (anyNA(selected_ids)) {
      stop("Each selected state must occur in state_levels.", call. = FALSE)
    }
    selected_ids <- unique(as.integer(selected_ids))
  }
  selected_ids
}

calculate_ppc_statistics <- function(y, state_id, zip_id, state_levels,
                                     selected_state_ids) {
  y <- as.numeric(y)
  state_id <- as.integer(state_id)
  zip_id <- as.integer(zip_id)
  zip_means <- tapply(y, zip_id, mean)

  overall_statistics <- data.frame(
    statistic = c(
      "overall_median_standardized_log_price",
      "overall_iqr_standardized_log_price",
      "overall_q90_standardized_log_price",
      "zip_mean_variance_standardized_log_price"
    ),
    scope = rep("overall", 4L),
    value = c(
      stats::median(y),
      stats::IQR(y),
      stats::quantile(y, probs = 0.90, names = FALSE),
      stats::var(zip_means)
    ),
    stringsAsFactors = FALSE
  )

  state_statistics <- data.frame(
    statistic = rep("state_median_standardized_log_price", length(selected_state_ids)),
    scope = state_levels[selected_state_ids],
    value = vapply(
      selected_state_ids,
      function(state) stats::median(y[state_id == state]),
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
  rbind(overall_statistics, state_statistics)
}

simulate_replicated_outcome <- function(draw, X, state_id, zip_id, model) {
  model <- match.arg(model, c("M0", "M1", "M2"))
  draw_names <- names(draw)
  draw <- as.numeric(draw)
  if (is.null(draw_names)) {
    stop("draw must be a named vector of posterior parameters.", call. = FALSE)
  }
  names(draw) <- draw_names

  if (identical(model, "M2")) {
    return(posterior_predict_m2(
      draw = draw,
      X = X,
      state_id = state_id,
      zip_id = zip_id,
      include_observed_zip = TRUE
    ))
  }

  beta_names <- required_indexed_parameter_names("beta", ncol(X))
  mu <- as.vector(draw["alpha"] + as.matrix(X) %*% draw[beta_names])
  if (identical(model, "M1")) {
    state_names <- required_indexed_parameter_names("u", max(state_id))
    mu <- mu + draw[state_names][state_id]
  }
  stats::rnorm(nrow(X), mean = mu, sd = draw["sigma"])
}

generate_replications <- function(draw_ids, draw_matrix, X, state, zip,
                                  model = "M2") {
  model <- match.arg(model, c("M0", "M1", "M2"))
  lapply(draw_ids, function(index) {
    draw <- as.numeric(draw_matrix[index, , drop = TRUE])
    names(draw) <- colnames(draw_matrix)
    simulate_replicated_outcome(
      draw = draw,
      X = X,
      state_id = state,
      zip_id = zip,
      model = model
    )
  })
}

summarise_ppc_statistics <- function(observed_statistics, replicated_statistics) {
  check_key <- function(statistics) {
    paste(statistics$statistic, statistics$scope, sep = "\r")
  }
  observed_keys <- check_key(observed_statistics)
  replicated_keys <- check_key(replicated_statistics)

  summary_rows <- lapply(seq_len(nrow(observed_statistics)), function(index) {
    values <- replicated_statistics$value[replicated_keys == observed_keys[index]]
    observed_value <- observed_statistics$value[index]
    data.frame(
      statistic = observed_statistics$statistic[index],
      scope = observed_statistics$scope[index],
      observed = observed_value,
      replicated_mean = mean(values),
      replicated_median = stats::median(values),
      replicated_lower_95 = stats::quantile(values, 0.025, names = FALSE),
      replicated_upper_95 = stats::quantile(values, 0.975, names = FALSE),
      posterior_predictive_tail_probability = mean(values >= observed_value),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, summary_rows)
}

ppc_check_label <- function(statistic, scope) {
  labels <- c(
    overall_median_standardized_log_price = "Overall median",
    overall_iqr_standardized_log_price = "Overall IQR",
    overall_q90_standardized_log_price = "Overall 90th percentile",
    zip_mean_variance_standardized_log_price = "ZIP-mean variance",
    state_median_standardized_log_price = "State median"
  )
  label <- unname(labels[statistic])
  ifelse(statistic == "state_median_standardized_log_price",
         paste(label, "-", scope), label)
}

save_ppc_plot <- function(replicated_statistics, observed_statistics,
                          model, figure_path) {
  replicated_statistics$check <- ppc_check_label(
    replicated_statistics$statistic,
    replicated_statistics$scope
  )
  observed_statistics$check <- ppc_check_label(
    observed_statistics$statistic,
    observed_statistics$scope
  )

  check_plot <- ggplot2::ggplot(
    replicated_statistics,
    ggplot2::aes(x = value)
  ) +
    ggplot2::geom_histogram(bins = 30L, fill = "#9ECAE1", colour = "white") +
    ggplot2::geom_vline(
      data = observed_statistics,
      ggplot2::aes(xintercept = value),
      colour = "#CB181D",
      linewidth = 0.9
    ) +
    ggplot2::facet_wrap(~check, scales = "free", ncol = 2L) +
    ggplot2::labs(
      title = paste0(model, " posterior predictive checks"),
      subtitle = "Histograms show replicated statistics; red lines show observed statistics.",
      x = "Statistic value (standardized log-price scale)",
      y = "Replicated data sets"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  ggplot2::ggsave(
    filename = figure_path,
    plot = check_plot,
    width = 11,
    height = max(5, 2.6 * ceiling(length(unique(replicated_statistics$check)) / 2)),
    dpi = 180
  )
}

run_model_ppc <- function(model, draw_matrix, analysis_objects, draw_ids,
                          selected_state_ids) {
  validate_ppc_inputs(draw_matrix, analysis_objects, model)
  jags_train <- analysis_objects$jags_train
  observed_statistics <- calculate_ppc_statistics(
    y = jags_train$y,
    state_id = jags_train$state,
    zip_id = jags_train$zip,
    state_levels = analysis_objects$state_levels,
    selected_state_ids = selected_state_ids
  )
  y_rep_list <- generate_replications(
    draw_ids = draw_ids,
    draw_matrix = draw_matrix,
    X = analysis_objects$X_train,
    state = jags_train$state,
    zip = jags_train$zip,
    model = model
  )
  replicated_statistics <- do.call(rbind, lapply(seq_along(y_rep_list), function(index) {
    statistics <- calculate_ppc_statistics(
      y = y_rep_list[[index]],
      state_id = jags_train$state,
      zip_id = jags_train$zip,
      state_levels = analysis_objects$state_levels,
      selected_state_ids = selected_state_ids
    )
    statistics$draw_id <- draw_ids[index]
    statistics
  }))
  replicated_statistics$model <- model
  observed_statistics$model <- model
  summary <- summarise_ppc_statistics(observed_statistics, replicated_statistics)
  summary$model <- model

  list(
    summary = summary,
    replicated_statistics = replicated_statistics,
    observed_statistics = observed_statistics,
    draw_ids = draw_ids
  )
}

run_posterior_predictive_checks <- function(
    draws_paths = list(
      M0 = here::here("output", "mcmc", "m0_draws.rds"),
      M1 = here::here("output", "mcmc", "m1_draws.rds"),
      M2 = here::here("output", "mcmc", "m2_draws.rds")
    ),
    analysis_objects_path = here::here("data", "processed", "analysis_objects.rds"),
    tables_dir = here::here("output", "tables"),
    figures_dir = here::here("output", "figures"),
    max_ppc_draws = 400L,
    selected_states = NULL,
    seed = 20260621L) {
  max_ppc_draws <- validate_positive_integer(max_ppc_draws, "max_ppc_draws")
  seed <- validate_positive_integer(seed, "seed")
  if (is.null(names(draws_paths)) ||
      length(draws_paths) == 0L ||
      anyDuplicated(names(draws_paths)) ||
      !all(names(draws_paths) %in% c("M0", "M1", "M2"))) {
    stop("draws_paths must be a named list with unique M0, M1, and/or M2 entries.", call. = FALSE)
  }
  if (!file.exists(analysis_objects_path)) {
    stop(
      "Prepared analysis objects not found: ", analysis_objects_path,
      ". Run scripts/03_prepare_jags_data.R first.",
      call. = FALSE
    )
  }
  draw_files <- unlist(draws_paths, use.names = TRUE)
  missing_draws <- draw_files[!file.exists(draw_files)]
  if (length(missing_draws) > 0L) {
    stop(
      "Posterior draws not found for ",
      paste(names(missing_draws), collapse = ", "),
      ". Run scripts/05_fit_models.R after convergence checks first.",
      call. = FALSE
    )
  }

  analysis_objects <- readRDS(analysis_objects_path)
  selected_state_ids <- select_ppc_state_ids(
    analysis_objects$state_levels,
    selected_states = selected_states
  )
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  results <- lapply(seq_along(draws_paths), function(index) {
    model <- names(draws_paths)[index]
    draw_matrix <- posterior::as_draws_matrix(readRDS(draws_paths[[index]]))
    set.seed(seed + index - 1L)
    draw_ids <- sample(
      seq_len(nrow(draw_matrix)),
      size = min(max_ppc_draws, nrow(draw_matrix)),
      replace = FALSE
    )
    run_model_ppc(
      model = model,
      draw_matrix = draw_matrix,
      analysis_objects = analysis_objects,
      draw_ids = draw_ids,
      selected_state_ids = selected_state_ids
    )
  })
  names(results) <- names(draws_paths)

  summary <- do.call(rbind, lapply(results, `[[`, "summary"))
  replicated_statistics <- do.call(
    rbind,
    lapply(results, `[[`, "replicated_statistics")
  )
  summary_path <- file.path(tables_dir, "posterior_predictive_check_summary.csv")
  replicated_path <- file.path(tables_dir, "posterior_predictive_check_replicates.csv")
  readr::write_csv(summary, summary_path)
  readr::write_csv(replicated_statistics, replicated_path)

  figure_paths <- vapply(names(results), function(model) {
    figure_path <- file.path(
      figures_dir,
      paste0(tolower(model), "_posterior_predictive_checks.png")
    )
    save_ppc_plot(
      replicated_statistics = results[[model]]$replicated_statistics,
      observed_statistics = results[[model]]$observed_statistics,
      model = model,
      figure_path = figure_path
    )
    figure_path
  }, character(1))

  message(
    "Saved posterior predictive checks for ", paste(names(results), collapse = ", "),
    ". Summary: ", summary_path,
    "; replicated statistics: ", replicated_path,
    "; figures: ", paste(figure_paths, collapse = ", ")
  )
  invisible(list(
    summary = summary,
    replicated_statistics = replicated_statistics,
    selected_states = analysis_objects$state_levels[selected_state_ids],
    draw_ids = lapply(results, `[[`, "draw_ids"),
    summary_path = summary_path,
    replicated_path = replicated_path,
    figure_paths = figure_paths
  ))
}

if (sys.nframe() == 0L) {
  run_posterior_predictive_checks()
}
