# Generate posterior predictive summaries for the converged M2 fit.
#
# This script implements point 15 of IMPLEMENTATION_GUIDE.md. It uses the
# saved parameter draws from M2 rather than monitoring a replicated outcome for
# every training listing in JAGS. The three report scenarios hold the property
# characteristics fixed at a representative in-range profile so that the
# difference between a low- and high-adjusted-price ZIP is interpretable as a
# geographic comparison. A third scenario predicts a new ZIP in the
# high-adjusted-price state by sampling its unobserved ZIP effect.

required_prediction_packages <- c("ggplot2", "here", "posterior", "readr")
missing_prediction_packages <- required_prediction_packages[!vapply(
  required_prediction_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_prediction_packages) > 0L) {
  stop(
    "Install the package(s) required for posterior prediction: ",
    paste(missing_prediction_packages, collapse = ", "),
    call. = FALSE
  )
}

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

is_heteroskedastic_m2_draws <- function(draw_matrix) {
  "gamma0" %in% colnames(draw_matrix)
}

validate_prediction_inputs <- function(draw_matrix, analysis_objects) {
  required_objects <- c(
    "train", "X_train", "predictor_names", "state_levels", "zip_levels",
    "zip_lookup", "y_mean", "y_sd"
  )
  missing_objects <- setdiff(required_objects, names(analysis_objects))
  if (length(missing_objects) > 0L) {
    stop(
      "Analysis objects are missing: ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  train <- analysis_objects$train
  required_train_columns <- c(
    "state_id", "zip_id_num", "status", "log_house_size", "log_acre_lot",
    "bed", "bath"
  )
  missing_train_columns <- setdiff(required_train_columns, names(train))
  if (length(missing_train_columns) > 0L) {
    stop(
      "Training data are missing: ",
      paste(missing_train_columns, collapse = ", "),
      call. = FALSE
    )
  }

  X_train <- as.matrix(analysis_objects$X_train)
  if (nrow(X_train) != nrow(train) ||
      ncol(X_train) != length(analysis_objects$predictor_names) ||
      !identical(colnames(X_train), analysis_objects$predictor_names)) {
    stop(
      "Saved training design matrix does not match the predictor metadata.",
      call. = FALSE
    )
  }
  if (any(!is.finite(X_train))) {
    stop("Saved training design matrix contains non-finite values.", call. = FALSE)
  }

  zip_lookup <- as.data.frame(analysis_objects$zip_lookup, stringsAsFactors = FALSE)
  required_lookup_columns <- c("zip_id", "state", "zip_id_num", "state_id")
  missing_lookup_columns <- setdiff(required_lookup_columns, names(zip_lookup))
  if (length(missing_lookup_columns) > 0L ||
      anyDuplicated(zip_lookup$zip_id_num) ||
      !setequal(zip_lookup$zip_id_num, seq_along(analysis_objects$zip_levels))) {
    stop("ZIP lookup is not a one-to-one JAGS ZIP index lookup.", call. = FALSE)
  }

  if (!is.matrix(draw_matrix) || nrow(draw_matrix) == 0L || ncol(draw_matrix) == 0L ||
      any(!is.finite(draw_matrix))) {
    stop("M2 parameter draws are empty or contain non-finite values.", call. = FALSE)
  }

  required_parameters <- c(
    "alpha", "sd_zip",
    required_indexed_parameter_names("beta", ncol(X_train)),
    required_indexed_parameter_names("u", length(analysis_objects$state_levels)),
    required_indexed_parameter_names("v", length(analysis_objects$zip_levels))
  )
  if (is_heteroskedastic_m2_draws(draw_matrix)) {
    required_parameters <- c(
      required_parameters,
      "gamma0", "sd_sigma_state",
      required_indexed_parameter_names("gamma", ncol(X_train)),
      required_indexed_parameter_names("h", length(analysis_objects$state_levels))
    )
  } else {
    required_parameters <- c(required_parameters, "sigma")
  }
  missing_parameters <- setdiff(required_parameters, colnames(draw_matrix))
  if (length(missing_parameters) > 0L) {
    stop(
      "M2 draws are missing required prediction parameter(s): ",
      paste(missing_parameters, collapse = ", "),
      call. = FALSE
    )
  }
  if (any(draw_matrix[, "sd_zip"] <= 0) ||
      (!is_heteroskedastic_m2_draws(draw_matrix) && any(draw_matrix[, "sigma"] <= 0)) ||
      (is_heteroskedastic_m2_draws(draw_matrix) && any(draw_matrix[, "sd_sigma_state"] <= 0))) {
    stop("M2 draws must have positive scale parameters.", call. = FALSE)
  }

  if (!is.numeric(analysis_objects$y_mean) || length(analysis_objects$y_mean) != 1L ||
      !is.finite(analysis_objects$y_mean) ||
      !is.numeric(analysis_objects$y_sd) || length(analysis_objects$y_sd) != 1L ||
      !is.finite(analysis_objects$y_sd) || analysis_objects$y_sd <= 0) {
    stop("Saved outcome scaling values are invalid.", call. = FALSE)
  }

  invisible(TRUE)
}

# Draw standardized-log-price predictive outcomes for one posterior draw. When
# `include_observed_zip` is FALSE, zip_id is intentionally ignored: a fresh
# ZIP random effect is drawn from Normal(0, sd_zip^2) for every listing.
posterior_predict_m2 <- function(draw, X, state_id, zip_id,
                                 include_observed_zip = TRUE) {
  X <- as.matrix(X)
  state_id <- as.integer(state_id)
  zip_id <- as.integer(zip_id)
  n_listings <- nrow(X)

  if (is.null(names(draw))) {
    stop("draw must be a named vector of M2 parameters.", call. = FALSE)
  }
  if (n_listings == 0L || length(state_id) != n_listings ||
      length(zip_id) != n_listings || anyNA(state_id) || any(state_id < 1L) ||
      (include_observed_zip && (anyNA(zip_id) || any(zip_id < 1L)))) {
    stop("X, state_id, and zip_id must describe the same non-empty listings.", call. = FALSE)
  }

  beta_names <- required_indexed_parameter_names("beta", ncol(X))
  state_names <- required_indexed_parameter_names("u", max(state_id))
  heteroskedastic <- "gamma0" %in% names(draw)
  required_names <- c("alpha", "sd_zip", beta_names, state_names)
  if (heteroskedastic) {
    required_names <- c(
      required_names,
      "gamma0",
      required_indexed_parameter_names("gamma", ncol(X)),
      required_indexed_parameter_names("h", max(state_id))
    )
  } else {
    required_names <- c(required_names, "sigma")
  }
  if (include_observed_zip) {
    required_names <- c(
      required_names,
      required_indexed_parameter_names("v", max(zip_id))
    )
  }
  missing_names <- setdiff(required_names, names(draw))
  if (length(missing_names) > 0L) {
    stop(
      "draw is missing required parameter(s): ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }
  if (any(!is.finite(X)) || any(!is.finite(draw[required_names])) ||
      draw["sd_zip"] <= 0 || (!heteroskedastic && draw["sigma"] <= 0)) {
    stop("Prediction inputs contain invalid covariates or scale parameters.", call. = FALSE)
  }

  beta <- unname(draw[beta_names])
  state_component <- unname(draw[state_names])[state_id]
  if (include_observed_zip) {
    zip_names <- required_indexed_parameter_names("v", max(zip_id))
    zip_component <- unname(draw[zip_names])[zip_id]
  } else {
    zip_component <- stats::rnorm(n_listings, mean = 0, sd = draw["sd_zip"])
  }

  mu <- as.vector(
    draw["alpha"] + X %*% beta + state_component + zip_component
  )
  if (heteroskedastic) {
    gamma_names <- required_indexed_parameter_names("gamma", ncol(X))
    state_sigma_names <- required_indexed_parameter_names("h", max(state_id))
    sigma <- exp(as.vector(
      draw["gamma0"] +
        X %*% unname(draw[gamma_names]) +
        unname(draw[state_sigma_names])[state_id]
    ))
  } else {
    sigma <- rep(draw["sigma"], n_listings)
  }
  stats::rnorm(n_listings, mean = mu, sd = sigma)
}

unscale_price <- function(y_standardised, y_mean, y_sd) {
  exp(y_standardised * y_sd + y_mean)
}

make_typical_property_profile <- function(analysis_objects) {
  train <- analysis_objects$train
  X_train <- as.matrix(analysis_objects$X_train)
  status_values <- as.character(train$status)
  status_counts <- sort(table(status_values), decreasing = TRUE)
  if (length(status_counts) == 0L) {
    stop("Cannot construct a profile without an observed listing status.", call. = FALSE)
  }
  common_status <- sort(names(status_counts)[status_counts == max(status_counts)])[1L]
  template_row <- which(status_values == common_status)[1L]

  profile_x <- X_train[template_row, ]
  continuous_predictors <- intersect(
    c("log_house_size_z", "log_acre_lot_z", "bed_z", "bath_z"),
    colnames(X_train)
  )
  profile_x[continuous_predictors] <- apply(
    X_train[, continuous_predictors, drop = FALSE],
    2L,
    stats::median
  )

  list(
    X = matrix(
      as.numeric(profile_x),
      nrow = 1L,
      dimnames = list(NULL, colnames(X_train))
    ),
    status = common_status,
    house_size_sqft = exp(stats::median(train$log_house_size)),
    acre_lot_acres = expm1(stats::median(train$log_acre_lot)),
    bedrooms = stats::median(train$bed),
    bathrooms = stats::median(train$bath)
  )
}

summarise_zip_location_effects <- function(draw_matrix, analysis_objects) {
  zip_lookup <- as.data.frame(analysis_objects$zip_lookup, stringsAsFactors = FALSE)
  zip_lookup <- zip_lookup[match(
    seq_along(analysis_objects$zip_levels),
    zip_lookup$zip_id_num
  ), ]

  location_effect_median <- vapply(seq_len(nrow(zip_lookup)), function(zip_index) {
    state_index <- zip_lookup$state_id[zip_index]
    stats::median(
      draw_matrix[, paste0("u[", state_index, "]")] +
        draw_matrix[, paste0("v[", zip_index, "]")]
    )
  }, numeric(1))

  zip_lookup$location_effect_median <- location_effect_median
  zip_lookup$zip_code <- sub("^.*__", "", zip_lookup$zip_id)
  zip_lookup
}

make_prediction_scenarios <- function(draw_matrix, analysis_objects) {
  profile <- make_typical_property_profile(analysis_objects)
  location_effects <- summarise_zip_location_effects(draw_matrix, analysis_objects)

  low_index <- order(
    location_effects$location_effect_median,
    location_effects$state,
    location_effects$zip_id
  )[1L]
  high_index <- order(
    -location_effects$location_effect_median,
    location_effects$state,
    location_effects$zip_id
  )[1L]
  low_location <- location_effects[low_index, ]
  high_location <- location_effects[high_index, ]

  scenario_metadata <- data.frame(
    scenario_id = c(
      "low_adjusted_price_observed_zip",
      "high_adjusted_price_observed_zip",
      "new_zip_in_high_price_state"
    ),
    scenario_label = c(
      "Representative home: low adjusted-price observed ZIP",
      "Representative home: high adjusted-price observed ZIP",
      "Representative home: new ZIP in high adjusted-price state"
    ),
    prediction_target = c(
      "New listing in an observed ZIP",
      "New listing in an observed ZIP",
      "New listing in a new ZIP within an observed state"
    ),
    state = c(low_location$state, high_location$state, high_location$state),
    zip_code = c(low_location$zip_code, high_location$zip_code, NA_character_),
    state_id = c(low_location$state_id, high_location$state_id, high_location$state_id),
    zip_id_num = c(low_location$zip_id_num, high_location$zip_id_num, NA_integer_),
    include_observed_zip = c(TRUE, TRUE, FALSE),
    location_effect_median = c(
      low_location$location_effect_median,
      high_location$location_effect_median,
      NA_real_
    ),
    status = rep(profile$status, 3L),
    house_size_sqft = rep(profile$house_size_sqft, 3L),
    acre_lot_acres = rep(profile$acre_lot_acres, 3L),
    bedrooms = rep(profile$bedrooms, 3L),
    bathrooms = rep(profile$bathrooms, 3L),
    stringsAsFactors = FALSE
  )

  list(
    X = profile$X[rep(1L, nrow(scenario_metadata)), , drop = FALSE],
    state_id = scenario_metadata$state_id,
    # The new-ZIP row does not use this placeholder because its ZIP effect is
    # newly drawn. Keeping a valid index makes the input shape explicit.
    zip_id = ifelse(
      is.na(scenario_metadata$zip_id_num),
      high_location$zip_id_num,
      scenario_metadata$zip_id_num
    ),
    include_observed_zip = scenario_metadata$include_observed_zip,
    metadata = scenario_metadata
  )
}

simulate_prediction_draws <- function(draw_matrix, scenarios, draw_ids) {
  n_scenarios <- nrow(scenarios$X)
  prediction_draws <- matrix(
    NA_real_,
    nrow = length(draw_ids),
    ncol = n_scenarios,
    dimnames = list(NULL, scenarios$metadata$scenario_id)
  )

  for (scenario_index in seq_len(n_scenarios)) {
    prediction_draws[, scenario_index] <- vapply(draw_ids, function(draw_id) {
      # `draws_matrix` preserves its class when subset to one iteration, so
      # explicitly restore a named numeric vector for the prediction helper.
      draw <- as.numeric(draw_matrix[draw_id, , drop = FALSE])
      names(draw) <- colnames(draw_matrix)
      posterior_predict_m2(
        draw = draw,
        X = scenarios$X[scenario_index, , drop = FALSE],
        state_id = scenarios$state_id[scenario_index],
        zip_id = scenarios$zip_id[scenario_index],
        include_observed_zip = scenarios$include_observed_zip[scenario_index]
      )
    }, numeric(1))
  }
  prediction_draws
}

summarise_predictive_draws <- function(prediction_draws, scenario_metadata,
                                       y_mean, y_sd, total_draws) {
  dollar_draws <- unscale_price(prediction_draws, y_mean = y_mean, y_sd = y_sd)
  summaries <- data.frame(
    scenario_metadata,
    n_parameter_draws = rep(total_draws, ncol(prediction_draws)),
    n_predictive_draws = rep(nrow(prediction_draws), ncol(prediction_draws)),
    median_standardised_log_price = apply(prediction_draws, 2L, stats::median),
    lower_50_standardised_log_price = apply(
      prediction_draws, 2L, stats::quantile, probs = 0.25, names = FALSE
    ),
    upper_50_standardised_log_price = apply(
      prediction_draws, 2L, stats::quantile, probs = 0.75, names = FALSE
    ),
    lower_95_standardised_log_price = apply(
      prediction_draws, 2L, stats::quantile, probs = 0.025, names = FALSE
    ),
    upper_95_standardised_log_price = apply(
      prediction_draws, 2L, stats::quantile, probs = 0.975, names = FALSE
    ),
    median_price_usd = apply(dollar_draws, 2L, stats::median),
    lower_50_price_usd = apply(
      dollar_draws, 2L, stats::quantile, probs = 0.25, names = FALSE
    ),
    upper_50_price_usd = apply(
      dollar_draws, 2L, stats::quantile, probs = 0.75, names = FALSE
    ),
    lower_95_price_usd = apply(
      dollar_draws, 2L, stats::quantile, probs = 0.025, names = FALSE
    ),
    upper_95_price_usd = apply(
      dollar_draws, 2L, stats::quantile, probs = 0.975, names = FALSE
    ),
    row.names = NULL,
    check.names = FALSE
  )
  summaries
}

save_prediction_interval_plot <- function(prediction_summary, figure_path) {
  plot_data <- prediction_summary[order(
    prediction_summary$median_price_usd,
    prediction_summary$scenario_label
  ), ]
  plot_data$scenario_label <- factor(
    plot_data$scenario_label,
    levels = plot_data$scenario_label
  )

  interval_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(y = scenario_label, x = median_price_usd)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower_95_price_usd, xmax = upper_95_price_usd),
      orientation = "y",
      width = 0,
      colour = "#9ECAE1",
      linewidth = 1.2
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower_50_price_usd, xmax = upper_50_price_usd),
      orientation = "y",
      width = 0,
      colour = "#2171B5",
      linewidth = 2.4
    ) +
    ggplot2::geom_point(colour = "#08306B", size = 2.5) +
    ggplot2::labs(
      title = "M2 heteroskedastic posterior predictive price intervals",
      subtitle = "Points are predictive medians; dark and light intervals are 50% and 95% predictive intervals.",
      x = "Predicted listing price (USD)",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11)

  ggplot2::ggsave(
    filename = figure_path,
    plot = interval_plot,
    width = 10,
    height = 4.5,
    dpi = 180
  )
}

run_posterior_predictions <- function(
    draws_path = here::here("output", "mcmc", "m2_heteroskedastic_xall_state_draws.rds"),
    analysis_objects_path = here::here("data", "processed", "analysis_objects.rds"),
    tables_dir = here::here("output", "tables"),
    figures_dir = here::here("output", "figures"),
    max_prediction_draws = 5000L,
    seed = 20260621L) {
  max_prediction_draws <- validate_positive_integer(
    max_prediction_draws,
    "max_prediction_draws"
  )
  seed <- validate_positive_integer(seed, "seed")

  if (!file.exists(draws_path)) {
    stop(
      "M2 production draws not found: ", draws_path,
      ". Run scripts/05_fit_models.R after convergence checks first.",
      call. = FALSE
    )
  }
  if (!file.exists(analysis_objects_path)) {
    stop(
      "Prepared analysis objects not found: ", analysis_objects_path,
      ". Run scripts/03_prepare_jags_data.R first.",
      call. = FALSE
    )
  }

  analysis_objects <- readRDS(analysis_objects_path)
  draw_matrix <- posterior::as_draws_matrix(readRDS(draws_path))
  validate_prediction_inputs(draw_matrix, analysis_objects)

  set.seed(seed)
  draw_ids <- sample(
    seq_len(nrow(draw_matrix)),
    size = min(max_prediction_draws, nrow(draw_matrix)),
    replace = FALSE
  )
  # Use all converged draws to identify the observed high/low adjusted-price
  # locations, then cap only the number of predictive outcome simulations.
  scenarios <- make_prediction_scenarios(draw_matrix, analysis_objects)
  prediction_draws <- simulate_prediction_draws(
    draw_matrix = draw_matrix,
    scenarios = scenarios,
    draw_ids = draw_ids
  )
  prediction_summary <- summarise_predictive_draws(
    prediction_draws = prediction_draws,
    scenario_metadata = scenarios$metadata,
    y_mean = analysis_objects$y_mean,
    y_sd = analysis_objects$y_sd,
    total_draws = nrow(draw_matrix)
  )

  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  table_path <- file.path(tables_dir, "m2_posterior_predictions.csv")
  figure_path <- file.path(figures_dir, "m2_posterior_prediction_intervals.png")
  readr::write_csv(prediction_summary, table_path)
  save_prediction_interval_plot(prediction_summary, figure_path)

  message(
    "Saved ", nrow(prediction_summary), " M2 posterior prediction scenarios using ",
    nrow(prediction_draws), " of ", nrow(draw_matrix), " parameter draws. Table: ",
    table_path, "; figure: ", figure_path
  )

  invisible(list(
    summary = prediction_summary,
    prediction_draws = prediction_draws,
    draw_ids = draw_ids,
    table_path = table_path,
    figure_path = figure_path
  ))
}

if (sys.nframe() == 0L) {
  run_posterior_predictions()
}
