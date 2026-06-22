# Summarise posterior inference from the production M2 fit.
#
# This script implements point 14 of IMPLEMENTATION_GUIDE.md. It creates the
# fixed-effect and hierarchical-variance summaries, complete labelled random-
# effect tables, and the two report-ready caterpillar plots. It does not make
# predictions, evaluate held-out performance, or run posterior predictive
# checks; those have separate implementation steps.

required_posterior_summary_packages <- c("ggplot2", "here", "posterior", "readr")
missing_posterior_summary_packages <- required_posterior_summary_packages[!vapply(
  required_posterior_summary_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_posterior_summary_packages) > 0L) {
  stop(
    "Install the package(s) required for posterior summaries: ",
    paste(missing_posterior_summary_packages, collapse = ", "),
    call. = FALSE
  )
}

draws_path <- here::here("output", "mcmc", "m2_draws.rds")
analysis_objects_path <- here::here("data", "processed", "analysis_objects.rds")
if (!file.exists(draws_path)) {
  stop(
    "M2 production draws not found: ", draws_path,
    ". Run scripts/05_fit_models.R first.",
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

tables_dir <- here::here("output", "tables")
figures_dir <- here::here("output", "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

analysis_objects <- readRDS(analysis_objects_path)
draw_matrix_m2 <- posterior::as_draws_matrix(readRDS(draws_path))
if (nrow(draw_matrix_m2) == 0L || ncol(draw_matrix_m2) == 0L) {
  stop("M2 production draws are empty.", call. = FALSE)
}
if (any(!is.finite(draw_matrix_m2))) {
  stop("M2 production draws contain non-finite values.", call. = FALSE)
}

required_analysis_objects <- c(
  "predictor_names", "state_levels", "zip_levels", "zip_lookup", "train"
)
missing_analysis_objects <- setdiff(required_analysis_objects, names(analysis_objects))
if (length(missing_analysis_objects) > 0L) {
  stop(
    "Analysis objects are missing: ",
    paste(missing_analysis_objects, collapse = ", "),
    call. = FALSE
  )
}

available_parameters <- colnames(draw_matrix_m2)
required_scalar_parameters <- c("sigma", "sd_state", "sd_zip")
missing_scalar_parameters <- setdiff(required_scalar_parameters, available_parameters)
if (length(missing_scalar_parameters) > 0L) {
  stop(
    "M2 draws are missing scalar parameter(s): ",
    paste(missing_scalar_parameters, collapse = ", "),
    call. = FALSE
  )
}

indexed_parameters <- function(prefix) {
  pattern <- paste0("^", prefix, "\\[([0-9]+)\\]$")
  parameters <- grep(pattern, available_parameters, value = TRUE)
  indices <- as.integer(sub(pattern, "\\1", parameters))
  parameters[order(indices)]
}

parameter_indices <- function(parameters) {
  as.integer(sub("^[^\\[]+\\[([0-9]+)\\]$", "\\1", parameters))
}

beta_parameters <- indexed_parameters("beta")
state_parameters <- indexed_parameters("u")
zip_parameters <- indexed_parameters("v")

if (length(beta_parameters) != length(analysis_objects$predictor_names)) {
  stop(
    "The number of beta draws does not match the saved predictor names.",
    call. = FALSE
  )
}
if (length(state_parameters) != length(analysis_objects$state_levels)) {
  stop(
    "The number of state-effect draws does not match the state lookup.",
    call. = FALSE
  )
}
if (length(zip_parameters) != length(analysis_objects$zip_levels)) {
  stop(
    "The number of ZIP-effect draws does not match the ZIP lookup.",
    call. = FALSE
  )
}

summarise_columns <- function(parameters) {
  values <- draw_matrix_m2[, parameters, drop = FALSE]
  data.frame(
    jags_parameter = parameters,
    posterior_median = apply(values, 2L, stats::median),
    lower_95 = apply(values, 2L, stats::quantile, probs = 0.025, names = FALSE),
    upper_95 = apply(values, 2L, stats::quantile, probs = 0.975, names = FALSE),
    prob_positive = colMeans(values > 0),
    row.names = NULL,
    check.names = FALSE
  )
}

fixed_effects <- summarise_columns(beta_parameters)
fixed_effects$parameter_index <- parameter_indices(fixed_effects$jags_parameter)
fixed_effects$predictor <- analysis_objects$predictor_names[fixed_effects$parameter_index]
fixed_effects <- fixed_effects[, c(
  "parameter_index", "predictor", "jags_parameter", "posterior_median",
  "lower_95", "upper_95", "prob_positive"
)]
readr::write_csv(fixed_effects, file.path(tables_dir, "m2_fixed_effects.csv"))

sd_difference <- as.numeric(
  draw_matrix_m2[, "sd_zip"] - draw_matrix_m2[, "sd_state"]
)
hierarchical_summary <- rbind(
  data.frame(
    quantity = c("Residual SD", "State SD", "ZIP SD", "ZIP SD - State SD"),
    posterior_median = c(
      stats::median(draw_matrix_m2[, "sigma"]),
      stats::median(draw_matrix_m2[, "sd_state"]),
      stats::median(draw_matrix_m2[, "sd_zip"]),
      stats::median(sd_difference)
    ),
    lower_95 = c(
      stats::quantile(draw_matrix_m2[, "sigma"], 0.025, names = FALSE),
      stats::quantile(draw_matrix_m2[, "sd_state"], 0.025, names = FALSE),
      stats::quantile(draw_matrix_m2[, "sd_zip"], 0.025, names = FALSE),
      stats::quantile(sd_difference, 0.025, names = FALSE)
    ),
    upper_95 = c(
      stats::quantile(draw_matrix_m2[, "sigma"], 0.975, names = FALSE),
      stats::quantile(draw_matrix_m2[, "sd_state"], 0.975, names = FALSE),
      stats::quantile(draw_matrix_m2[, "sd_zip"], 0.975, names = FALSE),
      stats::quantile(sd_difference, 0.975, names = FALSE)
    ),
    posterior_probability = c(NA_real_, NA_real_, NA_real_, mean(sd_difference > 0)),
    row.names = NULL
  ),
  data.frame(
    quantity = "P(ZIP SD > State SD)",
    posterior_median = NA_real_,
    lower_95 = NA_real_,
    upper_95 = NA_real_,
    posterior_probability = mean(sd_difference > 0),
    row.names = NULL
  )
)
readr::write_csv(
  hierarchical_summary,
  file.path(tables_dir, "m2_hierarchical_variance_summary.csv")
)

zip_minus_state_plot <- ggplot2::ggplot(
  data.frame(zip_minus_state_sd = sd_difference),
  ggplot2::aes(x = zip_minus_state_sd)
) +
  ggplot2::geom_density(fill = "#2C7FB8", alpha = 0.65, colour = "#1D4E89") +
  ggplot2::geom_vline(linetype = "dashed", colour = "#B2182B", xintercept = 0) +
  ggplot2::labs(
    title = "Posterior distribution of ZIP SD minus state SD",
    x = "ZIP SD - state SD (standardized log-price units)",
    y = "Posterior density"
  ) +
  ggplot2::theme_minimal(base_size = 12)
ggplot2::ggsave(
  filename = file.path(figures_dir, "m2_zip_minus_state_sd_posterior.png"),
  plot = zip_minus_state_plot,
  width = 8,
  height = 5,
  dpi = 180
)

state_effects <- summarise_columns(state_parameters)
state_effects$state_index <- parameter_indices(state_effects$jags_parameter)
state_effects$state <- analysis_objects$state_levels[state_effects$state_index]

state_counts <- tabulate(
  analysis_objects$train$state_id,
  nbins = length(analysis_objects$state_levels)
)
state_effects$training_listings <- state_counts[state_effects$state_index]
state_effects <- state_effects[, c(
  "state_index", "state", "jags_parameter", "training_listings",
  "posterior_median", "lower_95", "upper_95", "prob_positive"
)]
readr::write_csv(state_effects, file.path(tables_dir, "m2_state_effects.csv"))

zip_lookup <- as.data.frame(analysis_objects$zip_lookup, stringsAsFactors = FALSE)
required_zip_lookup_columns <- c("zip_id_num", "zip_id", "state")
missing_zip_lookup_columns <- setdiff(required_zip_lookup_columns, names(zip_lookup))
if (length(missing_zip_lookup_columns) > 0L ||
    anyDuplicated(zip_lookup$zip_id_num) ||
    !setequal(zip_lookup$zip_id_num, seq_along(analysis_objects$zip_levels))) {
  stop("The ZIP lookup is not a one-to-one JAGS ZIP index lookup.", call. = FALSE)
}
zip_lookup <- zip_lookup[match(seq_along(analysis_objects$zip_levels), zip_lookup$zip_id_num), ]
zip_code <- sub("^.*__", "", zip_lookup$zip_id)

zip_effects <- summarise_columns(zip_parameters)
zip_effects$zip_index <- parameter_indices(zip_effects$jags_parameter)
zip_effects$state <- zip_lookup$state[zip_effects$zip_index]
zip_effects$zip_id <- zip_lookup$zip_id[zip_effects$zip_index]
zip_effects$zip_code <- zip_code[zip_effects$zip_index]
zip_effects$label <- paste(zip_effects$state, zip_effects$zip_code, sep = " — ")

zip_counts <- tabulate(
  analysis_objects$train$zip_id_num,
  nbins = length(analysis_objects$zip_levels)
)
zip_effects$training_listings <- zip_counts[zip_effects$zip_index]
zip_effects <- zip_effects[, c(
  "zip_index", "state", "zip_code", "zip_id", "label", "jags_parameter",
  "training_listings", "posterior_median", "lower_95", "upper_95",
  "prob_positive"
)]
readr::write_csv(zip_effects, file.path(tables_dir, "m2_zip_effects.csv"))

# The main-report ZIP selection is deliberately made from training sample size
# only: take the 25 most-observed ZIPs, breaking ties by state and ZIP label.
# It therefore cannot be influenced by the fitted random-effect estimates.
selected_zip_count <- 25L
zip_selection_order <- order(
  -zip_effects$training_listings,
  zip_effects$state,
  zip_effects$zip_code
)
selected_zip_effects <- zip_effects[head(
  zip_selection_order,
  min(selected_zip_count, nrow(zip_effects))
), ]
selected_zip_effects$selection_rule <- paste0(
  "Top ", selected_zip_count,
  " ZIPs by training listings; ties by state then ZIP code."
)
readr::write_csv(
  selected_zip_effects,
  file.path(tables_dir, "m2_selected_zip_effects.csv")
)

make_caterpillar_plot <- function(effects, label_column, title) {
  plot_data <- effects[order(effects$posterior_median, effects[[label_column]]), ]
  plot_data$plot_label <- factor(
    plot_data[[label_column]],
    levels = plot_data[[label_column]]
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = plot_label, y = posterior_median)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "#666666") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower_95, ymax = upper_95),
      width = 0
    ) +
    ggplot2::geom_point(colour = "#2C7FB8", size = 2) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = "Conditional geographic effect (standardized log-price units)"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

state_effect_plot <- make_caterpillar_plot(
  state_effects,
  label_column = "state",
  title = "M2 state random effects with 95% credible intervals"
)
ggplot2::ggsave(
  filename = file.path(figures_dir, "m2_state_effects_caterpillar.png"),
  plot = state_effect_plot,
  width = 8,
  height = max(5, 0.28 * nrow(state_effects) + 1),
  dpi = 180
)

selected_zip_effect_plot <- make_caterpillar_plot(
  selected_zip_effects,
  label_column = "label",
  title = "M2 selected ZIP random effects with 95% credible intervals"
)
ggplot2::ggsave(
  filename = file.path(figures_dir, "m2_selected_zip_effects_caterpillar.png"),
  plot = selected_zip_effect_plot,
  width = 8,
  height = max(5, 0.28 * nrow(selected_zip_effects) + 1),
  dpi = 180
)

message(
  "Saved M2 posterior summaries: ", nrow(fixed_effects), " fixed effects, ",
  nrow(state_effects), " state effects, and ", nrow(zip_effects),
  " ZIP effects. Full ZIP results are in m2_zip_effects.csv; ",
  nrow(selected_zip_effects), " ZIPs use the pre-specified report rule."
)
