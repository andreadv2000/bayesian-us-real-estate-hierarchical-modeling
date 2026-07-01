# Produce report-ready figures and tables (IMPLEMENTATION_GUIDE.md point 19).
#
# This final artifact stage reads the saved products of the audit, sampling,
# simulation, fitting, diagnostics, comparison, posterior-predictive, and
# prediction workflows. It deliberately does not refit a model or recalculate
# a diagnostic. Alongside the existing M2 posterior summaries, it creates the
# remaining report figures, standardises the essential tables, and records a
# manifest with titles, units, interval definitions, and provenance.

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

make_m2_posterior_summary_assets <- function() {

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

invisible(list(
  fixed_effects = fixed_effects,
  hierarchical_summary = hierarchical_summary,
  state_effects = state_effects,
  zip_effects = zip_effects,
  selected_zip_effects = selected_zip_effects
))
}

read_report_csv <- function(path, required_columns = character()) {
  if (!file.exists(path)) {
    stop(
      "Required saved artifact not found: ", path,
      ". Run the upstream workflow that creates it before point 19.",
      call. = FALSE
    )
  }

  data <- as.data.frame(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Saved artifact has an unexpected schema: ", path,
      ". Missing column(s): ", paste(missing_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }
  data
}

write_report_csv <- function(data, tables_dir, file_name) {
  path <- file.path(tables_dir, file_name)
  readr::write_csv(data, path)
  path
}

save_fixed_effect_forest_plot <- function(fixed_effects, figure_path) {
  required_columns <- c("predictor", "posterior_median", "lower_95", "upper_95")
  missing_columns <- setdiff(required_columns, names(fixed_effects))
  if (length(missing_columns) > 0L) {
    stop(
      "M2 fixed-effect summary is missing: ",
      paste(missing_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }

  plot_data <- fixed_effects[order(
    fixed_effects$posterior_median,
    fixed_effects$predictor
  ), , drop = FALSE]
  plot_data$predictor <- factor(plot_data$predictor, levels = plot_data$predictor)

  forest_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = predictor, y = posterior_median)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "#666666") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower_95, ymax = upper_95),
      width = 0,
      colour = "#2171B5"
    ) +
    ggplot2::geom_point(size = 2.4, colour = "#08306B") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "M2 fixed-effect posterior medians and 95% credible intervals",
      subtitle = "Effects are conditional on the state and ZIP random intercepts.",
      x = NULL,
      y = "Coefficient (standardized log-price units)"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  ggplot2::ggsave(
    filename = figure_path,
    plot = forest_plot,
    width = 9,
    height = max(4.5, 0.52 * nrow(plot_data) + 1.2),
    dpi = 180
  )
}

save_variance_component_plot <- function(draw_matrix, figure_path) {
  parameters <- c("sigma", "sd_state", "sd_zip")
  missing_parameters <- setdiff(parameters, colnames(draw_matrix))
  if (length(missing_parameters) > 0L) {
    stop(
      "M2 draws are missing variance component(s): ",
      paste(missing_parameters, collapse = ", "), ".",
      call. = FALSE
    )
  }

  component_labels <- c(
    sigma = "Residual SD",
    sd_state = "State SD",
    sd_zip = "ZIP SD"
  )
  plot_data <- do.call(rbind, lapply(parameters, function(parameter) {
    data.frame(
      component = component_labels[[parameter]],
      value = as.numeric(draw_matrix[, parameter]),
      stringsAsFactors = FALSE
    )
  }))
  plot_data$component <- factor(
    plot_data$component,
    levels = unname(component_labels[parameters])
  )

  density_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value)) +
    ggplot2::geom_density(fill = "#6BAED6", colour = "#08519C", alpha = 0.7) +
    ggplot2::facet_wrap(~component, scales = "free_y", ncol = 1L) +
    ggplot2::labs(
      title = "M2 posterior distributions of residual and geographic standard deviations",
      x = "Standard deviation (standardized log-price units)",
      y = "Posterior density"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  ggplot2::ggsave(
    filename = figure_path,
    plot = density_plot,
    width = 9,
    height = 8,
    dpi = 180
  )
}

save_sampling_coverage_plots <- function(raw_state_counts, raw_zip_counts,
                                         sample_manifest, figures_dir) {
  required_state_columns <- c("state", "listings")
  required_zip_columns <- c("state", "zip_code", "listings")
  required_manifest_columns <- c("state", "zip_id", "split", "n")
  if (length(setdiff(required_state_columns, names(raw_state_counts))) > 0L ||
      length(setdiff(required_zip_columns, names(raw_zip_counts))) > 0L ||
      length(setdiff(required_manifest_columns, names(sample_manifest))) > 0L) {
    stop("Saved sampling inputs do not contain the expected state, ZIP, and count columns.", call. = FALSE)
  }

  normalise_zip_label <- function(values) {
    labels <- trimws(as.character(values))
    numeric_values <- suppressWarnings(as.numeric(labels))
    is_numeric_zip <- is.finite(numeric_values) & numeric_values >= 0 &
      numeric_values == floor(numeric_values)
    labels[is_numeric_zip] <- sprintf("%05d", as.integer(numeric_values[is_numeric_zip]))
    labels
  }
  raw_state_counts$state <- as.character(raw_state_counts$state)
  raw_zip_counts$state <- as.character(raw_zip_counts$state)
  raw_zip_counts$zip_code <- normalise_zip_label(raw_zip_counts$zip_code)
  sample_manifest$state <- as.character(sample_manifest$state)
  sample_manifest$zip_id <- as.character(sample_manifest$zip_id)

  selected_state_counts <- stats::aggregate(
    n ~ state,
    data = sample_manifest,
    FUN = sum
  )
  names(selected_state_counts)[2L] <- "sampled_listings"
  state_comparison <- merge(
    raw_state_counts,
    selected_state_counts,
    by = "state",
    all.y = TRUE,
    sort = FALSE
  )
  if (anyNA(state_comparison$listings) || any(state_comparison$listings <= 0) ||
      any(state_comparison$sampled_listings <= 0)) {
    stop("Each sampled state must have positive raw and sampled listing counts.", call. = FALSE)
  }
  state_comparison <- state_comparison[order(state_comparison$listings, state_comparison$state), ]
  state_long <- rbind(
    data.frame(
      state = state_comparison$state,
      stage = "Before sampling",
      listings = state_comparison$listings,
      stringsAsFactors = FALSE
    ),
    data.frame(
      state = state_comparison$state,
      stage = "Balanced analysis sample",
      listings = state_comparison$sampled_listings,
      stringsAsFactors = FALSE
    )
  )
  state_long$state <- factor(state_long$state, levels = state_comparison$state)
  state_plot <- ggplot2::ggplot(
    state_long,
    ggplot2::aes(x = state, y = listings, fill = stage)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_fill_manual(values = c(
      "Before sampling" = "#9ECAE1",
      "Balanced analysis sample" = "#2171B5"
    )) +
    ggplot2::labs(
      title = "Listings per selected state before and after sampling",
      subtitle = "The logarithmic axis preserves the raw-count comparison with the balanced sample.",
      x = NULL,
      y = "Listings (log10 scale)",
      fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11)
  state_path <- file.path(figures_dir, "sampled_listings_per_state_before_after.png")
  ggplot2::ggsave(
    filename = state_path,
    plot = state_plot,
    width = 10,
    height = max(6, 0.32 * nrow(state_comparison) + 1.4),
    dpi = 180
  )

  sample_manifest$zip_code <- normalise_zip_label(sub("^.*__", "", sample_manifest$zip_id))
  selected_zip_counts <- stats::aggregate(
    n ~ state + zip_code,
    data = sample_manifest,
    FUN = sum
  )
  names(selected_zip_counts)[3L] <- "sampled_listings"
  zip_comparison <- merge(
    raw_zip_counts,
    selected_zip_counts,
    by = c("state", "zip_code"),
    all.y = TRUE,
    sort = FALSE
  )
  if (anyNA(zip_comparison$listings) || any(zip_comparison$listings <= 0) ||
      any(zip_comparison$sampled_listings <= 0)) {
    stop("Each sampled state--ZIP group must have positive raw and sampled listing counts.", call. = FALSE)
  }
  zip_long <- rbind(
    data.frame(
      stage = "Before sampling",
      listings = zip_comparison$listings,
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "Balanced analysis sample",
      listings = zip_comparison$sampled_listings,
      stringsAsFactors = FALSE
    )
  )
  zip_long$stage <- factor(
    zip_long$stage,
    levels = c("Before sampling", "Balanced analysis sample")
  )
  zip_plot <- ggplot2::ggplot(zip_long, ggplot2::aes(x = stage, y = listings, fill = stage)) +
    ggplot2::geom_boxplot(width = 0.58, outlier.alpha = 0.25) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_fill_manual(values = c(
      "Before sampling" = "#9ECAE1",
      "Balanced analysis sample" = "#2171B5"
    )) +
    ggplot2::labs(
      title = "Listings per selected state--ZIP group before and after sampling",
      subtitle = "Each retained ZIP contributes exactly 30 listings before the train/test split.",
      x = NULL,
      y = "Listings per state--ZIP group (log10 scale)",
      fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "none")
  zip_path <- file.path(figures_dir, "sampled_listings_per_zip_before_after.png")
  ggplot2::ggsave(
    filename = zip_path,
    plot = zip_plot,
    width = 9,
    height = 5.5,
    dpi = 180
  )

  list(state_path = state_path, zip_path = zip_path)
}

save_simulation_recovery_plot <- function(recovery, figure_path) {
  required_columns <- c("label", "true_value", "posterior_mean", "q2_5", "q97_5")
  missing_columns <- setdiff(required_columns, names(recovery))
  if (length(missing_columns) > 0L) {
    stop(
      "Full-design simulation recovery results are missing: ",
      paste(missing_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (nrow(recovery) == 0L || any(!is.finite(as.matrix(recovery[, required_columns[-1L], drop = FALSE])))) {
    stop("Full-design simulation recovery results must be non-empty and finite.", call. = FALSE)
  }

  plot_data <- recovery[order(recovery$true_value, recovery$label), , drop = FALSE]
  plot_data$label <- factor(plot_data$label, levels = plot_data$label)
  recovery_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = label)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = q2_5, ymax = q97_5),
      width = 0,
      colour = "#2171B5"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = posterior_mean, shape = "Posterior mean"),
      colour = "#08519C",
      size = 2.3
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = true_value, shape = "Known simulated value"),
      colour = "#CB181D",
      size = 2.1
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_shape_manual(values = c("Posterior mean" = 16, "Known simulated value" = 4)) +
    ggplot2::labs(
      title = "M2 parameter recovery on synthetic data",
      subtitle = "Blue intervals are 95% credible intervals; red crosses are known simulation values.",
      x = NULL,
      y = "Standardized log-price units",
      shape = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(
    filename = figure_path,
    plot = recovery_plot,
    width = 9,
    height = max(5, 0.45 * nrow(plot_data) + 1.4),
    dpi = 180
  )
}

make_report_tables <- function(tables_dir, draws_path) {
  dictionary <- read_report_csv(
    file.path(tables_dir, "raw_data_dictionary.csv"),
    c("variable", "type", "units", "missing_n", "missing_prop")
  )
  missingness <- read_report_csv(
    file.path(tables_dir, "raw_missingness.csv"),
    c("variable", "missing_n", "missing_prop")
  )
  dictionary_match <- match(dictionary$variable, missingness$variable)
  dictionary$missing_n <- missingness$missing_n[dictionary_match]
  dictionary$missing_prop <- missingness$missing_prop[dictionary_match]
  dictionary$missing_percent <- 100 * dictionary$missing_prop
  data_dictionary_path <- write_report_csv(
    dictionary,
    tables_dir,
    "report_data_dictionary_missingness.csv"
  )

  row_accounting <- read_report_csv(
    file.path(tables_dir, "row_accounting.csv"),
    c("stage", "n")
  )
  sample_manifest <- read_report_csv(
    here::here("data", "processed", "sample_manifest.csv"),
    c("state", "zip_id", "split", "n")
  )
  sampling_accounting <- data.frame(
    stage = c(
      "Frozen balanced sample",
      "Frozen training sample",
      "Frozen held-out sample",
      "Selected states",
      "Selected state--ZIP groups"
    ),
    n = c(
      sum(sample_manifest$n),
      sum(sample_manifest$n[sample_manifest$split == "train"]),
      sum(sample_manifest$n[sample_manifest$split == "test"]),
      length(unique(sample_manifest$state)),
      length(unique(sample_manifest$zip_id))
    ),
    stringsAsFactors = FALSE
  )
  row_accounting_path <- write_report_csv(
    rbind(row_accounting, sampling_accounting),
    tables_dir,
    "report_cleaning_sampling_row_accounting.csv"
  )

  recovery_metrics <- read_report_csv(
    file.path(tables_dir, "simulation_recovery_repeated_parameter_metrics.csv"),
    c(
      "parameter", "label", "true_value", "n_successful_converged",
      "empirical_coverage_95", "posterior_mean_bias", "posterior_mean_rmse"
    )
  )
  recovery_metrics_path <- write_report_csv(
    recovery_metrics,
    tables_dir,
    "report_simulation_recovery_results.csv"
  )

  diagnostics <- read_report_csv(
    file.path(tables_dir, "m2_posterior_summary.csv"),
    c("variable", "rhat", "ess_bulk", "ess_tail", "mcse_mean", "converged")
  )
  draws_array <- posterior::as_draws_array(readRDS(draws_path))
  mcmc_settings <- data.frame(
    section = "MCMC setting",
    item = c(
      "Chains", "Adaptation iterations per chain", "Warm-up iterations per chain",
      "Sampling iterations per chain", "Thinning", "Retained posterior draws"
    ),
    value = c(
      as.character(posterior::nchains(draws_array)), "2000", "5000", "10000", "1",
      as.character(posterior::ndraws(draws_array))
    ),
    rhat = NA_real_,
    ess_bulk = NA_real_,
    ess_tail = NA_real_,
    mcse_mean = NA_real_,
    converged = NA,
    stringsAsFactors = FALSE
  )
  convergence_rows <- data.frame(
    section = "Convergence diagnostic",
    item = diagnostics$variable,
    value = "Posterior parameter",
    rhat = diagnostics$rhat,
    ess_bulk = diagnostics$ess_bulk,
    ess_tail = diagnostics$ess_tail,
    mcse_mean = diagnostics$mcse_mean,
    converged = diagnostics$converged,
    stringsAsFactors = FALSE
  )
  mcmc_table_path <- write_report_csv(
    rbind(mcmc_settings, convergence_rows),
    tables_dir,
    "report_mcmc_settings_and_convergence.csv"
  )

  fixed_effects <- read_report_csv(
    file.path(tables_dir, "m2_fixed_effects.csv"),
    c("predictor", "jags_parameter", "posterior_median", "lower_95", "upper_95", "prob_positive")
  )
  fixed_effects_report <- data.frame(
    component_type = "Fixed effect",
    component = fixed_effects$predictor,
    parameter = fixed_effects$jags_parameter,
    posterior_median = fixed_effects$posterior_median,
    lower_95 = fixed_effects$lower_95,
    upper_95 = fixed_effects$upper_95,
    posterior_probability = fixed_effects$prob_positive,
    units = "standardized log-price units",
    stringsAsFactors = FALSE
  )
  variance_summary <- read_report_csv(
    file.path(tables_dir, "m2_hierarchical_variance_summary.csv"),
    c("quantity", "posterior_median", "lower_95", "upper_95", "posterior_probability")
  )
  variance_report <- data.frame(
    component_type = "Variance component",
    component = variance_summary$quantity,
    parameter = NA_character_,
    posterior_median = variance_summary$posterior_median,
    lower_95 = variance_summary$lower_95,
    upper_95 = variance_summary$upper_95,
    posterior_probability = variance_summary$posterior_probability,
    units = "standardized log-price units",
    stringsAsFactors = FALSE
  )
  posterior_summary_path <- write_report_csv(
    rbind(fixed_effects_report, variance_report),
    tables_dir,
    "report_m2_fixed_effect_variance_summary.csv"
  )

  dic <- read_report_csv(
    file.path(tables_dir, "model_comparison_dic.csv"),
    c("model", "Dbar", "pD", "DIC", "delta_DIC")
  )
  dic$model_key <- sub(":.*$", "", dic$model)
  dic_decomposition_path <- write_report_csv(
    dic[, c("model", "Dbar", "pD", "DIC", "delta_DIC")],
    tables_dir,
    "report_dic_decomposition.csv"
  )
  heldout <- read_report_csv(
    file.path(tables_dir, "heldout_predictive_metrics.csv"),
    c("model", "n_test", "n_posterior_draws", "RMSE", "MAE", "coverage_95", "total_lppd", "MLPD")
  )
  heldout$model_key <- heldout$model
  heldout_metrics_path <- write_report_csv(
    heldout[, c("model", "n_test", "n_posterior_draws", "RMSE", "MAE", "coverage_95", "total_lppd", "MLPD")],
    tables_dir,
    "report_heldout_predictive_metrics.csv"
  )
  comparison <- merge(
    dic,
    heldout,
    by = "model_key",
    all = TRUE,
    sort = FALSE,
    suffixes = c("_dic", "_heldout")
  )
  comparison$model <- ifelse(
    !is.na(comparison$model_dic), comparison$model_dic, comparison$model_heldout
  )
  comparison <- comparison[, c(
    "model", "Dbar", "pD", "DIC", "delta_DIC", "n_test", "n_posterior_draws",
    "RMSE", "MAE", "coverage_95", "total_lppd", "MLPD"
  )]
  model_comparison_path <- write_report_csv(
    comparison,
    tables_dir,
    "report_dic_heldout_model_comparison.csv"
  )

  ppc_summary <- read_report_csv(
    file.path(tables_dir, "posterior_predictive_check_summary.csv"),
    c(
      "model", "statistic", "scope", "observed", "replicated_mean",
      "replicated_median", "replicated_lower_95", "replicated_upper_95",
      "posterior_predictive_tail_probability"
    )
  )
  ppc_table_path <- write_report_csv(
    ppc_summary,
    tables_dir,
    "report_posterior_predictive_test_statistics.csv"
  )

  list(
    data_dictionary_path = data_dictionary_path,
    row_accounting_path = row_accounting_path,
    recovery_metrics_path = recovery_metrics_path,
    mcmc_table_path = mcmc_table_path,
    posterior_summary_path = posterior_summary_path,
    dic_decomposition_path = dic_decomposition_path,
    heldout_metrics_path = heldout_metrics_path,
    model_comparison_path = model_comparison_path,
    ppc_table_path = ppc_table_path
  )
}

validate_upstream_report_figures <- function(figures_dir) {
  required_files <- c(
    "raw_price_histogram.png",
    "raw_log_price_histogram.png",
    "prior_predictive_price_ranges.png",
    "m2_traceplots.png",
    "m2_rank_histograms.png",
    "m2_acf.png",
    "m0_posterior_predictive_checks.png",
    "m1_posterior_predictive_checks.png",
    "m2_posterior_predictive_checks.png",
    "m2_posterior_prediction_intervals.png"
  )
  paths <- file.path(figures_dir, required_files)
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0L) {
    stop(
      "Point 19 requires saved upstream figure(s) that are missing: ",
      paste(missing_paths, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(paths)
}

make_report_artifact_manifest <- function(tables_dir, figures_dir) {
  artifacts <- data.frame(
    artifact_id = c(
      "raw_price_distributions", "row_accounting", "sampling_state_coverage",
      "sampling_zip_coverage", "prior_predictive", "parameter_recovery",
      "mcmc_trace", "mcmc_rank", "mcmc_acf", "fixed_effect_forest",
      "variance_component_posteriors", "state_caterpillar", "zip_caterpillar",
      "dic_heldout_comparison", "posterior_predictive_checks", "dollar_predictions",
      "data_dictionary_missingness", "simulation_recovery_results",
      "mcmc_settings_convergence", "m2_fixed_effect_variance_summary",
      "dic_decomposition", "heldout_metrics", "ppc_test_statistics"
    ),
    artifact_type = c(
      "figure", "table", rep("figure", 11L), "table", rep("figure", 2L),
      rep("table", 7L)
    ),
    file = c(
      "raw_price_histogram.png; raw_log_price_histogram.png",
      "report_cleaning_sampling_row_accounting.csv",
      "sampled_listings_per_state_before_after.png",
      "sampled_listings_per_zip_before_after.png",
      "prior_predictive_price_ranges.png",
      "simulation_recovery_scalar_parameter_recovery.png",
      "m2_traceplots.png", "m2_rank_histograms.png", "m2_acf.png",
      "m2_fixed_effect_forest.png", "m2_variance_component_posteriors.png",
      "m2_state_effects_caterpillar.png", "m2_selected_zip_effects_caterpillar.png",
      "report_dic_heldout_model_comparison.csv",
      "m0_posterior_predictive_checks.png; m1_posterior_predictive_checks.png; m2_posterior_predictive_checks.png",
      "m2_posterior_prediction_intervals.png",
      "report_data_dictionary_missingness.csv", "report_simulation_recovery_results.csv",
      "report_mcmc_settings_and_convergence.csv", "report_m2_fixed_effect_variance_summary.csv",
      "report_dic_decomposition.csv", "report_heldout_predictive_metrics.csv",
      "report_posterior_predictive_test_statistics.csv"
    ),
    title = c(
      "Raw and log listing-price distributions",
      "Cleaning and balanced-sampling row accounting",
      "Listings per selected state before and after sampling",
      "Listings per selected state--ZIP group before and after sampling",
      "Prior-predictive versus observed listing-price ranges",
      "Synthetic-data parameter recovery",
      "M2 trace diagnostics", "M2 rank diagnostics", "M2 autocorrelation diagnostics",
      "M2 fixed-effect forest plot",
      "M2 residual, state, and ZIP standard-deviation posteriors",
      "M2 state-effect caterpillar plot", "M2 selected ZIP-effect caterpillar plot",
      "DIC and held-out predictive model comparison",
      "Posterior-predictive observed-versus-replicated checks",
      "Dollar-scale posterior predictive distributions",
      "Data dictionary and missingness", "Simulation recovery results",
      "MCMC settings and convergence statistics",
      "M2 fixed-effect and variance-component posterior summaries",
      "DIC decomposition", "Held-out predictive metrics",
      "Posterior-predictive test statistics"
    ),
    units = c(
      "USD and log(USD)", "rows", "listings (log10 scale)", "listings per state--ZIP group (log10 scale)",
      "USD (log10 scale)", "standardized log-price units", rep("MCMC draws", 3L),
      "standardized log-price units", "standardized log-price units", rep("standardized log-price units", 2L),
      "DIC and standardized log-price predictive metrics", "standardized log-price units",
      "USD", "field-specific units", "standardized log-price units", "iterations and diagnostic statistics",
      "standardized log-price units", "DIC units", "standardized log-price units", "standardized log-price units"
    ),
    notes = c(
      "Generated by the raw-data audit.",
      "Combines saved audit row counts with the frozen sample manifest.",
      "Only states retained by the frozen sample are compared.",
      "Each retained ZIP has 30 listings before its train/test split.",
      "Generated by prior-predictive simulation before fitting.",
      "Blue intervals are 95% credible intervals; red crosses are known values.",
      "Pre-specified parameter selection; inspect all three diagnostic types together.",
      "Pre-specified parameter selection; inspect all three diagnostic types together.",
      "Pre-specified parameter selection; inspect all three diagnostic types together.",
      "Whiskers are central 95% posterior credible intervals.",
      "Densities use posterior draws from the converged M2 fit.",
      "Whiskers are central 95% posterior credible intervals.",
      "ZIPs are selected by training-sample size before inspecting effects.",
      "DIC is in-sample; held-out metrics evaluate the same frozen test listings.",
      "Red lines are observed statistics; histograms are replicated statistics.",
      "Dark and light intervals are 50% and 95% predictive intervals.",
      "Missingness is reported as a count and percentage.",
      "Empirical coverage is across converged repeated synthetic-data fits.",
      "R-hat, bulk/tail ESS, and MCSE use the pre-specified convergence criteria.",
      "Intervals are central 95% posterior credible intervals.",
      "Dbar + pD = DIC; delta DIC is relative to the smallest DIC.",
      "RMSE and MAE are on the standardized log-price scale; coverage is 95% predictive coverage.",
      "Replicated intervals are central 95% posterior-predictive intervals."
    ),
    stringsAsFactors = FALSE
  )
  expected_paths <- unlist(lapply(seq_len(nrow(artifacts)), function(index) {
    file_names <- trimws(strsplit(artifacts$file[index], ";", fixed = TRUE)[[1L]])
    base_dir <- if (artifacts$artifact_type[index] == "figure") figures_dir else tables_dir
    file.path(base_dir, file_names)
  }), use.names = FALSE)
  missing_paths <- expected_paths[!file.exists(expected_paths)]
  if (length(missing_paths) > 0L) {
    stop(
      "The report-artifact manifest cannot be completed because file(s) are missing: ",
      paste(missing_paths, collapse = ", "), ".",
      call. = FALSE
    )
  }
  artifacts$relative_directory <- ifelse(
    artifacts$artifact_type == "figure",
    file.path("output", "figures"),
    file.path("output", "tables")
  )
  write_report_csv(artifacts, tables_dir, "report_artifact_manifest.csv")
}

run_report_artifact_pipeline <- function() {
  tables_dir <- here::here("output", "tables")
  figures_dir <- here::here("output", "figures")
  draws_path <- here::here("output", "mcmc", "m2_draws.rds")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  posterior_assets <- make_m2_posterior_summary_assets()
  draw_matrix <- posterior::as_draws_matrix(readRDS(draws_path))
  save_fixed_effect_forest_plot(
    posterior_assets$fixed_effects,
    file.path(figures_dir, "m2_fixed_effect_forest.png")
  )
  save_variance_component_plot(
    draw_matrix,
    file.path(figures_dir, "m2_variance_component_posteriors.png")
  )

  raw_state_counts <- read_report_csv(
    file.path(tables_dir, "raw_listings_per_state.csv"),
    c("state", "listings")
  )
  raw_zip_counts <- read_report_csv(
    file.path(tables_dir, "raw_listings_per_state_zip.csv"),
    c("state", "zip_code", "listings")
  )
  sample_manifest <- read_report_csv(
    here::here("data", "processed", "sample_manifest.csv"),
    c("state", "zip_id", "split", "n")
  )
  sampling_assets <- save_sampling_coverage_plots(
    raw_state_counts,
    raw_zip_counts,
    sample_manifest,
    figures_dir
  )

  full_recovery <- read_report_csv(
    file.path(tables_dir, "simulation_recovery_full_scalar_parameters.csv"),
    c("label", "true_value", "posterior_mean", "q2_5", "q97_5")
  )
  save_simulation_recovery_plot(
    full_recovery,
    file.path(figures_dir, "simulation_recovery_scalar_parameter_recovery.png")
  )

  report_tables <- make_report_tables(tables_dir, draws_path)
  validate_upstream_report_figures(figures_dir)
  manifest_path <- make_report_artifact_manifest(tables_dir, figures_dir)

  message(
    "Point 19 complete. Saved report-ready tables and figures, including ",
    sampling_assets$state_path, ", ", sampling_assets$zip_path, ", and ",
    manifest_path, "."
  )
  invisible(list(
    posterior_assets = posterior_assets,
    sampling_assets = sampling_assets,
    report_tables = report_tables,
    manifest_path = manifest_path
  ))
}

if (sys.nframe() == 0L) {
  run_report_artifact_pipeline()
}
