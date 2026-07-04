# Diagnose convergence and Monte Carlo precision for the selected fitted M2
# model.
#
# This script evaluates every scalar parameter in the production M2 chain. It
# also creates the trace, rank-histogram, and autocorrelation plots required
# for the report. Run scripts/05_fit_models.R first.

required_diagnostic_packages <- c("bayesplot", "here", "posterior", "readr")
missing_diagnostic_packages <- required_diagnostic_packages[!vapply(
  required_diagnostic_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_diagnostic_packages) > 0L) {
  stop(
    "Install the package(s) required for MCMC diagnostics: ",
    paste(missing_diagnostic_packages, collapse = ", "),
    call. = FALSE
  )
}

draws_path <- here::here("output", "mcmc", "m2_heteroskedastic_xall_state_draws.rds")
if (!file.exists(draws_path)) {
  stop(
    "M2 production draws not found: ", draws_path,
    ". Run scripts/05_fit_models.R first.",
    call. = FALSE
  )
}

tables_dir <- here::here("output", "tables")
figures_dir <- here::here("output", "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

draws_m2 <- readRDS(draws_path)
draws_array_m2 <- posterior::as_draws_array(draws_m2)
if (posterior::ndraws(draws_array_m2) == 0L ||
    posterior::nchains(draws_array_m2) < 2L) {
  stop(
    "M2 diagnostics require non-empty draws from at least two chains.",
    call. = FALSE
  )
}

available_parameters <- posterior::variables(draws_array_m2)
required_parameters <- c("alpha", "gamma0", "sd_sigma_state", "sd_state", "sd_zip")
missing_parameters <- setdiff(required_parameters, available_parameters)
if (length(missing_parameters) > 0L) {
  stop(
    "M2 draws are missing required parameter(s): ",
    paste(missing_parameters, collapse = ", "),
    call. = FALSE
  )
}

summary_m2 <- posterior::summarise_draws(
  draws_array_m2,
  mean,
  median,
  sd,
  posterior::quantile2,
  posterior::rhat,
  posterior::ess_bulk,
  posterior::ess_tail,
  posterior::mcse_mean
)
summary_m2 <- as.data.frame(summary_m2)
names(summary_m2) <- sub("^posterior::", "", names(summary_m2))

# A mean MCSE below 5% of the posterior SD makes simulation error small
# relative to posterior uncertainty. For degenerate parameters, an undefined
# ratio is deliberately marked as failing instead of being silently accepted.
summary_m2$mcse_mean_sd_ratio <- summary_m2$mcse_mean / summary_m2$sd
summary_m2$rhat_pass <- is.finite(summary_m2$rhat) & summary_m2$rhat < 1.01
summary_m2$bulk_ess_pass <-
  is.finite(summary_m2$ess_bulk) & summary_m2$ess_bulk > 400
summary_m2$tail_ess_pass <-
  is.finite(summary_m2$ess_tail) & summary_m2$ess_tail > 400
summary_m2$mcse_pass <-
  is.finite(summary_m2$mcse_mean_sd_ratio) &
  summary_m2$mcse_mean_sd_ratio < 0.05
summary_m2$converged <-
  summary_m2$rhat_pass &
  summary_m2$bulk_ess_pass &
  summary_m2$tail_ess_pass &
  summary_m2$mcse_pass

readr::write_csv(
  summary_m2,
  file.path(tables_dir, "m2_posterior_summary.csv")
)

# The report selection is deterministic: the two first fixed-effect columns,
# then the first two states (alphabetical) and first three ZIPs (sorted) from
# the frozen training-data indices. This prevents choosing especially tidy
# chains after inspecting the diagnostics.
select_indexed_parameters <- function(prefix, count) {
  pattern <- paste0("^", prefix, "\\[([0-9]+)\\]$")
  matched <- grep(pattern, available_parameters, value = TRUE)
  indices <- as.integer(sub(pattern, "\\1", matched))
  matched[order(indices)][seq_len(min(count, length(matched)))]
}

beta_parameters <- select_indexed_parameters("beta", 2L)
state_parameters <- select_indexed_parameters("u", 2L)
zip_parameters <- select_indexed_parameters("v", 3L)
if (length(beta_parameters) < 2L || length(state_parameters) < 2L ||
    length(zip_parameters) < 3L) {
  stop(
    "M2 draws must contain at least two beta effects, two state effects, and ",
    "three ZIP effects for the pre-specified diagnostic plots.",
    call. = FALSE
  )
}

sigma_state_parameters <- select_indexed_parameters("h", 2L)
if (length(sigma_state_parameters) < 2L) {
  stop(
    "M2 draws must contain at least two residual-scale state effects for the ",
    "pre-specified diagnostic plots.",
    call. = FALSE
  )
}

selected_parameters <- c(
  "alpha", "gamma0", "sd_sigma_state", "sd_state", "sd_zip",
  beta_parameters, sigma_state_parameters, state_parameters, zip_parameters
)

analysis_objects_path <- here::here("data", "processed", "analysis_objects.rds")
selection_table <- data.frame(
  parameter = selected_parameters,
  diagnostic_group = c(
    rep("global", 5L),
    rep("fixed_effect", length(beta_parameters)),
    rep("residual_scale_state_effect", length(sigma_state_parameters)),
    rep("state_effect", length(state_parameters)),
    rep("zip_effect", length(zip_parameters))
  ),
  selection_rule = c(
    rep("Required global M2 parameter", 5L),
    rep("First two design-matrix columns", length(beta_parameters)),
    rep("First two alphabetically ordered states", length(sigma_state_parameters)),
    rep("First two alphabetically ordered states", length(state_parameters)),
    rep("First three sorted ZIP codes", length(zip_parameters))
  ),
  label = selected_parameters,
  stringsAsFactors = FALSE
)

if (file.exists(analysis_objects_path)) {
  analysis_objects <- readRDS(analysis_objects_path)
  index_for <- function(parameters) {
    as.integer(sub("^[^\\[]+\\[([0-9]+)\\]$", "\\1", parameters))
  }
  beta_indices <- index_for(beta_parameters)
  state_indices <- index_for(state_parameters)
  sigma_state_indices <- index_for(sigma_state_parameters)
  zip_indices <- index_for(zip_parameters)

  if (!is.null(analysis_objects$predictor_names) &&
      length(analysis_objects$predictor_names) >= max(beta_indices)) {
    selection_table$label[selection_table$diagnostic_group == "fixed_effect"] <-
      paste0(beta_parameters, " (", analysis_objects$predictor_names[beta_indices], ")")
  }
  if (!is.null(analysis_objects$state_levels) &&
      length(analysis_objects$state_levels) >= max(state_indices)) {
    selection_table$label[selection_table$diagnostic_group == "state_effect"] <-
      paste0(state_parameters, " (", analysis_objects$state_levels[state_indices], ")")
  }
  if (!is.null(analysis_objects$state_levels) &&
      length(analysis_objects$state_levels) >= max(sigma_state_indices)) {
    selection_table$label[selection_table$diagnostic_group == "residual_scale_state_effect"] <-
      paste0(
        sigma_state_parameters,
        " (", analysis_objects$state_levels[sigma_state_indices], ")"
      )
  }
  if (!is.null(analysis_objects$zip_levels) &&
      length(analysis_objects$zip_levels) >= max(zip_indices)) {
    selection_table$label[selection_table$diagnostic_group == "zip_effect"] <-
      paste0(zip_parameters, " (", analysis_objects$zip_levels[zip_indices], ")")
  }
}

selected_summary <- merge(
  selection_table,
  summary_m2,
  by.x = "parameter",
  by.y = "variable",
  all.x = TRUE,
  sort = FALSE
)
readr::write_csv(
  selected_summary,
  file.path(tables_dir, "m2_selected_diagnostics.csv")
)

save_diagnostic_plot <- function(path, height, plot_expression) {
  grDevices::png(path, width = 1800, height = height, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_expression)
}

save_diagnostic_plot(
  file.path(figures_dir, "m2_traceplots.png"),
  1400,
  bayesplot::mcmc_trace(draws_array_m2, pars = selected_parameters)
)
save_diagnostic_plot(
  file.path(figures_dir, "m2_rank_histograms.png"),
  1400,
  bayesplot::mcmc_rank_hist(draws_array_m2, pars = selected_parameters)
)
save_diagnostic_plot(
  file.path(figures_dir, "m2_acf.png"),
  1200,
  bayesplot::mcmc_acf(draws_array_m2, pars = selected_parameters)
)

failed_parameters <- summary_m2$variable[!summary_m2$converged]
message(
  "Saved M2 posterior diagnostics for ", nrow(summary_m2),
  " parameters. Convergence criteria passed for ",
  nrow(summary_m2) - length(failed_parameters), ".",
  if (length(failed_parameters) > 0L) {
    paste0(" Investigate: ", paste(failed_parameters, collapse = ", "), ".")
  } else {
    ""
  }
)
