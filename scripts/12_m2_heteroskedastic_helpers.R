# Shared helpers for the selected heteroskedastic M2 model.
#
# The current caveat-fix branch promotes the heteroskedastic state/ZIP model as
# the official M2 candidate.  These helpers are shared by fitting, DIC,
# diagnostics, held-out prediction, PPC, and report-artifact scripts so the
# pipeline does not depend on any removed experimental Student-t scripts.

required_hetero_helper_packages <- c(
  "coda", "ggplot2", "here", "posterior", "readr", "rjags"
)
missing_hetero_helper_packages <- required_hetero_helper_packages[!vapply(
  required_hetero_helper_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_hetero_helper_packages) > 0L) {
  stop(
    "Install the package(s) required for heteroskedastic M2 helpers: ",
    paste(missing_hetero_helper_packages, collapse = ", "),
    call. = FALSE
  )
}

source(here::here("scripts", "00_setup.R"))

hetero_read_positive_integer <- function(name, default) {
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

hetero_read_flag <- function(name, default = FALSE) {
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

hetero_indexed_parameter_names <- function(prefix, count) {
  paste0(prefix, "[", seq_len(count), "]")
}

validate_hetero_analysis_objects <- function(analysis_objects) {
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
  if (nrow(X_train) != analysis_objects$jags_train$N ||
      ncol(X_train) != length(analysis_objects$predictor_names) ||
      ncol(X_test) != length(analysis_objects$predictor_names) ||
      any(!is.finite(X_train)) ||
      any(!is.finite(X_test))) {
    stop("Saved design matrices are inconsistent with the JAGS objects.", call. = FALSE)
  }
  invisible(TRUE)
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

prepare_heteroskedastic_m2_data <- function(
    analysis_objects,
    sigma_prior_center = 0.5,
    gamma0_prior_scale = 0.5,
    gamma_prior_scale = 0.2,
    sigma_state_prior_scale = 0.2) {
  validate_hetero_analysis_objects(analysis_objects)

  if (any(c(
    sigma_prior_center, gamma0_prior_scale,
    gamma_prior_scale, sigma_state_prior_scale
  ) <= 0)) {
    stop("Heteroskedastic M2 prior scales must be positive.", call. = FALSE)
  }

  data_list <- augment_zip_centering_data(analysis_objects$jags_train)
  X_sigma <- as.matrix(analysis_objects$X_train)
  data_list$X_sigma <- X_sigma
  data_list$P_sigma <- as.integer(ncol(X_sigma))
  data_list$log_sigma_prior_mean <- log(sigma_prior_center)
  data_list$tau_gamma0_prior <- 1 / (gamma0_prior_scale^2)
  data_list$tau_gamma_prior <- 1 / (gamma_prior_scale^2)
  data_list$tau_sigma_state_prior <- 1 / (sigma_state_prior_scale^2)
  data_list
}

hetero_required_parameters <- function(n_predictors, n_states, n_zips) {
  c(
    "alpha", "gamma0", "sd_sigma_state", "sd_state", "sd_zip",
    hetero_indexed_parameter_names("beta", n_predictors),
    hetero_indexed_parameter_names("gamma", n_predictors),
    hetero_indexed_parameter_names("h", n_states),
    hetero_indexed_parameter_names("u", n_states),
    hetero_indexed_parameter_names("v", n_zips)
  )
}

hetero_chunk_sizes <- function(total_iter, chunk_iter) {
  starts <- seq.int(1L, total_iter, by = chunk_iter)
  as.integer(pmin(chunk_iter, total_iter - starts + 1L))
}

hetero_update_log_density_accumulator <- function(log_max, scaled_sum, log_density) {
  new_log_max <- pmax(log_max, log_density)
  list(
    log_max = new_log_max,
    scaled_sum =
      scaled_sum * exp(log_max - new_log_max) +
      exp(log_density - new_log_max)
  )
}

summarise_hetero_convergence <- function(draws) {
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

select_ppc_state_ids <- function(state_levels, selected_states = NULL,
                                 max_states = 4L) {
  state_levels <- as.character(state_levels)
  valid_max_states <- is.numeric(max_states) &&
    length(max_states) == 1L &&
    is.finite(max_states) &&
    max_states >= 1L &&
    max_states %% 1 == 0
  if (!valid_max_states) {
    stop("max_states must be one positive integer.", call. = FALSE)
  }
  max_states <- as.integer(max_states)
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

# Backward-compatible aliases used by the original heteroskedastic script.
student_t_read_positive_integer <- hetero_read_positive_integer
student_t_read_flag <- hetero_read_flag
student_t_indexed_parameter_names <- hetero_indexed_parameter_names
student_t_chunk_sizes <- hetero_chunk_sizes
student_t_update_log_density_accumulator <- hetero_update_log_density_accumulator
summarise_student_t_convergence <- summarise_hetero_convergence
validate_student_t_analysis_objects <- validate_hetero_analysis_objects
