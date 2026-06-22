# Prior-predictive simulation for the planned M2 state-and-ZIP model.
#
# This script is intentionally run before fitting any model.  It draws all
# parameters and geographic effects from their priors, transforms replicated
# standardized log prices back to dollars, and compares their ranges with the
# observed training-listing prices.

source(here::here("scripts", "00_setup.R"))

# The original weak priors generated billion-dollar homes for nearly every
# draw because a few standardized covariates are extreme.  These scales were
# tightened before any MCMC fit and are shared with all three JAGS models.
prior_scales <- c(
  alpha = 1,
  beta = 0.03,
  sigma = 0.5,
  sd_state = 0.3,
  sd_zip = 0.3
)

simulate_prior_y <- function(X, state, zip, S, Z, y_mean, y_sd) {
  alpha <- stats::rnorm(1, 0, prior_scales[["alpha"]])
  beta <- stats::rnorm(ncol(X), 0, prior_scales[["beta"]])
  sigma <- abs(stats::rnorm(1, 0, prior_scales[["sigma"]]))
  sd_state <- abs(stats::rnorm(1, 0, prior_scales[["sd_state"]]))
  sd_zip <- abs(stats::rnorm(1, 0, prior_scales[["sd_zip"]]))

  u <- stats::rnorm(S, 0, sd_state)
  v <- stats::rnorm(Z, 0, sd_zip)
  mu <- as.vector(alpha + X %*% beta + u[state] + v[zip])
  y_rep <- stats::rnorm(length(mu), mu, sigma)

  exp(y_rep * y_sd + y_mean)
}

validate_prior_inputs <- function(jags_data, y_mean, y_sd) {
  required_fields <- c("N", "P", "S", "Z", "X", "state", "zip")
  missing_fields <- setdiff(required_fields, names(jags_data))
  if (length(missing_fields) > 0) {
    stop(
      "JAGS data are missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.matrix(jags_data$X) ||
      nrow(jags_data$X) != jags_data$N ||
      ncol(jags_data$X) != jags_data$P ||
      any(!is.finite(jags_data$X))) {
    stop("`X` must be a finite N-by-P matrix.", call. = FALSE)
  }

  if (length(jags_data$state) != jags_data$N ||
      length(jags_data$zip) != jags_data$N ||
      anyNA(jags_data$state) || anyNA(jags_data$zip) ||
      any(jags_data$state < 1L | jags_data$state > jags_data$S) ||
      any(jags_data$zip < 1L | jags_data$zip > jags_data$Z)) {
    stop("State and ZIP indices must be one-based values within their declared ranges.", call. = FALSE)
  }

  if (!is.numeric(y_mean) || length(y_mean) != 1L || !is.finite(y_mean) ||
      !is.numeric(y_sd) || length(y_sd) != 1L || !is.finite(y_sd) || y_sd <= 0) {
    stop("Outcome scaling values must be finite, with a positive `y_sd`.", call. = FALSE)
  }
}

summarise_prices <- function(prices) {
  finite_prices <- prices[is.finite(prices) & prices > 0]
  summary_values <- c(
    n = length(prices),
    finite_fraction = length(finite_prices) / length(prices),
    min = NA_real_,
    p01 = NA_real_,
    median = NA_real_,
    p99 = NA_real_,
    max = NA_real_
  )

  if (length(finite_prices) > 0) {
    summary_values[c("min", "p01", "median", "p99", "max")] <- c(
      min(finite_prices),
      stats::quantile(finite_prices, 0.01, names = FALSE),
      stats::median(finite_prices),
      stats::quantile(finite_prices, 0.99, names = FALSE),
      max(finite_prices)
    )
  }

  summary_values
}

analysis_path <- here::here("data", "processed", "analysis_objects.rds")
if (!file.exists(analysis_path)) {
  stop(
    "Prepared analysis objects not found: ", analysis_path,
    ". Run scripts/03_prepare_jags_data.R first.",
    call. = FALSE
  )
}

analysis_objects <- readRDS(analysis_path)
required_analysis_fields <- c("jags_train", "train", "y_mean", "y_sd")
missing_analysis_fields <- setdiff(required_analysis_fields, names(analysis_objects))
if (length(missing_analysis_fields) > 0) {
  stop(
    "Analysis objects are missing required field(s): ",
    paste(missing_analysis_fields, collapse = ", "),
    call. = FALSE
  )
}
if (!"log_price" %in% names(analysis_objects$train) ||
    any(!is.finite(analysis_objects$train$log_price))) {
  stop("Training data must contain finite `log_price` values.", call. = FALSE)
}

validate_prior_inputs(
  jags_data = analysis_objects$jags_train,
  y_mean = analysis_objects$y_mean,
  y_sd = analysis_objects$y_sd
)

n_prior_draws <- 200L
prior_seed <- 20260621L
set.seed(prior_seed)

observed_prices <- exp(analysis_objects$train$log_price)
observed_summary <- summarise_prices(observed_prices)

prior_summary <- do.call(
  rbind,
  lapply(seq_len(n_prior_draws), function(draw) {
    replicated_prices <- simulate_prior_y(
      X = analysis_objects$jags_train$X,
      state = analysis_objects$jags_train$state,
      zip = analysis_objects$jags_train$zip,
      S = analysis_objects$jags_train$S,
      Z = analysis_objects$jags_train$Z,
      y_mean = analysis_objects$y_mean,
      y_sd = analysis_objects$y_sd
    )
    c(draw = draw, summarise_prices(replicated_prices))
  })
)
prior_summary <- as.data.frame(prior_summary, check.names = FALSE)

range_summary <- rbind(
  cbind(
    distribution = "observed",
    draw = NA_integer_,
    as.data.frame(as.list(observed_summary), check.names = FALSE)
  ),
  cbind(distribution = "prior_predictive", prior_summary)
)

tables_dir <- here::here("output", "tables")
figures_dir <- here::here("output", "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.csv(
  range_summary,
  file = file.path(tables_dir, "prior_predictive_range_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(
    n_prior_draws = n_prior_draws,
    seed = prior_seed,
    alpha_sd = prior_scales[["alpha"]],
    beta_sd = prior_scales[["beta"]],
    sigma_sd = prior_scales[["sigma"]],
    sd_state_sd = prior_scales[["sd_state"]],
    sd_zip_sd = prior_scales[["sd_zip"]],
    revision_reason = "Original scales generated pathological prior-predictive price ranges with extreme standardized covariates."
  ),
  file = file.path(tables_dir, "prior_predictive_settings.csv"),
  row.names = FALSE
)

plot_rows <- is.finite(prior_summary$min) & prior_summary$min > 0 &
  is.finite(prior_summary$median) & prior_summary$median > 0 &
  is.finite(prior_summary$max) & prior_summary$max > 0
if (!any(plot_rows)) {
  stop("No finite, positive prior-predictive price ranges are available to plot.", call. = FALSE)
}

plot_prior <- prior_summary[plot_rows, , drop = FALSE]
observed_log_range <- log10(observed_summary[c("min", "median", "max")])
prior_log_range <- log10(c(plot_prior$min, plot_prior$median, plot_prior$max))
y_limits <- range(c(observed_log_range, prior_log_range))
y_padding <- max(diff(y_limits) * 0.05, 0.05)

grDevices::png(
  filename = file.path(figures_dir, "prior_predictive_price_ranges.png"),
  width = 1600,
  height = 1000,
  res = 180
)
graphics::plot(
  x = plot_prior$draw,
  y = log10(plot_prior$median),
  type = "n",
  ylim = y_limits + c(-y_padding, y_padding),
  xlab = "Prior-predictive draw",
  ylab = expression(log[10]("listing price (USD)")),
  main = "Prior-predictive price ranges versus observed training prices"
)
graphics::segments(
  x0 = plot_prior$draw,
  y0 = log10(plot_prior$min),
  x1 = plot_prior$draw,
  y1 = log10(plot_prior$max),
  col = grDevices::adjustcolor("steelblue", alpha.f = 0.25)
)
graphics::points(
  x = plot_prior$draw,
  y = log10(plot_prior$median),
  pch = 16,
  cex = 0.35,
  col = "steelblue4"
)
graphics::abline(h = observed_log_range["min"], col = "firebrick", lty = 2)
graphics::abline(h = observed_log_range["median"], col = "firebrick", lty = 1, lwd = 2)
graphics::abline(h = observed_log_range["max"], col = "firebrick", lty = 2)
graphics::legend(
  "topleft",
  legend = c("Prior draw min-to-max range", "Prior draw median", "Observed min / median / max"),
  col = c("steelblue", "steelblue4", "firebrick"),
  lty = c(1, NA, 1),
  lwd = c(1, NA, 2),
  pch = c(NA, 16, NA),
  bty = "n"
)
grDevices::dev.off()

if (any(prior_summary$finite_fraction < 1)) {
  warning(
    "At least one prior draw produced non-finite dollar prices. Review the prior-predictive range summary before fitting MCMC.",
    call. = FALSE
  )
}

message(
  "Prior-predictive checks complete: ", n_prior_draws,
  " draws. Saved range summaries to output/tables and the comparison figure to output/figures."
)
