# Pilot JAGS fit for the state-and-ZIP hierarchical model.
#
# This is deliberately not the final MCMC analysis. It uses a small,
# deterministic subset to expose model-compilation, indexing, initial-value,
# and runtime problems before a full fit is attempted.

required_fit_packages <- c("coda", "here", "rjags")
missing_fit_packages <- required_fit_packages[!vapply(
  required_fit_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_fit_packages) > 0L) {
  stop(
    "Install the package(s) required for the pilot JAGS fit: ",
    paste(missing_fit_packages, collapse = ", "),
    call. = FALSE
  )
}

source(here::here("scripts", "00_setup.R"))

make_pilot_data <- function(analysis_objects,
                            n_states = 5L,
                            n_zips_per_state = 8L,
                            n_per_zip = 20L,
                            seed = 20260621L) {
  required_objects <- c("train", "X_train")
  missing_objects <- setdiff(required_objects, names(analysis_objects))
  if (length(missing_objects) > 0L) {
    stop(
      "Analysis objects are missing required item(s): ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  train <- analysis_objects$train
  X_train <- analysis_objects$X_train
  required_columns <- c("y", "state", "zip_id")
  missing_columns <- setdiff(required_columns, names(train))
  if (length(missing_columns) > 0L) {
    stop(
      "Training data are missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.matrix(X_train) || nrow(X_train) != nrow(train) ||
      any(!is.finite(X_train)) || any(!is.finite(train$y))) {
    stop(
      "Training outcome and design matrix must be finite and aligned.",
      call. = FALSE
    )
  }

  group_counts <- stats::aggregate(
    rep(1L, nrow(train)),
    by = list(state = as.character(train$state), zip_id = as.character(train$zip_id)),
    FUN = sum
  )
  names(group_counts)[3L] <- "n_available"
  eligible_groups <- group_counts[group_counts$n_available >= n_per_zip, , drop = FALSE]
  eligible_states <- names(which(table(eligible_groups$state) >= n_zips_per_state))

  if (length(eligible_states) < n_states) {
    stop(
      "Only ", length(eligible_states), " states have at least ",
      n_zips_per_state, " ZIP groups with ", n_per_zip,
      " training listings; the pilot requires ", n_states, ".",
      call. = FALSE
    )
  }

  set.seed(seed)
  selected_states <- sort(sample(eligible_states, n_states))
  selected_groups <- do.call(
    rbind,
    lapply(selected_states, function(selected_state) {
      state_groups <- eligible_groups[
        eligible_groups$state == selected_state,
        ,
        drop = FALSE
      ]
      state_groups[sample(seq_len(nrow(state_groups)), n_zips_per_state), , drop = FALSE]
    })
  )
  selected_groups <- selected_groups[order(selected_groups$state, selected_groups$zip_id), ]

  pilot_rows <- unlist(lapply(seq_len(nrow(selected_groups)), function(group_index) {
    group_rows <- which(
      as.character(train$state) == selected_groups$state[group_index] &
        as.character(train$zip_id) == selected_groups$zip_id[group_index]
    )
    sort(sample(group_rows, n_per_zip))
  }), use.names = FALSE)

  pilot_train <- train[pilot_rows, , drop = FALSE]
  state_levels <- sort(unique(as.character(pilot_train$state)))
  zip_levels <- sort(unique(as.character(pilot_train$zip_id)))
  state_index <- match(as.character(pilot_train$state), state_levels)
  zip_index <- match(as.character(pilot_train$zip_id), zip_levels)
  zip_state_counts <- stats::aggregate(
    state_index,
    by = list(zip_index = zip_index),
    FUN = function(values) length(unique(values))
  )

  if (anyNA(state_index) || anyNA(zip_index) ||
      any(zip_state_counts$x != 1L)) {
    stop("Pilot state and ZIP indices are not a valid nested hierarchy.", call. = FALSE)
  }

  pilot_data <- list(
    N = length(pilot_rows),
    P = ncol(X_train),
    S = length(state_levels),
    Z = length(zip_levels),
    y = as.numeric(pilot_train$y),
    X = unname(X_train[pilot_rows, , drop = FALSE]),
    state = as.integer(state_index),
    zip = as.integer(zip_index)
  )

  expected_n <- n_states * n_zips_per_state * n_per_zip
  if (pilot_data$N != expected_n || pilot_data$S != n_states ||
      pilot_data$Z != n_states * n_zips_per_state) {
    stop("Pilot subset does not have the requested balanced dimensions.", call. = FALSE)
  }

  list(
    data = pilot_data,
    selection = selected_groups[, c("state", "zip_id", "n_available"), drop = FALSE],
    settings = list(
      n_states = n_states,
      n_zips_per_state = n_zips_per_state,
      n_per_zip = n_per_zip,
      seed = seed
    )
  )
}

fit_jags_model <- function(model_file, data_list, parameters,
                           n_chains = 4L, n_adapt = 2000L,
                           n_burn = 5000L, n_iter = 10000L,
                           thin = 1L, seed = 20260621L) {
  if (!file.exists(model_file)) {
    stop("JAGS model file not found: ", model_file, call. = FALSE)
  }
  if (!is.character(parameters) || length(parameters) == 0L) {
    stop("Provide at least one parameter to monitor.", call. = FALSE)
  }
  if (any(c(n_chains, n_adapt, n_burn, n_iter, thin) < 1L)) {
    stop("MCMC settings must all be positive integers.", call. = FALSE)
  }

  set.seed(seed)

  make_inits <- function(chain_id) {
    list(
      alpha = stats::rnorm(1L, 0, 0.2),
      beta = stats::rnorm(data_list$P, 0, 0.1),
      sigma = stats::runif(1L, 0.5, 1.5),
      sd_state = stats::runif(1L, 0.1, 0.8),
      sd_zip = stats::runif(1L, 0.1, 0.8),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = as.integer(seed + 1000L * chain_id)
    )
  }

  inits <- lapply(seq_len(n_chains), make_inits)
  model_text <- paste(readLines(model_file, warn = FALSE), collapse = "\n")
  if (!grepl("sd_state", model_text, fixed = TRUE)) {
    inits <- lapply(inits, function(values) values[setdiff(names(values), "sd_state")])
  }
  if (!grepl("sd_zip", model_text, fixed = TRUE)) {
    inits <- lapply(inits, function(values) values[setdiff(names(values), "sd_zip")])
  }

  model <- rjags::jags.model(
    file = model_file,
    data = data_list,
    inits = inits,
    n.chains = n_chains,
    n.adapt = n_adapt,
    quiet = TRUE
  )
  stats::update(model, n.iter = n_burn, progress.bar = "none")
  draws <- rjags::coda.samples(
    model = model,
    variable.names = parameters,
    n.iter = n_iter,
    thin = thin,
    progress.bar = "none"
  )

  list(model = model, draws = draws, inits = inits)
}

validate_pilot_draws <- function(draws) {
  chain_matrices <- lapply(draws, as.matrix)
  if (length(chain_matrices) < 2L || any(vapply(
    chain_matrices,
    function(values) nrow(values) < 2L || ncol(values) == 0L,
    FUN.VALUE = logical(1)
  ))) {
    stop("Pilot fit did not retain usable draws from every chain.", call. = FALSE)
  }

  all_draws <- do.call(rbind, chain_matrices)
  if (any(!is.finite(all_draws))) {
    stop("Pilot fit produced non-finite monitored values.", call. = FALSE)
  }

  parameter_sd <- apply(all_draws, 2L, stats::sd)
  if (any(!is.finite(parameter_sd)) || any(parameter_sd == 0)) {
    stuck <- names(parameter_sd)[!is.finite(parameter_sd) | parameter_sd == 0]
    stop(
      "Pilot fit has monitored value(s) that do not move: ",
      paste(stuck, collapse = ", "),
      call. = FALSE
    )
  }

  sd_parameters <- grep("^(sigma|sd_state|sd_zip)$", colnames(all_draws), value = TRUE)
  if (length(sd_parameters) > 0L) {
    sd_draws <- all_draws[, sd_parameters, drop = FALSE]
    if (any(sd_draws <= 0) || any(apply(sd_draws, 2L, max) < 1e-6)) {
      stop("Pilot fit has a standard deviation at its lower boundary.", call. = FALSE)
    }
  }

  coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf
}

save_pilot_outputs <- function(pilot_fit, pilot_data, elapsed_seconds, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  diagnostics <- validate_pilot_draws(pilot_fit$draws)
  summary_values <- summary(pilot_fit$draws)$statistics
  monitored <- c("alpha", "sigma", "sd_state", "sd_zip")
  monitored <- intersect(monitored, colnames(as.matrix(pilot_fit$draws[[1L]])))

  saveRDS(
    list(
      draws = pilot_fit$draws,
      selection = pilot_data$selection,
      settings = pilot_data$settings,
      elapsed_seconds = elapsed_seconds,
      diagnostics = diagnostics
    ),
    file.path(output_dir, "pilot_m2_fit.rds")
  )
  utils::write.csv(
    data.frame(
      parameter = rownames(summary_values),
      mean = summary_values[, "Mean"],
      sd = summary_values[, "SD"],
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "pilot_m2_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(
      n_states = pilot_data$settings$n_states,
      n_zips_per_state = pilot_data$settings$n_zips_per_state,
      n_per_zip = pilot_data$settings$n_per_zip,
      n_chains = length(pilot_fit$draws),
      elapsed_seconds = elapsed_seconds,
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "pilot_m2_runtime.csv"),
    row.names = FALSE
  )

  grDevices::png(
    file.path(output_dir, "pilot_m2_traceplots.png"),
    width = 1800,
    height = 1200,
    res = 160
  )
  if (length(monitored) > 0L) {
    plot(pilot_fit$draws[, monitored, drop = FALSE])
  }
  grDevices::dev.off()

  invisible(diagnostics)
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
pilot_data <- make_pilot_data(analysis_objects)
pilot_parameters <- c("alpha", "beta", "sigma", "sd_state", "sd_zip", "u", "v")
pilot_model <- here::here("models", "m2_state_zip.jags")

started_at <- Sys.time()
pilot_fit <- fit_jags_model(
  model_file = pilot_model,
  data_list = pilot_data$data,
  parameters = pilot_parameters,
  n_adapt = 1000L,
  n_burn = 2000L,
  n_iter = 4000L
)
elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))

diagnostics <- save_pilot_outputs(
  pilot_fit = pilot_fit,
  pilot_data = pilot_data,
  elapsed_seconds = elapsed_seconds,
  output_dir = here::here("output", "mcmc")
)

message(
  "M2 pilot fit passed compilation and draw checks: ",
  pilot_data$data$N, " listings across ", pilot_data$data$S,
  " states and ", pilot_data$data$Z, " ZIP groups; runtime ",
  format(round(elapsed_seconds, 1), nsmall = 1), " seconds. ",
  "Pilot artefacts saved under output/mcmc/."
)
