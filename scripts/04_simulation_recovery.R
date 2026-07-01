# Simulation-based parameter recovery for the M2 state-and-ZIP model.
#
# The full-design recovery uses the frozen training design matrix exactly as it
# will be used in the observed-data fit.  The repeated experiment keeps the
# same M2 code, priors, initialisation, and MCMC settings, but uses a balanced
# subset of that design so that 20 replications are practical to run.

required_recovery_packages <- c("coda", "here", "posterior", "rjags")
missing_recovery_packages <- required_recovery_packages[!vapply(
  required_recovery_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_recovery_packages) > 0L) {
  stop(
    "Install the package(s) required for simulation recovery: ",
    paste(missing_recovery_packages, collapse = ", "),
    call. = FALSE
  )
}

source(here::here("scripts", "00_setup.R"))

recovery_seed <- 20260621L
model_path <- here::here("models", "m2_state_zip.jags")
analysis_path <- here::here("data", "processed", "analysis_objects.rds")
tables_dir <- here::here("output", "tables")
figures_dir <- here::here("output", "figures")

# These are the production MCMC settings specified for M2.  Environment
# variables are deliberately available for a short smoke test only; do not use
# reduced iterations for the substantive recovery evidence.
read_positive_integer <- function(name, default) {
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

read_flag <- function(name, default = FALSE) {
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

mcmc_settings <- list(
  n_chains = read_positive_integer("RECOVERY_N_CHAINS", 4L),
  n_adapt = read_positive_integer("RECOVERY_N_ADAPT", 2000L),
  n_burn = read_positive_integer("RECOVERY_N_BURN", 5000L),
  n_iter = read_positive_integer("RECOVERY_N_ITER", 10000L),
  thin = read_positive_integer("RECOVERY_THIN", 1L)
)
n_replications <- read_positive_integer("RECOVERY_REPLICATIONS", 20L)
reduced_design_settings <- list(
  n_states = read_positive_integer("RECOVERY_REDUCED_N_STATES", 5L),
  n_zips_per_state = read_positive_integer("RECOVERY_REDUCED_N_ZIPS_PER_STATE", 5L),
  n_per_zip = read_positive_integer("RECOVERY_REDUCED_N_PER_ZIP", 12L)
)
default_reduced_design_settings <- list(
  n_states = 5L,
  n_zips_per_state = 5L,
  n_per_zip = 12L
)
make_safe_label <- function(value) {
  value <- gsub("[^A-Za-z0-9_-]+", "_", value)
  value <- gsub("_+", "_", value)
  gsub("^_|_$", "", value)
}
custom_reduced_design <- !identical(
  unlist(reduced_design_settings, use.names = TRUE),
  unlist(default_reduced_design_settings, use.names = TRUE)
)
reduced_design_label <- Sys.getenv("RECOVERY_REDUCED_LABEL", unset = "")
if (!nzchar(reduced_design_label) && custom_reduced_design) {
  reduced_design_label <- sprintf(
    "%dstates_%dzips_%dperzip",
    reduced_design_settings$n_states,
    reduced_design_settings$n_zips_per_state,
    reduced_design_settings$n_per_zip
  )
}
reduced_design_label <- make_safe_label(reduced_design_label)
reduced_output_suffix <- if (nzchar(reduced_design_label)) {
  paste0("_", reduced_design_label)
} else {
  ""
}
# The full frozen design has 14,400 observations and 600 ZIP effects, so its
# production recovery fit is intentionally opt-in rather than an accidental
# multi-hour interactive job. It remains mandatory before real-data claims.
run_full_recovery <- read_flag("RECOVERY_RUN_FULL", default = FALSE)
# By default, run the practical repeated reduced-design check when the
# full-design check is not requested.  If RECOVERY_RUN_FULL=true, the script
# now runs only that full-design check unless RECOVERY_RUN_REPEATED=true is
# also set.  This avoids accidentally launching 21 production-setting JAGS
# fits from a command that was intended to refresh the full-design evidence.
run_repeated_recovery <- read_flag(
  "RECOVERY_RUN_REPEATED",
  default = !run_full_recovery
)
force_full_recovery <- read_flag("RECOVERY_FORCE_FULL", default = FALSE)
force_repeated_recovery <- read_flag("RECOVERY_FORCE_REPEATED", default = FALSE)
resume_repeated_recovery <- read_flag("RECOVERY_RESUME_REPEATED", default = TRUE)
if (!run_full_recovery && !run_repeated_recovery) {
  stop(
    "Nothing to run: set RECOVERY_RUN_FULL=true or RECOVERY_RUN_REPEATED=true.",
    call. = FALSE
  )
}

validate_jags_data <- function(jags_data) {
  required_fields <- c("N", "P", "S", "Z", "y", "X", "state", "zip")
  missing_fields <- setdiff(required_fields, names(jags_data))
  if (length(missing_fields) > 0L) {
    stop(
      "JAGS data are missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  dimensions <- unlist(jags_data[c("N", "P", "S", "Z")], use.names = TRUE)
  if (any(!is.finite(dimensions)) || any(dimensions < 1L) ||
      any(dimensions != as.integer(dimensions))) {
    stop("N, P, S, and Z must be positive integers.", call. = FALSE)
  }
  if (!is.matrix(jags_data$X) || nrow(jags_data$X) != jags_data$N ||
      ncol(jags_data$X) != jags_data$P || any(!is.finite(jags_data$X))) {
    stop("`X` must be a finite N-by-P matrix.", call. = FALSE)
  }
  if (length(jags_data$y) != jags_data$N || any(!is.finite(jags_data$y))) {
    stop("`y` must contain one finite value per row.", call. = FALSE)
  }
  for (index_name in c("state", "zip")) {
    upper_bound <- if (index_name == "state") jags_data$S else jags_data$Z
    index <- jags_data[[index_name]]
    if (length(index) != jags_data$N || anyNA(index) ||
        any(index != as.integer(index)) || any(index < 1L | index > upper_bound)) {
      stop(
        "`", index_name, "` must contain one-based integer indices within range.",
        call. = FALSE
      )
    }
  }

  zip_state <- tapply(jags_data$state, jags_data$zip, function(values) length(unique(values)))
  if (length(zip_state) != jags_data$Z || any(zip_state != 1L)) {
    stop("Each ZIP must belong to exactly one state in the M2 hierarchy.", call. = FALSE)
  }

  invisible(TRUE)
}

zip_state_lookup <- function(jags_data) {
  validate_jags_data(jags_data)
  vapply(seq_len(jags_data$Z), function(zip_id) {
    states <- unique(as.integer(jags_data$state[jags_data$zip == zip_id]))
    states <- states[!is.na(states)]
    if (length(states) != 1L) {
      stop("Each ZIP must belong to exactly one state.", call. = FALSE)
    }
    states
  }, integer(1L))
}

augment_zip_centering_data <- function(jags_data) {
  zip_state <- zip_state_lookup(jags_data)
  n_zips_by_state <- tabulate(zip_state, nbins = jags_data$S)
  if (any(n_zips_by_state < 1L)) {
    stop("Every state must have at least one ZIP for M2 ZIP centering.", call. = FALSE)
  }

  k_max <- max(n_zips_by_state)
  zip_in_state <- matrix(1L, nrow = jags_data$S, ncol = k_max)
  zip_in_state_weight <- matrix(0, nrow = jags_data$S, ncol = k_max)
  for (state_id in seq_len(jags_data$S)) {
    state_zips <- which(zip_state == state_id)
    zip_in_state[state_id, seq_along(state_zips)] <- state_zips
    zip_in_state_weight[state_id, seq_along(state_zips)] <- 1
  }

  jags_data$K <- as.integer(k_max)
  jags_data$zip_state <- as.integer(zip_state)
  jags_data$n_zips_by_state <- as.integer(n_zips_by_state)
  jags_data$zip_in_state <- zip_in_state
  jags_data$zip_in_state_weight <- zip_in_state_weight
  jags_data
}

make_truth <- function(P) {
  beta <- c(0.45, 0.08, 0.06, 0.16, rep(0.05, P - 4L))
  if (length(beta) != P) {
    stop("The M2 truth does not match the design-matrix column count.", call. = FALSE)
  }
  list(
    alpha = 0.00,
    beta = beta,
    sigma = 0.55,
    sd_state = 0.30,
    sd_zip = 0.25
  )
}

truth_vector <- function(truth) {
  c(
    alpha = truth$alpha,
    stats::setNames(truth$beta, paste0("beta[", seq_along(truth$beta), "]")),
    sigma = truth$sigma,
    sd_state = truth$sd_state,
    sd_zip = truth$sd_zip
  )
}

simulate_m2_data <- function(jags_data, truth) {
  validate_jags_data(jags_data)
  jags_data <- augment_zip_centering_data(jags_data)
  scalar_truth <- truth_vector(truth)
  if (length(truth$beta) != jags_data$P ||
      any(!is.finite(c(scalar_truth, truth$sd_state, truth$sd_zip))) ||
      truth$sigma <= 0 || truth$sd_state <= 0 || truth$sd_zip <= 0) {
    stop("Simulation truth is incompatible with the M2 design.", call. = FALSE)
  }

  u_raw_true <- stats::rnorm(jags_data$S, 0, truth$sd_state)
  u_true <- u_raw_true - mean(u_raw_true)
  v_raw_true <- stats::rnorm(jags_data$Z, 0, truth$sd_zip)
  v_state_mean <- vapply(seq_len(jags_data$S), function(state_id) {
    mean(v_raw_true[jags_data$zip_state == state_id])
  }, numeric(1L))
  v_true <- v_raw_true - v_state_mean[jags_data$zip_state]
  mu <- as.vector(
    truth$alpha +
      jags_data$X %*% truth$beta +
      u_true[jags_data$state] +
      v_true[jags_data$zip]
  )
  y <- stats::rnorm(jags_data$N, mu, truth$sigma)

  list(
    data = utils::modifyList(jags_data, list(y = y)),
    truth = scalar_truth,
    u_true = u_true,
    v_true = v_true
  )
}

# This is intentionally identical in structure to the production M2
# initialisation in scripts/05_fit_models.R, with separate seeds per recovery
# replication so all chains remain reproducible and dispersed.
make_m2_inits <- function(data_list, n_chains, seed) {
  set.seed(seed)
  lapply(seq_len(n_chains), function(chain_id) {
    list(
      alpha = stats::rnorm(1L, 0, 0.2),
      beta = stats::rnorm(data_list$P, 0, 0.1),
      sigma = stats::runif(1L, 0.5, 1.5),
      sd_state = stats::runif(1L, 0.1, 0.8),
      sd_zip = stats::runif(1L, 0.1, 0.8),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = as.integer(seed + 1000L * chain_id)
    )
  })
}

fit_m2_recovery <- function(data_list, settings, seed) {
  validate_jags_data(data_list)
  data_list <- augment_zip_centering_data(data_list)
  if (!file.exists(model_path)) {
    stop("M2 model file not found: ", model_path, call. = FALSE)
  }

  model <- rjags::jags.model(
    file = model_path,
    data = data_list,
    inits = make_m2_inits(data_list, settings$n_chains, seed),
    n.chains = settings$n_chains,
    n.adapt = settings$n_adapt,
    quiet = TRUE
  )
  stats::update(model, n.iter = settings$n_burn, progress.bar = "none")
  draws <- rjags::coda.samples(
    model = model,
    variable.names = c("alpha", "beta", "sigma", "sd_state", "sd_zip", "u", "v"),
    n.iter = settings$n_iter,
    thin = settings$thin,
    progress.bar = "none"
  )

  list(model = model, draws = draws)
}

draw_matrix <- function(draws) {
  values <- do.call(rbind, lapply(draws, as.matrix))
  if (nrow(values) == 0L || ncol(values) == 0L || any(!is.finite(values))) {
    stop("MCMC output contains no usable finite draws.", call. = FALSE)
  }
  values
}

diagnose_draws <- function(draws) {
  values <- draw_matrix(draws)
  draws_array <- posterior::as_draws_array(draws)
  diagnostics <- as.data.frame(posterior::summarise_draws(
    draws_array,
    posterior::rhat,
    posterior::ess_bulk,
    posterior::ess_tail,
    posterior::mcse_mean
  ))
  names(diagnostics) <- sub("^posterior::", "", names(diagnostics))
  diagnostics <- diagnostics[match(colnames(values), diagnostics$variable), , drop = FALSE]
  if (anyNA(diagnostics$variable)) {
    stop("Could not align posterior diagnostics with sampled parameters.", call. = FALSE)
  }

  diagnostics$converged <-
    is.finite(diagnostics$rhat) & diagnostics$rhat < 1.01 &
    is.finite(diagnostics$ess_bulk) & diagnostics$ess_bulk > 400 &
    is.finite(diagnostics$ess_tail) & diagnostics$ess_tail > 400
  diagnostics
}

summarise_scalar_recovery <- function(draws, truth, predictor_names) {
  values <- draw_matrix(draws)
  truth <- truth[intersect(names(truth), colnames(values))]
  if (length(truth) == 0L || !identical(names(truth), names(truth_vector(make_truth(length(predictor_names)))))) {
    stop("Scalar truth names do not match the M2 parameterisation.", call. = FALSE)
  }

  intervals <- t(vapply(names(truth), function(parameter) {
    stats::quantile(values[, parameter], probs = c(0.025, 0.975), names = FALSE)
  }, numeric(2L)))
  label <- names(truth)
  beta_rows <- grepl("^beta\\[", label)
  label[beta_rows] <- paste0(
    label[beta_rows], " (", predictor_names[as.integer(sub("^beta\\[([0-9]+)\\]$", "\\1", label[beta_rows]))], ")"
  )

  data.frame(
    parameter = names(truth),
    label = label,
    true_value = as.numeric(truth),
    posterior_mean = colMeans(values[, names(truth), drop = FALSE]),
    q2_5 = intervals[, 1L],
    q97_5 = intervals[, 2L],
    covers_truth = intervals[, 1L] <= truth & truth <= intervals[, 2L],
    bias = colMeans(values[, names(truth), drop = FALSE]) - as.numeric(truth),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

summarise_random_effects <- function(draws, simulated, jags_data, state_labels, zip_labels) {
  values <- draw_matrix(draws)
  effect_specs <- list(
    state = list(prefix = "u", truth = simulated$u_true, index = jags_data$state, labels = state_labels),
    zip = list(prefix = "v", truth = simulated$v_true, index = jags_data$zip, labels = zip_labels)
  )

  effect_results <- lapply(names(effect_specs), function(effect_name) {
    spec <- effect_specs[[effect_name]]
    parameters <- paste0(spec$prefix, "[", seq_along(spec$truth), "]")
    posterior_mean <- colMeans(values[, parameters, drop = FALSE])
    n_listings <- tabulate(spec$index, nbins = length(spec$truth))
    effect_data <- data.frame(
      effect = effect_name,
      group_id = seq_along(spec$truth),
      group_label = as.character(spec$labels),
      listings = n_listings,
      true_value = as.numeric(spec$truth),
      posterior_mean = posterior_mean,
      error = posterior_mean - as.numeric(spec$truth),
      stringsAsFactors = FALSE
    )
    metrics <- data.frame(
      effect = effect_name,
      n_groups = nrow(effect_data),
      correlation = stats::cor(effect_data$true_value, effect_data$posterior_mean),
      rmse = sqrt(mean(effect_data$error^2)),
      mean_bias = mean(effect_data$error),
      stringsAsFactors = FALSE
    )
    list(data = effect_data, metrics = metrics)
  })
  names(effect_results) <- names(effect_specs)

  list(
    effects = do.call(rbind, lapply(effect_results, `[[`, "data")),
    metrics = do.call(rbind, lapply(effect_results, `[[`, "metrics"))
  )
}

plot_random_effect_recovery <- function(effect_data, output_path) {
  max_abs <- max(abs(c(effect_data$true_value, effect_data$posterior_mean)))
  limit <- max(max_abs * 1.08, 0.05)
  point_size <- 0.7 + 1.8 * effect_data$listings / max(effect_data$listings)

  grDevices::png(output_path, width = 1500, height = 1200, res = 180)
  graphics::plot(
    effect_data$true_value,
    effect_data$posterior_mean,
    xlim = c(-limit, limit),
    ylim = c(-limit, limit),
    pch = 16,
    cex = point_size,
    col = grDevices::adjustcolor("steelblue4", alpha.f = 0.65),
    xlab = "Known simulated effect",
    ylab = "Posterior mean effect",
    main = paste0("M2 ", effect_data$effect[1L], "-effect parameter recovery")
  )
  graphics::abline(a = 0, b = 1, lty = 2, lwd = 2, col = "firebrick")
  graphics::mtext("Point area is proportional to listings in the group.", side = 3, line = 0.2, cex = 0.8)
  grDevices::dev.off()
}

full_recovery_paths <- function(tables_dir, figures_dir) {
  c(
    file.path(tables_dir, "simulation_recovery_full_scalar_parameters.csv"),
    file.path(tables_dir, "simulation_recovery_full_diagnostics.csv"),
    file.path(tables_dir, "simulation_recovery_full_random_effects.csv"),
    file.path(tables_dir, "simulation_recovery_full_random_effect_metrics.csv"),
    file.path(figures_dir, "simulation_recovery_full_state_effects.png"),
    file.path(figures_dir, "simulation_recovery_full_zip_effects.png")
  )
}

all_paths_exist <- function(paths) {
  all(file.exists(paths))
}

all_paths_current_for_model <- function(paths, model_path) {
  if (!all_paths_exist(paths) || !file.exists(model_path)) {
    return(FALSE)
  }
  model_mtime <- file.info(model_path)$mtime
  all(file.info(paths)$mtime >= model_mtime)
}

read_convergence_status <- function(path) {
  if (!file.exists(path)) {
    return(NA)
  }
  diagnostics <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!"converged" %in% names(diagnostics) || nrow(diagnostics) == 0L) {
    return(NA)
  }
  all(as.logical(diagnostics$converged))
}

repeated_status_path <- function(tables_dir, suffix = "") {
  file.path(tables_dir, paste0("simulation_recovery_repeated_fit_status", suffix, ".csv"))
}

repeated_metrics_path <- function(tables_dir, suffix = "") {
  file.path(tables_dir, paste0("simulation_recovery_repeated_parameter_metrics", suffix, ".csv"))
}

repeated_scalar_path <- function(tables_dir, replication, suffix = "") {
  file.path(
    tables_dir,
    sprintf(
      "simulation_recovery_repeated_scalar_replication_%02d%s.csv",
      replication,
      suffix
    )
  )
}

normalise_repeated_status <- function(status) {
  if (is.null(status) || nrow(status) == 0L) {
    return(status)
  }
  status$replication <- as.integer(status$replication)
  status$completed <- as.logical(status$completed)
  status$converged <- as.logical(status$converged)
  status
}

read_repeated_status <- function(path) {
  if (!file.exists(path)) {
    return(data.frame())
  }
  normalise_repeated_status(utils::read.csv(path, stringsAsFactors = FALSE))
}

write_repeated_status <- function(status_list, path) {
  available_status <- Filter(Negate(is.null), status_list)
  if (length(available_status) == 0L) {
    return(invisible(NULL))
  }
  status <- do.call(rbind, available_status)
  status <- status[order(status$replication), , drop = FALSE]
  utils::write.csv(status, path, row.names = FALSE)
  invisible(status)
}

existing_repeated_success <- function(existing_status, replication, scalar_path,
                                      model_path = NULL) {
  if (nrow(existing_status) == 0L || !file.exists(scalar_path)) {
    return(FALSE)
  }
  if (!is.null(model_path) && file.exists(model_path) &&
      file.info(scalar_path)$mtime < file.info(model_path)$mtime) {
    return(FALSE)
  }
  row <- existing_status[existing_status$replication == replication, , drop = FALSE]
  nrow(row) == 1L && isTRUE(row$completed) && isTRUE(row$converged)
}

empty_repeated_metrics <- function(truth) {
  data.frame(
    parameter = names(truth_vector(truth)),
    label = names(truth_vector(truth)),
    true_value = as.numeric(truth_vector(truth)),
    n_successful_converged = 0L,
    empirical_coverage_95 = NA_real_,
    posterior_mean_bias = NA_real_,
    posterior_mean_rmse = NA_real_,
    stringsAsFactors = FALSE
  )
}

# Build a fixed, balanced reduced design from the real M2 design.  Only the
# outcome is resimulated across replications, making the coverage experiment
# directly comparable from one synthetic data set to the next.
make_reduced_design <- function(jags_data, n_states = 5L, n_zips_per_state = 5L,
                                n_per_zip = 12L, seed = recovery_seed) {
  validate_jags_data(jags_data)
  zip_state <- vapply(seq_len(jags_data$Z), function(zip_id) {
    unique(jags_data$state[jags_data$zip == zip_id])
  }, integer(1L))
  zip_n <- tabulate(jags_data$zip, nbins = jags_data$Z)
  eligible_zips <- which(zip_n >= n_per_zip)
  eligible_states <- which(tabulate(zip_state[eligible_zips], nbins = jags_data$S) >= n_zips_per_state)
  if (length(eligible_states) < n_states) {
    stop(
      "The real design has only ", length(eligible_states), " states with at least ",
      n_zips_per_state, " ZIPs containing ", n_per_zip, " rows.",
      call. = FALSE
    )
  }

  set.seed(seed)
  selected_states <- sort(sample(eligible_states, n_states))
  selected_zips <- unlist(lapply(selected_states, function(state_id) {
    state_zips <- eligible_zips[zip_state[eligible_zips] == state_id]
    sort(sample(state_zips, n_zips_per_state))
  }), use.names = FALSE)
  selected_rows <- unlist(lapply(selected_zips, function(zip_id) {
    sample(which(jags_data$zip == zip_id), n_per_zip)
  }), use.names = FALSE)
  selected_rows <- sort(selected_rows)

  old_state <- jags_data$state[selected_rows]
  old_zip <- jags_data$zip[selected_rows]
  state_map <- match(old_state, selected_states)
  zip_map <- match(old_zip, selected_zips)
  reduced_data <- list(
    N = length(selected_rows),
    P = jags_data$P,
    S = length(selected_states),
    Z = length(selected_zips),
    y = as.numeric(jags_data$y[selected_rows]),
    X = unname(jags_data$X[selected_rows, , drop = FALSE]),
    state = as.integer(state_map),
    zip = as.integer(zip_map)
  )
  validate_jags_data(reduced_data)

  list(
    data = reduced_data,
    original_state_ids = selected_states,
    original_zip_ids = selected_zips,
    settings = list(
      n_states = n_states,
      n_zips_per_state = n_zips_per_state,
      n_per_zip = n_per_zip,
      seed = seed
    )
  )
}

run_one_recovery <- function(jags_data, truth, settings, seed, predictor_names,
                             state_labels, zip_labels) {
  started_at <- Sys.time()
  set.seed(seed)
  simulated <- simulate_m2_data(jags_data, truth)
  fit <- fit_m2_recovery(simulated$data, settings, seed + 500L)
  diagnostics <- diagnose_draws(fit$draws)
  scalar <- summarise_scalar_recovery(fit$draws, simulated$truth, predictor_names)
  random_effects <- summarise_random_effects(
    fit$draws, simulated, simulated$data, state_labels, zip_labels
  )

  list(
    scalar = scalar,
    diagnostics = diagnostics,
    random_effects = random_effects,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started_at, units = "secs")),
    converged = all(diagnostics$converged),
    max_rhat = max(diagnostics$rhat, na.rm = TRUE),
    min_bulk_ess = min(diagnostics$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(diagnostics$ess_tail, na.rm = TRUE)
  )
}

if (!file.exists(analysis_path)) {
  stop(
    "Prepared analysis objects not found: ", analysis_path,
    ". Run scripts/03_prepare_jags_data.R first.",
    call. = FALSE
  )
}
if (!file.exists(model_path)) {
  stop("M2 model file not found: ", model_path, call. = FALSE)
}

analysis_objects <- readRDS(analysis_path)
required_analysis_objects <- c("jags_train", "predictor_names", "state_levels", "zip_levels")
missing_analysis_objects <- setdiff(required_analysis_objects, names(analysis_objects))
if (length(missing_analysis_objects) > 0L) {
  stop(
    "Analysis objects are missing required item(s): ",
    paste(missing_analysis_objects, collapse = ", "),
    call. = FALSE
  )
}

full_data <- analysis_objects$jags_train
validate_jags_data(full_data)
truth <- make_truth(full_data$P)
if (length(analysis_objects$predictor_names) != full_data$P) {
  stop("Predictor names do not align with the M2 design matrix.", call. = FALSE)
}

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

settings_table <- data.frame(
  component = c("full_design", "repeated_reduced_design"),
  design_description = c(
    "Frozen M2 training design: all observations, states, and ZIPs.",
    "Balanced subset of the frozen M2 design; only outcomes vary across replications."
  ),
  n_chains = mcmc_settings$n_chains,
  n_adapt = mcmc_settings$n_adapt,
  n_burn = mcmc_settings$n_burn,
  n_iter = mcmc_settings$n_iter,
  thin = mcmc_settings$thin,
  n_replications = c(1L, n_replications),
  n_states = c(full_data$S, reduced_design_settings$n_states),
  n_zips_per_state = c(NA_integer_, reduced_design_settings$n_zips_per_state),
  n_per_zip = c(NA_integer_, reduced_design_settings$n_per_zip),
  reduced_design_label = c(NA_character_, reduced_design_label),
  run_full_recovery = run_full_recovery,
  run_repeated_recovery = run_repeated_recovery,
  force_full_recovery = force_full_recovery,
  force_repeated_recovery = force_repeated_recovery,
  resume_repeated_recovery = resume_repeated_recovery,
  stringsAsFactors = FALSE
)
utils::write.csv(
  settings_table,
  file.path(tables_dir, "simulation_recovery_settings.csv"),
  row.names = FALSE
)

full_result <- NULL
full_result_converged <- NA
if (run_full_recovery) {
  existing_full_paths <- full_recovery_paths(tables_dir, figures_dir)
  existing_full_complete <- all_paths_current_for_model(existing_full_paths, model_path)
  if (existing_full_complete && !force_full_recovery) {
    message(
      "Skipping full-design recovery because all expected outputs already exist. ",
      "Set RECOVERY_FORCE_FULL=true to refit it."
    )
    full_result_converged <- read_convergence_status(
      file.path(tables_dir, "simulation_recovery_full_diagnostics.csv")
    )
  } else {
    full_result <- run_one_recovery(
      jags_data = full_data,
      truth = truth,
      settings = mcmc_settings,
      seed = recovery_seed,
      predictor_names = analysis_objects$predictor_names,
      state_labels = analysis_objects$state_levels,
      zip_labels = analysis_objects$zip_levels
    )
    full_result$scalar$converged <- full_result$converged
    full_result_converged <- full_result$converged
    utils::write.csv(
      full_result$scalar,
      file.path(tables_dir, "simulation_recovery_full_scalar_parameters.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      full_result$diagnostics,
      file.path(tables_dir, "simulation_recovery_full_diagnostics.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      full_result$random_effects$effects,
      file.path(tables_dir, "simulation_recovery_full_random_effects.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      full_result$random_effects$metrics,
      file.path(tables_dir, "simulation_recovery_full_random_effect_metrics.csv"),
      row.names = FALSE
    )
    plot_random_effect_recovery(
      subset(full_result$random_effects$effects, effect == "state"),
      file.path(figures_dir, "simulation_recovery_full_state_effects.png")
    )
    plot_random_effect_recovery(
      subset(full_result$random_effects$effects, effect == "zip"),
      file.path(figures_dir, "simulation_recovery_full_zip_effects.png")
    )
  }
}

reduced_truth <- make_truth(full_data$P)
replication_status <- data.frame()
successful_replications <- integer(0L)
if (run_repeated_recovery) {
  reduced_design <- make_reduced_design(
    full_data,
    n_states = reduced_design_settings$n_states,
    n_zips_per_state = reduced_design_settings$n_zips_per_state,
    n_per_zip = reduced_design_settings$n_per_zip,
    seed = recovery_seed
  )
  reduced_truth <- make_truth(reduced_design$data$P)
  status_path <- repeated_status_path(tables_dir, reduced_output_suffix)
  existing_status <- if (resume_repeated_recovery && !force_repeated_recovery) {
    read_repeated_status(status_path)
  } else {
    data.frame()
  }
  replication_status_list <- vector("list", n_replications)
  if (nrow(existing_status) > 0L && !force_repeated_recovery) {
    for (replication in seq_len(n_replications)) {
      scalar_path <- repeated_scalar_path(tables_dir, replication, reduced_output_suffix)
      if (existing_repeated_success(existing_status, replication, scalar_path, model_path)) {
        replication_status_list[[replication]] <- existing_status[
          existing_status$replication == replication, ,
          drop = FALSE
        ]
      }
    }
  }

  for (replication in seq_len(n_replications)) {
    scalar_path <- repeated_scalar_path(tables_dir, replication, reduced_output_suffix)
    if (!force_repeated_recovery &&
        existing_repeated_success(existing_status, replication, scalar_path, model_path)) {
      replication_status_list[[replication]] <- existing_status[
        existing_status$replication == replication, ,
        drop = FALSE
      ]
      message(
        "Skipping reduced-design recovery replication ", replication, "/",
        n_replications, " because a converged checkpoint exists",
        if (nzchar(reduced_design_label)) paste0(" for ", reduced_design_label) else "",
        "."
      )
      write_repeated_status(replication_status_list, status_path)
      next
    }

    message(
      "Running reduced-design recovery replication ", replication, "/",
      n_replications,
      " (", reduced_design$data$S, " states, ", reduced_design$data$Z,
      " ZIPs, ", reduced_design$data$N, " rows)."
    )
    result <- tryCatch(
      run_one_recovery(
        jags_data = reduced_design$data,
        truth = reduced_truth,
        settings = mcmc_settings,
        seed = recovery_seed + 10000L * replication,
        predictor_names = analysis_objects$predictor_names,
        state_labels = analysis_objects$state_levels[reduced_design$original_state_ids],
        zip_labels = analysis_objects$zip_levels[reduced_design$original_zip_ids]
      ),
      error = function(error) error
    )

    if (inherits(result, "error")) {
      replication_status_list[[replication]] <- data.frame(
        replication = replication,
        completed = FALSE,
        converged = FALSE,
        elapsed_seconds = NA_real_,
        max_rhat = NA_real_,
        min_bulk_ess = NA_real_,
        min_tail_ess = NA_real_,
        error = conditionMessage(result),
        stringsAsFactors = FALSE
      )
    } else {
      scalar <- result$scalar
      scalar$replication <- replication
      utils::write.csv(scalar, scalar_path, row.names = FALSE)
      replication_status_list[[replication]] <- data.frame(
        replication = replication,
        completed = TRUE,
        converged = result$converged,
        elapsed_seconds = result$elapsed_seconds,
        max_rhat = result$max_rhat,
        min_bulk_ess = result$min_bulk_ess,
        min_tail_ess = result$min_tail_ess,
        error = NA_character_,
        stringsAsFactors = FALSE
      )
    }
    write_repeated_status(replication_status_list, status_path)
    status <- replication_status_list[[replication]]
    message(
      "Finished reduced-design recovery replication ", replication, "/",
      n_replications, ": completed=", status$completed,
      ", converged=", status$converged,
      ", elapsed_seconds=", round(status$elapsed_seconds, 1), "."
    )
  }
  replication_status <- write_repeated_status(replication_status_list, status_path)
  replication_status <- normalise_repeated_status(replication_status)

  successful_replications <- replication_status$replication[
    replication_status$completed & replication_status$converged
  ]
  successful_scalar_paths <- repeated_scalar_path(
    tables_dir,
    successful_replications,
    reduced_output_suffix
  )
  successful_replications <- successful_replications[file.exists(successful_scalar_paths)]
  if (length(successful_replications) > 0L) {
    repeated_scalar <- do.call(rbind, lapply(successful_replications, function(replication) {
      utils::read.csv(
        repeated_scalar_path(tables_dir, replication, reduced_output_suffix),
        stringsAsFactors = FALSE
      )
    }))
    repeated_metrics <- do.call(rbind, lapply(split(repeated_scalar, repeated_scalar$parameter), function(values) {
      data.frame(
        parameter = values$parameter[1L],
        label = values$label[1L],
        true_value = values$true_value[1L],
        n_successful_converged = nrow(values),
        empirical_coverage_95 = mean(values$covers_truth),
        posterior_mean_bias = mean(values$bias),
        posterior_mean_rmse = sqrt(mean(values$bias^2)),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    repeated_metrics <- empty_repeated_metrics(reduced_truth)
  }
  utils::write.csv(
    repeated_metrics,
    repeated_metrics_path(tables_dir, reduced_output_suffix),
    row.names = FALSE
  )
} else {
  message(
    "Skipping repeated reduced-design recovery because ",
    "RECOVERY_RUN_REPEATED=false."
  )
  if (file.exists(repeated_status_path(tables_dir, reduced_output_suffix))) {
    replication_status <- normalise_repeated_status(
      read_repeated_status(repeated_status_path(tables_dir, reduced_output_suffix))
    )
    successful_replications <- replication_status$replication[
      replication_status$completed & replication_status$converged
    ]
  }
}

if (!run_full_recovery) {
  warning(
    "The mandatory full-design recovery was not run. Set RECOVERY_RUN_FULL=true ",
    "for the final validation run before making any real-data claims.",
    call. = FALSE
  )
}
if (!is.na(full_result_converged) && !full_result_converged) {
  warning(
    "The full-design simulation recovery did not meet the pre-specified R-hat/ESS criteria. ",
    "Inspect simulation_recovery_full_diagnostics.csv before interpreting recovery estimates.",
    call. = FALSE
  )
}
if (run_repeated_recovery && length(successful_replications) < n_replications) {
  warning(
    length(successful_replications), " of ", n_replications,
    " reduced-design replications completed and met the pre-specified convergence criteria. ",
    "Inspect ",
    basename(repeated_status_path(tables_dir, reduced_output_suffix)),
    " before interpreting coverage.",
    call. = FALSE
  )
}

message(
  "Simulation recovery complete. ",
  if (run_full_recovery) "Full-design recovery was checked. " else "Full-design recovery was not run; set RECOVERY_RUN_FULL=true for final validation. ",
  if (run_repeated_recovery) {
    paste0(
      length(successful_replications), "/", n_replications,
      " reduced-design replications completed and converged",
      if (nzchar(reduced_design_label)) paste0(" for ", reduced_design_label) else "",
      ". "
    )
  } else {
    "Repeated reduced-design recovery was skipped. "
  },
  "Results are under output/tables and output/figures."
)
