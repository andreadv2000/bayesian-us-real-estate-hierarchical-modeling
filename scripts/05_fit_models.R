# Fit the M0, M1, and M2 models on the frozen training sample.
#
# Each model uses the identical `jags_train` object.  The only difference
# between fits is the geographic structure specified by the JAGS model file.

required_fit_packages <- c("coda", "here", "rjags")
missing_fit_packages <- required_fit_packages[!vapply(
  required_fit_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_fit_packages) > 0L) {
  stop(
    "Install the package(s) required for model fitting: ",
    paste(missing_fit_packages, collapse = ", "),
    call. = FALSE
  )
}

source(here::here("scripts", "00_setup.R"))

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
  if (any(c(n_chains, n_adapt, n_burn, n_iter, thin) < 1L) ||
      any(c(n_chains, n_adapt, n_burn, n_iter, thin) !=
        as.integer(c(n_chains, n_adapt, n_burn, n_iter, thin)))) {
    stop("MCMC settings must be positive integers.", call. = FALSE)
  }

  required_data <- c("N", "P", "y", "X")
  missing_data <- setdiff(required_data, names(data_list))
  if (length(missing_data) > 0L) {
    stop(
      "JAGS data are missing required item(s): ",
      paste(missing_data, collapse = ", "),
      call. = FALSE
    )
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

analysis_path <- here::here("data", "processed", "analysis_objects.rds")
if (!file.exists(analysis_path)) {
  stop(
    "Prepared analysis objects not found: ", analysis_path,
    ". Run scripts/03_prepare_jags_data.R first.",
    call. = FALSE
  )
}

analysis_objects <- readRDS(analysis_path)
if (!"jags_train" %in% names(analysis_objects)) {
  stop("Analysis objects do not contain the frozen JAGS training data.", call. = FALSE)
}

parameters_m0 <- c("alpha", "beta", "sigma")
parameters_m1 <- c("alpha", "beta", "sigma", "sd_state", "u")
parameters_m2 <- c("alpha", "beta", "sigma", "sd_state", "sd_zip", "u", "v")

fit_m0 <- fit_jags_model(
  model_file = here::here("models", "m0_pooled.jags"),
  data_list = analysis_objects$jags_train,
  parameters = parameters_m0
)
fit_m1 <- fit_jags_model(
  model_file = here::here("models", "m1_state.jags"),
  data_list = analysis_objects$jags_train,
  parameters = parameters_m1
)
fit_m2 <- fit_jags_model(
  model_file = here::here("models", "m2_state_zip.jags"),
  data_list = analysis_objects$jags_train,
  parameters = parameters_m2
)

mcmc_output_dir <- here::here("output", "mcmc")
dir.create(mcmc_output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(fit_m0$draws, file.path(mcmc_output_dir, "m0_draws.rds"))
saveRDS(fit_m1$draws, file.path(mcmc_output_dir, "m1_draws.rds"))
saveRDS(fit_m2$draws, file.path(mcmc_output_dir, "m2_draws.rds"))

source(here::here("scripts", "07_model_comparison.R"))
dic_results <- calculate_and_write_dic(fit_m0, fit_m1, fit_m2)

message(
  "Saved four-chain production draws to ", mcmc_output_dir,
  " and the DIC comparison to ", here::here("output", "tables", "model_comparison_dic.csv"),
  "."
)
