# Prepare the frozen balanced sample for the JAGS models.
#
# Scaling parameters, factor levels, and geographic indices are learned from
# training observations only.  The resulting object is the single source of
# truth for the later fitting, diagnostics, and prediction scripts.

required_preparation_packages <- c("arrow", "dplyr", "forcats", "here")
missing_preparation_packages <- required_preparation_packages[!vapply(
  required_preparation_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_preparation_packages) > 0) {
  stop(
    "Install the package(s) required for JAGS data preparation: ",
    paste(missing_preparation_packages, collapse = ", ")
  )
}

source(here::here("scripts", "00_setup.R"))

sample_path <- here::here("data", "processed", "analysis_sample.parquet")
if (!file.exists(sample_path)) {
  stop(
    "Frozen analysis sample not found: ", sample_path,
    ". Run scripts/02_make_sample.R first."
  )
}

analysis_sample <- arrow::read_parquet(sample_path)

required_columns <- c(
  "split", "log_price", "log_house_size", "log_acre_lot", "bed", "bath",
  "status", "state", "zip_id"
)
missing_columns <- setdiff(required_columns, names(analysis_sample))
if (length(missing_columns) > 0) {
  stop(
    "Analysis sample is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

for (variable in c("status", "state", "zip_id")) {
  values <- as.character(analysis_sample[[variable]])
  if (anyNA(values) || any(!nzchar(trimws(values)))) {
    stop("`", variable, "` must be non-missing and non-blank in the frozen sample.")
  }
}

unexpected_splits <- setdiff(unique(analysis_sample$split), c("train", "test"))
if (length(unexpected_splits) > 0 || anyNA(analysis_sample$split)) {
  stop("`split` must contain only non-missing 'train' and 'test' labels.")
}

train <- analysis_sample |>
  dplyr::filter(split == "train") |>
  dplyr::mutate(status = forcats::fct_drop(factor(status)))
test <- analysis_sample |>
  dplyr::filter(split == "test")

if (nrow(train) == 0 || nrow(test) == 0) {
  stop("The frozen sample must contain both training and test observations.")
}

# An unseen categorical level has no fitted coefficient.  Stop rather than
# silently convert it to NA and drop a held-out observation from evaluation.
unseen_test_statuses <- setdiff(
  unique(as.character(test$status)),
  levels(train$status)
)
if (length(unseen_test_statuses) > 0 || anyNA(test$status)) {
  stop(
    "Test data contain status level(s) absent from training: ",
    paste(unseen_test_statuses, collapse = ", "),
    ". Recreate the frozen sample with each status represented in training."
  )
}
test$status <- factor(as.character(test$status), levels = levels(train$status))

continuous_vars <- c("log_house_size", "log_acre_lot", "bed", "bath")
variables_to_scale <- c("log_price", continuous_vars)
for (variable in variables_to_scale) {
  if (!is.numeric(train[[variable]]) || !is.numeric(test[[variable]]) ||
      any(!is.finite(train[[variable]])) || any(!is.finite(test[[variable]]))) {
    stop("`", variable, "` must be finite and numeric in both data splits.")
  }
}

# Compute centring/scaling constants once from training data and apply those
# same constants to held-out data.  Test observations never affect a prior or
# a fitted coefficient through this preprocessing step.
scaling_values <- train |>
  dplyr::summarise(dplyr::across(
    dplyr::all_of(variables_to_scale),
    list(mean = mean, sd = stats::sd)
  ))

scaling_sd_names <- paste0(variables_to_scale, "_sd")
if (any(!is.finite(unlist(scaling_values[scaling_sd_names]))) ||
    any(unlist(scaling_values[scaling_sd_names]) <= 0)) {
  stop("Every continuous training variable must have a finite, positive SD.")
}

scale_from_train <- function(data, variable, centre, scale) {
  (data[[variable]] - centre) / scale
}

for (variable in continuous_vars) {
  variable_mean <- scaling_values[[paste0(variable, "_mean")]]
  variable_sd <- scaling_values[[paste0(variable, "_sd")]]
  standardised_name <- paste0(variable, "_z")

  train[[standardised_name]] <- scale_from_train(
    train, variable, variable_mean, variable_sd
  )
  test[[standardised_name]] <- scale_from_train(
    test, variable, variable_mean, variable_sd
  )
}

y_mean <- scaling_values$log_price_mean
y_sd <- scaling_values$log_price_sd
train$y <- scale_from_train(train, "log_price", y_mean, y_sd)
test$y <- scale_from_train(test, "log_price", y_mean, y_sd)

# Use one terms object to guarantee that every model receives the same columns
# for training and test data.  The intercept belongs to alpha in the JAGS
# models, so it is deliberately removed from X.
predictor_formula <- ~
  log_house_size_z +
  log_acre_lot_z +
  bed_z +
  bath_z +
  status

terms_x <- stats::terms(predictor_formula, data = train)
X_train_full <- stats::model.matrix(terms_x, data = train)
X_test_full <- stats::model.matrix(terms_x, data = test)

X_train <- X_train_full[, colnames(X_train_full) != "(Intercept)", drop = FALSE]
X_test <- X_test_full[, colnames(X_test_full) != "(Intercept)", drop = FALSE]
if (!identical(colnames(X_train), colnames(X_test))) {
  stop("Training and test design matrices do not have identical columns.")
}
P <- ncol(X_train)
if (P == 0) {
  stop("The predictor design matrix has no non-intercept columns.")
}

# JAGS receives one-based integer indices.  State and ZIP levels come only
# from training so prediction is explicitly limited to observed locations.
state_levels <- sort(unique(train$state))
zip_levels <- sort(unique(train$zip_id))

train$state_id <- match(train$state, state_levels)
train$zip_id_num <- match(train$zip_id, zip_levels)
test$state_id <- match(test$state, state_levels)
test$zip_id_num <- match(test$zip_id, zip_levels)

if (anyNA(train$state_id) || anyNA(train$zip_id_num) ||
    anyNA(test$state_id) || anyNA(test$zip_id_num)) {
  stop(
    "Every test state and ZIP must be present in training. ",
    "Recreate the frozen sample with training observations per ZIP."
  )
}

zip_lookup <- train |>
  dplyr::distinct(zip_id, state) |>
  dplyr::mutate(
    zip_id_num = match(zip_id, zip_levels),
    state_id = match(state, state_levels)
  ) |>
  dplyr::arrange(zip_id_num)

if (nrow(zip_lookup) != length(zip_levels) || anyDuplicated(zip_lookup$zip_id_num)) {
  stop("Each training ZIP must map to exactly one state.")
}

jags_train <- list(
  N = nrow(train),
  P = P,
  S = length(state_levels),
  Z = length(zip_levels),
  y = train$y,
  X = X_train,
  state = train$state_id,
  zip = train$zip_id_num
)

analysis_objects <- list(
  train = train,
  test = test,
  X_train = X_train,
  X_test = X_test,
  jags_train = jags_train,
  state_levels = state_levels,
  zip_levels = zip_levels,
  zip_lookup = zip_lookup,
  y_mean = y_mean,
  y_sd = y_sd,
  predictor_names = colnames(X_train),
  sampling_seed = 20260621
)

stopifnot(
  jags_train$N == nrow(train),
  nrow(jags_train$X) == jags_train$N,
  ncol(jags_train$X) == jags_train$P,
  length(jags_train$y) == jags_train$N,
  length(jags_train$state) == jags_train$N,
  length(jags_train$zip) == jags_train$N
)

output_path <- here::here("data", "processed", "analysis_objects.rds")
saveRDS(analysis_objects, output_path)

message(
  "JAGS data preparation complete: ",
  jags_train$N, " training rows, ", nrow(test), " test rows, ",
  jags_train$P, " predictors, ", jags_train$S, " states, and ",
  jags_train$Z, " ZIP codes. Saved: ", output_path
)
