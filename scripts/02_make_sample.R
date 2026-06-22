# Build the clean modelling frame used by the later sampling and JAGS scripts.
# This script intentionally retains only model-relevant variables; broker and
# street are high-cardinality identifiers that do not answer the location
# question.

required_cleaning_packages <- c(
  "data.table", "janitor", "dplyr", "stringr", "forcats",
  "tibble", "readr", "here"
)
missing_cleaning_packages <- required_cleaning_packages[!vapply(
  required_cleaning_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_cleaning_packages) > 0) {
  stop(
    "Install the package(s) required for the clean modelling frame: ",
    paste(missing_cleaning_packages, collapse = ", ")
  )
}

source(here::here("scripts", "00_setup.R"))

raw_path <- here::here("data", "raw", "realtor-data.csv")
if (!file.exists(raw_path)) {
  stop("Raw data file not found: ", raw_path)
}

raw <- data.table::fread(
  raw_path,
  na.strings = c("", "NA", "N/A", "null", "NULL")
) |>
  janitor::clean_names()

required_core <- c(
  "price", "bed", "bath", "house_size",
  "state", "zip_code", "status"
)
stopifnot(all(required_core %in% names(raw)))

cleaned <- raw |>
  dplyr::transmute(
    price = as.numeric(price),
    bed = as.numeric(bed),
    bath = as.numeric(bath),
    house_size = as.numeric(house_size),
    acre_lot = as.numeric(acre_lot),
    status = as.character(status),
    state = stringr::str_squish(as.character(state)),
    city = stringr::str_squish(as.character(city)),
    zip_code_raw = stringr::str_squish(as.character(zip_code)),
    prev_sold_date = as.character(prev_sold_date)
  ) |>
  dplyr::filter(
    is.finite(price), price > 0,
    is.finite(house_size), house_size > 0,
    is.finite(bed), bed > 0,
    is.finite(bath), bath > 0,
    !is.na(state), state != "",
    !is.na(zip_code_raw), zip_code_raw != "",
    !is.na(status), status != ""
  ) |>
  dplyr::mutate(
    zip_code_label = stringr::str_pad(
      stringr::str_remove(zip_code_raw, "\\.0$"),
      width = 5, side = "left", pad = "0"
    ),
    zip_id = paste(state, zip_code_label, sep = "__"),
    log_price = log(price),
    log_house_size = log(house_size),
    log_acre_lot = dplyr::if_else(
      !is.na(acre_lot) & acre_lot >= 0,
      log1p(acre_lot),
      NA_real_
    )
  )

# Lot size is retained in the core models.  The recorded missingness lets the
# next implementation step assess whether this complete-case choice preserves
# enough state and ZIP coverage for its frozen MCMC sample.
lot_size_missing_fraction <- mean(is.na(cleaned$log_acre_lot))
analysis_frame <- cleaned |>
  dplyr::filter(!is.na(log_acre_lot)) |>
  dplyr::mutate(status = forcats::fct_infreq(factor(status)))

row_accounting <- tibble::tibble(
  stage = c("Raw rows", "Core valid rows", "Lot-size complete rows"),
  n = c(nrow(raw), nrow(cleaned), nrow(analysis_frame))
)

dir.create(here::here("output", "tables"), recursive = TRUE,
           showWarnings = FALSE)
readr::write_csv(
  row_accounting,
  here::here("output", "tables", "row_accounting.csv")
)

print(row_accounting)
message(
  "Missing log-acre-lot fraction among core-valid rows: ",
  format(lot_size_missing_fraction, digits = 4)
)

# Freeze a balanced, computationally tractable MCMC sample.  The state--ZIP
# hierarchy is sampled first so every retained ZIP contributes the same number
# of listings and has observations in both the training and test partitions.
required_sampling_packages <- "arrow"
missing_sampling_packages <- required_sampling_packages[!vapply(
  required_sampling_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_sampling_packages) > 0) {
  stop(
    "Install the package(s) required for stratified sampling: ",
    paste(missing_sampling_packages, collapse = ", ")
  )
}

set.seed(20260621)

N_STATES <- 25L
N_ZIPS_PER_STATE <- 24L
N_PER_ZIP <- 30L
N_TRAIN_PER_ZIP <- 24L

eligible_zips <- analysis_frame |>
  dplyr::count(state, zip_id, name = "n_available") |>
  dplyr::filter(n_available >= N_PER_ZIP)

eligible_states <- eligible_zips |>
  dplyr::count(state, name = "n_eligible_zips") |>
  dplyr::filter(n_eligible_zips >= N_ZIPS_PER_STATE)

if (nrow(eligible_states) < N_STATES) {
  stop(
    "The modelling frame has only ", nrow(eligible_states),
    " states with at least ", N_ZIPS_PER_STATE,
    " ZIP codes containing ", N_PER_ZIP, " listings; ",
    N_STATES, " are required."
  )
}

selected_states <- eligible_states |>
  dplyr::slice_sample(n = N_STATES) |>
  dplyr::select(state)

selected_zips <- eligible_zips |>
  dplyr::semi_join(selected_states, by = "state") |>
  dplyr::group_by(state) |>
  dplyr::slice_sample(n = N_ZIPS_PER_STATE) |>
  dplyr::ungroup()

analysis_sample <- analysis_frame |>
  dplyr::semi_join(selected_zips, by = c("state", "zip_id")) |>
  dplyr::group_by(state, zip_id) |>
  dplyr::slice_sample(n = N_PER_ZIP) |>
  dplyr::mutate(
    split = sample(
      c(
        rep("train", N_TRAIN_PER_ZIP),
        rep("test", N_PER_ZIP - N_TRAIN_PER_ZIP)
      ),
      size = N_PER_ZIP,
      replace = FALSE
    )
  ) |>
  dplyr::ungroup()

sample_manifest <- analysis_sample |>
  dplyr::count(state, zip_id, split, name = "n") |>
  dplyr::arrange(state, zip_id, split)

train_per_zip <- analysis_sample |>
  dplyr::filter(split == "train") |>
  dplyr::count(state, zip_id, name = "n")
test_per_zip <- analysis_sample |>
  dplyr::filter(split == "test") |>
  dplyr::count(state, zip_id, name = "n")

stopifnot(
  nrow(analysis_sample) == N_STATES * N_ZIPS_PER_STATE * N_PER_ZIP,
  dplyr::n_distinct(analysis_sample$state) == N_STATES,
  dplyr::n_distinct(analysis_sample$zip_id) == N_STATES * N_ZIPS_PER_STATE,
  all(train_per_zip$n == N_TRAIN_PER_ZIP),
  all(test_per_zip$n == N_PER_ZIP - N_TRAIN_PER_ZIP)
)

processed_dir <- here::here("data", "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  sample_manifest,
  file.path(processed_dir, "sample_manifest.csv")
)
arrow::write_parquet(
  analysis_sample,
  file.path(processed_dir, "analysis_sample.parquet")
)

sample_validation <- list(
  n_rows = nrow(analysis_sample),
  n_train = sum(analysis_sample$split == "train"),
  n_test = sum(analysis_sample$split == "test"),
  n_states = dplyr::n_distinct(analysis_sample$state),
  n_zips = dplyr::n_distinct(analysis_sample$zip_id),
  train_per_zip = sort(unique(train_per_zip$n)),
  test_per_zip = sort(unique(test_per_zip$n))
)
print(sample_validation)
message(
  "Balanced MCMC sample complete. Manifest and Parquet data: ", processed_dir
)
