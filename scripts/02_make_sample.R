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
