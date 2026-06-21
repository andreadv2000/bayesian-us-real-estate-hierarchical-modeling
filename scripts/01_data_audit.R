# Audit the preserved raw listings without applying modelling filters.
source(here::here("scripts", "00_setup.R"))

raw_path <- here::here("data", "raw", "realtor-data.csv")
tables_dir <- here::here("output", "tables")
figures_dir <- here::here("output", "figures")

stopifnot(file.exists(raw_path))
fs::dir_create(c(tables_dir, figures_dir), recurse = TRUE)

raw <- data.table::fread(
  raw_path,
  na.strings = c("", "NA", "N/A", "null", "NULL")
)
raw <- janitor::clean_names(raw)

required_columns <- c(
  "status", "price", "bed", "bath", "acre_lot", "city", "state",
  "zip_code", "house_size"
)
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns) > 0) {
  stop("Raw data are missing expected columns: ", paste(missing_columns, collapse = ", "))
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else median(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

units <- c(
  brokered_by = "broker identifier",
  status = "listing-status category",
  price = "USD",
  bed = "bedrooms",
  bath = "bathrooms",
  acre_lot = "acres",
  street = "street address",
  city = "city name",
  state = "state name",
  zip_code = "ZIP-code label",
  house_size = "square feet",
  prev_sold_date = "calendar date"
)

# Basic import record: this gives a quick, reproducible check of data coverage.
dataset_overview <- tibble::tibble(
  metric = c("rows", "columns", "unique_states", "unique_cities", "unique_zip_codes"),
  value = c(
    nrow(raw),
    ncol(raw),
    dplyr::n_distinct(raw$state, na.rm = TRUE),
    dplyr::n_distinct(raw$city, na.rm = TRUE),
    dplyr::n_distinct(raw$zip_code, na.rm = TRUE)
  )
)
readr::write_csv(dataset_overview, fs::path(tables_dir, "raw_dataset_overview.csv"))

missingness <- tibble::tibble(
  variable = names(raw),
  missing_n = vapply(raw, function(x) sum(is.na(x)), numeric(1)),
  missing_prop = vapply(raw, function(x) mean(is.na(x)), numeric(1))
) |>
  dplyr::arrange(dplyr::desc(missing_prop), variable)
readr::write_csv(missingness, fs::path(tables_dir, "raw_missingness.csv"))

# Describe every raw field. Numerical summaries are intentionally blank for
# non-numerical fields so identifiers are not accidentally interpreted as values.
data_dictionary <- tibble::tibble(
  variable = names(raw),
  type = vapply(raw, function(x) paste(class(x), collapse = "/"), character(1)),
  units = unname(units[names(raw)]),
  missing_n = vapply(raw, function(x) sum(is.na(x)), numeric(1)),
  missing_prop = vapply(raw, function(x) mean(is.na(x)), numeric(1)),
  distinct_n = vapply(raw, function(x) dplyr::n_distinct(x, na.rm = TRUE), numeric(1)),
  min = vapply(raw, function(x) if (is.numeric(x)) safe_min(x) else NA_real_, numeric(1)),
  median = vapply(raw, function(x) if (is.numeric(x)) safe_median(x) else NA_real_, numeric(1)),
  max = vapply(raw, function(x) if (is.numeric(x)) safe_max(x) else NA_real_, numeric(1))
)
readr::write_csv(data_dictionary, fs::path(tables_dir, "raw_data_dictionary.csv"))

status_counts <- raw |>
  dplyr::count(status, name = "listings", sort = TRUE) |>
  dplyr::mutate(status = as.character(status))
readr::write_csv(status_counts, fs::path(tables_dir, "raw_status_counts.csv"))

audit_variables <- c("price", "bed", "bath", "house_size", "acre_lot")
value_quality_counts <- purrr::map_dfr(audit_variables, function(variable) {
  values <- raw[[variable]]
  tibble::tibble(
    variable = variable,
    total_n = length(values),
    missing_n = sum(is.na(values)),
    zero_n = sum(values == 0, na.rm = TRUE),
    negative_n = sum(values < 0, na.rm = TRUE),
    non_positive_n = sum(values <= 0, na.rm = TRUE)
  )
})
readr::write_csv(value_quality_counts, fs::path(tables_dir, "raw_value_quality_counts.csv"))

# Freeze the modelling decisions that the raw-data audit can support. These do
# not alter `raw`; the filtering itself belongs to the next implementation step.
quality_value <- function(variable, column) {
  value_quality_counts[[column]][match(variable, value_quality_counts$variable)]
}

exact_duplicate_rows <- sum(duplicated(raw))
audit_decisions <- tibble::tibble(
  topic = c(
    "listing_status",
    "lot_size",
    "exact_duplicates",
    "impossible_price_and_house_size",
    "high_price_outliers"
  ),
  decision = c(
    "Use status as a categorical predictor; do not restrict the audit to one listing status.",
    "Retain acre_lot as a core covariate using log1p on observed nonnegative values; use a complete-case analysis rather than median imputation.",
    "Retain all records: the audit found no exact duplicate rows.",
    "Exclude non-positive price and non-positive house_size when constructing the modelling frame.",
    "Do not remove expensive homes solely for being extreme; model price on the log scale."
  ),
  audit_evidence = c(
    sprintf("%d observed status levels are present and status is a planned covariate.", nrow(status_counts)),
    sprintf(
      "acre_lot is a planned covariate; %.2f%% of raw values are missing and %s values are zero, which log1p can represent.",
      100 * quality_value("acre_lot", "missing_n") / nrow(raw),
      format(quality_value("acre_lot", "zero_n"), big.mark = ",")
    ),
    sprintf(
      "%s exact duplicates among %s raw rows.",
      format(exact_duplicate_rows, big.mark = ","),
      format(nrow(raw), big.mark = ",")
    ),
    sprintf(
      "%s prices are zero and %s are negative; %s house-size values are non-positive.",
      format(quality_value("price", "zero_n"), big.mark = ","),
      format(quality_value("price", "negative_n"), big.mark = ","),
      format(quality_value("house_size", "non_positive_n"), big.mark = ",")
    ),
    "The audit records distributions without defining an upper price cutoff."
  )
)
readr::write_csv(audit_decisions, fs::path(tables_dir, "raw_data_decisions.csv"))

state_listing_counts <- raw |>
  dplyr::count(state, name = "listings", sort = TRUE) |>
  dplyr::mutate(state = as.character(state))
readr::write_csv(state_listing_counts, fs::path(tables_dir, "raw_listings_per_state.csv"))

# Pair state and ZIP because the ZIP itself is a categorical label and may not
# be globally unique in imperfect listing data.
zip_listing_counts <- raw |>
  dplyr::mutate(
    state = as.character(state),
    zip_code = as.character(zip_code)
  ) |>
  dplyr::count(state, zip_code, name = "listings", sort = TRUE)
readr::write_csv(zip_listing_counts, fs::path(tables_dir, "raw_listings_per_state_zip.csv"))

duplicate_summary <- tibble::tibble(
  total_rows = nrow(raw),
  exact_duplicate_rows = exact_duplicate_rows,
  unique_rows = nrow(raw) - exact_duplicate_rows
)
readr::write_csv(duplicate_summary, fs::path(tables_dir, "raw_exact_duplicate_summary.csv"))

save_histogram <- function(data, x, x_label, file_name) {
  # Retain all finite raw values, including zero and negative values, for the
  # audit. Non-finite values cannot be represented on a continuous histogram.
  plot_data <- dplyr::filter(data, is.finite(.data[[x]]))
  if (nrow(plot_data) == 0) {
    stop("Cannot create histogram for ", x_label, ": no finite values found.")
  }

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x]])) +
    ggplot2::geom_histogram(bins = 60, colour = "white", fill = "#2C7FB8") +
    ggplot2::labs(
      title = paste("Raw-data distribution:", x_label),
      x = x_label,
      y = "Listings"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  ggplot2::ggsave(
    filename = fs::path(figures_dir, file_name),
    plot = plot,
    width = 8,
    height = 5,
    dpi = 300
  )
}

save_histogram(raw, "price", "Listing price (USD)", "raw_price_histogram.png")
save_histogram(
  dplyr::filter(raw, price > 0) |>
    dplyr::transmute(log_price = log(price)),
  "log_price",
  "Log listing price",
  "raw_log_price_histogram.png"
)
save_histogram(raw, "house_size", "House size (square feet)", "raw_house_size_histogram.png")
save_histogram(
  dplyr::filter(raw, house_size > 0) |>
    dplyr::transmute(log_house_size = log(house_size)),
  "log_house_size",
  "Log house size",
  "raw_log_house_size_histogram.png"
)
save_histogram(raw, "acre_lot", "Lot size (acres)", "raw_acre_lot_histogram.png")

zip_sample_size_plot <- ggplot2::ggplot(
  zip_listing_counts,
  ggplot2::aes(x = listings)
) +
  ggplot2::geom_histogram(bins = 60, colour = "white", fill = "#41AB5D") +
  ggplot2::labs(
    title = "Raw listings per state-ZIP group",
    x = "Listings per state-ZIP group",
    y = "Number of state-ZIP groups"
  ) +
  ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(
  filename = fs::path(figures_dir, "raw_zip_sample_size_histogram.png"),
  plot = zip_sample_size_plot,
  width = 8,
  height = 5,
  dpi = 300
)

message("Raw-data audit complete. Tables: ", tables_dir, "; figures: ", figures_dir)
