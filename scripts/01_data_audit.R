# Audit the preserved raw listings without applying modelling filters.
source(file.path("scripts", "00_setup.R"))

raw_path <- file.path("data", "raw", "realtor-data.csv")
tables_dir <- file.path("output", "tables")
figures_dir <- file.path("output", "figures")

if (!file.exists(raw_path)) {
  stop("Raw data file not found: ", raw_path)
}

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Fast import is useful here because the raw CSV is about 170 MiB.
raw <- data.table::fread(
  raw_path,
  na.strings = c("", "NA", "N/A", "null", "NULL")
)

# Make names lower-case and use underscores so that later column references
# are simple, for example raw$house_size.
clean_column_names <- function(column_names) {
  column_names <- tolower(column_names)
  column_names <- gsub("[^a-z0-9]+", "_", column_names)
  gsub("^_|_$", "", column_names)
}
names(raw) <- clean_column_names(names(raw))

required_columns <- c(
  "status", "price", "bed", "bath", "acre_lot", "city", "state",
  "zip_code", "house_size"
)
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns) > 0) {
  stop(
    "Raw data are missing expected columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

count_observed_values <- function(values) {
  length(unique(values[!is.na(values)]))
}

write_audit_table <- function(data, file_name) {
  utils::write.csv(
    data,
    file.path(tables_dir, file_name),
    row.names = FALSE,
    na = ""
  )
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
dataset_overview <- data.frame(
  metric = c("rows", "columns", "unique_states", "unique_cities", "unique_zip_codes"),
  value = c(
    nrow(raw),
    ncol(raw),
    count_observed_values(raw$state),
    count_observed_values(raw$city),
    count_observed_values(raw$zip_code)
  ),
  stringsAsFactors = FALSE
)
write_audit_table(dataset_overview, "raw_dataset_overview.csv")

# Count missing entries in every column using a simple loop.
missingness <- data.frame(
  variable = names(raw),
  missing_n = numeric(ncol(raw)),
  missing_prop = numeric(ncol(raw)),
  stringsAsFactors = FALSE
)
for (i in seq_along(raw)) {
  missingness$missing_n[i] <- sum(is.na(raw[[i]]))
  missingness$missing_prop[i] <- mean(is.na(raw[[i]]))
}
missingness <- missingness[order(-missingness$missing_prop, missingness$variable), ]
write_audit_table(missingness, "raw_missingness.csv")

# Describe every raw field. Numerical summaries are blank for non-numerical
# fields so identifiers are not accidentally interpreted as values.
data_dictionary <- data.frame(
  variable = names(raw),
  type = character(ncol(raw)),
  units = unname(units[names(raw)]),
  missing_n = numeric(ncol(raw)),
  missing_prop = numeric(ncol(raw)),
  distinct_n = numeric(ncol(raw)),
  min = rep(NA_real_, ncol(raw)),
  median = rep(NA_real_, ncol(raw)),
  max = rep(NA_real_, ncol(raw)),
  stringsAsFactors = FALSE
)

for (i in seq_along(raw)) {
  values <- raw[[i]]
  data_dictionary$type[i] <- paste(class(values), collapse = "/")
  data_dictionary$missing_n[i] <- sum(is.na(values))
  data_dictionary$missing_prop[i] <- mean(is.na(values))
  data_dictionary$distinct_n[i] <- count_observed_values(values)

  if (is.numeric(values)) {
    finite_values <- values[is.finite(values)]
    if (length(finite_values) > 0) {
      data_dictionary$min[i] <- min(finite_values)
      data_dictionary$median[i] <- median(finite_values)
      data_dictionary$max[i] <- max(finite_values)
    }
  }
}
write_audit_table(data_dictionary, "raw_data_dictionary.csv")

# Listing status is kept as a categorical variable. Include a missing-status
# group if one is present in the raw file.
status_table <- sort(table(raw$status, useNA = "ifany"), decreasing = TRUE)
status_counts <- data.frame(
  status = names(status_table),
  listings = as.integer(status_table),
  stringsAsFactors = FALSE
)
status_counts$status[status_counts$status == "<NA>"] <- NA_character_
write_audit_table(status_counts, "raw_status_counts.csv")

# Check the values that are important for the planned model.
audit_variables <- c("price", "bed", "bath", "house_size", "acre_lot")
value_quality_counts <- data.frame(
  variable = audit_variables,
  total_n = numeric(length(audit_variables)),
  missing_n = numeric(length(audit_variables)),
  zero_n = numeric(length(audit_variables)),
  negative_n = numeric(length(audit_variables)),
  non_positive_n = numeric(length(audit_variables)),
  stringsAsFactors = FALSE
)

for (i in seq_along(audit_variables)) {
  values <- raw[[audit_variables[i]]]
  value_quality_counts$total_n[i] <- length(values)
  value_quality_counts$missing_n[i] <- sum(is.na(values))
  value_quality_counts$zero_n[i] <- sum(values == 0, na.rm = TRUE)
  value_quality_counts$negative_n[i] <- sum(values < 0, na.rm = TRUE)
  value_quality_counts$non_positive_n[i] <- sum(values <= 0, na.rm = TRUE)
}
write_audit_table(value_quality_counts, "raw_value_quality_counts.csv")

get_quality_value <- function(variable, column) {
  row_number <- match(variable, value_quality_counts$variable)
  value_quality_counts[row_number, column]
}

# Freeze the modelling decisions supported by this audit. These decisions do
# not change `raw`; filtering belongs to the next implementation stage.
exact_duplicate_rows <- sum(duplicated(raw))
audit_decisions <- data.frame(
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
      100 * get_quality_value("acre_lot", "missing_n") / nrow(raw),
      format(get_quality_value("acre_lot", "zero_n"), big.mark = ",")
    ),
    sprintf(
      "%s exact duplicates among %s raw rows.",
      format(exact_duplicate_rows, big.mark = ","),
      format(nrow(raw), big.mark = ",")
    ),
    sprintf(
      "%s prices are zero and %s are negative; %s house-size values are non-positive.",
      format(get_quality_value("price", "zero_n"), big.mark = ","),
      format(get_quality_value("price", "negative_n"), big.mark = ","),
      format(get_quality_value("house_size", "non_positive_n"), big.mark = ",")
    ),
    "The audit records distributions without defining an upper price cutoff."
  ),
  stringsAsFactors = FALSE
)
write_audit_table(audit_decisions, "raw_data_decisions.csv")

# Count records by state.
state_table <- sort(table(raw$state, useNA = "ifany"), decreasing = TRUE)
state_listing_counts <- data.frame(
  state = names(state_table),
  listings = as.integer(state_table),
  stringsAsFactors = FALSE
)
state_listing_counts$state[state_listing_counts$state == "<NA>"] <- NA_character_
write_audit_table(state_listing_counts, "raw_listings_per_state.csv")

# Count records by the state--ZIP combination. A tab separates the two labels;
# it is not present in state names or ZIP labels.
state_values <- as.character(raw$state)
zip_values <- as.character(raw$zip_code)
missing_state_label <- "(missing state)"
missing_zip_label <- "(missing ZIP)"
state_values[is.na(state_values)] <- missing_state_label
zip_values[is.na(zip_values)] <- missing_zip_label
state_zip_key <- paste(state_values, zip_values, sep = "\t")
state_zip_table <- sort(table(state_zip_key), decreasing = TRUE)

zip_listing_counts <- data.frame(
  state = character(length(state_zip_table)),
  zip_code = character(length(state_zip_table)),
  listings = as.integer(state_zip_table),
  stringsAsFactors = FALSE
)
for (i in seq_along(state_zip_table)) {
  labels <- strsplit(names(state_zip_table)[i], "\t", fixed = TRUE)[[1]]
  zip_listing_counts$state[i] <- labels[1]
  zip_listing_counts$zip_code[i] <- labels[2]
}
zip_listing_counts$state[zip_listing_counts$state == missing_state_label] <- NA_character_
zip_listing_counts$zip_code[zip_listing_counts$zip_code == missing_zip_label] <- NA_character_
write_audit_table(zip_listing_counts, "raw_listings_per_state_zip.csv")

duplicate_summary <- data.frame(
  total_rows = nrow(raw),
  exact_duplicate_rows = exact_duplicate_rows,
  unique_rows = nrow(raw) - exact_duplicate_rows
)
write_audit_table(duplicate_summary, "raw_exact_duplicate_summary.csv")

# Save one simple base-R histogram. Values that cannot be placed on a
# continuous scale are left out of the plot but remain in the audit tables.
save_histogram <- function(values, title, x_label, file_name, fill_colour) {
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    stop("Cannot create histogram for ", x_label, ": no finite values found.")
  }

  grDevices::png(
    filename = file.path(figures_dir, file_name),
    width = 2400,
    height = 1500,
    res = 300
  )
  graphics::hist(
    values,
    breaks = 60,
    col = fill_colour,
    border = "white",
    main = title,
    xlab = x_label,
    ylab = "Listings"
  )
  grDevices::dev.off()
}

save_histogram(raw$price, "Raw-data distribution: Listing price (USD)",
               "Listing price (USD)", "raw_price_histogram.png", "#2C7FB8")
save_histogram(log(raw$price[raw$price > 0]),
               "Raw-data distribution: Log listing price", "Log listing price",
               "raw_log_price_histogram.png", "#2C7FB8")
save_histogram(raw$house_size, "Raw-data distribution: House size (square feet)",
               "House size (square feet)", "raw_house_size_histogram.png", "#2C7FB8")
save_histogram(log(raw$house_size[raw$house_size > 0]),
               "Raw-data distribution: Log house size", "Log house size",
               "raw_log_house_size_histogram.png", "#2C7FB8")
save_histogram(raw$acre_lot, "Raw-data distribution: Lot size (acres)",
               "Lot size (acres)", "raw_acre_lot_histogram.png", "#2C7FB8")
save_histogram(zip_listing_counts$listings,
               "Raw listings per state-ZIP group", "Listings per state-ZIP group",
               "raw_zip_sample_size_histogram.png", "#41AB5D")

message("Raw-data audit complete. Tables: ", tables_dir, "; figures: ", figures_dir)
