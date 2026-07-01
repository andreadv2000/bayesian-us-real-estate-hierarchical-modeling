# Settings used by the scripts that are implemented so far.
options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  stringsAsFactors = FALSE,
  scipen = 999
)

# `fread()` is retained because the raw CSV is large.  All other work in the
# completed audit uses base R functions.
required_packages <- "data.table"
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]

if (length(missing_packages) > 0) {
  stop(
    "Install the required package before running the audit: ",
    paste(missing_packages, collapse = ", ")
  )
}

set.seed(20260621)
