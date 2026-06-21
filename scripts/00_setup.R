options(repos = c(CRAN = "https://cloud.r-project.org"))

required_packages <- c(
  "data.table", "dplyr", "tidyr", "stringr", "forcats",
  "ggplot2", "readr", "janitor", "here", "fs",
  "rjags", "coda", "posterior", "bayesplot",
  "mcmcse", "broom", "knitr", "kableExtra",
  "arrow", "sessioninfo", "purrr", "tibble"
)

missing <- required_packages[
  !required_packages %in% rownames(installed.packages())
]
if (length(missing) > 0) {
  install.packages(missing)
}

invisible(lapply(required_packages, library, character.only = TRUE))

jags_runtime_version <- rjags::jags.version()
message("JAGS version: ", as.character(jags_runtime_version))

set.seed(20260621)
options(stringsAsFactors = FALSE, scipen = 999)
