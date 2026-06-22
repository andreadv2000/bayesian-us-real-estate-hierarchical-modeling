# Calculate and record DIC for the three fitted JAGS models.
#
# This file defines functions only. Source it from scripts/05_fit_models.R while
# the live JAGS model objects are available: JAGS model objects cannot be
# restored from an RDS file and then used with rjags::dic.samples().

required_comparison_packages <- c("dplyr", "here", "readr", "rjags", "tibble")
missing_comparison_packages <- required_comparison_packages[!vapply(
  required_comparison_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing_comparison_packages) > 0L) {
  stop(
    "Install the package(s) required for DIC comparison: ",
    paste(missing_comparison_packages, collapse = ", "),
    call. = FALSE
  )
}

validate_dic_iterations <- function(n_iter) {
  valid <- is.numeric(n_iter) &&
    length(n_iter) == 1L &&
    is.finite(n_iter) &&
    n_iter >= 1L &&
    n_iter %% 1 == 0
  if (!valid) {
    stop("n_iter must be one positive integer.", call. = FALSE)
  }
  as.integer(n_iter)
}

sample_dic <- function(jags_model, n_iter = 5000L) {
  if (!inherits(jags_model, "jags")) {
    stop("jags_model must be a live rjags model object.", call. = FALSE)
  }

  rjags::dic.samples(
    model = jags_model,
    n.iter = validate_dic_iterations(n_iter),
    type = "pD",
    progress.bar = "none"
  )
}

summarise_dic <- function(jags_model, n_iter = 5000L) {
  dic_object <- sample_dic(jags_model, n_iter = n_iter)
  d_bar <- sum(as.numeric(dic_object$deviance))
  p_d <- sum(as.numeric(dic_object$penalty))

  if (!all(is.finite(c(d_bar, p_d)))) {
    stop("DIC calculation returned a non-finite deviance or penalty.", call. = FALSE)
  }

  tibble::tibble(Dbar = d_bar, pD = p_d, DIC = d_bar + p_d)
}

calculate_and_write_dic <- function(fit_m0, fit_m1, fit_m2,
                                    n_iter = 5000L,
                                    table_path = here::here(
                                      "output", "tables", "model_comparison_dic.csv"
                                    ),
                                    log_path = here::here(
                                      "output", "logs", "model_comparison_dic.log"
                                    )) {
  fit_objects <- list(
    "M0: pooled" = fit_m0,
    "M1: state" = fit_m1,
    "M2: state + ZIP" = fit_m2
  )
  invalid_fits <- names(fit_objects)[!vapply(
    fit_objects,
    function(fit) is.list(fit) && inherits(fit$model, "jags"),
    FUN.VALUE = logical(1)
  )]
  if (length(invalid_fits) > 0L) {
    stop(
      "DIC requires live JAGS fits for: ",
      paste(invalid_fits, collapse = ", "),
      ". Run scripts/05_fit_models.R in this R session.",
      call. = FALSE
    )
  }

  n_iter <- validate_dic_iterations(n_iter)
  dic_objects <- lapply(fit_objects, function(fit) {
    sample_dic(fit$model, n_iter = n_iter)
  })
  dic_table <- dplyr::bind_rows(Map(
    function(model_name, dic_object) {
      d_bar <- sum(as.numeric(dic_object$deviance))
      p_d <- sum(as.numeric(dic_object$penalty))
      tibble::tibble(model = model_name, Dbar = d_bar, pD = p_d, DIC = d_bar + p_d)
    },
    names(dic_objects),
    dic_objects
  )) |>
    dplyr::arrange(DIC) |>
    dplyr::mutate(delta_DIC = DIC - min(DIC))

  dic_components <- unlist(dic_table[c("Dbar", "pD", "DIC")], use.names = FALSE)
  if (any(!is.finite(dic_components))) {
    stop("DIC comparison returned invalid values.", call. = FALSE)
  }

  dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(dic_table, table_path)

  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  log_lines <- c(
    "DIC calculation from rjags::dic.samples(type = \"pD\").",
    paste("rjags version:", as.character(utils::packageVersion("rjags"))),
    paste("DIC iterations per model:", n_iter),
    "Dbar is posterior mean deviance; pD is the effective parameter count; DIC = Dbar + pD.",
    "All models were fitted on the same frozen training observations in scripts/05_fit_models.R.",
    "",
    "Printed DIC object for M0: pooled:",
    utils::capture.output(print(dic_objects[["M0: pooled"]]))
  )
  writeLines(log_lines, con = log_path, useBytes = TRUE)

  list(table = dic_table, dic_objects = dic_objects)
}
