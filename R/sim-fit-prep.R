#' Prepare Simulated Outcomes for Model Fitting
#'
#' Build a fitting data set and formula metadata from an `mlsim_data` object
#' that contains an outcome generator.
#'
#' @param sim An `mlsim_data` object returned by [simulate_data()].
#' @param outcome Optional name of the outcome generator to prepare. Required
#'   when `sim` contains more than one outcome generator.
#' @param target Fitting target. `"generic"` returns ordinary data and formulas;
#'   `"brmcoda"` prepares compositional outcomes with
#'   [complr()]; `"auto"` chooses `"brmcoda"` for
#'   compositional outcomes and `"generic"` otherwise.
#' @param drop_incomplete Logical; if `TRUE`, remove rows with missing values in
#'   the generated fitting formulas.
#' @param ... Additional arguments passed to [complr()] when
#'   `target = "brmcoda"`.
#'
#' @return An `mlsim_fit_prep` object containing fitting data, formulas,
#'   helper-column names, term maps, residual-correlation metadata, and target
#'   specific objects such as `complr` for `brmcoda`.
#'
#' @details
#' Outcome formulas can contain simulation helpers such as `lag1()`,
#' `within()`, `between()`, and `ar1()`. `prepare_outcome_fit()` maps those
#' helpers to concrete columns suitable for fitting and records a term map from
#' simulation terms to fitting terms.
#'
#' @examples
#' sim <- simulate_data(
#'   n_groups = 2,
#'   n_per_group = 3,
#'   seed = 20,
#'   generators = list(
#'     x = gen_normal("x"),
#'     y = gen_outcome(y ~ lag1(x), residual_cov = matrix(0, 1, 1, dimnames = list("y", "y")))
#'   )
#' )
#' prep <- prepare_outcome_fit(sim)
#' prep$formula
#' prep$helper_columns
#'
#' @export
prepare_outcome_fit <- function(sim, outcome = NULL,
                                target = c("auto", "generic", "brmcoda"),
                                drop_incomplete = TRUE, ...) {
  if (!inherits(sim, "mlsim_data")) {
    .mlsim_stop("`sim` must be an object returned by `simulate_data()`.")
  }
  target <- match.arg(target)
  drop_incomplete <- .mlsim_check_scalar_logical(drop_incomplete, "drop_incomplete")
  outcome <- .mlsim_resolve_outcome_name(sim, outcome)
  metadata <- sim$generator_metadata[[outcome]]

  if (identical(target, "auto")) {
    target <- if (isTRUE(metadata$compositional)) "brmcoda" else "generic"
  }

  if (identical(target, "generic")) {
    return(.mlsim_prepare_generic_fit(sim, outcome, metadata, drop_incomplete, ...))
  }
  .mlsim_prepare_brmcoda_fit(sim, outcome, metadata, drop_incomplete, ...)
}

#' Internal Outcome Fit Preparation Helpers
#'
#' Resolve outcome generator names, transform formulas and helper columns, and
#' build the internal `mlsim_fit_prep` object used by [prepare_outcome_fit()].
#'
#' @param sim An `mlsim_data` object.
#' @param outcome Outcome generator name.
#' @param metadata Outcome generator metadata.
#' @param drop_incomplete Logical flag controlling row filtering.
#' @param ... Additional target-specific arguments.
#' @param response Outcome response names.
#' @param response_map Mapping between simulated and target response columns.
#' @param data Data table or data frame.
#' @param source Source column name.
#' @param idvar Optional grouping identifier.
#' @param formula Formula to transform or inspect.
#' @param replacements Named character vector of symbol replacements.
#' @param term Model term string.
#'
#' @return Internal helper outputs used by [prepare_outcome_fit()].
#'
#' @examples
#' multilevelcoda:::.mlsim_map_term_string("x:y", c(x = "x_fit"))
#' multilevelcoda:::.mlsim_complete_formula_rows(
#'   data.table::data.table(y = c(1, NA), x = c(2, 3)),
#'   y ~ x
#' )
#'
#' @keywords internal
#' @name multilevelcoda-internal-fit-prep
NULL

#' @rdname multilevelcoda-internal-fit-prep
.mlsim_resolve_outcome_name <- function(sim, outcome) {
  specs <- sim$generator_specs
  outcome_names <- names(specs)[vapply(specs, function(spec) {
    identical(spec$type, "outcome")
  }, logical(1))]

  if (!is.null(outcome)) {
    outcome <- .mlsim_check_vars(outcome, 1L, "outcome")
    if (!(outcome %anyin% outcome_names)) {
      .mlsim_stop("`outcome` must name an outcome generator.")
    }
    return(outcome)
  }

  if (length(outcome_names) == 0L) {
    .mlsim_stop("`sim` does not contain an outcome generator.")
  }
  if (length(outcome_names) > 1L) {
    .mlsim_stop(
      "`sim` contains multiple outcome generators; specify `outcome` as one of: %s.",
      paste(outcome_names, collapse = ", ")
    )
  }
  outcome_names[[1L]]
}

#' @rdname multilevelcoda-internal-fit-prep
.mlsim_prepare_generic_fit <- function(sim, outcome, metadata, drop_incomplete, ...) {
  if (length(list(...)) > 0L) {
    .mlsim_stop("Additional arguments in `...` are only used for `target = \"brmcoda\"`.")
  }
  if (isTRUE(metadata$compositional)) {
    .mlsim_stop("Compositional outcomes default to `target = \"brmcoda\"`; generic fitting prep is not supported yet.")
  }

  data <- data.table::copy(sim$data)
  formula <- metadata$fit_formula
  fixed_formula <- metadata$fit_fixed_formula
  complete_rows <- .mlsim_complete_formula_rows(data, formula)
  if (!is.null(metadata$scale_fit_formula)) {
    complete_rows <- complete_rows & .mlsim_complete_formula_rows(data, metadata$scale_fit_formula)
  }
  if (isTRUE(drop_incomplete)) {
    data <- data[which(complete_rows), ]
  }

  .mlsim_fit_prep(
    target = "generic",
    outcome = outcome,
    data = data,
    formula = formula,
    fixed_formula = fixed_formula,
    term_map = metadata$fit_term_map,
    helper_columns = metadata$helper_columns,
    scale_formula = metadata$scale_fit_formula,
    scale_term_map = metadata$scale_fit_term_map,
    scale_helper_columns = metadata$scale_helper_columns,
    residual_cor = metadata$residual_cor,
    complete_rows = complete_rows
  )
}

#' @rdname multilevelcoda-internal-fit-prep
.mlsim_prepare_brmcoda_fit <- function(sim, outcome, metadata, drop_incomplete, ...) {
  if (!isTRUE(metadata$compositional)) {
    .mlsim_stop("`target = \"brmcoda\"` requires a compositional outcome.")
  }

  data <- data.table::copy(sim$data)
  idvar <- sim$metadata$group_id
  if (is.null(sim$metadata$n_groups)) {
    idvar <- NULL
  }

  complr_args <- c(
    list(
      data = data,
      parts = metadata$parts,
      sbp = metadata$sbp,
      total = metadata$total,
      idvar = idvar,
      transform = "ilr"
    ),
    list(...)
  )
  complr <- do.call("complr", complr_args)
  dataout <- data.table::as.data.table(complr$dataout)

  order_cols <- c(idvar, sim$metadata$time_id %||% sim$metadata$obs_id)
  order_cols <- (order_cols[!is.na(order_cols)]) %sin% names(dataout)
  if (length(order_cols) > 0L) {
    data.table::setorderv(dataout, order_cols)
  }

  response_map <- .mlsim_brmcoda_response_map(metadata$response)
  helper_map <- .mlsim_brmcoda_helper_map(metadata, response_map)
  if (nrow(helper_map) > 0L) {
    for (i in seq_len(nrow(helper_map))) {
      if (isTRUE(helper_map$response_derived[[i]])) {
        source <- helper_map$brmcoda_source[[i]]
        column <- helper_map$brmcoda_column[[i]]
        if (column %anyin% names(dataout)) {
          .mlsim_stop("Fit helper column `%s` would overwrite an existing column.", column)
        }
        data.table::set(dataout, j = column, value = .mlsim_group_lag_values(dataout, source, idvar))
      }
    }
  }

  symbol_map <- c(
    stats::setNames(response_map$brmcoda_response, response_map$response),
    stats::setNames(helper_map$brmcoda_column, helper_map$column)
  )
  formula <- .mlsim_formula_replace_symbols(metadata$fit_formula, symbol_map)
  term_map <- metadata$fit_term_map
  if (!is.null(term_map) && nrow(term_map) > 0L) {
    term_map$fit_term <- vapply(
      term_map$fit_term,
      .mlsim_map_term_string,
      character(1),
      replacements = symbol_map
    )
  }
  scale_formula <- metadata$scale_fit_formula
  scale_term_map <- metadata$scale_fit_term_map
  if (!is.null(scale_formula)) {
    scale_formula <- .mlsim_formula_replace_symbols(scale_formula, symbol_map)
  }
  if (!is.null(scale_term_map) && nrow(scale_term_map) > 0L) {
    scale_term_map$fit_term <- vapply(
      scale_term_map$fit_term,
      .mlsim_map_term_string,
      character(1),
      replacements = symbol_map
    )
  }

  complete_rows <- .mlsim_complete_formula_rows(dataout, formula)
  if (!is.null(scale_formula)) {
    complete_rows <- complete_rows & .mlsim_complete_formula_rows(dataout, scale_formula)
  }
  if (isTRUE(drop_incomplete)) {
    dataout <- dataout[which(complete_rows), ]
  }
  complr$dataout <- dataout

  .mlsim_fit_prep(
    target = "brmcoda",
    outcome = outcome,
    data = dataout,
    complr = complr,
    formula = formula,
    term_map = term_map,
    helper_columns = unique(helper_map$brmcoda_column),
    scale_formula = scale_formula,
    scale_term_map = scale_term_map,
    scale_helper_columns = .mlsim_brmcoda_helper_columns(
      metadata$scale_helper_columns,
      helper_map
    ),
    residual_cor = metadata$residual_cor,
    complete_rows = complete_rows,
    parts = metadata$parts,
    sbp = metadata$sbp
  )
}

#' @rdname multilevelcoda-internal-fit-prep
.mlsim_brmcoda_response_map <- function(response) {
  data.frame(
    response = response,
    brmcoda_response = paste0("z", seq_along(response), "_1"),
    stringsAsFactors = FALSE
  )
}

#' @rdname multilevelcoda-internal-fit-prep
.mlsim_brmcoda_helper_map <- function(metadata, response_map) {
  helpers <- metadata$fit_helpers
  if (is.null(helpers) || nrow(helpers) == 0L) {
    return(data.frame(
      column = character(),
      brmcoda_column = character(),
      brmcoda_source = character(),
      response_derived = logical(),
      stringsAsFactors = FALSE
    ))
  }
  helpers$brmcoda_column <- helpers$column
  helpers$brmcoda_source <- helpers$source

  for (i in seq_len(nrow(response_map))) {
    source <- response_map$response[[i]]
    z <- response_map$brmcoda_response[[i]]
    match <- helpers$response_derived & helpers$source == source
    helpers$brmcoda_source[match] <- z
    helpers$brmcoda_column[match] <- paste0("lag_", z)
  }

  helpers
}

#' @rdname multilevelcoda-internal-fit-prep
.mlsim_brmcoda_helper_columns <- function(helper_columns, helper_map) {
  helper_columns <- helper_columns %||% character()
  if (length(helper_columns) == 0L) {
    return(character())
  }
  if (is.null(helper_map) || nrow(helper_map) == 0L) {
    .mlsim_stop(
      "Fit helper metadata is inconsistent; unknown helper columns: %s.",
      paste(unique(helper_columns), collapse = ", ")
    )
  }

  unmatched <- helper_columns[!(helper_columns %in% helper_map$column)]
  if (length(unmatched) > 0L) {
    .mlsim_stop(
      "Fit helper metadata is inconsistent; unknown helper columns: %s.",
      paste(unique(unmatched), collapse = ", ")
    )
  }

  unique(helper_map$brmcoda_column[match(helper_columns, helper_map$column)])
}

#' @rdname multilevelcoda-internal-fit-prep
.mlsim_group_lag_values <- function(data, source, idvar) {
  if (source %any!in% names(data)) {
    .mlsim_stop("Cannot create lag helper because `%s` is not available.", source)
  }
  if (is.null(idvar)) {
    return(data.table::shift(data[[source]]))
  }
  if (idvar %any!in% names(data)) {
    .mlsim_stop("Cannot create lag helper because `%s` is not available.", idvar)
  }

  out <- rep(NA, nrow(data))
  for (idx in split(seq_len(nrow(data)), data[[idvar]])) {
    out[idx] <- data.table::shift(data[[source]][idx])
  }
  out
}

#' @rdname multilevelcoda-internal-fit-prep
.mlsim_formula_replace_symbols <- function(formula, replacements) {
  stats::as.formula(
    .mlsim_replace_symbols(formula, replacements),
    env = environment(formula)
  )
}

#' @rdname multilevelcoda-internal-fit-prep
.mlsim_map_term_string <- function(term, replacements) {
  if (identical(term, "(Intercept)")) {
    return(term)
  }
  pieces <- strsplit(term, ":", fixed = TRUE)[[1L]]
  pieces <- vapply(pieces, function(piece) {
    if (piece %anyin% names(replacements)) {
      replacements[[piece]]
    } else {
      piece
    }
  }, character(1L))
  paste(pieces, collapse = ":")
}

#' @rdname multilevelcoda-internal-fit-prep
.mlsim_complete_formula_rows <- function(data, formula) {
  vars <- all.vars(formula)
  vars <- vars %sin% names(data)
  if (length(vars) == 0L) {
    return(rep(TRUE, nrow(data)))
  }
  stats::complete.cases(as.data.frame(data[, vars, with = FALSE]))
}

#' @rdname multilevelcoda-internal-fit-prep
.mlsim_fit_prep <- function(...) {
  structure(list(...), class = "mlsim_fit_prep")
}

#' Print an Outcome Fit Preparation Object
#'
#' @param x An `mlsim_fit_prep` object.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @method print mlsim_fit_prep
#' @export
print.mlsim_fit_prep <- function(x, ...) {
  cat("<mlsim_fit_prep>\n")
  cat("  target: ", x$target, "\n", sep = "")
  cat("  outcome: ", x$outcome, "\n", sep = "")
  if (!is.null(x$data)) {
    cat("  rows: ", nrow(x$data), "\n", sep = "")
  }
  if (!is.null(x$helper_columns) && length(x$helper_columns) > 0L) {
    cat("  helper columns: ", paste(x$helper_columns, collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$formula)) {
    cat("  formula: ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  }
  invisible(x)
}
