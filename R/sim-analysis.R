#' Prepare Simulated Outcome Data for Analysis
#'
#' Convert an [simulate_data()] result with a [gen_outcome()] generator into an
#' analysis-ready data set and inferred `brms` formula. The helper recomputes
#' observed between- and within-person predictors from the raw simulated
#' columns, creates within-person centered lag columns for `ar1()` terms, and
#' creates a [complr()] object when the outcome is compositional.
#'
#' The inferred analysis model is a deliberately *pragmatic default estimator*
#' built from observed data, not the matched model for the simulated
#' data-generating process. See the "Pragmatic default estimator" section
#' before interpreting parameter recovery results.
#'
#' @param sim An `mlsim_data` object returned by [simulate_data()].
#' @param outcome Optional character scalar naming the outcome generator. When
#'   `NULL`, the helper uses the only generator with
#'   `distribution = "outcome"`.
#' @param drop_lag_na Logical scalar. When `FALSE`, the default, first rows in
#'   each series are retained and lag-derived columns are left as `NA`. When
#'   `TRUE`, rows with missing lag-derived predictors are removed.
#'
#' @return An `mlsim_analysis` object, a list with:
#'   \item{`data`}{Analysis data with derived between/within and lag columns.}
#'   \item{`formula`}{An inferred `brms` formula.}
#'   \item{`complr`}{A `complr` object for compositional outcomes, otherwise
#'   `NULL`.}
#'   \item{`metadata`}{Preparation metadata, including derived column names and
#'   formula mappings.}
#'
#' @details
#' The analysis formula is inferred from the stored `gen_outcome()` formula.
#' Simulation terms `between(x)` and `within(x)` become observed-data columns
#' named `x_between` and `x_within`, computed from `x` by the simulation group
#' identifier. These columns are recomputed even when columns with the same
#' names already exist in `sim$data`.
#'
#' For dynamic formulas, `ar1()` is translated to lagged observed response
#' columns centered at each person's observed mean of all their values (not
#' the mean of the lagged values). For compositional outcomes, the helper
#' rebuilds the ILR coordinates through [complr()] using the simulator's parts
#' and SBP metadata, then lags the generated `z` coordinates used by
#' [brmcoda()].
#'
#' @section Pragmatic default estimator:
#' [gen_outcome()] simulates latent residual AR/VAR dynamics around the
#' model-implied mean and resolves `between()`/`within()` from latent
#' generating components supplied by upstream predictor generators.
#' `prep_sim_analysis()` instead constructs the model applied analysts
#' commonly fit to observed data:
#' \itemize{
#'   \item `between(x)` and `within(x)` become `x_between` and `x_within`,
#'     recomputed from realised person means of the observed `x` (manifest
#'     centering), not from the latent generating components.
#'   \item `ar1()` becomes person-mean-centered lagged \emph{observed}
#'     response predictors (`lag_<response>_within`), not the latent residual
#'     state.
#' }
#' These observed-data constructions target different estimands from the
#' simulation truth. Manifest person-mean centering and observed-score lagged
#' regression are known to yield biased estimates of between-person effects
#' and of inertia and cross-lag parameters relative to the latent generating
#' values, especially with short series (Ludtke et al. 2008; Hamaker &
#' Grasman 2014). This mismatch is intentional: it lets simulation studies
#' quantify the bias of the pragmatic estimator. Do not interpret systematic
#' discrepancies between estimates from this default analysis model and the
#' `gen_outcome()` truth parameters as errors in the simulator. For
#' matched-model recovery studies, construct the analysis model by hand
#' instead of using this helper.
#'
#' @examples
#' params <- list(
#'   location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
#'   scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y")))
#' )
#' sim <- simulate_data(
#'   n = 5,
#'   seed = 1,
#'   generators = list(
#'     outcome = gen_outcome(
#'       y ~ 1,
#'       scale = sigma ~ 1,
#'       params = params,
#'       burnin = 0
#'     )
#'   )
#' )
#' analysis <- prep_sim_analysis(sim)
#' analysis$formula
#'
#' @export
prep_sim_analysis <- function(sim, outcome = NULL, drop_lag_na = FALSE) {
  if (!inherits(sim, "mlsim_data")) {
    .mlsim_stop("`sim` must be an `mlsim_data` object returned by `simulate_data()`.")
  }
  drop_lag_na <- .mlsim_check_scalar_logical(drop_lag_na, "drop_lag_na")
  outcome <- .mlsim_analysis_resolve_outcome(sim, outcome)
  outcome_metadata <- sim$generator_metadata[[outcome]]
  parsed <- outcome_metadata$parsed %||% list()
  outcomes <- parsed$outcomes
  if (length(outcomes) == 0L) {
    .mlsim_stop("Outcome generator `%s` does not contain parsed outcome metadata.", outcome)
  }

  data <- data.table::copy(sim$data)
  group_id <- sim$metadata$group_id
  grouped <- !is.null(sim$metadata$n_groups)
  time_id <- sim$metadata$time_id
  has_ar <- length(parsed$ar_terms %||% character()) > 0L

  if (isTRUE(grouped) && !group_id %in% names(data)) {
    .mlsim_stop("The active grouping column `%s` is missing from `sim$data`.", group_id)
  }
  if (isTRUE(has_ar)) {
    if (is.null(time_id) || !time_id %in% names(data)) {
      .mlsim_stop("Preparing `ar1()` analysis terms requires the simulation `time_id` column.")
    }
    order_cols <- if (isTRUE(grouped)) c(group_id, time_id) else time_id
    data.table::setorderv(data, order_cols)
  }

  derived_roles <- .mlsim_analysis_add_role_columns(data, outcome_metadata, grouped, group_id)
  compositional <- isTRUE(outcome_metadata$compositional)
  response_map <- NULL
  complr_object <- NULL

  if (isTRUE(compositional)) {
    if (is.null(outcome_metadata$parts) || is.null(outcome_metadata$sbp)) {
      .mlsim_stop("Compositional outcome `%s` is missing parts or SBP metadata.", outcome)
    }
    complr_object <- complr(
      data = data,
      parts = outcome_metadata$parts,
      sbp = outcome_metadata$sbp,
      idvar = if (isTRUE(grouped)) group_id else NULL,
      total = outcome_metadata$total %||% 1
    )
    data <- data.table::copy(complr_object$dataout)
    response_names <- colnames(complr_object$output[[1L]]$Z)
    response_map <- stats::setNames(response_names, outcomes)
  } else {
    response_map <- stats::setNames(outcomes, outcomes)
  }

  lag_info <- .mlsim_analysis_add_lag_columns(
    data = data,
    response_names = unname(response_map),
    has_ar = has_ar,
    grouped = grouped,
    group_id = group_id,
    drop_lag_na = drop_lag_na
  )
  data <- lag_info$data

  if (!is.null(complr_object)) {
    complr_object$dataout <- data.table::copy(data)
    if (isTRUE(lag_info$dropped_rows > 0L)) {
      keep_rows <- lag_info$keep_rows
      complr_object$datain <- complr_object$datain[keep_rows]
      complr_object$output <- lapply(complr_object$output, function(out) {
        for (element in c("X", "bX", "wX", "Z", "bZ", "wZ", "dataout")) {
          if (!is.null(out[[element]])) {
            out[[element]] <- out[[element]][keep_rows, , drop = FALSE]
          }
        }
        out
      })
    }
  }

  formula <- .mlsim_analysis_build_brms_formula(
    formula = outcome_metadata$formula,
    scale = outcome_metadata$scale,
    response_map = response_map,
    lag_columns = lag_info$lag_columns
  )

  metadata <- list(
    outcome_generator = outcome,
    compositional = compositional,
    group_id = if (isTRUE(grouped)) group_id else NULL,
    time_id = time_id,
    drop_lag_na = drop_lag_na,
    original_formula = outcome_metadata$formula,
    original_scale = outcome_metadata$scale,
    response_map = response_map,
    derived_roles = derived_roles,
    lag_columns = lag_info$lag_columns,
    dropped_rows = lag_info$dropped_rows,
    formula = formula
  )

  structure(
    list(
      data = data,
      formula = formula,
      complr = complr_object,
      metadata = metadata
    ),
    class = "mlsim_analysis"
  )
}

#' @export
print.mlsim_analysis <- function(x, ...) {
  cat("<mlsim_analysis>\n")
  cat("  outcome: ", x$metadata$outcome_generator, "\n", sep = "")
  cat("  rows: ", nrow(x$data), "\n", sep = "")
  cat("  compositional: ", if (isTRUE(x$metadata$compositional)) "yes" else "no", "\n", sep = "")
  if (length(x$metadata$derived_roles$column %||% character()) > 0L) {
    cat("  derived predictors: ", paste(unique(x$metadata$derived_roles$column), collapse = ", "), "\n", sep = "")
  }
  if (length(x$metadata$lag_columns) > 0L) {
    cat("  lag predictors: ", paste(x$metadata$lag_columns, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

.mlsim_analysis_resolve_outcome <- function(sim, outcome) {
  metadata <- sim$generator_metadata
  if (length(metadata) == 0L) {
    .mlsim_stop("Could not find an outcome generator in `sim$generator_metadata`.")
  }
  outcome_generators <- names(metadata)[vapply(metadata, function(x) {
    identical(x$distribution, "outcome")
  }, logical(1))]
  if (is.null(outcome)) {
    if (length(outcome_generators) != 1L) {
      .mlsim_stop(
        "Could not infer the outcome generator; found %d outcome generators. Supply `outcome`.",
        length(outcome_generators)
      )
    }
    return(outcome_generators[[1L]])
  }
  .mlsim_check_vars(outcome, 1L, "outcome")
  if (!outcome %in% names(metadata)) {
    .mlsim_stop("`outcome` must name a generator in `sim$generator_metadata`.")
  }
  if (!identical(metadata[[outcome]]$distribution, "outcome")) {
    .mlsim_stop("Generator `%s` is not a `gen_outcome()` outcome generator.", outcome)
  }
  outcome
}

.mlsim_analysis_add_role_columns <- function(data, outcome_metadata, grouped, group_id) {
  roles <- outcome_metadata$selected_column_roles
  if (is.null(roles) || nrow(roles) == 0L) {
    return(data.table::data.table(
      variable = character(),
      component = character(),
      column = character(),
      overwritten = logical()
    ))
  }
  roles <- unique(roles[, c("variable", "component"), with = FALSE])
  roles <- roles[roles$component %in% c("between", "within")]
  if (nrow(roles) == 0L) {
    return(data.table::data.table(
      variable = character(),
      component = character(),
      column = character(),
      overwritten = logical()
    ))
  }
  if (!isTRUE(grouped)) {
    .mlsim_stop("Preparing `between()` or `within()` analysis terms requires a grouped simulation design.")
  }
  if (!group_id %in% names(data)) {
    .mlsim_stop("The active grouping column `%s` is missing from `sim$data`.", group_id)
  }

  out <- data.table::data.table(
    variable = roles$variable,
    component = roles$component,
    column = paste(roles$variable, roles$component, sep = "_"),
    overwritten = paste(roles$variable, roles$component, sep = "_") %in% names(data)
  )
  for (variable in unique(roles$variable)) {
    if (!variable %in% names(data)) {
      .mlsim_stop("Cannot prepare analysis terms because observed predictor `%s` is missing.", variable)
    }
    if (!is.numeric(data[[variable]])) {
      .mlsim_stop("Observed predictor `%s` must be numeric to compute between/within analysis terms.", variable)
    }
    needed <- roles$component[roles$variable == variable]
    mean_col <- paste(variable, "between", sep = "_")
    within_col <- paste(variable, "within", sep = "_")
    mean_tmp <- paste0(".__mlsim_analysis_mean_", make.names(variable), "__")
    data[, (mean_tmp) := mean(get(variable), na.rm = TRUE), by = group_id]
    if ("between" %in% needed) {
      data[, (mean_col) := get(mean_tmp)]
    }
    if ("within" %in% needed) {
      data[, (within_col) := get(variable) - get(mean_tmp)]
    }
    data[, (mean_tmp) := NULL]
  }
  out
}

.mlsim_analysis_add_lag_columns <- function(data, response_names, has_ar, grouped,
                                            group_id, drop_lag_na) {
  if (!isTRUE(has_ar)) {
    return(list(
      data = data,
      lag_columns = character(),
      keep_rows = rep(TRUE, nrow(data)),
      dropped_rows = 0L
    ))
  }
  missing_responses <- response_names[!response_names %in% names(data)]
  if (length(missing_responses) > 0L) {
    .mlsim_stop(
      "Cannot prepare `ar1()` analysis terms because response columns are missing: %s.",
      paste(missing_responses, collapse = ", ")
    )
  }
  lag_columns <- paste0("lag_", response_names, "_within")
  for (i in seq_along(response_names)) {
    response <- response_names[[i]]
    lag_col <- lag_columns[[i]]
    raw_col <- paste0(".__mlsim_lag_", response, "__")
    if (isTRUE(grouped)) {
      data[, (raw_col) := data.table::shift(get(response)), by = group_id]
      data[, (lag_col) := get(raw_col) - mean(get(response), na.rm = TRUE), by = group_id]
    } else {
      data[, (raw_col) := data.table::shift(get(response))]
      data[, (lag_col) := get(raw_col) - mean(get(response), na.rm = TRUE)]
    }
    data[, (raw_col) := NULL]
  }
  keep_rows <- rep(TRUE, nrow(data))
  if (isTRUE(drop_lag_na)) {
    keep_rows <- stats::complete.cases(data[, lag_columns, with = FALSE])
    data <- data[keep_rows]
  }
  list(
    data = data,
    lag_columns = lag_columns,
    keep_rows = keep_rows,
    dropped_rows = sum(!keep_rows)
  )
}

.mlsim_analysis_build_brms_formula <- function(formula, scale, response_map, lag_columns) {
  rhs <- .mlsim_analysis_transform_expr(formula[[3L]], lag_columns)
  scale_rhs <- .mlsim_analysis_transform_expr(scale[[3L]], character())
  response_names <- unname(response_map)

  formulas <- lapply(response_names, function(response) {
    brms::bf(
      .mlsim_analysis_formula(as.name(response), rhs),
      .mlsim_analysis_formula(as.name("sigma"), scale_rhs),
      family = stats::gaussian()
    )
  })
  if (length(formulas) == 1L) {
    return(formulas[[1L]])
  }
  Reduce(`+`, formulas) + brms::set_rescor(TRUE)
}

.mlsim_analysis_formula <- function(lhs, rhs) {
  out <- eval(call("~", lhs, rhs), envir = parent.frame())
  environment(out) <- parent.frame()
  out
}

.mlsim_analysis_transform_expr <- function(expr, lag_columns) {
  if (is.symbol(expr) || is.atomic(expr)) {
    return(expr)
  }
  if (!is.call(expr)) {
    return(expr)
  }
  op <- as.character(expr[[1L]])
  if (identical(op, "(")) {
    return(call("(", .mlsim_analysis_transform_expr(expr[[2L]], lag_columns)))
  }
  if (identical(op, "between") || identical(op, "within")) {
    if (length(expr) != 2L || !is.symbol(expr[[2L]])) {
      .mlsim_stop("`%s()` terms must contain one variable name.", op)
    }
    return(as.name(paste(as.character(expr[[2L]]), op, sep = "_")))
  }
  if (identical(op, "ar1")) {
    if (length(expr) != 1L) {
      .mlsim_stop("`ar1()` does not take arguments.")
    }
    if (length(lag_columns) == 0L) {
      .mlsim_stop("Cannot translate `ar1()` because no lag columns were prepared.")
    }
    return(.mlsim_analysis_sum_expr(lapply(lag_columns, as.name)))
  }
  if (identical(op, "+")) {
    return(.mlsim_analysis_sum_expr(list(
      .mlsim_analysis_transform_expr(expr[[2L]], lag_columns),
      .mlsim_analysis_transform_expr(expr[[3L]], lag_columns)
    )))
  }
  if (identical(op, ":")) {
    left <- .mlsim_analysis_split_sum(.mlsim_analysis_transform_expr(expr[[2L]], lag_columns))
    right <- .mlsim_analysis_split_sum(.mlsim_analysis_transform_expr(expr[[3L]], lag_columns))
    expanded <- unlist(lapply(left, function(l) {
      lapply(right, function(r) call(":", l, r))
    }), recursive = FALSE)
    return(.mlsim_analysis_sum_expr(expanded))
  }
  if (identical(op, "|")) {
    return(call("|", .mlsim_analysis_transform_expr(expr[[2L]], lag_columns), expr[[3L]]))
  }
  args <- lapply(as.list(expr[-1L]), .mlsim_analysis_transform_expr, lag_columns = lag_columns)
  as.call(c(list(expr[[1L]]), args))
}

.mlsim_analysis_split_sum <- function(expr) {
  if (is.call(expr) && identical(as.character(expr[[1L]]), "+")) {
    return(c(.mlsim_analysis_split_sum(expr[[2L]]), .mlsim_analysis_split_sum(expr[[3L]])))
  }
  list(expr)
}

.mlsim_analysis_sum_expr <- function(exprs) {
  exprs <- unlist(lapply(exprs, .mlsim_analysis_split_sum), recursive = FALSE)
  if (length(exprs) == 0L) {
    return(1)
  }
  Reduce(function(x, y) call("+", x, y), exprs)
}
