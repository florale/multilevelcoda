# Outcome formula parsing, helper evaluation, and model-matrix internals.

#' @rdname multilevelcoda-internal-outcome
.mlsim_check_outcome_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    .mlsim_stop("`formula` must be a two-sided formula.")
  }
  formula
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_check_scale_formula <- function(scale_formula) {
  if (is.null(scale_formula)) {
    return(NULL)
  }
  if (!inherits(scale_formula, "formula") || !(length(scale_formula) %anyin% c(2L, 3L))) {
    .mlsim_stop("`scale_formula` must be a one-sided formula or `sigma ~ ...`.")
  }
  if (length(scale_formula) == 3L) {
    lhs <- scale_formula[[2L]]
    if (!is.symbol(lhs) || !identical(as.character(lhs), "sigma")) {
      .mlsim_stop("The left-hand side of `scale_formula` must be `sigma`.")
    }
  }
  scale_formula
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_outcome_response <- function(formula) {
  lhs <- formula[[2L]]
  if (is.symbol(lhs)) {
    return(.mlsim_check_vars(as.character(lhs), 1L, "response"))
  }
  if (is.call(lhs) && .mlsim_call_name(lhs) %anyin% c("cbind", "mvbind")) {
    responses <- vapply(as.list(lhs)[-1L], function(x) {
      if (!is.symbol(x)) {
        .mlsim_stop("Multivariate outcome responses must be bare column names.")
      }
      as.character(x)
    }, character(1L))
    return(.mlsim_check_vars(responses, arg = "response"))
  }
  .mlsim_stop("Outcome response must be a name, `cbind(...)`, or `mvbind(...)`.")
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_outcome_composition <- function(response, compositional, parts, sbp, total, keep_ilr) {
  d <- length(response)
  compositional <- .mlsim_check_scalar_logical(compositional, "compositional")
  keep_ilr <- .mlsim_check_scalar_logical(keep_ilr, "keep_ilr")
  if (compositional) {
    composition <- .mlsim_resolve_composition(parts, sbp, d)
    total <- .mlsim_check_positive(total, "total")
    if (length(total) != 1L) {
      .mlsim_stop("`total` must be a scalar positive value.")
    }
    .mlsim_check_composition_output_names(response, composition$parts, keep_ilr)
    output_names <- if (keep_ilr) c(response, composition$parts) else composition$parts
    return(list(
      compositional = TRUE,
      response = response,
      output_names = output_names,
      parts = composition$parts,
      sbp = composition$sbp,
      total = total,
      keep_ilr = keep_ilr
    ))
  }

  if (!is.null(parts) || !is.null(sbp)) {
    .mlsim_stop("`parts` and `sbp` require `compositional = TRUE`.")
  }
  list(
    compositional = FALSE,
    response = response,
    output_names = response,
    parts = NULL,
    sbp = NULL,
    total = NULL,
    keep_ilr = FALSE
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_check_burnin <- function(burnin) {
  if (identical(burnin, "auto")) {
    return(burnin)
  }
  burnin <- .mlsim_recycle_integer(burnin, 1L, "burnin")
  burnin
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_outcome_state <- function(data, context, response) {
  state <- new.env(parent = emptyenv())
  state$data <- as.data.frame(data)
  state$context <- context
  state$response <- response
  state$name_map <- character()
  state$brms_map <- character()
  state$fit_map <- character()
  state$fit_helpers <- list()
  state$ar_internal <- character()
  state$ar_names <- character()
  state
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_transform_outcome_expr <- function(expr, state, allow_ar1 = TRUE) {
  if (!is.call(expr)) {
    return(list(expr = expr, state = state))
  }

  name <- .mlsim_call_name(expr)
  if (identical(name, "(")) {
    inner <- .mlsim_transform_outcome_expr(expr[[2L]], state, allow_ar1 = allow_ar1)$expr
    return(list(expr = as.call(list(as.name("("), inner)), state = state))
  }
  if (name %anyin% c("lag1", "within", "between", "ar1")) {
    if (identical(name, "ar1") && !isTRUE(allow_ar1)) {
      .mlsim_stop("`ar1()` is not supported in `scale_formula`.")
    }
    return(list(expr = .mlsim_outcome_helper_expr(expr, state), state = state))
  }

  args <- as.list(expr)
  args[-1L] <- lapply(args[-1L], function(arg) {
    .mlsim_transform_outcome_expr(arg, state, allow_ar1 = allow_ar1)$expr
  })
  list(expr = as.call(args), state = state)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_outcome_helper_expr <- function(expr, state) {
  name <- .mlsim_call_name(expr)
  args <- as.list(expr)[-1L]
  if (identical(name, "ar1")) {
    if (length(args) != 0L) {
      .mlsim_stop("`ar1()` does not accept arguments.")
    }
    internals <- vapply(state$response, function(response) {
      internal <- .mlsim_internal_name("ar1", response)
      if (internal %any!in% names(state$data)) {
        state$data[[internal]] <- 0
        state$name_map[[internal]] <- sprintf("ar1(%s)", response)
        state$brms_map[[internal]] <- .mlsim_internal_name("lag1", response)
        fit_column <- .mlsim_fit_column_name("ar1", response)
        .mlsim_register_fit_helper(
          state = state,
          type = "ar1",
          source = response,
          internal = internal,
          column = fit_column,
          sim_term = sprintf("ar1(%s)", response),
          response_derived = TRUE
        )
        state$ar_internal <- c(state$ar_internal, internal)
        state$ar_names <- c(state$ar_names, response)
      }
      internal
    }, character(1L))
    return(.mlsim_sum_symbols(internals))
  }

  if (length(args) != 1L || !is.symbol(args[[1L]])) {
    .mlsim_stop("`%s()` requires one bare column name.", name)
  }
  var <- as.character(args[[1L]])
  if (var %any!in% names(state$data)) {
    .mlsim_stop("`%s` is not available for outcome simulation.", var)
  }
  internal <- .mlsim_internal_name(name, var)
  if (internal %any!in% names(state$data)) {
    values <- switch(
      name,
      lag1 = .mlsim_lag1(state$data[[var]], state$context),
      within = .mlsim_within_group(state$data[[var]], state$context, var),
      between = .mlsim_between_group(state$data[[var]], state$context, var)
    )
    state$data[[internal]] <- values
    state$name_map[[internal]] <- sprintf("%s(%s)", name, var)
    state$brms_map[[internal]] <- internal
    .mlsim_register_fit_helper(
      state = state,
      type = name,
      source = var,
      internal = internal,
      column = .mlsim_fit_column_name(name, var),
      sim_term = sprintf("%s(%s)", name, var),
      response_derived = FALSE
    )
  }
  as.name(internal)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_model_matrix_row <- function(setup, row_data) {
  .mlsim_model_matrix_data(setup, row_data)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_model_matrix_data <- function(setup, data) {
  x <- .mlsim_model_matrix(
    setup$fixed_formula,
    data = data,
    contrasts.arg = setup$fixed_contrasts,
    xlev = setup$fixed_xlevels
  )
  .mlsim_align_model_matrix(x, colnames(setup$fixed_matrix))
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_random_model_matrix <- function(spec, data) {
  z <- .mlsim_model_matrix(
    spec$z_formula,
    data = data,
    contrasts.arg = spec$z_contrasts,
    xlev = spec$z_xlevels
  )
  .mlsim_align_model_matrix(z, colnames(spec$z))
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_align_model_matrix <- function(x, columns) {
  missing <- setdiff(columns, colnames(x))
  if (length(missing) > 0L) {
    x <- cbind(x, matrix(0, nrow = nrow(x), ncol = length(missing), dimnames = list(NULL, missing)))
  }
  x[, columns, drop = FALSE]
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_model_xlevels <- function(formula, data) {
  frame <- .mlsim_model_frame(formula, data)
  stats::.getXlevels(attr(frame, "terms"), frame)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_model_matrix <- function(formula, data, contrasts.arg = NULL, xlev = NULL) {
  frame <- .mlsim_model_frame(formula, data, xlev = xlev)
  stats::model.matrix(attr(frame, "terms"), frame, contrasts.arg = contrasts.arg)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_model_frame <- function(formula, data, xlev = NULL) {
  stats::model.frame(
    formula,
    data = data,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE,
    xlev = xlev
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_lag1 <- function(x, context) {
  out <- rep(NA, length(x))
  groups <- if (is.null(context$n_groups)) rep(1L, length(x)) else context$group_index
  for (idx in split(seq_along(x), groups)) {
    if (length(idx) > 1L) {
      out[idx[-1L]] <- x[idx[-length(idx)]]
    }
  }
  out
}


#' @rdname multilevelcoda-internal-outcome
#' @importFrom stats ave
.mlsim_within_group <- function(x, context, var) {
  if (is.null(context$n_groups)) {
    .mlsim_stop("`within(%s)` requires grouped data.", var)
  }
  if (!is.numeric(x)) {
    .mlsim_stop("`within(%s)` requires a numeric column.", var)
  }
  means <- ave(x, context$group_index, FUN = function(z) mean(z, na.rm = TRUE))
  x - means
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_between_group <- function(x, context, var) {
  if (is.null(context$n_groups)) {
    .mlsim_stop("`between(%s)` requires grouped data.", var)
  }
  if (!is.numeric(x)) {
    .mlsim_stop("`between(%s)` requires a numeric column.", var)
  }
  ave(x, context$group_index, FUN = function(z) mean(z, na.rm = TRUE))
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_replace_symbols <- function(expr, replacements) {
  if (is.symbol(expr)) {
    name <- as.character(expr)
    if (name %anyin% names(replacements)) {
      return(as.name(replacements[[name]]))
    }
    return(expr)
  }
  if (!is.call(expr)) {
    return(expr)
  }
  args <- as.list(expr)
  args[-1L] <- lapply(args[-1L], .mlsim_replace_symbols, replacements = replacements)
  as.call(args)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_pretty_design_names <- function(names, state) {
  pretty <- names
  for (internal in names(state$name_map)) {
    pretty <- gsub(internal, state$name_map[[internal]], pretty, fixed = TRUE)
  }
  pretty
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_internal_name <- function(prefix, name) {
  paste0(".mlsim_", prefix, "_", make.names(name))
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_rhs_formula <- function(expr, env) {
  stats::as.formula(as.call(list(as.name("~"), expr)), env = env)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_plus_formula <- function(left, right) {
  .mlsim_rhs_formula(.mlsim_plus_expr(left[[2L]], right[[2L]]), environment(left))
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_sum_symbols <- function(names) {
  expr <- as.name(names[[1L]])
  if (length(names) == 1L) {
    return(expr)
  }
  for (name in names[-1L]) {
    expr <- as.call(list(as.name("+"), expr, as.name(name)))
  }
  as.call(list(as.name("("), expr))
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_plus_expr <- function(left, right) {
  if (is.null(left)) {
    return(right)
  }
  if (is.null(right)) {
    return(left)
  }
  as.call(list(as.name("+"), left, right))
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_is_random_call <- function(expr) {
  is.call(expr) && identical(.mlsim_call_name(expr), "|")
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_unwrap_parens <- function(expr) {
  while (is.call(expr) && identical(.mlsim_call_name(expr), "(")) {
    expr <- expr[[2L]]
  }
  expr
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_call_name <- function(expr) {
  if (!is.call(expr)) {
    return(NULL)
  }
  as.character(expr[[1L]])
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_deparse_one <- function(expr) {
  paste(deparse(expr, width.cutoff = 500L), collapse = "")
}
