# Outcome coefficient and covariance validation internals.

#' @rdname multilevelcoda-internal-outcome
.mlsim_default_coefficients <- function(setup) {
  values <- matrix(0, nrow = length(setup$fixed_names), ncol = length(setup$response))
  rownames(values) <- setup$fixed_names
  colnames(values) <- setup$response
  values
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_default_scale_coefficients <- function(setup) {
  if (!isTRUE(setup$scale$enabled)) {
    return(NULL)
  }
  values <- matrix(0, nrow = length(setup$scale$fixed_names), ncol = length(setup$response))
  rownames(values) <- setup$scale$fixed_names
  colnames(values) <- setup$response
  values
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_default_random_cov <- function(setup) {
  lapply(setup$joint_random_specs, function(spec) {
    cov <- matrix(0, nrow = length(spec$coef_names), ncol = length(spec$coef_names))
    rownames(cov) <- colnames(cov) <- spec$coef_names
    cov
  })
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_check_outcome_coefficients <- function(coefficients, setup) {
  default <- .mlsim_default_coefficients(setup)
  if (is.null(coefficients)) {
    return(default)
  }
  if (is.vector(coefficients) && length(setup$response) == 1L) {
    if (is.null(names(coefficients))) {
      .mlsim_stop("`coefficients` must be named.")
    }
    coefficients <- matrix(coefficients, ncol = 1L, dimnames = list(names(coefficients), setup$response))
  } else {
    coefficients <- as.matrix(coefficients)
  }
  if (!is.numeric(coefficients) || anyNA(coefficients) || any(!is.finite(coefficients))) {
    .mlsim_stop("`coefficients` must contain finite numbers.")
  }
  if (!identical(rownames(coefficients), rownames(default)) ||
      !identical(colnames(coefficients), colnames(default))) {
    .mlsim_stop("`coefficients` names must match `outcome_template()` exactly.")
  }
  coefficients
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_check_scale_coefficients <- function(scale_coefficients, setup) {
  default <- .mlsim_default_scale_coefficients(setup)
  if (!isTRUE(setup$scale$enabled)) {
    if (!is.null(scale_coefficients)) {
      .mlsim_stop("`scale_coefficients` requires `scale_formula`.")
    }
    return(NULL)
  }
  if (is.null(scale_coefficients)) {
    return(default)
  }
  if (is.vector(scale_coefficients) && length(setup$response) == 1L) {
    if (is.null(names(scale_coefficients))) {
      .mlsim_stop("`scale_coefficients` must be named.")
    }
    scale_coefficients <- matrix(scale_coefficients, ncol = 1L, dimnames = list(names(scale_coefficients), setup$response))
  } else {
    scale_coefficients <- as.matrix(scale_coefficients)
  }
  if (!is.numeric(scale_coefficients) || anyNA(scale_coefficients) || any(!is.finite(scale_coefficients))) {
    .mlsim_stop("`scale_coefficients` must contain finite numbers.")
  }
  if (!identical(rownames(scale_coefficients), rownames(default)) ||
      !identical(colnames(scale_coefficients), colnames(default))) {
    .mlsim_stop("`scale_coefficients` names must match `outcome_template()` exactly.")
  }
  scale_coefficients
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_check_named_cov <- function(cov, names, arg) {
  cov <- .mlsim_as_matrix(cov, length(names), arg)
  if (is.null(rownames(cov)) && is.null(colnames(cov)) && length(names) > 1L) {
    .mlsim_stop("`%s` row and column names must match the outcome names.", arg)
  }
  if (!is.null(rownames(cov)) || !is.null(colnames(cov))) {
    row_names <- rownames(cov)
    col_names <- colnames(cov)
    if (is.null(row_names) || is.null(col_names) ||
        anyNA(row_names) || anyNA(col_names) ||
        row_names %any==% "" || col_names %any==% "" ||
        anyDuplicated(row_names) || anyDuplicated(col_names) ||
        !setequal(row_names, names) || !setequal(col_names, names)) {
      .mlsim_stop("`%s` row and column names must match the outcome names.", arg)
    }
    cov <- cov[names, names, drop = FALSE]
  }
  dimnames(cov) <- list(names, names)
  .mlsim_check_cov(cov, length(names), arg)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_check_residual_cor <- function(cor, names) {
  cor <- .mlsim_check_named_cov(cor, names, "residual_cor")
  if (any(abs(diag(cor) - 1) > 1e-8)) {
    .mlsim_stop("`residual_cor` must have ones on the diagonal.")
  }
  if (any(cor < -1 - 1e-8 | cor > 1 + 1e-8)) {
    .mlsim_stop("`residual_cor` values must be correlations in [-1, 1].")
  }
  cor
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_check_random_cov <- function(random_cov, setup) {
  default <- .mlsim_default_random_cov(setup)
  if (length(default) == 0L) {
    return(list())
  }
  if (is.null(random_cov)) {
    .mlsim_stop(
      "`random_cov` is required when outcome or scale formulas contain random-effect terms; use `outcome_template()` to create a named covariance template."
    )
  }
  if (is.matrix(random_cov) && length(default) == 1L) {
    random_cov <- setNames(list(random_cov), names(default))
  }
  random_cov_names <- names(random_cov)
  if (!is.list(random_cov) ||
      is.null(random_cov_names) ||
      anyNA(random_cov_names) ||
      random_cov_names %any==% "" ||
      anyDuplicated(random_cov_names) ||
      !setequal(random_cov_names, names(default))) {
    .mlsim_stop("`random_cov` must be a named list matching `outcome_template()`.")
  }
  random_cov <- random_cov[names(default)]
  for (key in names(default)) {
    expected <- rownames(default[[key]])
    block <- .mlsim_as_matrix(random_cov[[key]], length(expected), "random_cov")
    if (is.null(rownames(block)) && is.null(colnames(block)) && length(expected) > 1L) {
      .mlsim_stop("`random_cov` block `%s` names must match `outcome_template()` exactly.", key)
    }
    if (!is.null(rownames(block)) || !is.null(colnames(block))) {
      row_names <- rownames(block)
      col_names <- colnames(block)
      if (is.null(row_names) || is.null(col_names) ||
          anyNA(row_names) || anyNA(col_names) ||
          row_names %any==% "" || col_names %any==% "" ||
          anyDuplicated(row_names) || anyDuplicated(col_names) ||
          !setequal(row_names, expected) || !setequal(col_names, expected)) {
        .mlsim_stop("`random_cov` block `%s` names must match `outcome_template()` exactly.", key)
      }
      block <- block[expected, expected, drop = FALSE]
    }
    dimnames(block) <- list(expected, expected)
    random_cov[[key]] <- .mlsim_check_cov(block, length(expected), "random_cov")
  }
  random_cov
}
