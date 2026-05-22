# Outcome residual and scale-model internals.

#' @rdname multilevelcoda-internal-outcome
.mlsim_scale_setup <- function(scale_formula, state, response) {
  if (is.null(scale_formula)) {
    return(list(
      enabled = FALSE,
      formula = NULL,
      fit_formula = NULL,
      fit_term_map = NULL,
      helper_columns = character(),
      fixed_matrix = NULL,
      fixed_names = character(),
      random_specs = list()
    ))
  }

  rhs <- .mlsim_scale_rhs(scale_formula)
  parsed <- .mlsim_extract_random_terms(rhs)
  fixed_expr <- .mlsim_transform_outcome_expr(parsed$fixed_expr %||% 1, state, allow_ar1 = FALSE)$expr
  fixed_formula <- .mlsim_rhs_formula(fixed_expr, environment(scale_formula))
  fixed_matrix <- .mlsim_model_matrix(fixed_formula, data = state$data)
  fixed_xlevels <- .mlsim_model_xlevels(fixed_formula, state$data)
  fixed_contrasts <- attr(fixed_matrix, "contrasts")
  fixed_names <- .mlsim_pretty_design_names(colnames(fixed_matrix), state)
  random_specs <- .mlsim_scale_random_specs(parsed$random, state, environment(scale_formula), response)
  fit_formula <- .mlsim_scale_fit_formula(fixed_expr, parsed$random, state, environment(scale_formula))
  fit_term_map <- .mlsim_scale_fit_term_map(fixed_names, fixed_expr, state, environment(scale_formula))
  helper_columns <- .mlsim_scale_helper_columns(fit_formula, state)

  list(
    enabled = TRUE,
    formula = scale_formula,
    fit_formula = fit_formula,
    fit_term_map = fit_term_map,
    helper_columns = helper_columns,
    fixed_formula = fixed_formula,
    fixed_matrix = fixed_matrix,
    fixed_xlevels = fixed_xlevels,
    fixed_contrasts = fixed_contrasts,
    fixed_names = fixed_names,
    random_specs = random_specs
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_scale_rhs <- function(scale_formula) {
  if (length(scale_formula) == 2L) {
    scale_formula[[2L]]
  } else {
    scale_formula[[3L]]
  }
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_residual_components <- function(setup, residual_cov, residual_cor,
                                       scale_coefficients, random) {
  if (!isTRUE(setup$scale$enabled)) {
    if (!is.null(scale_coefficients)) {
      .mlsim_stop("`scale_coefficients` requires `scale_formula`.")
    }
    if (!is.null(residual_cor)) {
      .mlsim_stop("`residual_cor` requires `scale_formula`.")
    }
    residual_cov <- .mlsim_check_named_cov(
      residual_cov %||% .mlsim_default_residual_cov(setup),
      setup$response,
      "residual_cov"
    )
    residual <- .mlsim_rmvnorm(nrow(setup$data), rep(0, length(setup$response)), residual_cov)
    colnames(residual) <- setup$response
    return(list(
      type = "global",
      residual = residual,
      residual_cov = residual_cov,
      residual_cor = NULL,
      scale_coefficients = NULL,
      scale_linear_predictor = NULL,
      residual_sd = NULL
    ))
  }

  if (!is.null(residual_cov)) {
    .mlsim_stop("`residual_cov` cannot be supplied with `scale_formula`; use `residual_cor` and `scale_coefficients`.")
  }
  scale_coefficients <- .mlsim_check_scale_coefficients(scale_coefficients, setup)
  residual_cor <- .mlsim_check_residual_cor(
    residual_cor %||% .mlsim_default_residual_cor(setup),
    setup$response
  )
  scale_eta <- setup$scale$fixed_matrix %*% scale_coefficients + random$scale_contribution
  colnames(scale_eta) <- setup$response
  residual_sd <- exp(scale_eta)
  colnames(residual_sd) <- setup$response
  residual <- .mlsim_rmvnorm_with_row_sd(residual_sd, residual_cor)
  colnames(residual) <- setup$response

  list(
    type = "scale",
    residual = residual,
    residual_cov = NULL,
    residual_cor = residual_cor,
    scale_coefficients = scale_coefficients,
    scale_linear_predictor = scale_eta,
    residual_sd = residual_sd
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_burnin_residual_draw <- function(residual_info, row) {
  d <- ncol(residual_info$residual)
  if (identical(residual_info$type, "scale")) {
    draw <- as.numeric(.mlsim_rmvnorm(1L, rep(0, d), residual_info$residual_cor))
    return(draw * residual_info$residual_sd[row, ])
  }
  as.numeric(.mlsim_rmvnorm(1L, rep(0, d), residual_info$residual_cov))
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_rmvnorm_with_row_sd <- function(residual_sd, residual_cor) {
  residual_sd <- as.matrix(residual_sd)
  base <- .mlsim_rmvnorm(nrow(residual_sd), rep(0, ncol(residual_sd)), residual_cor)
  base * residual_sd
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_default_residual_cov <- function(setup) {
  cov <- diag(length(setup$response))
  rownames(cov) <- colnames(cov) <- setup$response
  cov
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_default_residual_cor <- function(setup) {
  cor <- diag(length(setup$response))
  rownames(cor) <- colnames(cor) <- setup$response
  cor
}


