# Outcome setup and simulation orchestration internals.

#' @rdname multilevelcoda-internal-outcome
.mlsim_outcome_setup <- function(formula, context, composition, burnin, scale_formula = NULL) {
  if (!is.null(context$n_groups) && context$group_id %any!in% names(context$data)) {
    .mlsim_stop("Grouped outcome simulation requires `%s` in `context$data`.", context$group_id)
  }
  parsed <- .mlsim_extract_random_terms(formula[[3L]])
  state <- .mlsim_outcome_state(context$data, context, composition$response)
  fixed_expr <- .mlsim_transform_outcome_expr(parsed$fixed_expr %||% 1, state)$expr
  fixed_formula <- .mlsim_rhs_formula(fixed_expr, environment(formula))
  fixed_matrix <- .mlsim_model_matrix(fixed_formula, data = state$data)
  fixed_xlevels <- .mlsim_model_xlevels(fixed_formula, state$data)
  fixed_contrasts <- attr(fixed_matrix, "contrasts")
  fixed_names <- .mlsim_pretty_design_names(colnames(fixed_matrix), state)

  random_specs <- .mlsim_random_specs(parsed$random, state, environment(formula), composition$response)
  scale <- .mlsim_scale_setup(scale_formula, state, composition$response)
  brms_formula <- .mlsim_brms_formula(formula, parsed$fixed_expr %||% 1, parsed$random, state)
  fit_formula <- .mlsim_fit_formula(formula, parsed$fixed_expr %||% 1, parsed$random, state)
  fit_fixed_formula <- .mlsim_fit_fixed_formula(formula, parsed$fixed_expr %||% 1, state)
  fit_term_map <- .mlsim_fit_term_map(fixed_names, parsed$fixed_expr %||% 1, state, environment(formula))
  fit_helpers <- .mlsim_fit_helpers_dataframe(state)
  helper_columns <- fit_helpers$column %||% character()
  duplicated_helpers <- helper_columns[duplicated(helper_columns)]
  if (length(duplicated_helpers) > 0L) {
    .mlsim_stop(
      "Fit helper columns must be unique; duplicated names: %s.",
      paste(unique(duplicated_helpers), collapse = ", ")
    )
  }
  complete_rows <- stats::complete.cases(fixed_matrix)
  for (spec in random_specs) {
    complete_rows <- complete_rows & stats::complete.cases(spec$z)
  }
  if (isTRUE(scale$enabled)) {
    complete_rows <- complete_rows & stats::complete.cases(scale$fixed_matrix)
    for (spec in scale$random_specs) {
      complete_rows <- complete_rows & stats::complete.cases(spec$z)
    }
  }

  setup <- list(
    formula = formula,
    brms_formula = brms_formula,
    fit_formula = fit_formula,
    fit_fixed_formula = fit_fixed_formula,
    fit_term_map = fit_term_map,
    fit_helpers = fit_helpers,
    helper_columns = helper_columns,
    context = context,
    data = state$data,
    response = composition$response,
    output_names = composition$output_names,
    fixed_formula = fixed_formula,
    fixed_matrix = fixed_matrix,
    fixed_xlevels = fixed_xlevels,
    fixed_contrasts = fixed_contrasts,
    fixed_names = fixed_names,
    complete_rows = complete_rows,
    random_specs = random_specs,
    joint_random_specs = .mlsim_joint_random_specs(random_specs, scale$random_specs),
    ar_internal = state$ar_internal,
    ar_names = state$ar_names,
    scale = scale,
    composition = composition,
    burnin = burnin
  )
  setup$ar_info <- .mlsim_ar_info(setup)
  setup
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_simulate_outcome <- function(setup, coefficients, residual_cov, random_cov,
                                    scale_coefficients = NULL, residual_cor = NULL) {
  coefficients <- .mlsim_check_outcome_coefficients(coefficients, setup)
  random_cov <- .mlsim_check_random_cov(random_cov, setup)
  random <- .mlsim_random_contribution(setup, random_cov)
  residual_info <- .mlsim_residual_components(
    setup = setup,
    residual_cov = residual_cov,
    residual_cor = residual_cor,
    scale_coefficients = scale_coefficients,
    random = random
  )
  ar <- NULL
  if (length(setup$ar_internal) > 0L) {
    ar <- .mlsim_ar_components(setup, coefficients, random)
  }
  .mlsim_validate_ar_stability(setup, ar)

  residual <- residual_info$residual

  if (length(setup$ar_internal) > 0L) {
    values <- .mlsim_simulate_ar_outcome(setup, coefficients, random, residual, residual_info, ar)
  } else {
    eta <- setup$fixed_matrix %*% coefficients + random$contribution
    values <- eta + residual
    values[!setup$complete_rows, ] <- NA_real_
  }
  colnames(values) <- setup$response

  metadata <- .mlsim_outcome_template_metadata(setup)
  metadata$coefficients <- coefficients
  metadata$residual_cov <- residual_info$residual_cov
  metadata$residual_cor <- residual_info$residual_cor
  metadata$scale_coefficients <- residual_info$scale_coefficients
  metadata$scale_linear_predictor <- residual_info$scale_linear_predictor
  metadata$residual_sd <- residual_info$residual_sd
  metadata$random_cov <- random_cov
  metadata$random_effects <- random$effects
  metadata$complete_rows <- setup$complete_rows

  out_values <- values
  if (isTRUE(setup$composition$compositional)) {
    comp <- .mlsim_ilr_inverse(
      values,
      parts = setup$composition$parts,
      total = setup$composition$total,
      coordinate_names = setup$response,
      sbp = setup$composition$sbp
    )
    out_values <- if (isTRUE(setup$composition$keep_ilr)) cbind(values, comp$values) else comp$values
    metadata$basis <- comp$basis
    metadata$ilr_coordinate_map <- comp$coordinate_map
    metadata$total <- comp$total
  }
  fit_values <- .mlsim_fit_helper_values(setup, values)
  if (ncol(fit_values) > 0L) {
    repeated_helpers <- intersect(colnames(fit_values), colnames(out_values))
    if (length(repeated_helpers) > 0L) {
      .mlsim_stop(
        "Fit helper columns would overwrite generated outcome columns: %s.",
        paste(repeated_helpers, collapse = ", ")
      )
    }
    out_values <- cbind(out_values, fit_values)
  }
  .mlsim_result(out_values, colnames(out_values), metadata)
}


