#' Generate Gaussian Outcomes from a Formula
#'
#' Create outcome generator specifications and parameter templates for
#' Gaussian, multivariate Gaussian, and compositional outcomes.
#'
#' @param formula Two-sided outcome formula. The left-hand side is a response
#'   name or `cbind()`/`mvbind()` response vector. The right-hand side may use
#'   ordinary model terms, random-effect terms such as `(1 | group_id)` or
#'   `(1 + x | id | group_id)`, and simulation helpers `lag1(x)`,
#'   `within(x)`, `between(x)`, and `ar1()`.
#' @param coefficients Optional coefficient vector or matrix matching the
#'   template returned by `outcome_template()`. If omitted, coefficients default
#'   to zero.
#' @param residual_cov Residual covariance matrix for outcomes without
#'   `scale_formula`.
#' @param random_cov Named random-effect covariance block or list of blocks
#'   matching `outcome_template()`. Required when `formula` or `scale_formula`
#'   contains random-effect terms.
#' @param compositional Logical; if `TRUE`, treat multivariate responses as ILR
#'   coordinates and emit closed composition parts.
#' @param parts Character vector naming composition parts. Must have one more
#'   entry than the number of ILR coordinates unless supplied by `sbp`.
#' @param sbp Optional sequential binary partition matrix.
#' @param total Positive scalar total for closed compositions.
#' @param keep_ilr Logical; if `TRUE`, emit both ILR coordinates and parts. If
#'   `FALSE`, emit only parts.
#' @param burnin `"auto"` or a non-negative integer. Used to initialize
#'   autoregressive outcome terms before the first observed row in each unit.
#' @param scale_formula Optional one-sided formula, or `sigma ~ ...`, defining
#'   a log residual standard-deviation model.
#' @param scale_coefficients Optional coefficient vector or matrix for
#'   `scale_formula`, matching `outcome_template()`.
#' @param residual_cor Residual correlation matrix used with `scale_formula`.
#'
#' @return `gen_outcome()` returns an `mlsim_generator_spec` for use in
#'   [simulate_data()]. `outcome_template()` returns a generator specification
#'   that emits no columns but records default coefficient, covariance, formula,
#'   helper-column, and fitting metadata in `generator_metadata`.
#'
#' @details
#' `outcome_template()` is the safest way to discover the exact coefficient and
#' covariance names required by `gen_outcome()`. Outcome helper terms are
#' evaluated during simulation and are mapped to concrete helper columns for
#' fitting by [prepare_outcome_fit()].
#'
#' @examples
#' template <- simulate_data(
#'   n_groups = 2,
#'   n_per_group = 3,
#'   generators = list(
#'     x = gen_normal("x"),
#'     y_template = outcome_template(y ~ x + (1 | group_id))
#'   )
#' )$generator_metadata$y_template
#'
#' coefficients <- template$coefficients
#' coefficients["(Intercept)", "y"] <- 1
#' coefficients["x", "y"] <- 0.5
#'
#' sim <- simulate_data(
#'   n_groups = 2,
#'   n_per_group = 3,
#'   seed = 30,
#'   generators = list(
#'     x = gen_normal("x"),
#'     y = gen_outcome(
#'       y ~ x + (1 | group_id),
#'       coefficients = coefficients,
#'       residual_cov = template$residual_cov,
#'       random_cov = template$random_cov
#'     )
#'   )
#' )
#' sim$data
#'
#' @family outcome generators
#' @name outcome-generators
NULL

#' @rdname outcome-generators
#' @export
gen_outcome <- function(formula, coefficients = NULL, residual_cov = NULL,
                        random_cov = NULL, compositional = FALSE,
                        parts = NULL, sbp = NULL, total = 1,
                        keep_ilr = TRUE, burnin = "auto",
                        scale_formula = NULL, scale_coefficients = NULL,
                        residual_cor = NULL) {
  formula <- .mlsim_check_outcome_formula(formula)
  response <- .mlsim_outcome_response(formula)
  composition <- .mlsim_outcome_composition(response, compositional, parts, sbp, total, keep_ilr)
  burnin <- .mlsim_check_burnin(burnin)
  scale_formula <- .mlsim_check_scale_formula(scale_formula)

  simulate <- function(context) {
    setup <- .mlsim_outcome_setup(formula, context, composition, burnin, scale_formula)
    .mlsim_simulate_outcome(
      setup, coefficients, residual_cov, random_cov,
      scale_coefficients = scale_coefficients,
      residual_cor = residual_cor
    )
  }

  .mlsim_spec(
    "outcome", composition$output_names, "single", simulate,
    formula = formula,
    compositional = composition$compositional,
    parts = composition$parts,
    sbp = composition$sbp,
    total = total,
    keep_ilr = composition$keep_ilr,
    burnin = burnin,
    scale_formula = scale_formula
  )
}

#' @rdname outcome-generators
#' @export
outcome_template <- function(formula, compositional = FALSE, parts = NULL,
                             sbp = NULL, total = 1, keep_ilr = TRUE,
                             burnin = "auto", scale_formula = NULL,
                             residual_cor = NULL) {
  formula <- .mlsim_check_outcome_formula(formula)
  response <- .mlsim_outcome_response(formula)
  composition <- .mlsim_outcome_composition(response, compositional, parts, sbp, total, keep_ilr)
  burnin <- .mlsim_check_burnin(burnin)
  scale_formula <- .mlsim_check_scale_formula(scale_formula)

  simulate <- function(context) {
    setup <- .mlsim_outcome_setup(formula, context, composition, burnin, scale_formula)
    metadata <- .mlsim_outcome_template_metadata(setup)
    if (isTRUE(setup$scale$enabled)) {
      metadata$residual_cor <- .mlsim_check_residual_cor(
        residual_cor %||% .mlsim_default_residual_cor(setup),
        setup$response
      )
    } else if (!is.null(residual_cor)) {
      .mlsim_stop("`residual_cor` requires `scale_formula`.")
    }
    list(data = data.table::data.table(), names = character(), metadata = metadata)
  }

  .mlsim_spec(
    "outcome_template", character(), "single", simulate,
    formula = formula,
    compositional = composition$compositional,
    parts = composition$parts,
    sbp = composition$sbp,
    total = total,
    keep_ilr = composition$keep_ilr,
    burnin = burnin,
    scale_formula = scale_formula
  )
}

#' Internal Outcome Simulation Helpers
#'
#' Parse outcome formulas, transform helper terms, build model matrices,
#' validate coefficient and covariance inputs, simulate residuals and random
#' effects, handle autoregressive state, and create fitting metadata.
#'
#' @param formula,scale_formula Formula objects used for simulation or fitting.
#' @param coefficients,scale_coefficients Coefficient inputs or validated
#'   coefficient matrices.
#' @param residual_cov,residual_cor Residual covariance or correlation inputs.
#' @param random_cov Random-effect covariance input.
#' @param compositional Logical flag for composition outputs.
#' @param parts,response Character vectors naming composition parts or response
#'   coordinates.
#' @param sbp Sequential binary partition matrix.
#' @param total Composition total.
#' @param keep_ilr Logical flag for emitting ILR coordinates.
#' @param burnin Autoregressive burn-in setting.
#' @param context Simulation context supplied by [simulate_data()].
#' @param composition Resolved composition metadata.
#' @param setup Outcome setup object.
#' @param state Mutable outcome parsing state.
#' @param expr,left,right,term Formula expressions or terms.
#' @param allow_ar1 Logical flag controlling `ar1()` support.
#' @param random_terms,location_specs,scale_specs Random-effect specifications.
#' @param specs,spec Internal random-effect specification lists.
#' @param env Formula environment.
#' @param random,effects,z_by_key Random-effect draws and model matrices.
#' @param residual,residual_info Residual draws and residual metadata.
#' @param ar Autoregressive transition metadata.
#' @param base,transition,zero_data Intermediate autoregressive linearity-check
#'   objects.
#' @param data,row_data Data frames used for model-matrix evaluation.
#' @param fixed_expr,fixed_names,fit_formula Fixed-effect expressions, names,
#'   and transformed fitting formulas.
#' @param replacements Named symbol replacement vector.
#' @param names,name,columns,arg,type,source,internal,column,prefix Helper names
#'   and diagnostic labels.
#' @param actual,expected Matrices compared for internal linearity checks.
#' @param rows Row-selection vector.
#' @param tolerance Numeric comparison tolerance.
#' @param row,j,n,d,value,values,first_rows Indices, dimensions, values, or row
#'   selections used by simulation helpers.
#' @param x,cov,cor,contrasts.arg,xlev Matrix, formula, vector, expression, or
#'   model-matrix controls, depending on helper.
#' @param var Variable name used in within/between transformations.
#' @param sim_term,response_derived Fit-helper metadata.
#'
#' @return Internal helper outputs used by [gen_outcome()],
#'   [outcome_template()], and [prepare_outcome_fit()].
#'
#' @examples
#' multilevelcoda:::.mlsim_outcome_response(mvbind(y1, y2) ~ x)
#' multilevelcoda:::.mlsim_check_scale_formula(sigma ~ x)
#' multilevelcoda:::.mlsim_deparse_one(quote(x + y))
#'
#' @keywords internal
#' @name multilevelcoda-internal-outcome
NULL
