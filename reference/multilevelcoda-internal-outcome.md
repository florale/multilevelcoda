# Internal Outcome Simulation Helpers

Parse outcome formulas, transform helper terms, build model matrices,
validate coefficient and covariance inputs, simulate residuals and
random effects, handle autoregressive state, and create fitting
metadata.

## Usage

``` r
.mlsim_simulate_ar_outcome(
  setup,
  coefficients,
  random,
  residual,
  residual_info,
  ar
)

.mlsim_simulate_ar_outcome_rowwise(
  setup,
  coefficients,
  random,
  residual,
  residual_info
)

.mlsim_burnin_prior_ar(setup, ar, residual_info, first_rows)

.mlsim_burnin_prior(setup, coefficients, random, residual_info, first_rows)

.mlsim_ar_info(setup)

.mlsim_ar_components(setup, coefficients, random)

.mlsim_ar_eta(setup, coefficients, random, data)

.mlsim_ar_data(setup, value)

.mlsim_ar_components_are_linear(
  setup,
  coefficients,
  random,
  zero_data,
  base,
  transition
)

.mlsim_ar_transition_delta(transition, j, n, d)

.mlsim_same_matrix(actual, expected, rows, tolerance = 1e-08)

.mlsim_validate_ar_stability(setup, ar = NULL)

.mlsim_ar_transition_row(ar, row, d, response)

.mlsim_outcome_setup(
  formula,
  context,
  composition,
  burnin,
  scale_formula = NULL
)

.mlsim_simulate_outcome(
  setup,
  coefficients,
  residual_cov,
  random_cov,
  scale_coefficients = NULL,
  residual_cor = NULL
)

.mlsim_outcome_template_metadata(setup)

.mlsim_scale_fit_formula(fixed_expr, random_terms, state, env)

.mlsim_scale_fit_term_map(fixed_names, fixed_expr, state, env)

.mlsim_scale_helper_columns(fit_formula, state)

.mlsim_register_fit_helper(
  state,
  type,
  source,
  internal,
  column,
  sim_term,
  response_derived
)

.mlsim_fit_column_name(type, source)

.mlsim_fit_helpers_dataframe(state)

.mlsim_fit_data_from_state(state)

.mlsim_fit_helper_values(setup, values)

.mlsim_brms_formula(formula, fixed_expr, random_terms, state)

.mlsim_transform_for_brms(expr, state)

.mlsim_fit_formula(formula, fixed_expr, random_terms, state)

.mlsim_fit_fixed_formula(formula, fixed_expr, state)

.mlsim_transform_for_fit(expr, state)

.mlsim_fit_random_term_expr(term, state)

.mlsim_fit_term_map(fixed_names, fixed_expr, state, env)

.mlsim_random_term_expr(term, state)

.mlsim_check_outcome_formula(formula)

.mlsim_check_scale_formula(scale_formula)

.mlsim_outcome_response(formula)

.mlsim_outcome_composition(
  response,
  compositional,
  parts,
  sbp,
  total,
  keep_ilr
)

.mlsim_check_burnin(burnin)

.mlsim_outcome_state(data, context, response)

.mlsim_transform_outcome_expr(expr, state, allow_ar1 = TRUE)

.mlsim_outcome_helper_expr(expr, state)

.mlsim_model_matrix_row(setup, row_data)

.mlsim_model_matrix_data(setup, data)

.mlsim_random_model_matrix(spec, data)

.mlsim_align_model_matrix(x, columns)

.mlsim_model_xlevels(formula, data)

.mlsim_model_matrix(formula, data, contrasts.arg = NULL, xlev = NULL)

.mlsim_model_frame(formula, data, xlev = NULL)

.mlsim_lag1(x, context)

.mlsim_within_group(x, context, var)

.mlsim_between_group(x, context, var)

.mlsim_replace_symbols(expr, replacements)

.mlsim_pretty_design_names(names, state)

.mlsim_internal_name(prefix, name)

.mlsim_rhs_formula(expr, env)

.mlsim_plus_formula(left, right)

.mlsim_sum_symbols(names)

.mlsim_plus_expr(left, right)

.mlsim_is_random_call(expr)

.mlsim_unwrap_parens(expr)

.mlsim_call_name(expr)

.mlsim_deparse_one(expr)

.mlsim_extract_random_terms(expr)

.mlsim_extract_random_terms_rec(expr)

.mlsim_parse_random_term(expr)

.mlsim_random_specs(random_terms, state, env, response)

.mlsim_merge_random_specs(specs)

.mlsim_scale_random_specs(random_terms, state, env, response)

.mlsim_joint_random_specs(location_specs, scale_specs)

.mlsim_random_contribution(setup, random_cov)

.mlsim_random_contribution_data(setup, effects, data)

.mlsim_random_contribution_from_z(setup, effects, z_by_key)

.mlsim_random_contribution_row(setup, effects, row, row_data)

.mlsim_scale_random_contribution(setup, effects)

.mlsim_random_model_matrix_row(spec, row_data)

.mlsim_scale_setup(scale_formula, state, response)

.mlsim_scale_rhs(scale_formula)

.mlsim_residual_components(
  setup,
  residual_cov,
  residual_cor,
  scale_coefficients,
  random
)

.mlsim_burnin_residual_draw(residual_info, row)

.mlsim_rmvnorm_with_row_sd(residual_sd, residual_cor)

.mlsim_default_residual_cov(setup)

.mlsim_default_residual_cor(setup)

.mlsim_default_coefficients(setup)

.mlsim_default_scale_coefficients(setup)

.mlsim_default_random_cov(setup)

.mlsim_check_outcome_coefficients(coefficients, setup)

.mlsim_check_scale_coefficients(scale_coefficients, setup)

.mlsim_check_named_cov(cov, names, arg)

.mlsim_check_residual_cor(cor, names)

.mlsim_check_random_cov(random_cov, setup)
```

## Arguments

- setup:

  Outcome setup object.

- coefficients, scale_coefficients:

  Coefficient inputs or validated coefficient matrices.

- random, effects, z_by_key:

  Random-effect draws and model matrices.

- residual, residual_info:

  Residual draws and residual metadata.

- ar:

  Autoregressive transition metadata.

- data, row_data:

  Data frames used for model-matrix evaluation.

- base, transition, zero_data:

  Intermediate autoregressive linearity-check objects.

- actual, expected:

  Matrices compared for internal linearity checks.

- rows:

  Row-selection vector.

- tolerance:

  Numeric comparison tolerance.

- row, j, n, d, value, values, first_rows:

  Indices, dimensions, values, or row selections used by simulation
  helpers.

- formula, scale_formula:

  Formula objects used for simulation or fitting.

- context:

  Simulation context supplied by
  [`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

- composition:

  Resolved composition metadata.

- burnin:

  Autoregressive burn-in setting.

- residual_cov, residual_cor:

  Residual covariance or correlation inputs.

- random_cov:

  Random-effect covariance input.

- fixed_expr, fixed_names, fit_formula:

  Fixed-effect expressions, names, and transformed fitting formulas.

- random_terms, location_specs, scale_specs:

  Random-effect specifications.

- state:

  Mutable outcome parsing state.

- env:

  Formula environment.

- sim_term, response_derived:

  Fit-helper metadata.

- expr, left, right, term:

  Formula expressions or terms.

- compositional:

  Logical flag for composition outputs.

- parts, response:

  Character vectors naming composition parts or response coordinates.

- sbp:

  Sequential binary partition matrix.

- total:

  Composition total.

- keep_ilr:

  Logical flag for emitting ILR coordinates.

- allow_ar1:

  Logical flag controlling `ar1()` support.

- x, cov, cor, contrasts.arg, xlev:

  Matrix, formula, vector, expression, or model-matrix controls,
  depending on helper.

- var:

  Variable name used in within/between transformations.

- replacements:

  Named symbol replacement vector.

- names, name, columns, arg, type, source, internal, column, prefix:

  Helper names and diagnostic labels.

- specs, spec:

  Internal random-effect specification lists.

## Value

Internal helper outputs used by
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/outcome-generators.md),
[`outcome_template()`](https://florale.github.io/multilevelcoda/reference/outcome-generators.md),
and
[`prepare_outcome_fit()`](https://florale.github.io/multilevelcoda/reference/prepare_outcome_fit.md).

## Examples

``` r
multilevelcoda:::.mlsim_outcome_response(mvbind(y1, y2) ~ x)
#> [1] "y1" "y2"
multilevelcoda:::.mlsim_check_scale_formula(sigma ~ x)
#> sigma ~ x
#> <environment: 0x55663f735cd0>
multilevelcoda:::.mlsim_deparse_one(quote(x + y))
#> [1] "x + y"
```
