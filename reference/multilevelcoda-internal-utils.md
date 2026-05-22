# Internal Validation, Recycling, and Random-Parameter Helpers

Validate user inputs, recycle scalar parameters, draw
covariance-structured random values, preserve random seeds, and assemble
metadata shared by generator implementations.

## Usage

``` r
.mlsim_match_level(level)

.mlsim_stop(..., call. = FALSE)

.mlsim_format_args(args)

.mlsim_warn_ignored_args(generator, used, ignored)

.mlsim_check_multilevel_unused(level, supplied, generator, use)

.mlsim_check_multilevel_required(level, supplied, generator)

.mlsim_check_vars(vars, n = NULL, arg = "vars", allow_empty = FALSE)

.mlsim_check_new_names_compatible(new, existing, arg = "generated names")

.mlsim_check_finite_numeric(x, arg)

.mlsim_check_prob(x, arg = "prob")

.mlsim_check_open_prob(x, arg = "prob")

.mlsim_check_categories(categories)

.mlsim_reference_index(reference, labels)

.mlsim_check_categorical_ordered(ordered, output)

.mlsim_check_categorical_logits(logits, nonreference, arg)

.mlsim_check_categorical_random_cov(random_cov, nonreference)

.mlsim_uniform_category_prob(n, labels)

.mlsim_resolve_category_prob(prob, n, labels, arg = "prob")

.mlsim_check_category_prob_matrix(prob, arg = "prob")

.mlsim_baseline_category_prob(eta, reference_index, labels)

.mlsim_sample_categories(prob)

.mlsim_encode_categories(codes, labels, output, ordered)

.mlsim_category_code_map(labels)

.mlsim_check_positive(x, arg)

.mlsim_check_nonnegative(x, arg)

.mlsim_check_scalar_number(x, arg)

.mlsim_check_scalar_logical(x, arg)

.mlsim_check_whole_nonnegative(x, arg)

.mlsim_recycle(x, n, arg)

.mlsim_align_named_vector(x, expected_names, arg)

.mlsim_recycle_integer(x, n, arg)

.mlsim_as_matrix(x, n, arg)

.mlsim_align_square_matrix_dimnames(x, expected_names, arg)

.mlsim_check_named_square_cov(
  cov,
  expected_names,
  arg,
  default = c("identity", "zero"),
  allow_scalar = TRUE,
  allow_zero = TRUE
)

.mlsim_check_cov(x, n, arg, allow_zero = TRUE)

.mlsim_rmvnorm(n, mean, cov)

.mlsim_with_seed(seed, expr)

.mlsim_n_for_level(level, context)

.mlsim_expand_level2(values, level, context)

.mlsim_expand_level2_vector(values, level, context)

.mlsim_check_multilevel_only_args(level, supplied)

.mlsim_random_intercept_names(vars, scale = FALSE)

.mlsim_check_position_cov(cov, names, arg = "random_cov")

.mlsim_joint_random_intercepts(vars, random_cov, context, scale = FALSE)

.mlsim_check_scale_fixed_intercept(
  scale_fixed_intercept,
  n = 1L,
  expected_names = NULL
)

.mlsim_random_eta(
  type,
  fixed_intercept,
  random_cov,
  context,
  vars,
  scale_fixed_intercept = NULL
)

.mlsim_check_simple_residual_cor(residual_cor, vars)

.mlsim_mvn_row_residual(residual_sd, residual_cor)

.mlsim_mvn_row_covariance_array(residual_sd, residual_cor)

.mlsim_check_scale_or_arg(arg, arg_name, scale_fixed_intercept, generator)

.mlsim_check_positive_or_null(x, arg)

.mlsim_simple_random_metadata(
  fixed_intercept,
  random,
  scale_fixed_intercept = NULL,
  scale_linear_predictor = NULL
)

.mlsim_validate_scale_multilevel_args(
  level,
  scale_fixed_intercept,
  random_cov = NULL,
  residual_cor = NULL
)

.mlsim_mean_from_eta(type, eta)

.mlsim_parameter_metadata(level, context, ...)

.mlsim_resolve_size(size, n, context)

.mlsim_result(values, names, metadata)

.mlsim_custom_result(result, vars, level, context, wrapper_metadata)

.mlsim_draw_count_distribution(spec, n)

x %||% y
```

## Arguments

- level:

  Simulation level.

- ...:

  Values passed to formatting or error helpers.

- call.:

  Logical flag passed to [`stop()`](https://rdrr.io/r/base/stop.html).

- args, used, ignored, supplied, use:

  Character vectors or named logicals used to build diagnostics.

- generator:

  Generator label used in diagnostic messages.

- vars, n, names, arg, arg_name:

  Variable names, target lengths, or argument labels.

- allow_empty:

  Logical flag allowing an empty name vector.

- new, existing:

  Candidate and existing generated names.

- x, y:

  Generic vectors, matrices, formulas, or fallback values.

- categories, reference, labels:

  Category values and reference labels.

- ordered, output:

  Categorical output controls.

- logits, nonreference:

  Baseline-category logit parameters.

- random_cov, cov, residual_cor, residual_sd:

  Covariance, correlation, or standard-deviation inputs.

- prob, eta, reference_index, codes:

  Category probabilities, linear predictors, and sampled category codes.

- expected_names:

  Optional names used to align named vector and matrix inputs.

- allow_zero:

  Logical flag allowing positive semi-definite covariance matrices.

- mean:

  Mean vector for multivariate normal draws.

- seed:

  Optional random seed.

- expr:

  Expression evaluated under a temporary seed.

- context:

  Simulation context supplied by
  [`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

- values, result:

  Generated values or raw custom-generator result.

- scale:

  Logical flag for joint location/scale random effects.

- type:

  Distribution type.

- fixed_intercept, scale_fixed_intercept:

  Location and scale intercepts.

- random:

  Random-effect draw metadata.

- scale_linear_predictor:

  Row-level scale linear predictor.

- size:

  Count-distribution size specification.

- metadata, wrapper_metadata:

  Metadata lists.

- spec:

  Count-distribution specification list.

## Value

Internal helper outputs used by generator implementations.

## Examples

``` r
multilevelcoda:::.mlsim_recycle(1, 3, "x")
#> [1] 1 1 1
multilevelcoda:::.mlsim_check_cov(diag(2), 2, "cov")
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
multilevelcoda:::.mlsim_mean_from_eta("poisson", log(2))
#> [1] 2
```
