# Generate Gaussian Outcomes from a Formula

Create outcome generator specifications and parameter templates for
Gaussian, multivariate Gaussian, and compositional outcomes.

## Usage

``` r
gen_outcome(
  formula,
  coefficients = NULL,
  residual_cov = NULL,
  random_cov = NULL,
  compositional = FALSE,
  parts = NULL,
  sbp = NULL,
  total = 1,
  keep_ilr = TRUE,
  burnin = "auto",
  scale_formula = NULL,
  scale_coefficients = NULL,
  residual_cor = NULL
)

outcome_template(
  formula,
  compositional = FALSE,
  parts = NULL,
  sbp = NULL,
  total = 1,
  keep_ilr = TRUE,
  burnin = "auto",
  scale_formula = NULL,
  residual_cor = NULL
)
```

## Arguments

- formula:

  Two-sided outcome formula. The left-hand side is a response name or
  [`cbind()`](https://rdrr.io/r/base/cbind.html)/`mvbind()` response
  vector. The right-hand side may use ordinary model terms,
  random-effect terms such as `(1 | group_id)` or
  `(1 + x | id | group_id)`, and simulation helpers `lag1(x)`,
  `within(x)`, `between(x)`, and `ar1()`.

- coefficients:

  Optional coefficient vector or matrix matching the template returned
  by `outcome_template()`. If omitted, coefficients default to zero.

- residual_cov:

  Residual covariance matrix for outcomes without `scale_formula`.

- random_cov:

  Named random-effect covariance block or list of blocks matching
  `outcome_template()`. Required when `formula` or `scale_formula`
  contains random-effect terms.

- compositional:

  Logical; if `TRUE`, treat multivariate responses as ILR coordinates
  and emit closed composition parts.

- parts:

  Character vector naming composition parts. Must have one more entry
  than the number of ILR coordinates unless supplied by `sbp`.

- sbp:

  Optional sequential binary partition matrix.

- total:

  Positive scalar total for closed compositions.

- keep_ilr:

  Logical; if `TRUE`, emit both ILR coordinates and parts. If `FALSE`,
  emit only parts.

- burnin:

  `"auto"` or a non-negative integer. Used to initialize autoregressive
  outcome terms before the first observed row in each unit.

- scale_formula:

  Optional one-sided formula, or `sigma ~ ...`, defining a log residual
  standard-deviation model.

- scale_coefficients:

  Optional coefficient vector or matrix for `scale_formula`, matching
  `outcome_template()`.

- residual_cor:

  Residual correlation matrix used with `scale_formula`.

## Value

`gen_outcome()` returns an `mlsim_generator_spec` for use in
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).
`outcome_template()` returns a generator specification that emits no
columns but records default coefficient, covariance, formula,
helper-column, and fitting metadata in `generator_metadata`.

## Details

`outcome_template()` is the safest way to discover the exact coefficient
and covariance names required by `gen_outcome()`. Outcome helper terms
are evaluated during simulation and are mapped to concrete helper
columns for fitting by
[`prepare_outcome_fit()`](https://florale.github.io/multilevelcoda/reference/prepare_outcome_fit.md).

## Examples

``` r
template <- simulate_data(
  n_groups = 2,
  n_per_group = 3,
  generators = list(
    x = gen_normal("x"),
    y_template = outcome_template(y ~ x + (1 | group_id))
  )
)$generator_metadata$y_template

coefficients <- template$coefficients
coefficients["(Intercept)", "y"] <- 1
coefficients["x", "y"] <- 0.5

sim <- simulate_data(
  n_groups = 2,
  n_per_group = 3,
  seed = 30,
  generators = list(
    x = gen_normal("x"),
    y = gen_outcome(
      y ~ x + (1 | group_id),
      coefficients = coefficients,
      residual_cov = template$residual_cov,
      random_cov = template$random_cov
    )
  )
)
sim$data
#>    group_id obs_id          x          y
#>       <int>  <int>      <num>      <num>
#> 1:        1      1 -1.2885182 -0.3141561
#> 2:        1      2 -0.3476894  1.1006750
#> 3:        1      3 -0.5216288 -0.2840864
#> 4:        2      1  1.2734732 -0.1826613
#> 5:        2      2  1.8245206  1.2444705
#> 6:        2      3 -1.5113079  0.1850480
```
