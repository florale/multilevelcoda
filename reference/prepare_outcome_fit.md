# Prepare Simulated Outcomes for Model Fitting

Build a fitting data set and formula metadata from an `mlsim_data`
object that contains an outcome generator.

## Usage

``` r
prepare_outcome_fit(
  sim,
  outcome = NULL,
  target = c("auto", "generic", "brmcoda"),
  drop_incomplete = TRUE,
  ...
)
```

## Arguments

- sim:

  An `mlsim_data` object returned by
  [`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

- outcome:

  Optional name of the outcome generator to prepare. Required when `sim`
  contains more than one outcome generator.

- target:

  Fitting target. `"generic"` returns ordinary data and formulas;
  `"brmcoda"` prepares compositional outcomes with
  [`complr()`](https://florale.github.io/multilevelcoda/reference/complr.md);
  `"auto"` chooses `"brmcoda"` for compositional outcomes and
  `"generic"` otherwise.

- drop_incomplete:

  Logical; if `TRUE`, remove rows with missing values in the generated
  fitting formulas.

- ...:

  Additional arguments passed to
  [`complr()`](https://florale.github.io/multilevelcoda/reference/complr.md)
  when `target = "brmcoda"`.

## Value

An `mlsim_fit_prep` object containing fitting data, formulas,
helper-column names, term maps, residual-correlation metadata, and
target specific objects such as `complr` for `brmcoda`.

## Details

Outcome formulas can contain simulation helpers such as `lag1()`,
[`within()`](https://rdrr.io/r/base/with.html), `between()`, and
`ar1()`. `prepare_outcome_fit()` maps those helpers to concrete columns
suitable for fitting and records a term map from simulation terms to
fitting terms.

## Examples

``` r
sim <- simulate_data(
  n_groups = 2,
  n_per_group = 3,
  seed = 20,
  generators = list(
    x = gen_normal("x"),
    y = gen_outcome(y ~ lag1(x), residual_cov = matrix(0, 1, 1, dimnames = list("y", "y")))
  )
)
prep <- prepare_outcome_fit(sim)
prep$formula
#> y ~ lag_x
#> <environment: 0x55663ede2518>
prep$helper_columns
#> [1] "lag_x"
```
