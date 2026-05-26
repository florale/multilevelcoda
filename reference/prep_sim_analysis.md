# Prepare Simulated Outcome Data for Analysis

Convert an
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md)
result with a
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
generator into an analysis-ready data set and inferred `brms` formula.
The helper recomputes observed between- and within-person predictors
from the raw simulated columns, creates within-person centered lag
columns for `ar1()` terms, and creates a
[`complr()`](https://florale.github.io/multilevelcoda/reference/complr.md)
object when the outcome is compositional.

## Usage

``` r
prep_sim_analysis(sim, outcome = NULL, drop_lag_na = FALSE)
```

## Arguments

- sim:

  An `mlsim_data` object returned by
  [`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

- outcome:

  Optional character scalar naming the outcome generator. When `NULL`,
  the helper uses the only generator with `distribution = "outcome"`.

- drop_lag_na:

  Logical scalar. When `FALSE`, the default, first rows in each series
  are retained and lag-derived columns are left as `NA`. When `TRUE`,
  rows with missing lag-derived predictors are removed.

## Value

An `mlsim_analysis` object, a list with:

- `data`:

  Analysis data with derived between/within and lag columns.

- `formula`:

  An inferred `brms` formula.

- `complr`:

  A `complr` object for compositional outcomes, otherwise `NULL`.

- `metadata`:

  Preparation metadata, including derived column names and formula
  mappings.

## Details

The analysis formula is inferred from the stored
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
formula. Simulation terms `between(x)` and `within(x)` become
observed-data columns named `x_between` and `x_within`, computed from
`x` by the simulation group identifier. These columns are recomputed
even when columns with the same names already exist in `sim$data`.

For dynamic formulas, `ar1()` is translated to within-person centered
lag predictors. For compositional outcomes, the helper rebuilds the ILR
coordinates through
[`complr()`](https://florale.github.io/multilevelcoda/reference/complr.md)
using the simulator's parts and SBP metadata, then lags the generated
`z` coordinates used by
[`brmcoda()`](https://florale.github.io/multilevelcoda/reference/brmcoda.md).

## Examples

``` r
params <- list(
  location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
  scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y")))
)
sim <- simulate_data(
  n = 5,
  seed = 1,
  generators = list(
    outcome = gen_outcome(
      y ~ 1,
      scale = sigma ~ 1,
      params = params,
      burnin = 0
    )
  )
)
analysis <- prep_sim_analysis(sim)
analysis$formula
#> y ~ 1 
#> sigma ~ 1
```
