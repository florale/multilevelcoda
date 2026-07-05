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

The inferred analysis model is a deliberately *pragmatic default
estimator* built from observed data, not the matched model for the
simulated data-generating process. See the "Pragmatic default estimator"
section before interpreting parameter recovery results.

The analysis formula is inferred from the stored
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
formula. Simulation terms `between(x)` and `within(x)` become
observed-data columns named `x_between` and `x_within`, computed from
`x` by the simulation group identifier. These columns are recomputed
even when columns with the same names already exist in `sim$data`.

For dynamic formulas, `ar1()` is translated to lagged observed response
columns centered at each person's observed mean of all their values (not
the mean of the lagged values). For compositional outcomes, the helper
rebuilds the ILR coordinates through
[`complr()`](https://florale.github.io/multilevelcoda/reference/complr.md)
using the simulator's parts and SBP metadata, then lags the generated
`z` coordinates used by
[`brmcoda()`](https://florale.github.io/multilevelcoda/reference/brmcoda.md).

## Pragmatic default estimator

[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
simulates latent residual AR/VAR dynamics around the model-implied mean
and resolves `between()`/[`within()`](https://rdrr.io/r/base/with.html)
from latent generating components supplied by upstream predictor
generators. `prep_sim_analysis()` instead constructs the model applied
analysts commonly fit to observed data:

- `between(x)` and `within(x)` become `x_between` and `x_within`,
  recomputed from realised person means of the observed `x` (manifest
  centering), not from the latent generating components.

- `ar1()` becomes person-mean-centered lagged *observed* response
  predictors (`lag_<response>_within`), not the latent residual state.

These observed-data constructions target different estimands from the
simulation truth. Manifest person-mean centering and observed-score
lagged regression are known to yield biased estimates of between-person
effects and of inertia and cross-lag parameters relative to the latent
generating values, especially with short series (Ludtke et al. 2008;
Hamaker & Grasman 2014). This mismatch is intentional: it lets
simulation studies quantify the bias of the pragmatic estimator. Do not
interpret systematic discrepancies between estimates from this default
analysis model and the
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
truth parameters as errors in the simulator. For matched-model recovery
studies, construct the analysis model by hand instead of using this
helper.

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
