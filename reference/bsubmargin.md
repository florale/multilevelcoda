# Between-person Average Substitution

This function is an alias of
[`substitution`](https://florale.github.io/multilevelcoda/reference/substitution.md)
to estimates the difference in an outcome when compositional parts are
substituted for specific unit(s) at *between* level using cluster mean
(e.g., compositional mean at individual level) as reference composition.
It is recommended that users run substitution model using the
[`substitution`](https://florale.github.io/multilevelcoda/reference/substitution.md)
function.

## Usage

``` r
bsubmargin(
  object,
  delta,
  ref = "clustermean",
  level = "between",
  summary = TRUE,
  at = NULL,
  parts = 1,
  base,
  type = "one-to-one",
  weight = "proportional",
  scale = c("response", "linear"),
  cores = NULL,
  ...
)
```

## Arguments

- object:

  A fitted
  [`brmcoda`](https://florale.github.io/multilevelcoda/reference/brmcoda.md)
  object.

- delta:

  A integer, numeric value or vector indicating the amount of
  substituted change between compositional parts.

- ref:

  Either a character value or vector or a dataset. Can be `"grandmean"`
  and/or `"clustermean"`, or a `data.frame` or `data.table` of user's
  specified reference grid consisting of combinations of covariates over
  which predictions are made. User's specified reference grid is only
  possible for simple substitution. Single level models are default to
  `"grandmean"`.

- level:

  A character string or vector. Should the estimate of multilevel models
  focus on the `"between"` and/or `"within"` or `"aggregate"` variance?
  Single-level models are default to `"aggregate"`.

- summary:

  A logical value to obtain summary statistics instead of the raw
  values. Default is `TRUE`. Currently only support outputing raw values
  for model using grandmean as reference composition.

- at:

  An optional named list of levels for the corresponding variables in
  the reference grid.

- parts:

  A optional character string specifying names of compositional parts
  that should be considered in the substitution analysis. This should
  correspond to a single set of names of compositional parts specified
  in the `complr` object. Default to the first composition in the
  `complr` object.

- base:

  An optional base substitution. Can be a `data.frame` or `data.table`
  of the base possible substitution of compositional parts, which can be
  computed using function
  [`build.base`](https://florale.github.io/multilevelcoda/reference/build.base.md).

- type:

  A character string to indicate the type of substitution. If
  `"one-to-all"`, all possible one-to-remaining reallocations are
  estimated. If `"one-to-one"`, all possible one-to-one reallocations
  are estimated. If `"equal"`, give equal weight to units (e.g.,
  individuals). If `"proportional"`, weights in proportion to the
  frequencies of units being averaged (e.g., observations across
  individuals). Default to `"equal"` for `ref = "grandmean"` and
  `"proportional"` for `ref = "clustermean"`.

- weight:

  A character value specifying the weight to use in calculation of the
  reference composition.

- scale:

  Either `"response"` or `"linear"`. If `"response"`, results are
  returned on the scale of the response variable. If `"linear"`, results
  are returned on the scale of the linear predictor term, that is
  without applying the inverse link function or other transformations.

- cores:

  Number of cores to use when executing the chains in parallel, we
  recommend setting the `mc.cores` option to be as many processors as
  the hardware and RAM allow (up to the number of compositional parts).
  For non-Windows OS in non-interactive R sessions, forking is used
  instead of PSOCK clusters. Default to `"one-to-one"`.

- ...:

  Further arguments passed to
  [`posterior_summary`](https://paulbuerkner.com/brms/reference/posterior_summary.html).

## Value

A list containing the results of multilevel compositional substitution
model. The first six lists contain the results of the substitution
estimation for a compositional part.

- `Mean`:

  Posterior means.

- `CI_low` and `CI_high`:

  95% credible intervals.

- `Delta`:

  Amount substituted across compositional parts.

- `From`:

  Compositional part that is substituted from.

- `To`:

  Compositional parts that is substituted to.

- `Level`:

  Level where changes in composition takes place.

- `Reference`:

  Either `grandmean`, `clustermean`, or `users`.

## See also

[`substitution`](https://florale.github.io/multilevelcoda/reference/substitution.md)

## Examples

``` r
# \donttest{
if(requireNamespace("cmdstanr")){
cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp, 
               parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)

m <- brmcoda(complr = cilr, 
             formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 + 
                                wz1_1 + wz2_1 + wz3_1 + wz4_1 + 
                                Female + (1 | ID), 
             chains = 1, iter = 500,
             backend = "cmdstanr")
             
subm <- bsubmargin(object = m, base = psub, delta = 5)
}# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
