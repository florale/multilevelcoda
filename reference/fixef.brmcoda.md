# Population-Level Estimates

Extract the population-level ('fixed') effects from the `brmsfit` object
in a `brmcoda` object.

## Usage

``` r
# S3 method for class 'brmcoda'
fixef(object, ...)
```

## Arguments

- object:

  An object of class `brmcoda`.

- ...:

  Further arguments passed to
  [`fixef.brmsfit`](https://paulbuerkner.com/brms/reference/fixef.brmsfit.html).

## Value

If `summary` is `TRUE`, a matrix returned by
[`posterior_summary`](https://paulbuerkner.com/brms/reference/posterior_summary.html)
for the population-level effects. If `summary` is `FALSE`, a matrix with
one row per posterior draw and one column per population-level effect.

## See also

[`fixef.brmsfit`](https://paulbuerkner.com/brms/reference/fixef.brmsfit.html)

## Examples

``` r
# \donttest{
## fit a model
if(requireNamespace("cmdstanr")){
  ## fit a model
  m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
                               parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
                               idvar = "ID", total = 1440),
               formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                                  wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
                                  chain = 1, iter = 500,
                                  backend = "cmdstanr")

  ## extract population-Level coefficients
  fixef(m)
}# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
