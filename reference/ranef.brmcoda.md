# Group-Level Estimates

Extract the group-level ('random') effects of each level of the
`brmsfit` object in a `brmcoda` object.

## Usage

``` r
# S3 method for class 'brmcoda'
ranef(object, ...)
```

## Arguments

- object:

  An object of class `brmcoda`.

- ...:

  Further arguments passed to
  [`ranef.brmsfit`](https://paulbuerkner.com/brms/reference/ranef.brmsfit.html).

## Value

A list of 3D arrays (one per grouping factor). If `summary` is `TRUE`,
the 1st dimension contains the factor levels, the 2nd dimension contains
the summary statistics (see
[`posterior_summary`](https://paulbuerkner.com/brms/reference/posterior_summary.html)),
and the 3rd dimension contains the group-level effects. If `summary` is
`FALSE`, the 1st dimension contains the posterior draws, the 2nd
dimension contains the factor levels, and the 3rd dimension contains the
group-level effects.

## See also

[`ranef.brmsfit`](https://paulbuerkner.com/brms/reference/ranef.brmsfit.html)

## Examples

``` r
# \donttest{
## fit a model
if(requireNamespace("cmdstanr")){
  m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
                               parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
                               idvar = "ID", total = 1440),
               formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                                  wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
                                  chain = 1, iter = 500,
                                  backend = "cmdstanr")

  ## extract group-level coefficients
  ranef(m)
}# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
