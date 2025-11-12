# Create a Summary of a fitted `brmsfit` model in a `brmcoda` object

Create a Summary of a fitted `brmsfit` model in a `brmcoda` object

## Usage

``` r
# S3 method for class 'brmcoda'
summary(object, ...)
```

## Arguments

- object:

  An object of class `brmcoda`.

- ...:

  Other arguments passed to
  [`summary.brmsfit`](https://paulbuerkner.com/brms/reference/summary.brmsfit.html).

## Examples

``` r
# \donttest{
if(requireNamespace("cmdstanr")){
  m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
                                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
                                 idvar = "ID", total = 1440),
  formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                     wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
  chain = 1, iter = 500,
  backend = "cmdstanr")

  summary(m)
}# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
