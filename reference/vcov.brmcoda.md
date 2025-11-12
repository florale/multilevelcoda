# Covariance and Correlation Matrix of Population-Level Effects

Get a point estimate of the covariance or correlation matrix of
population-level parameters of the `brmsfit` object in a `brmcoda`
object.

## Usage

``` r
# S3 method for class 'brmcoda'
vcov(object, ...)
```

## Arguments

- object:

  An object of class `brmcoda`.

- ...:

  Further arguments passed to
  [`vcov.brmsfit`](https://paulbuerkner.com/brms/reference/vcov.brmsfit.html).

## Value

covariance or correlation matrix of population-level parameters

## See also

[`vcov.brmsfit`](https://paulbuerkner.com/brms/reference/vcov.brmsfit.html)

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

  vcov(m)
}# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
