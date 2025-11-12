# Posterior Draws of Residuals/Predictive Errors

Compute posterior draws of residuals/predictive errors

## Usage

``` r
# S3 method for class 'brmcoda'
residuals(object, ...)
```

## Arguments

- object:

  An object of class `brmcoda`.

- ...:

  Further arguments passed to
  [`residuals.brmsfit`](https://paulbuerkner.com/brms/reference/residuals.brmsfit.html).

## Value

An `array` of predictive error/residual draws. If `summary = FALSE` the
output resembles those of
[`predictive_error.brmsfit`](https://paulbuerkner.com/brms/reference/predictive_error.brmsfit.html).
If `summary = TRUE` the output is an N x E matrix, where N is the number
of observations and E denotes the summary statistics computed from the
draws.

## See also

[`residuals.brmsfit`](https://paulbuerkner.com/brms/reference/residuals.brmsfit.html)

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

  ## extract residuals
  res <- residuals(m)
  head(res)
}# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
