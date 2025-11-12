# Draws from the Posterior Predictive Distribution

Compute posterior draws of the posterior predictive distribution of a
`brmsfit` model in the `brmcoda` object. Can be performed for the data
used to fit the model (posterior predictive checks) or for new data. By
definition, these draws have higher variance than draws of the expected
value of the posterior predictive distribution computed by
[`fitted.brmcoda`](https://florale.github.io/multilevelcoda/reference/fitted.brmcoda.md).
This is because the residual error is incorporated in
`posterior_predict`. However, the estimated means of both methods
averaged across draws should be very similar.

## Usage

``` r
# S3 method for class 'brmcoda'
predict(object, scale = c("linear", "response"), parts = 1, ...)
```

## Arguments

- object:

  An object of class `brmcoda`.

- scale:

  Specifically for models with compositional responses, either
  `"response"` or `"linear"`. If `"linear"`, results are returned on the
  log-ratio scale. If `"response"`, results are returned on the
  compositional scale of the response variable.

- parts:

  Only for models with compositional response A optional character
  string specifying names of compositional parts that make up the
  response in `brmcoda` model. This should correspond to a single set of
  names of compositional parts specified in the `complr` object. Default
  to the first composition in the `complr` object.

- ...:

  Further arguments passed to
  [`predict.brmsfit`](https://paulbuerkner.com/brms/reference/predict.brmsfit.html)
  that control additional aspects of prediction.

## Value

An `array` of predicted response values. If `summary = FALSE` the output
resembles those of
[`posterior_predict.brmsfit`](https://paulbuerkner.com/brms/reference/posterior_predict.brmsfit.html).

If `summary = TRUE` the output depends on the family: For categorical
and ordinal families, the output is an N x C matrix, where N is the
number of observations, C is the number of categories, and the values
are predicted category probabilities. For all other families, the output
is a N x E matrix where E = `2 + length(probs)` is the number of summary
statistics: The `Estimate` column contains point estimates (either mean
or median depending on argument `robust`), while the `Est.Error` column
contains uncertainty estimates (either standard deviation or median
absolute deviation depending on argument `robust`). The remaining
columns starting with `Q` contain quantile estimates as specified via
argument `probs`.

## See also

[`predict.brmsfit`](https://paulbuerkner.com/brms/reference/predict.brmsfit.html)

## Examples

``` r
# \donttest{
if(requireNamespace("cmdstanr")){
  ## fit a model
  x <- complr(data = mcompd, sbp = sbp,
                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
                 idvar = "ID", total = 1440)

  m1 <- brmcoda(complr = x,
                formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
                chain = 1, iter = 500,
                backend = "cmdstanr")

  ## predicted responses
  pred <- predict(m1)
  head(pred)

  ## fit a model with compositional outcome
  m2 <- brmcoda(complr = x,
                formula = mvbind(z1_1, z2_1, z3_1, z4_1) ~ 
                          bz1_1 + bz2_1 + bz3_1 + bz4_1 + Female + (1 | ID),
                chain = 1, iter = 500,
                backend = "cmdstanr")

  ## predicted responses on ilr scale
  predilr <- predict(m2, scale = "linear")
  head(predilr)

  ## predicted responses on compositional scale
  predcomp <- predict(m2, scale = "response")
  head(predcomp)
}# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
