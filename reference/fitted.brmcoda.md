# Expected Values of the Posterior Predictive Distribution

Compute posterior draws of the expected value of the posterior
predictive distribution of a `brmsfit` model in the `brmcoda` object.
Can be performed for the data used to fit the model (posterior
predictive checks) or for new data. By definition, these predictions
have smaller variance than the posterior predictions performed by the
[`predict.brmcoda`](https://florale.github.io/multilevelcoda/reference/predict.brmcoda.md)
method. This is because only the uncertainty in the expected value of
the posterior predictive distribution is incorporated in the draws
computed by `fitted` while the residual error is ignored there. However,
the estimated means of both methods averaged across draws should be very
similar.

## Usage

``` r
# S3 method for class 'brmcoda'
fitted(object, scale = c("linear", "response"), parts = 1, summary = TRUE, ...)
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

- summary:

  Should summary statistics be returned instead of the raw values?
  Default is `TRUE`.

- ...:

  Further arguments passed to
  [`fitted.brmsfit`](https://paulbuerkner.com/brms/reference/fitted.brmsfit.html)
  that control additional aspects of prediction.

## Value

An `array` of predicted *mean* response values. If `summary = FALSE` the
output resembles those of
[`posterior_epred.brmsfit`](https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.html).

If `summary = TRUE` the output depends on the family: For categorical
and ordinal families, the output is an N x E x C array, where N is the
number of observations, E is the number of summary statistics, and C is
the number of categories. For all other families, the output is an N x E
matrix. The number of summary statistics E is equal to
`2 + length(probs)`: The `Estimate` column contains point estimates
(either mean or median depending on argument `robust`), while the
`Est.Error` column contains uncertainty estimates (either standard
deviation or median absolute deviation depending on argument `robust`).
The remaining columns starting with `Q` contain quantile estimates as
specified via argument `probs`.

In multivariate models, an additional dimension is added to the output
which indexes along the different response variables.

## See also

[`fitted.brmsfit`](https://paulbuerkner.com/brms/reference/fitted.brmsfit.html)

## Examples

``` r
# \donttest{
## fit a model
if(requireNamespace("cmdstanr")){
  ## compute composition and ilr coordinates
  x <- complr(data = mcompd, sbp = sbp,
                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
                 idvar = "ID", total = 1440)

  ## fit a model
  m1 <- brmcoda(complr = x,
                formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
                chain = 1, iter = 500,
                backend = "cmdstanr")

  ## compute expected predictions
  epred <- fitted(m1)
  head(epred)

  ## fit a model with compositional outcome
  m2 <- brmcoda(complr = x,
                formula = mvbind(z1_1, z2_1, z3_1, z4_1) ~ Stress + Female + (1 | ID),
                chain = 1, iter = 500,
                backend = "cmdstanr")

  ## expected predictions on compositional scale
  epredcomp <- fitted(m2, scale = "response")
  head(epredcomp)
}# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
