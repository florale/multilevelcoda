# MCMC Plots Implemented in bayesplot

Call MCMC plotting functions implemented in the bayesplot package.

## Usage

``` r
# S3 method for class 'brmcoda'
mcmc_plot(object, ...)
```

## Arguments

- object:

  A `brmcoda` class object.

- ...:

  Further arguments passed to
  [`mcmc_plot.brmsfit`](https://paulbuerkner.com/brms/reference/mcmc_plot.brmsfit.html).

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html) object
that can be further customized using the ggplot2 package.

## See also

[`mcmc_plot.brmsfit`](https://paulbuerkner.com/brms/reference/mcmc_plot.brmsfit.html)

## Examples

``` r
if (FALSE) { # \dontrun{
cilr <- complr(data = mcompd, sbp = sbp,
        parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")

# model with compositional predictor at between and within-person levels
fit <- brmcoda(complr = cilr,
                formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
               chain = 1, iter = 500)
mcmc_plot(fit)
} # }
```
