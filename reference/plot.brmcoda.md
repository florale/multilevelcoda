# Trace and Density Plots for MCMC Draws plot

Make a plot of `brmcoda` model results.

## Usage

``` r
# S3 method for class 'brmcoda'
plot(x, ...)
```

## Arguments

- x:

  A
  [`brmcoda`](https://florale.github.io/multilevelcoda/reference/brmcoda.md)
  class object.

- ...:

  Further arguments passed to
  [`plot.brmsfit`](https://paulbuerkner.com/brms/reference/plot.brmsfit.html).

## Value

An invisible list of
[`gtable`](https://gtable.r-lib.org/reference/gtable.html) objects.

## See also

[`plot.brmsfit`](https://paulbuerkner.com/brms/reference/plot.brmsfit.html)

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
plot(fit)
} # }
```
