# Create a matrix of output plots from a [`brmcoda`](https://florale.github.io/multilevelcoda/reference/brmcoda.md)'s [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html) object

A [`pairs`](https://rdrr.io/r/graphics/pairs.html) method that is
customized for MCMC output.

## Usage

``` r
# S3 method for class 'brmcoda'
pairs(x, ...)
```

## Arguments

- x:

  A `brmcoda` class object.

- ...:

  Further arguments passed to
  [`pairs.brmsfit`](https://paulbuerkner.com/brms/reference/pairs.brmsfit.html).

## See also

[`pairs.brmsfit`](https://paulbuerkner.com/brms/reference/pairs.brmsfit.html)

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
pairs(fit)
} # }
```
