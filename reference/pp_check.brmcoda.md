# Posterior Predictive Checks for `brmcoda` Objects

Perform posterior predictive checks with the help of the bayesplot
package.

## Usage

``` r
# S3 method for class 'brmcoda'
pp_check(object, ...)
```

## Arguments

- object:

  An object of class `brmcoda`.

- ...:

  Further arguments passed to
  [`predict.brmsfit`](https://paulbuerkner.com/brms/reference/predict.brmsfit.html)
  as well as to the PPC function specified in `type`.

## See also

[`pp_check.brmsfit`](https://paulbuerkner.com/brms/reference/pp_check.brmsfit.html)
