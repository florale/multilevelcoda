# Efficient approximate leave-one-out cross-validation (LOO)

Perform approximate leave-one-out cross-validation based on the
posterior likelihood using the loo package. For more details see
[`loo`](https://mc-stan.org/loo/reference/loo.html).

## Usage

``` r
# S3 method for class 'brmcoda'
loo(x, ...)
```

## Arguments

- x:

  A `brmcoda` object.

- ...:

  More `brmsfit` objects or further arguments passed to the underlying
  post-processing functions. In particular, see
  [`prepare_predictions`](https://paulbuerkner.com/brms/reference/prepare_predictions.html)
  for further supported arguments.

## Value

If just one object is provided, an object of class `loo`. If multiple
objects are provided, an object of class `loolist`.

## See also

[`loo.brmsfit`](https://paulbuerkner.com/brms/reference/loo.brmsfit.html)
