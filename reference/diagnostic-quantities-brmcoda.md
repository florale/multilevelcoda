# Extract Diagnostic Quantities from `brmsfit` Models in `brmcoda`

Extract Diagnostic Quantities from `brmsfit` Models in `brmcoda`

## Usage

``` r
# S3 method for class 'brmcoda'
log_posterior(object, ...)

# S3 method for class 'brmcoda'
nuts_params(object, ...)

# S3 method for class 'brmcoda'
rhat(x, ...)

# S3 method for class 'brmcoda'
neff_ratio(object, ...)
```

## Arguments

- ...:

  Arguments passed to individual methods (if applicable).

- x, object:

  A `brmcoda` object or another R object for which the methods are
  defined.

## Value

The exact form of the output depends on the method.

## See also

[`log_posterior.brmsfit`](https://paulbuerkner.com/brms/reference/diagnostic-quantities.html)

[`nuts_params.brmsfit`](https://paulbuerkner.com/brms/reference/diagnostic-quantities.html)

[`rhat.brmsfit`](https://paulbuerkner.com/brms/reference/diagnostic-quantities.html)

[`neff_ratio.brmsfit`](https://paulbuerkner.com/brms/reference/diagnostic-quantities.html)
