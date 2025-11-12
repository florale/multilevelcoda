# Extract Variance and Correlation Components

Calculates the estimated standard deviations, correlations and
covariances of the group-level terms of the `brmsfit` object in a
`brmcoda` object.

## Usage

``` r
# S3 method for class 'brmcoda'
VarCorr(x, ...)
```

## Arguments

- x:

  An object of class `brmcoda`.

- ...:

  Further arguments passed to
  [`VarCorr.brmsfit`](https://paulbuerkner.com/brms/reference/VarCorr.brmsfit.html).

## Value

A list of lists (one per grouping factor), each with three elements: a
matrix containing the standard deviations, an array containing the
correlation matrix, and an array containing the covariance matrix with
variances on the diagonal.

## See also

[`VarCorr.brmsfit`](https://paulbuerkner.com/brms/reference/VarCorr.brmsfit.html)

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

  VarCorr(m)
}# }
#> Loading required namespace: cmdstanr
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
