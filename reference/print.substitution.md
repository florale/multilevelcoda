# Print a Summary for a `substitution` object

Print a Summary for a `substitution` object

## Usage

``` r
# S3 method for class 'substitution'
print(x, ...)
```

## Arguments

- x:

  A `substitution` object.

- ...:

  Additional arguments to be passed to to method `summary` of
  `substitution`.

## See also

[`summary.substitution`](https://florale.github.io/multilevelcoda/reference/summary.substitution.md)

## Examples

``` r
# \donttest{
if(requireNamespace("cmdstanr")){
  ## fit a model with compositional predictor at between and between-person levels
  m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
                               parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
                               idvar = "ID", total = 1440),
  formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                     wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
  chain = 1, iter = 500,
  backend = "cmdstanr")

  subm <- substitution(object = m, delta = 5)
  print(subm)
}# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
