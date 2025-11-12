# Create a Summary of a fitted `brmsfit` model from a `pivot_coord` object

Create a Summary of a fitted `brmsfit` model from a `pivot_coord` object

## Usage

``` r
# S3 method for class 'pivot_coord'
summary(object, digits = 2, ...)
```

## Arguments

- object:

  An object of class `pivot_coord`.

- digits:

  A integer value used for number formatting. Default is `2`.

- ...:

  currently ignored.

## Value

A data table of results.

## Examples

``` r
# \donttest{
if(requireNamespace("cmdstanr")){
  m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
                               parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
                               idvar = "ID", total = 1440),
  formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                     wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
  chain = 1, iter = 500,
  backend = "cmdstanr")

  m_pc <- pivot_coord(m, method = "refit")
  summary(m_pc)
}# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
