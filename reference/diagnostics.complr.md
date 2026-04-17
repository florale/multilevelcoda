# Robust multivariate normal diagnostics for `complr` coordinates

Uses robust minimum covariance determinant estimates to calculate
squared Mahalanobis distances for the total, between, or within
log-ratio coordinates stored in a
[`complr`](https://florale.github.io/multilevelcoda/reference/complr.md)
object. This can be used to identify extreme values (outliers).

## Usage

``` r
# S3 method for class 'complr'
diagnostics(
  x,
  level = c("total", "between", "within"),
  parts = 1,
  ev.perc = 0.001,
  extremevalues = c("theoretical", "empirical"),
  ...
)
```

## Arguments

- x:

  A
  [`complr`](https://florale.github.io/multilevelcoda/reference/complr.md)
  object.

- level:

  A character string indicating which coordinates to diagnose:
  `"total"`, `"between"`, or `"within"`.

- parts:

  Optional character vector or integer specifying which set of
  compositional parts should be used. Defaults to the first composition
  in the `complr` object.

- ev.perc:

  Proportion in the upper tail used to flag extreme values. Defaults to
  `.001`.

- extremevalues:

  Character string indicating whether extreme values should be flagged
  using a `"theoretical"` chi-squared cutoff or an `"empirical"` cutoff
  from the observed robust distances.

- ...:

  Additional arguments passed to
  [`covMcd`](https://rdrr.io/pkg/robustbase/man/covMcd.html).

## Value

A `diagnostics` object with the raw composition matrix in `x`, a
`data.table` of robust distances in `distance`, the extreme-value cutoff
in `cutoff`, the requested upper-tail proportion in `ev.perc`, the
cutoff type in `extremevalues`, and the diagnosed coordinate level in
`levels`. The `distance` element includes an `is_extremevalue` column
indicating whether each fitted observation is flagged as an extreme
value. For `level = "between"`, diagnostics are fitted using one
observation per `idvar`.

## Examples

``` r
data(mcompd)
data(sbp)

ids <- unique(mcompd$ID)[1:20]
cilr <- complr(
  data = mcompd[ID %in% ids, .SD[1:3], by = ID],
  sbp = sbp,
  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
  idvar = "ID",
  total = 1440
)

# One diagnostic object per coordinate level
total_diag <- diagnostics(cilr, level = "total")
between_diag <- diagnostics(cilr, level = "between")
within_diag <- diagnostics(cilr, level = "within")

# Use an empirical cutoff and a larger tail proportion
empirical_diag <- diagnostics(
  cilr,
  level = "between",
  ev.perc = .05,
  extremevalues = "empirical"
)

is.diagnostics(between_diag)
#> [1] TRUE
head(between_diag$distance)
#>     distance expected_distance   deviates is_extremevalue
#>        <num>             <num>      <num>          <lgcl>
#> 1: 14.364835          5.672230  8.6926050           FALSE
#> 2:  1.417544          1.218762  0.1987823           FALSE
#> 3:  3.403716          3.518969 -0.1152530           FALSE
#> 4: 26.231626          6.342329 19.8892973            TRUE
#> 5:  2.545905          2.898220 -0.3523146           FALSE
#> 6: 57.552119         11.143287 46.4088323            TRUE
```
