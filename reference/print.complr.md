# Print a Summary for a `complr` object

Print a Summary for a `complr` object

## Usage

``` r
# S3 method for class 'complr'
print(x, ...)
```

## Arguments

- x:

  An object of class `complr`.

- ...:

  Other arguments passed to
  [`summary.complr`](https://florale.github.io/multilevelcoda/reference/summary.complr.md).

## See also

[`summary.complr`](https://florale.github.io/multilevelcoda/reference/summary.complr.md)

## Examples

``` r
cilr <- complr(data = mcompd, sbp = sbp,
                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
                idvar = "ID")
print(cilr)
#> Summary of Composition and Logratio Variables
#> =============================================
#> 
#>  Number of observations: 3540
#>  Number of groups (IDs): 266
#> 
#>  Transformed variables:
#>  Composition geometry: acomp, logratio ransformation: ilr
#> 
#>                       Composition 1
#> X     tTST, tWAKE, tMVPA, tLPA, tSB
#> bX    bTST, bWAKE, bMVPA, bLPA, bSB
#> wX    wTST, wWAKE, wMVPA, wLPA, wSB
#> Z            z1_1, z2_1, z3_1, z4_1
#> bZ       bz1_1, bz2_1, bz3_1, bz4_1
#> wZ       wz1_1, wz2_1, wz3_1, wz4_1
#> total                             1
```
