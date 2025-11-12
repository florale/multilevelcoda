# Create a Summary of a `complr` object

Create a Summary of a `complr` object

## Usage

``` r
# S3 method for class 'complr'
summary(object, ...)
```

## Arguments

- object:

  An object of class `complr`.

- ...:

  generic argument, not in use.

## Examples

``` r
x1 <- complr(data = mcompd, sbp = sbp,
               parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
               idvar = "ID")
summary(x1)
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
x2 <- complr(data = mcompd, sbp = sbp,
               parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
summary(x2)
#> Summary of Composition and Logratio Variables
#> =============================================
#> 
#>  Number of observations: 3540
#>  Transformed variables:
#>  Composition geometry: acomp, logratio ransformation: ilr
#> 
#>                       Composition 1
#> X     tTST, tWAKE, tMVPA, tLPA, tSB
#> Z            z1_1, z2_1, z3_1, z4_1
#> total                             1
```
