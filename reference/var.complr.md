# Variance of compositions presented in a `complr` object.

Variance of compositions presented in a `complr` object.

## Usage

``` r
# S3 method for class 'complr'
var(x, weight = c("equal", "proportional"), parts = 1, ...)
```

## Arguments

- x:

  An object of class `complr`.

- weight:

  A character value specifying the weight to use in calculation of the
  reference composition. If `"equal"`, give equal weight to units (e.g.,
  individuals). If `"proportional"`, weights in proportion to the
  frequencies of units being averaged (e.g., observations across
  individuals) Default is `equal`.

- parts:

  A optional character string specifying names of compositional parts
  that should be considered in the substitution analysis. This should
  correspond to a single set of names of compositional parts specified
  in the `complr` object. Default to the first composition in the
  `complr` object.

- ...:

  generic argument, not in use.

## Examples

``` r
x <- complr(data = mcompd, sbp = sbp,
                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
                idvar = "ID")
multilevelcoda:::var.complr(x)
#> $X
#>             tTST     tWAKE     tMVPA       tLPA       tSB
#> tTST  0.00000000 0.2407016 0.2147058 0.06574185 0.1677480
#> tWAKE 0.24070156 0.0000000 0.4651644 0.27222191 0.4910045
#> tMVPA 0.21470581 0.4651644 0.0000000 0.19851113 0.4266598
#> tLPA  0.06574185 0.2722219 0.1985111 0.00000000 0.2176891
#> tSB   0.16774797 0.4910045 0.4266598 0.21768908 0.0000000
#> 
#> $bX
#>             bTST     bWAKE     bMVPA       bLPA       bSB
#> bTST  0.00000000 0.2407008 0.2146939 0.06573977 0.1677550
#> bWAKE 0.24070082 0.0000000 0.4651665 0.27222554 0.4910121
#> bMVPA 0.21469390 0.4651665 0.0000000 0.19850554 0.4266513
#> bLPA  0.06573977 0.2722255 0.1985055 0.00000000 0.2176890
#> bSB   0.16775502 0.4910121 0.4266513 0.21768902 0.0000000
#> 
#> $wX
#>               wTST        wWAKE        wMVPA         wLPA          wSB
#> wTST  0.0000000000 0.0003468755 0.0007531461 0.0001656332 0.0006786266
#> wWAKE 0.0003468755 0.0000000000 0.0008223263 0.0005175806 0.0009265183
#> wMVPA 0.0007531461 0.0008223263 0.0000000000 0.0007049798 0.0018683307
#> wLPA  0.0001656332 0.0005175806 0.0007049798 0.0000000000 0.0009627149
#> wSB   0.0006786266 0.0009265183 0.0018683307 0.0009627149 0.0000000000
#> 
```
