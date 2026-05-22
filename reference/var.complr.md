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
## ensure dispatch to the s3 generic var function
## defined in the compositions package
compositions::var(x)
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
#> bTST  0.00000000 0.2414272 0.2151267 0.06580703 0.1724031
#> bWAKE 0.24142719 0.0000000 0.4691787 0.27169303 0.4901455
#> bMVPA 0.21512670 0.4691787 0.0000000 0.20059492 0.4367313
#> bLPA  0.06580703 0.2716930 0.2005949 0.00000000 0.2198574
#> bSB   0.17240306 0.4901455 0.4367313 0.21985738 0.0000000
#> 
#> $wX
#>               wTST        wWAKE        wMVPA         wLPA          wSB
#> wTST  0.0000000000 0.0003499222 0.0007538029 0.0001665439 0.0006894323
#> wWAKE 0.0003499222 0.0000000000 0.0008186745 0.0005160081 0.0009567599
#> wMVPA 0.0007538029 0.0008186745 0.0000000000 0.0007023513 0.0018927582
#> wLPA  0.0001665439 0.0005160081 0.0007023513 0.0000000000 0.0009834585
#> wSB   0.0006894323 0.0009567599 0.0018927582 0.0009834585 0.0000000000
#> 
```
