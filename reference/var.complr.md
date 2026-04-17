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
#> tTST  0.00000000 0.2414276 0.2151381 0.06580904 0.1723962
#> tWAKE 0.24142765 0.0000000 0.4691750 0.27168933 0.4901382
#> tMVPA 0.21513810 0.4691750 0.0000000 0.20059960 0.4367386
#> tLPA  0.06580904 0.2716893 0.2005996 0.00000000 0.2198577
#> tSB   0.17239620 0.4901382 0.4367386 0.21985775 0.0000000
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
#> wTST  0.0000000000 0.0003499222 0.0007538029 0.0001665439 0.0006894323
#> wWAKE 0.0003499222 0.0000000000 0.0008186745 0.0005160081 0.0009567599
#> wMVPA 0.0007538029 0.0008186745 0.0000000000 0.0007023513 0.0018927582
#> wLPA  0.0001665439 0.0005160081 0.0007023513 0.0000000000 0.0009834585
#> wSB   0.0006894323 0.0009567599 0.0018927582 0.0009834585 0.0000000000
#> 
```
