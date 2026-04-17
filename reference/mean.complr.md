# Mean amounts and mean compositions presented in a `complr` object.

Mean amounts and mean compositions presented in a `complr` object.

## Usage

``` r
# S3 method for class 'complr'
mean(x, weight = c("equal", "proportional"), parts = 1, ...)
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
mean(x)
#> $X
#>         tTST        tWAKE        tMVPA         tLPA          tSB 
#> "0.31090985" "0.04219846" "0.09859420" "0.38164618" "0.16665131" 
#> attr(,"class")
#> [1] "acomp"
#> 
#> $bX
#>         bTST        bWAKE        bMVPA         bLPA          bSB 
#> "0.31123086" "0.04213887" "0.09890581" "0.38162842" "0.16609604" 
#> attr(,"class")
#> [1] "acomp"
#> 
#> $wX
#>        wTST       wWAKE       wMVPA        wLPA         wSB 
#> "0.2024762" "0.1957899" "0.1958865" "0.2045439" "0.2013036" 
#> attr(,"class")
#> [1] "acomp"
#> 
#> $Z
#>       z1_1       z2_1       z3_1       z4_1 
#> -0.4658468  1.4872448 -0.8680194  0.6631910 
#> 
#> $bZ
#>      bz1_1      bz2_1      bz3_1      bz4_1 
#> -0.5212779  1.4167793 -0.7802634  0.5796594 
#> 
#> $wZ
#>       wz1_1       wz2_1       wz3_1       wz4_1 
#>  0.05543108  0.07046554 -0.08775604  0.08353160 
#> 
```
