# Build Sequential Binary Partition

Build a default sequential binary partition for `complr` object. The
default sequential binary partition is a pivot balance that allows the
effect of this first balance coordinate to be interpreted as the change
in the prediction for the dependent variable when that given part
increases while all remaining parts decrease by a common proportion.

## Usage

``` r
build.sbp(parts)
```

## Arguments

- parts:

  A character vector specifying the names of compositional variables to
  be used.

## Value

A matrix sequential binary partition.

## Examples

``` r
sbp1 <- build.sbp(parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
print(sbp1)
#>      TST WAKE MVPA LPA SB
#> [1,]   1   -1   -1  -1 -1
#> [2,]   0    1   -1  -1 -1
#> [3,]   0    0    1  -1 -1
#> [4,]   0    0    0   1 -1

sbp2 <- build.sbp(c("WAKE", "MVPA", "LPA", "SB"))
print(sbp2)
#>      WAKE MVPA LPA SB
#> [1,]    1   -1  -1 -1
#> [2,]    0    1  -1 -1
#> [3,]    0    0   1 -1
```
