# Build Base Pairwise Substitution

Make a data set of all possible pairwise substitution of a composition
which can be used as the base for substitution models.

## Usage

``` r
build.base(parts, type = NULL)
```

## Arguments

- parts:

  A character vector specifying the names of compositional variables to
  be used.

- type:

  Either `"one-to-one"` or `"one-to-all"`. Default is `"one-to-one"`.

## Value

A data table of all possible pairwise substitution.

## Examples

``` r
ps1 <- build.base(parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
print(ps1)
#>       TST  WAKE  MVPA   LPA    SB
#>     <num> <num> <num> <num> <num>
#>  1:     1    -1     0     0     0
#>  2:     1     0    -1     0     0
#>  3:     1     0     0    -1     0
#>  4:     1     0     0     0    -1
#>  5:    -1     1     0     0     0
#>  6:     0     1    -1     0     0
#>  7:     0     1     0    -1     0
#>  8:     0     1     0     0    -1
#>  9:    -1     0     1     0     0
#> 10:     0    -1     1     0     0
#> 11:     0     0     1    -1     0
#> 12:     0     0     1     0    -1
#> 13:    -1     0     0     1     0
#> 14:     0    -1     0     1     0
#> 15:     0     0    -1     1     0
#> 16:     0     0     0     1    -1
#> 17:    -1     0     0     0     1
#> 18:     0    -1     0     0     1
#> 19:     0     0    -1     0     1
#> 20:     0     0     0    -1     1
#>       TST  WAKE  MVPA   LPA    SB

ps2 <- build.base(c("WAKE", "MVPA", "LPA", "SB"), type = "one-to-all")
print(ps2)
#>          WAKE       MVPA        LPA         SB
#>         <num>      <num>      <num>      <num>
#> 1:  1.0000000 -0.3333333 -0.3333333 -0.3333333
#> 2: -0.3333333  1.0000000 -0.3333333 -0.3333333
#> 3: -0.3333333 -0.3333333  1.0000000 -0.3333333
#> 4: -0.3333333 -0.3333333 -0.3333333  1.0000000
#> 5: -1.0000000  0.3333333  0.3333333  0.3333333
#> 6:  0.3333333 -1.0000000  0.3333333  0.3333333
#> 7:  0.3333333  0.3333333 -1.0000000  0.3333333
#> 8:  0.3333333  0.3333333  0.3333333 -1.0000000
```
