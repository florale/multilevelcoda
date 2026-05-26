# Print Simulated Data

Print a compact overview of an `mlsim_data` object.

## Usage

``` r
# S3 method for class 'mlsim_data'
print(x, ...)
```

## Arguments

- x:

  An `mlsim_data` object returned by
  [`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

- ...:

  Ignored.

## Value

`x`, invisibly.

## Examples

``` r
sim <- simulate_data(
  n = 3,
  generators = list(
    x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1)
  )
)
print(sim)
#> <mlsim_data>
#>   rows: 3
#>   columns: 2 (generated: 1)
#>   seed: <none>
#>   grouping: none
#>   generators: 1
#> 
#> Generators:
#>  generator distribution  level vars
#>          x          mvn single    x
```
