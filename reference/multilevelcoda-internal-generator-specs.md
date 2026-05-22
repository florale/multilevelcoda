# Internal Generator Specification Helpers

Create and normalize generator specification objects used by
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).
These helpers are internal implementation details for the public
`gen_*()` API.

## Usage

``` r
.mlsim_spec(type, vars, level, simulate, ...)
```

## Arguments

- type:

  Character scalar naming the generator type.

- vars:

  Character vector of generated variable names.

- level:

  Simulation level: `"single"`, `"level2"`, or `"multilevel"`.

- simulate:

  Function called by
  [`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md)
  with a simulation context.

- ...:

  Additional metadata stored on the generator specification.

## Value

Internal helper outputs used by the simulation engine.

## Examples

``` r
spec <- multilevelcoda:::.mlsim_spec(
  "custom", "x", "single",
  function(context) {
    multilevelcoda:::.mlsim_result(rep(0, context$n_rows), "x", list())
  }
)
inherits(spec, "mlsim_generator_spec")
#> [1] TRUE
```
