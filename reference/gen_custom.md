# Generate Variables with a User-Supplied Function

Wrap a custom generator function so it can participate in
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

## Usage

``` r
gen_custom(
  vars,
  generator,
  level = c("single", "level2", "multilevel"),
  ...,
  metadata = list()
)
```

## Arguments

- vars:

  Character vector naming the generated variables.

- generator:

  Function called as
  `generator(context = context, vars = vars, level = level, ...)`. It
  must return a list with `data` and `names` elements, and may include
  `metadata`.

- level:

  Simulation level. A level-2 custom generator may return either one row
  per group or one row per simulated observation.

- ...:

  Additional arguments passed to `generator`.

- metadata:

  List of wrapper metadata stored on the generator specification.

## Value

An `mlsim_generator_spec` object for use in
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

## See also

Other predictor generators:
[`continuous-generators`](https://florale.github.io/multilevelcoda/reference/continuous-generators.md),
[`count-generators`](https://florale.github.io/multilevelcoda/reference/count-generators.md),
[`gen_categorical()`](https://florale.github.io/multilevelcoda/reference/gen_categorical.md),
[`gen_mvn()`](https://florale.github.io/multilevelcoda/reference/gen_mvn.md),
[`gen_normal()`](https://florale.github.io/multilevelcoda/reference/gen_normal.md)

## Examples

``` r
constant_generator <- function(context, vars, level, value = 1) {
  list(data = data.frame(value = rep(value, context$n_rows)), names = vars)
}

sim <- simulate_data(
  n = 3,
  generators = list(x = gen_custom("x", constant_generator, value = 7))
)
sim$data
#>    obs_id     x
#>     <int> <num>
#> 1:      1     7
#> 2:      2     7
#> 3:      3     7
```
