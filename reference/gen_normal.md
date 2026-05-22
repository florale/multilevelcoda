# Generate Normal Variables

Create a generator specification for univariate Gaussian variables.

## Usage

``` r
gen_normal(
  vars,
  level = c("single", "level2", "multilevel"),
  mean = 0,
  sd = NULL,
  fixed_intercept = NULL,
  random_cov = NULL,
  residual_var = NULL,
  scale_fixed_intercept = NULL
)
```

## Arguments

- vars:

  Character scalar naming the generated variable.

- level:

  Simulation level. `"single"` generates one row-level value per
  observation, `"level2"` generates one group-level value and expands it
  to each row in the group, and `"multilevel"` generates row-level
  values from a random-intercept model.

- mean, sd:

  Direct mean and standard deviation for `"single"` and `"level2"`
  generators. Values are recycled to the number of generated units.

- fixed_intercept:

  Location intercept for multilevel generators. For normal variables
  this is on the response scale.

- random_cov:

  Group-level random-intercept covariance. A scalar is accepted for
  univariate location-only generators; a matrix is required when
  location and scale random effects are modeled jointly.

- residual_var:

  Residual variance for multilevel normal generators when no scale model
  is used.

- scale_fixed_intercept:

  Optional log residual standard-deviation intercept for multilevel
  generators. When supplied, residual standard deviations are
  `exp(scale_fixed_intercept + group_random_effect)`.

## Value

An `mlsim_generator_spec` object for use in
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

## Details

For non-multilevel normal variables, `mean` and `sd` define independent
draws directly. For multilevel variables, `fixed_intercept`,
`random_cov`, and either `residual_var` or `scale_fixed_intercept`
define a random-intercept data-generating model.

## See also

Other predictor generators:
[`continuous-generators`](https://florale.github.io/multilevelcoda/reference/continuous-generators.md),
[`count-generators`](https://florale.github.io/multilevelcoda/reference/count-generators.md),
[`gen_categorical()`](https://florale.github.io/multilevelcoda/reference/gen_categorical.md),
[`gen_custom()`](https://florale.github.io/multilevelcoda/reference/gen_custom.md),
[`gen_mvn()`](https://florale.github.io/multilevelcoda/reference/gen_mvn.md)

## Examples

``` r
sim <- simulate_data(
  n_groups = 2,
  n_per_group = 3,
  seed = 1,
  generators = list(
    x = gen_normal("x", mean = 0, sd = 1),
    u = gen_normal(
      "u",
      level = "multilevel",
      fixed_intercept = 0,
      random_cov = 0.2,
      residual_var = 1
    )
  )
)
sim$data
#>    group_id obs_id          x           u
#>       <int>  <int>      <num>       <num>
#> 1:        1      1 -0.6264538  0.79376625
#> 2:        1      2  0.1836433 -0.08740349
#> 3:        1      3 -0.8356286  1.72976607
#> 4:        2      1  1.5952808  0.72003208
#> 5:        2      2  0.3295078 -0.29105173
#> 6:        2      3 -0.8204684 -1.88451104
```
