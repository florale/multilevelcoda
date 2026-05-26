# Generate Positive and Bounded Continuous Variables

Create generator specifications for gamma and beta variables.

## Usage

``` r
gen_gamma(
  vars,
  level = c("single", "level2", "multilevel"),
  fixed_intercept = NULL,
  ...,
  scale_fixed_intercept = NULL,
  random_cov = NULL
)

gen_beta(
  vars,
  level = c("single", "level2", "multilevel"),
  fixed_intercept = NULL,
  ...,
  scale_fixed_intercept = NULL,
  random_cov = NULL
)
```

## Arguments

- vars:

  Character scalar naming the generated variable.

- level:

  Simulation level. `"single"` generates row-level values, `"level2"`
  generates group-level values, and `"multilevel"` uses a
  random-intercept model.

- fixed_intercept:

  Link-scale intercept. Gamma uses the log link and beta uses the logit
  link.

- ...:

  Removed direct distribution parameters are rejected.

- scale_fixed_intercept:

  Log shape intercept for gamma variables or log precision intercept for
  beta variables.

- random_cov:

  Group-level random-intercept covariance for multilevel generators.

## Value

An `mlsim_generator_spec` object for use in
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

## See also

Other predictor generators:
[`count-generators`](https://florale.github.io/multilevelcoda/reference/count-generators.md),
[`gen_categorical()`](https://florale.github.io/multilevelcoda/reference/gen_categorical.md),
[`gen_custom()`](https://florale.github.io/multilevelcoda/reference/gen_custom.md),
[`gen_mvn()`](https://florale.github.io/multilevelcoda/reference/gen_mvn.md),
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md),
[`gen_template()`](https://florale.github.io/multilevelcoda/reference/gen_template.md)

## Examples

``` r
sim <- simulate_data(
  n = 5,
  seed = 5,
  generators = list(
    positive = gen_gamma(
      "positive",
      fixed_intercept = log(1.5),
      scale_fixed_intercept = log(2)
    ),
    proportion = gen_beta(
      "proportion",
      fixed_intercept = stats::qlogis(0.4),
      scale_fixed_intercept = log(10)
    )
  )
)
sim$data
#>    obs_id  positive proportion
#>     <int>     <num>      <num>
#> 1:      1 0.4851945  0.4264494
#> 2:      2 0.6624663  0.2618874
#> 3:      3 1.1903527  0.6327267
#> 4:      4 3.2462519  0.5897335
#> 5:      5 0.6393495  0.5073695
```
