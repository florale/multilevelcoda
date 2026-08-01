# Generate Count Variables

Create generator specifications for binomial, Poisson, and negative
binomial count variables.

## Usage

``` r
gen_binomial(
  vars,
  level = c("single", "level2", "multilevel"),
  size = NULL,
  fixed_intercept = NULL,
  random_cov = NULL
)

gen_poisson(
  vars,
  level = c("single", "level2", "multilevel"),
  fixed_intercept = NULL,
  random_cov = NULL
)

gen_negbin(
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

- size:

  Binomial trial size. May be a scalar, vector, function of `n`, or
  count-distribution list. Count-distribution `minimum`/`maximum` bounds
  clamp draws to the bounds (censoring, not truncation).

- fixed_intercept:

  Link-scale intercept. Binomial uses the logit link; Poisson and
  negative binomial use the log link.

- random_cov:

  Group-level random-intercept covariance for multilevel generators.

- ...:

  Removed direct distribution parameters are rejected.

- scale_fixed_intercept:

  Log negative-binomial size intercept.

## Value

An `mlsim_generator_spec` object for use in
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

## Details

Count generators use link-scale intercepts at every level.
Negative-binomial size is `exp(scale_fixed_intercept)`.

## See also

Other predictor generators:
[`continuous-generators`](https://florale.github.io/multilevelcoda/reference/continuous-generators.md),
[`gen_categorical()`](https://florale.github.io/multilevelcoda/reference/gen_categorical.md),
[`gen_custom()`](https://florale.github.io/multilevelcoda/reference/gen_custom.md),
[`gen_mvn()`](https://florale.github.io/multilevelcoda/reference/gen_mvn.md),
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md),
[`gen_template()`](https://florale.github.io/multilevelcoda/reference/gen_template.md)

## Examples

``` r
sim <- simulate_data(
  n = 5,
  seed = 4,
  generators = list(
    events = gen_poisson("events", fixed_intercept = log(2)),
    successes = gen_binomial(
      "successes",
      size = 4,
      fixed_intercept = stats::qlogis(0.5)
    ),
    overdispersed = gen_negbin(
      "overdispersed",
      fixed_intercept = log(2),
      scale_fixed_intercept = log(3)
    )
  )
)
sim$data
#>    obs_id events successes overdispersed
#>     <int>  <int>     <int>         <num>
#> 1:      1      2         1             1
#> 2:      2      0         3             4
#> 3:      3      1         3             8
#> 4:      4      1         4             8
#> 5:      5      3         1             2
```
