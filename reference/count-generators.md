# Generate Count Variables

Create generator specifications for binomial, Poisson, and negative
binomial count variables.

## Usage

``` r
gen_binomial(
  vars,
  level = c("single", "level2", "multilevel"),
  size = 1,
  prob = NULL,
  fixed_intercept = NULL,
  random_cov = NULL
)

gen_poisson(
  vars,
  level = c("single", "level2", "multilevel"),
  lambda = NULL,
  fixed_intercept = NULL,
  random_cov = NULL
)

gen_negbin(
  vars,
  level = c("single", "level2", "multilevel"),
  size = NULL,
  mu = NULL,
  fixed_intercept = NULL,
  random_cov = NULL,
  scale_fixed_intercept = NULL
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

  Binomial trial size or negative-binomial size/dispersion parameter.
  May be a scalar, vector, function of `n`, or count-distribution list
  where supported.

- prob:

  Binomial success probability for `"single"` and `"level2"` generators.

- fixed_intercept:

  Link-scale intercept for direct-parameter defaults and multilevel
  models. Binomial uses the logit link; Poisson and negative binomial
  use the log link.

- random_cov:

  Group-level random-intercept covariance for multilevel generators.

- lambda:

  Poisson mean for `"single"` and `"level2"` generators.

- mu:

  Negative-binomial mean for `"single"` and `"level2"` generators.

- scale_fixed_intercept:

  Optional log size intercept for multilevel negative-binomial
  generators.

## Value

An `mlsim_generator_spec` object for use in
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

## Details

For multilevel count generators, `fixed_intercept` and `random_cov`
define the link-scale mean. Non-multilevel generators can use direct
distribution parameters or a link-scale `fixed_intercept`.

## See also

Other predictor generators:
[`continuous-generators`](https://florale.github.io/multilevelcoda/reference/continuous-generators.md),
[`gen_categorical()`](https://florale.github.io/multilevelcoda/reference/gen_categorical.md),
[`gen_custom()`](https://florale.github.io/multilevelcoda/reference/gen_custom.md),
[`gen_mvn()`](https://florale.github.io/multilevelcoda/reference/gen_mvn.md),
[`gen_normal()`](https://florale.github.io/multilevelcoda/reference/gen_normal.md)

## Examples

``` r
sim <- simulate_data(
  n = 5,
  seed = 4,
  generators = list(
    events = gen_poisson("events", lambda = 2),
    successes = gen_binomial("successes", size = 4, prob = 0.5),
    overdispersed = gen_negbin("overdispersed", size = 3, mu = 2)
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
