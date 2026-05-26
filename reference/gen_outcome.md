# Generate Dynamic Gaussian and Compositional Outcomes

Create an outcome generator for
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).
Outcomes are simulated on a Gaussian scale, optionally as ILR
coordinates that are back-transformed to strictly positive compositions.
If `ar1()` appears in the location formula, it defines a residual VAR(1)
process, not an observed lagged-outcome predictor.

## Usage

``` r
gen_outcome(
  formula,
  scale,
  params,
  burnin,
  composition = list(parts = NULL, total = 24, sbp = NULL, keep_ilr = TRUE),
  ar_stability = c("resample", "shrink", "error"),
  max_stability_attempts = 1000,
  shrink_target_radius = 0.98
)
```

## Arguments

- formula:

  Outcome location formula. The left-hand side may be a single outcome
  (`y`) or `mvbind(y1, y2, ...)`. The right-hand side may include
  ordinary model terms, `between(x)`, `within(x)`, `ar1()`,
  interactions, and one brms/lme4-style grouping term.

- scale:

  Required scale formula with `sigma` on the left-hand side, for example
  `sigma ~ 1` or `sigma ~ treatment + (1 | ID)`. The scale model is on
  the log conditional standard-deviation scale.

- params:

  List of true parameters. Required components are
  `params$location$beta`, `params$scale$beta`, and, for multivariate
  outcomes, `params$scale$correlation`. When AR terms are present,
  `params$ar$beta` is required. When grouping terms are present,
  `params$random[[group]]$covariance` is required.

- burnin:

  Fixed non-negative integer burn-in length used when AR is active.
  Ignored when no AR terms are present.

- composition:

  List controlling optional ILR back-transformation. Use `parts` or
  `sbp` to request compositional output, `total` for the closure total,
  and `keep_ilr` to keep ILR coordinates alongside parts.

- ar_stability:

  Handling for unstable AR matrices: `"resample"`, `"shrink"`, or
  `"error"`.

- max_stability_attempts:

  Positive integer maximum number of resampling or shrinkage attempts.

- shrink_target_radius:

  Target spectral radius used by `ar_stability = "shrink"`.

## Value

An `mlsim_generator_spec` object for use in
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

## See also

Other predictor generators:
[`continuous-generators`](https://florale.github.io/multilevelcoda/reference/continuous-generators.md),
[`count-generators`](https://florale.github.io/multilevelcoda/reference/count-generators.md),
[`gen_categorical()`](https://florale.github.io/multilevelcoda/reference/gen_categorical.md),
[`gen_custom()`](https://florale.github.io/multilevelcoda/reference/gen_custom.md),
[`gen_mvn()`](https://florale.github.io/multilevelcoda/reference/gen_mvn.md),
[`gen_template()`](https://florale.github.io/multilevelcoda/reference/gen_template.md)

## Examples

``` r
beta_location <- matrix(
  c(0, 0, 0.2, -0.1),
  nrow = 2,
  dimnames = list(c("(Intercept)", "treatmenttreatment"), c("ilr1", "ilr2"))
)
beta_scale <- matrix(
  log(c(0.4, 0.35)),
  nrow = 1,
  dimnames = list("(Intercept)", c("ilr1", "ilr2"))
)
beta_ar <- array(
  c(0.25, 0.02, -0.01, 0.2),
  dim = c(1, 2, 2),
  dimnames = list("ar1()", c("ilr1", "ilr2"), c("ilr1", "ilr2"))
)
corr <- diag(2)
dimnames(corr) <- list(c("ilr1", "ilr2"), c("ilr1", "ilr2"))

sim <- simulate_data(
  n_groups = 4,
  n_per_group = 4,
  group_id = "ID",
  time_id = "day",
  seed = 1,
  generators = list(
    treatment = gen_categorical(
      "treatment",
      level = "level2",
      categories = c("control", "treatment"),
      fixed_intercept = stats::qlogis(0.5)
    ),
    outcome = gen_outcome(
      mvbind(ilr1, ilr2) ~ treatment + ar1(),
      scale = sigma ~ 1,
      params = list(
        location = list(beta = beta_location),
        scale = list(beta = beta_scale, correlation = corr),
        ar = list(beta = beta_ar)
      ),
      burnin = 10,
      composition = list(parts = c("sleep", "active", "sedentary"), total = 24)
    )
  )
)
head(sim$data)
#>       ID obs_id   day treatment        ilr1      ilr2    sleep    active
#>    <int>  <int> <int>    <fctr>       <num>     <num>    <num>     <num>
#> 1:     1      1     1   control -0.13966660 0.8202144 6.343236 13.442538
#> 2:     1      2     2   control -0.28904909 0.3408946 6.102631 11.064925
#> 3:     1      3     3   control -0.01135302 0.7371613 7.251181 12.382951
#> 4:     1      4     4   control  0.18305016 0.1609256 9.200034  8.238424
#> 5:     2      1     1   control -0.14991849 0.3436463 6.907676 10.582903
#> 6:     2      2     2   control -0.17536396 0.6210080 6.448046 12.399700
#>    sedentary
#>        <num>
#> 1:  4.214225
#> 2:  6.832444
#> 3:  4.365868
#> 4:  6.561543
#> 5:  6.509421
#> 6:  5.152254
```
