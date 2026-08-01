# Create a Parameter Template for `gen_outcome()`

Build a generator that parses a
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
specification and records a complete parameter template without
simulating outcome values. Use this inside
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md)
after any generators that define predictors used by the outcome formula,
especially predictors referenced by `between()` or
[`within()`](https://rdrr.io/r/base/with.html).

## Usage

``` r
gen_template(
  formula,
  scale = NULL,
  burnin = NULL,
  family = "gaussian",
  composition = list(parts = NULL, total = 24, sbp = NULL, keep_ilr = TRUE),
  ar_stability = c("resample", "shrink", "error"),
  max_stability_attempts = 1000,
  shrink_target_radius = 0.98
)
```

## Arguments

- formula:

  Outcome location formula. The left-hand side may be a single outcome
  (`y`) or, for `family = "gaussian"` only, `mvbind(y1, y2, ...)`. The
  right-hand side may include ordinary model terms, `between(x)`,
  `within(x)`, interactions, one brms/lme4-style grouping term, and, for
  `family = "gaussian"` only, `ar1()`.
  `between()`/[`within()`](https://rdrr.io/r/base/with.html) resolve
  against the column roles emitted by earlier predictor generators; this
  includes the ILR coordinates of multilevel compositional
  [`gen_mvn()`](https://florale.github.io/multilevelcoda/reference/gen_mvn.md)
  generators (for example `between(ilr1)`), enabling simulation of
  models where a scalar outcome depends on the between- and
  within-person parts of a composition. Roles apply to ILR coordinates,
  not to composition parts.

- scale:

  Scale formula with the family's brms distributional-parameter name on
  the left-hand side: `sigma` for `"gaussian"` (log conditional SD, for
  example `sigma ~ 1` or `sigma ~ treatment + (1 | ID)`), `shape` for
  `"negbin"` (log size) and `"gamma"` (log shape), and `phi` for
  `"beta"` (log precision). Required for these families; must be omitted
  for `"poisson"` and `"binomial"`, which have no auxiliary scale
  parameter.

- burnin:

  Fixed integer burn-in length for the residual AR/VAR process. Required
  when `ar1()` appears in the location formula and ignored otherwise
  (non-AR Gaussian models and all non-Gaussian families have no burn-in
  phase). Each series is initialized at zero residuals and iterated
  `burnin` steps before the first observed row, so `burnin = 0` starts
  the residual process exactly at zero and the first observations are
  under-dispersed relative to the stationary distribution. Do not use
  `burnin = 0` or other small values for simulation studies: choose a
  burn-in long enough for the process to forget its start, for example
  at least `log(0.01) / log(rho)` steps, where `rho` is the largest
  spectral radius of the AR matrices (about 90 steps at `rho = 0.95`).

- family:

  Character scalar naming the outcome family: `"gaussian"` (default),
  `"poisson"` (log link), `"binomial"` (logit link, requires `trials`),
  `"negbin"` (log link, `scale` is log size), `"gamma"` (log link,
  `scale` is log shape), or `"beta"` (logit link, `scale` is log
  precision). Non-Gaussian families are univariate only and do not
  support `ar1()` or `composition`.

- composition:

  List controlling optional ILR back-transformation
  (`family = "gaussian"` only). Use `parts` or `sbp` to request
  compositional output, `total` for the closure total, and `keep_ilr` to
  keep ILR coordinates alongside parts.

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
It emits no data columns and stores the parameter template in generator
metadata.

## Details

The generated template is stored at
`sim$generator_metadata[[name]]$params`, where `name` is the name given
to this generator in the `generators` list. The template must be created
with the same simulation design, previous generators, factor levels,
formula, scale formula, and composition settings that will be used for
the later
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
call.

Template values are initialized to zero for regression, AR, and
group-level covariance parameters, and to identity for residual ILR
correlations. The object is ready to edit and pass as `params` to
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md).

Templates are family-aware: `params$scale` is only included for families
with an auxiliary scale model (`"gaussian"`, `"negbin"`, `"gamma"`,
`"beta"`), `params$scale$correlation` only for `"gaussian"`, and
`params$ar` only when `ar1()` appears (which requires
`family = "gaussian"`). `trials` is not needed to build a binomial
template because no parameter block depends on it.

## See also

Other predictor generators:
[`continuous-generators`](https://florale.github.io/multilevelcoda/reference/continuous-generators.md),
[`count-generators`](https://florale.github.io/multilevelcoda/reference/count-generators.md),
[`gen_categorical()`](https://florale.github.io/multilevelcoda/reference/gen_categorical.md),
[`gen_custom()`](https://florale.github.io/multilevelcoda/reference/gen_custom.md),
[`gen_mvn()`](https://florale.github.io/multilevelcoda/reference/gen_mvn.md),
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)

## Examples

``` r
template_sim <- simulate_data(
  n_groups = 3,
  n_per_group = 2,
  group_id = "ID",
  seed = 1,
  generators = list(
    treatment = gen_categorical(
      "treatment",
      level = "level2",
      categories = c("control", "treatment"),
      fixed_intercept = stats::qlogis(0.5)
    ),
    outcome_template = gen_template(
      mvbind(ilr1, ilr2) ~ treatment,
      scale = sigma ~ 1,
      composition = list(parts = c("sleep", "active", "sedentary"), total = 24)
    )
  )
)

params <- template_sim$generator_metadata$outcome_template$params
params$location$beta["treatmenttreatment", "ilr1"] <- 0.2
```
