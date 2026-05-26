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
  scale,
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
      burnin = 0,
      composition = list(parts = c("sleep", "active", "sedentary"), total = 24)
    )
  )
)

params <- template_sim$generator_metadata$outcome_template$params
params$location$beta["treatmenttreatment", "ilr1"] <- 0.2
```
