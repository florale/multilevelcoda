# Generate Dynamic Gaussian, Compositional, and GLM-Family Outcomes

Create an outcome generator for
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).
With the default `family = "gaussian"`, outcomes are simulated on a
Gaussian scale, optionally as ILR coordinates that are back-transformed
to strictly positive compositions, and `ar1()` in the location formula
defines a residual VAR(1) process, not an observed lagged-outcome
predictor. Non-Gaussian families (`"poisson"`, `"binomial"`, `"negbin"`,
`"gamma"`, `"beta"`) simulate a single outcome from the family
distribution with the location linear predictor mapped through the
family link.

## Usage

``` r
gen_outcome(
  formula,
  scale = NULL,
  params,
  burnin = NULL,
  family = "gaussian",
  trials = NULL,
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

- params:

  List of true parameters. Required components are
  `params$location$beta` and, for families with a scale model,
  `params$scale$beta`. For multivariate Gaussian outcomes,
  `params$scale$correlation` is required. When AR terms are present,
  `params$ar$beta` is required. When grouping terms are present,
  `params$random[[group]]$covariance` is required.

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

- trials:

  Number of binomial trials, required for and only allowed with
  `family = "binomial"`. Accepts a scalar, a vector recycled to the
  number of rows, a function of `n`, or a count-distribution list as in
  [`gen_binomial()`](https://florale.github.io/multilevelcoda/reference/count-generators.md)
  (count-distribution `minimum`/`maximum` bounds clamp draws to the
  bounds; censoring, not truncation). Resolved trial sizes are written
  into the simulated data as a `<outcome>_trials` column for use with
  `y | trials(y_trials)` in brms.

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

## Non-Gaussian families

Non-Gaussian outcomes are drawn independently per row from the family
distribution after mapping the location linear predictor through the
inverse link, with the auxiliary parameter (size, shape, or precision)
taken as [`exp()`](https://rdrr.io/r/base/Log.html) of the scale linear
predictor. `ar1()` is currently only supported for
`family = "gaussian"`. Group-level effects on the location and scale
linear predictors are drawn jointly from
`params$random[[group]]$covariance` exactly as for Gaussian outcomes.
Because the links are nonlinear, fixed effects are conditional
(subject-specific) effects: the marginal mean of the outcome is not the
inverse link of the fixed-effect predictor alone. Extreme
linear-predictor values that produce non-finite means or auxiliary
parameters are an error; draws that land on a boundary the family
likelihood cannot support (beta outcomes numerically at 0 or 1, gamma
outcomes at 0) produce a warning.

## Time spacing

When `ar1()` is present, time must be complete and equally spaced within
each participant (or within the single series). Participants may have
different start times, end times, and numbers of observations, and the
simulator does not check or enforce equal spacing between participants:
different participants may also use different step sizes. AR and VAR
coefficients are defined per observation step, not per unit of real
time, so dynamic parameters are only comparable across participants in
real-time units when all participants share the same step size.

## AR stability and realized moderator draws

Stability is enforced through the row-wise spectral radius of the
assembled AR coefficient matrices for every observed row. When `ar1()`
interacts with predictors (for example `within(stress):ar1()`), the
row-specific AR matrix depends on the moderator values realized earlier
in the generator pipeline. Stability is therefore a property of the AR
parameters jointly with the realized predictor data, not of the
parameters alone: the same AR parameter values may be accepted in one
simulated data set and rejected in another with more extreme moderator
draws, and the chance of an unstable row grows with the number of rows.
`ar_stability = "resample"` redraws only group-level effects, so it
cannot repair instability caused by the population-level part of a
moderated AR term; that case errors instead.

Two further diagnostics guard the validity of accepted draws. First,
when the AR matrices vary over time within a series (moderated `ar1()`
with more than one outcome), a spectral radius below 1 for every row is
necessary but not sufficient for stability of the time-varying process:
products of individually stable matrices can still diverge. The
simulator therefore also computes the spectral norm (largest singular
value) of each row-specific AR matrix; a largest norm below 1 guarantees
stability. When the matrices vary within a series and the largest
spectral norm is 1 or more, stability cannot be guaranteed and a warning
is issued; the result is recorded in the stability metadata as
`time_varying_stability_guaranteed`. Second, `ar_stability = "resample"`
truncates the group-level effect distribution: each group's entire
random-effect vector (location, scale, and AR effects jointly) is
redrawn until the AR part is stable, so when the acceptance rate is
below 1 the realized random-effect covariance is smaller than
`params$random[[group]]$covariance` requested. The realized covariance
is always recorded in the generator metadata as
`realised_group_level_effect_covariance`, and a message is emitted when
the acceptance rate falls below 0.95.

## Performance

The simulator is written to scale to large designs without changing the
data-generating model. Innovations are drawn in one vectorized step from
the fixed conditional correlation matrix and scaled by the row-wise
conditional SDs (an exact draw from the row-specific innovation
distribution, since its covariance is the SD-scaled correlation matrix).
Row-specific AR matrices are assembled with vectorized array operations,
spectral radii are computed once per unique AR matrix so row-constant AR
designs cost one eigendecomposition per participant, and stability
resampling or shrinkage re-evaluates only the affected participant's
rows. Because the order in which random numbers are consumed is part of
the implementation, a given seed maps to a particular realization only
within a package version; the distribution of simulated data is
unaffected.

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
#>       ID obs_id   day treatment       ilr1        ilr2     sleep    active
#>    <int>  <int> <int>    <fctr>      <num>       <num>     <num>     <num>
#> 1:     1      1     1   control -0.2064913 -0.03609880  6.710760  8.424008
#> 2:     1      2     2   control -0.2868224  0.70699869  5.708368 13.371705
#> 3:     1      3     3   control -0.4443665  0.41099101  5.225170 12.041240
#> 4:     1      4     4   control -0.4260561 -0.05385306  5.488733  8.903350
#> 5:     2      1     1   control -0.1319125  0.28382912  7.062593 10.145892
#> 6:     2      2     2   control  0.7619243  0.47254122 13.108039  7.200856
#>    sedentary
#>        <num>
#> 1:  8.865232
#> 2:  4.919927
#> 3:  6.733590
#> 4:  9.607916
#> 5:  6.791514
#> 6:  3.691105
```
