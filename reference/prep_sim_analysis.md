# Prepare Simulated Outcome Data for Analysis

Convert an
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md)
result with a
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
generator into an analysis-ready data set and inferred `brms` formula.
The helper recomputes observed between- and within-person predictors
from the raw simulated columns, creates lagged outcome predictor columns
for `ar1()` terms (person-mean-centered by default; see `lag_center`),
and creates a
[`complr()`](https://florale.github.io/multilevelcoda/reference/complr.md)
object when the outcome is compositional or when
`between()`/[`within()`](https://rdrr.io/r/base/with.html) terms
reference ILR coordinates of a compositional predictor generator.

## Usage

``` r
prep_sim_analysis(
  sim,
  outcome = NULL,
  drop_lag_na = FALSE,
  time_step = NULL,
  lag_center = c("within", "none"),
  centering = c("manifest", "latent"),
  link_random = TRUE
)
```

## Arguments

- sim:

  An `mlsim_data` object returned by
  [`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

- outcome:

  Optional character scalar naming the outcome generator. When `NULL`,
  the helper uses the only generator with `distribution = "outcome"`.

- drop_lag_na:

  Logical scalar. When `FALSE`, the default, first rows in each series
  are retained and lag-derived columns are left as `NA`. When `TRUE`,
  rows with missing lag-derived predictors are removed. If that removes
  every row of a group – for example because the group's time spacing
  does not match the (inferred) `time_step` – a warning lists the
  removed groups, which are also recorded in `metadata$dropped_groups`.

- time_step:

  Optional positive scalar giving the spacing between consecutive time
  points of the simulation `time_id`, used to build `ar1()` lag columns
  with
  [`lag_by_time()`](https://florale.github.io/multilevelcoda/reference/lag_by_time.md).
  When `NULL`, the default, the step is inferred as the smallest
  positive within-group time difference. If groups are spaced
  differently, the smallest step wins: every lag of a more widely spaced
  group is `NA` because that group has no observation one inferred step
  earlier. A warning then lists the affected groups; supply `time_step`
  explicitly to silence it when the chosen step is intended. For `Date`
  time the unit is days; for `POSIXct` time the unit is seconds. Ignored
  when the outcome formula has no `ar1()` term.

- lag_center:

  Character scalar controlling *whether* `ar1()` lag columns are
  centered; `centering` controls *which values* that centering uses.
  With `"within"`, the default, the lag is centered: under the default
  `centering = "manifest"` at the person's observed mean of the
  response, producing `lag_<response>_within` columns (cluster-mean
  centering, CMC), and under `centering = "latent"` at the latent
  conditional mean, which is the recorded latent residual state,
  producing `lag_<response>_latent` columns. With `"none"`, the raw
  lagged response is used as-is under either `centering`, producing
  `lag_<response>` columns (no centering, NC) – the non-centered
  parametrization of Hamaker & Grasman (2014). Under `"none"` the model
  intercept estimates `(1 - phi_i) * mu_i` (multivariate:
  `(I - Phi_i) mu_i`) rather than the person mean `mu_i`, so location
  parameters lose their direct correspondence with the simulation truth;
  see the "Pragmatic default estimator" section. Ignored when the
  outcome formula has no `ar1()` term.

- centering:

  Character scalar controlling which values the analysis centers with,
  for `between()` and [`within()`](https://rdrr.io/r/base/with.html)
  columns and for a centered `ar1()` lag alike. With `"manifest"`, the
  default, `between()`/[`within()`](https://rdrr.io/r/base/with.html)
  columns are recomputed from the observed data as applied analysts
  would – person means for scalars and the closed arithmetic-mean
  composition for compositional predictors – which are observed proxies
  for the generating components rather than the components themselves.
  With `"latent"`, they are taken from the generating components the
  simulator recorded, so they match the simulated quantity exactly.
  `"latent"` requires a grouped design and generators that emitted
  latent components, and is available only in simulated data; see the
  "Latent centering" section. When the formula contains `ar1()` and the
  lag is centered (`lag_center = "within"`, the default), `"latent"`
  also centers the lag with latent values, which makes it the recorded
  latent residual state; see the "Latent centering and `ar1()`" section
  for why even that oracle leaves the autoregressive parameters biased
  in models with group-level effects. With `lag_center = "none"` there
  is no lag centering for `centering` to act on, so the two settings
  compose freely. Ignored when the outcome formula has no `between()`,
  [`within()`](https://rdrr.io/r/base/with.html), or `ar1()` terms.

- link_random:

  Logical scalar. When `TRUE`, the default, grouping factors that appear
  in both the mean and the scale formulas emit brms ID-linked random
  effects (for example `(1 | p1 | ID)`), so the analysis model can
  estimate the cross-parameter random-effect correlations that the
  simulator generates from one joint covariance. Set to `FALSE` to emit
  separate, uncorrelated random-effect blocks instead; this drops those
  correlations from the analysis model but keeps large models tractable
  for approximate algorithms (variational inference and pathfinder often
  fail on large linked random-effect blocks).

## Value

An `mlsim_analysis` object, a list with:

- `data`:

  Analysis data with derived between/within and lag columns.

- `formula`:

  An inferred `brms` formula.

- `complr`:

  A `complr` object covering every composition used by the analysis
  model (compositional predictors referenced through
  `between()`/[`within()`](https://rdrr.io/r/base/with.html) and, when
  applicable, the compositional outcome), otherwise `NULL`. This object
  can be passed directly to
  [`brmcoda()`](https://florale.github.io/multilevelcoda/reference/brmcoda.md)
  and used with
  [`substitution()`](https://florale.github.io/multilevelcoda/reference/substitution.md).

- `truth`:

  A
  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  of the true generating parameter values labeled with the parameter
  names of the inferred analysis model, ready to be compared against the
  fitted model with
  [`sim_recovery()`](https://florale.github.io/multilevelcoda/reference/sim_recovery.md).
  `NULL` when the analysis model cannot be aligned with the simulator
  parameters (the reason is recorded in
  `metadata$truth_unavailable_reason`).

- `metadata`:

  Preparation metadata, including derived column names and formula
  mappings.

## Details

The inferred analysis model is a deliberately *pragmatic default
estimator* built from observed data, not the matched model for the
simulated data-generating process. See the "Pragmatic default estimator"
section before interpreting parameter recovery results.

The analysis formula is inferred from the stored
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
formula. Simulation terms `between(x)` and `within(x)` become
observed-data columns named `x_between` and `x_within`, computed from
`x` by the simulation group identifier. These columns are recomputed
even when columns with the same names already exist in `sim$data`.

When `between()`/[`within()`](https://rdrr.io/r/base/with.html)
reference ILR coordinates of a compositional
[`gen_mvn()`](https://florale.github.io/multilevelcoda/reference/gen_mvn.md)
predictor, the helper instead builds a
[`complr()`](https://florale.github.io/multilevelcoda/reference/complr.md)
object for that predictor composition (using the simulator's parts, SBP,
and total) and maps the terms to the complr between/within coordinates
(`bz<k>_<j>` and `wz<k>_<j>`, where `k` is the coordinate and `j` the
composition index), exactly as in the standard `multilevelcoda`
workflow. When the outcome is also compositional, one multi-composition
`complr` object covers both the predictor and outcome compositions.

The returned `brms` formula carries the family recorded by the outcome
generator ([`gaussian()`](https://rdrr.io/r/stats/family.html),
[`poisson()`](https://rdrr.io/r/stats/family.html),
[`binomial()`](https://rdrr.io/r/stats/family.html),
[`brms::negbinomial()`](https://paulbuerkner.com/brms/reference/brmsfamily.html),
`Gamma(link = "log")`, or
[`brms::Beta()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)).
When the generator has a scale model, it becomes the matching
distributional formula (`sigma`, `shape`, or `phi`). Binomial outcomes
use `y | trials(<outcome>_trials)` with the trials column generated by
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md).
When the same grouping factor appears in both the mean and the scale
formulas, the random effects are emitted with brms ID-linked syntax (for
example `(1 | p1 | ID)`) by default, because the simulator draws these
group-level effects from one joint covariance and the linked syntax lets
the analysis model estimate their correlation; see `link_random` to opt
out for large models fitted with approximate algorithms.

For dynamic formulas, `ar1()` is translated to lagged response columns
whose construction and name follow the two centering arguments:

- `lag_center = "within"`, `centering = "manifest"` (both defaults): the
  lagged observed response centered at each person's observed mean of
  all their values (not the mean of the lagged values),
  `lag_<response>_within`.

- `lag_center = "within"`, `centering = "latent"`: the lagged response
  centered at the latent conditional mean, which is the recorded latent
  residual state, `lag_<response>_latent`.

- `lag_center = "none"`, either `centering`: the raw lagged observed
  response, `lag_<response>`.

For compositional outcomes, the helper rebuilds the ILR coordinates
through
[`complr()`](https://florale.github.io/multilevelcoda/reference/complr.md)
using the simulator's parts and SBP metadata, then lags the generated
`z` coordinates used by
[`brmcoda()`](https://florale.github.io/multilevelcoda/reference/brmcoda.md).

The lags are built by
[`lag_by_time()`](https://florale.github.io/multilevelcoda/reference/lag_by_time.md)
as *time-based* lags on the simulation `time_id`: a row's lag comes from
the observation at `time - time_step`, not from the previous row. When a
time point is missing – as can happen when analysing real data with
skipped observations – the following row's lag is `NA` rather than the
value from before the gap. The autoregressive coefficient therefore
refers to one `time_step`. For data simulated with `ar1()` the time grid
is complete and equally spaced, so the result is identical to a
positional shift. The number of gap-affected rows is recorded in
`metadata$lag_gaps`.

Because a single `time_step` is applied to all groups, groups whose own
spacing is wider than the inferred step end up with only `NA` lags, and
`drop_lag_na = TRUE` would then remove them from the analysis data
entirely. `prep_sim_analysis()` warns in both situations and records the
diagnostics in the metadata: `metadata$time_step_by_group` (each group's
smallest positive spacing, `NA` for single-observation groups),
`metadata$time_step_heterogeneous` (whether those spacings differ),
`metadata$lag_na_groups` (groups without a single usable lagged row),
and `metadata$dropped_groups` (groups removed by `drop_lag_na`).

## Pragmatic default estimator

[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
simulates latent residual AR/VAR dynamics around the model-implied mean
and resolves `between()`/[`within()`](https://rdrr.io/r/base/with.html)
from latent generating components supplied by upstream predictor
generators. By default (`centering = "manifest"`) `prep_sim_analysis()`
instead constructs the model applied analysts commonly fit to observed
data:

- `between(x)` and `within(x)` become `x_between` and `x_within`,
  recomputed from realised person means of the observed `x` (manifest
  centering), not from the latent generating components.

- `ar1()` becomes, by default, person-mean-centered lagged *observed*
  response predictors (`lag_<response>_within`), not the latent residual
  state.

- For ILR coordinates of compositional predictors, `between(ilr)` and
  `within(ilr)` become the
  [`complr()`](https://florale.github.io/multilevelcoda/reference/complr.md)
  coordinates `bz*`/`wz*`, where the between composition is the ILR of
  each person's closed *arithmetic-mean* composition of the observed
  parts. The simulator's `between(ilr)`, in contrast, is the latent
  group-level ILR mean, whose back-transform is the closed *geometric*
  (Aitchison) center. These are different estimands, and the gap has two
  parts: a structural Jensen gap between the arithmetic and geometric
  centers, which does *not* shrink as the number of observations per
  group grows, plus ordinary sampling error, which does. Manifest
  centering is therefore biased for latent between-person effects even
  in long series, unlike the scalar case where only sampling error is
  involved.

These observed-data constructions target different estimands from the
simulation truth. Manifest person-mean centering and observed-score
lagged regression are known to yield biased estimates of between-person
effects and of inertia and cross-lag parameters relative to the latent
generating values, especially with short series (Ludtke et al. 2008;
Hamaker & Grasman 2014). This mismatch is intentional: it lets
simulation studies quantify the bias of the pragmatic estimator. Do not
interpret systematic discrepancies between estimates from this default
analysis model and the
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
truth parameters as errors in the simulator.

## Latent centering

`centering = "latent"` replaces the manifest construction with the
generating components themselves, so for the affected
`between()`/[`within()`](https://rdrr.io/r/base/with.html) terms the
fitted model targets exactly the quantity the simulator recorded:

- Scalar predictors take `x_between` from the column the generator
  labelled with component `"between"` – for `level = "level2"`
  predictors that is the variable itself – and set
  `x_within <- x - x_between`.

- Compositional predictors keep the returned
  [`complr()`](https://florale.github.io/multilevelcoda/reference/complr.md)
  object, but its between and within blocks (`bX`, `wX`, `bZ`, `wZ`, and
  the matching `b<part>`/`w<part>`/`bz*`/`wz*` columns) are rebuilt from
  the latent between composition emitted by the generator. Column names
  are unchanged, so
  [`brmcoda()`](https://florale.github.io/multilevelcoda/reference/brmcoda.md)
  and
  [`substitution()`](https://florale.github.io/multilevelcoda/reference/substitution.md)
  work as usual.

This is an oracle: the latent components exist only in simulated data,
so `"latent"` has no counterpart in a real analysis. Use it to separate
estimator bias from simulator error – if a parameter is recovered under
`"latent"` but not under `"manifest"`, the gap is the estimator's, not
the simulator's. It requires a grouped design and a generator that
emitted latent components;
[`gen_custom()`](https://florale.github.io/multilevelcoda/reference/gen_custom.md)
generators that label roles without emitting them can only be analysed
with `"manifest"`.

## Latent centering and `ar1()`

When the outcome formula contains `ar1()` and the lag is centered
(`lag_center = "within"`, the default), `centering = "latent"` centers
the lag with latent values, which replaces the lag predictor entirely:
instead of a lagged *observed* response centered at the person's
observed mean, the analysis uses the latent residual state the simulator
recorded, in columns named `lag_<response>_latent`. The two
constructions are the same operation applied with different values,
because the recorded residual is `y_t - mu_t`, the response centered at
its latent conditional mean. The latent state needs no further
centering: person-mean centering it on top would reintroduce exactly the
attenuation this construction avoids.

With `lag_center = "none"` there is no lag centering for `centering` to
act on, so the raw lagged observed response is used under `"latent"`
just as it is under `"manifest"`. That combination is the point of
keeping the two arguments separate: a model can carry
`between()`/[`within()`](https://rdrr.io/r/base/with.html) predictors
that have nothing to do with the lag, and those can be built from the
latent generating components while the AR part keeps the non-centered
parametrization discussed below.

A latent lag makes the AR, scale and `rescor` parameters estimate
exactly what the simulator recorded **only when the analysis model has
no group-level effects**. The reason is not the proxy problem that
`lag_center` addresses. Writing the generating model as
`z_t = X_t beta + Phi e_{t-1} + eps_t`, the regressor `e_{t-1}` is
*predetermined* but not *strictly exogenous*: `Cov(eps_t, e_t)` equals
the innovation variance, and `e_t` is the next row's regressor. Any
estimator that mixes rows within a group – which is what the
random-effects GLS transformation implied by a `(1 | group)` term does –
picks up that cross-row covariance and incurs an O(1/T) dynamic-panel
(Nickell) bias. Substituting the true latent state removes the
measurement problem but not this one. So with any group-level effect, no
parameter of a dynamic model is guaranteed to keep its generating
estimand, under either centering.

For a matched-model recovery of a *multilevel* dynamic model, the
appropriate analysis is not lag regression at all but a residual
autocorrelation structure, for example
`brms::bf(y ~ 1 + (1 | ID)) + brms::ar(time = day, gr = ID, p = 1, cov = TRUE)`.
That form is outside what this helper builds – brms
[`ar()`](https://rdrr.io/r/stats/ar.html) supports only a scalar
unmoderated AR per response and parametrizes `sigma` as the marginal
rather than the innovation standard deviation – so construct it by hand.

Note that under `"manifest"` the latent columns of the same name are
overwritten in `analysis$data`; `sim$data` always keeps the originals.

Whether the `ar1()` lag columns are centered at all is controlled by
`lag_center`, independently of `centering`. Hamaker & Grasman (2014)
show that in multilevel autoregressive models, cluster-mean centering
the lagged outcome (the `"within"` default) attenuates the average
autoregressive coefficient, whereas the non-centered lagged outcome
(`"none"`) recovers it nearly unbiased. The trade-off is that under
`"none"` the whole location mean structure is reparametrized: the
intercept estimates `(I - Phi_i) mu_i` rather than the person mean, and
any location covariate effects absorb omitted `-Phi x_{t-1} beta` terms.
Centering does not fully resolve this either: the centered observed lag
still carries the lagged mean structure, so when a predictor varies
within series its current coefficient can absorb omitted
lagged-covariate terms under `"within"` as well. Whenever the formula
contains `ar1()` and the lag is an observed response, then, no parameter
is guaranteed to keep its generating estimand, under either `lag_center`
setting. Use `"none"` when the average autoregressive (inertia and
cross-lag) coefficients are the estimands of interest, and the default
`"within"` when person means and their predictors must remain
interpretable; the data-simulation vignette compares parameter recovery
under both.

## References

Ludtke, O., et al. (2008). The multilevel latent covariate model.
*Psychological Methods*, 13(3), 203-229.

Hamaker, E. L., & Grasman, R. P. P. P. (2014). To center or not to
center? Investigating inertia with a multilevel autoregressive model.
*Frontiers in Psychology*, 5, 1492.

## Examples

``` r
params <- list(
  location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
  scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y")))
)
sim <- simulate_data(
  n = 5,
  seed = 1,
  generators = list(
    outcome = gen_outcome(
      y ~ 1,
      scale = sigma ~ 1,
      params = params
    )
  )
)
analysis <- prep_sim_analysis(sim)
analysis$formula
#> y ~ 1 
#> sigma ~ 1
```
