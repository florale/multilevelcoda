# Align Fitted Estimates with Simulation Truth

`sim_recovery()` joins the posterior summaries of a model fitted to
simulated data with the true generating parameter values recorded by
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md),
using the brms parameter names on both sides. It replaces the
error-prone manual name matching otherwise needed to compare
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
truth with
[`brms::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html) and
[`brms::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) output.

## Usage

``` r
sim_recovery(fit, analysis, probs = c(0.025, 0.975))
```

## Arguments

- fit:

  A
  [`brmcoda()`](https://florale.github.io/multilevelcoda/reference/brmcoda.md)
  model object or
  [brms::brmsfit](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  fitted to `analysis$data` (or `analysis$complr`) using
  `analysis$formula`.

- analysis:

  An `mlsim_analysis` object from
  [`prep_sim_analysis()`](https://florale.github.io/multilevelcoda/reference/prep_sim_analysis.md)
  with a non-`NULL` `$truth` element.

- probs:

  A length-two numeric vector of lower and upper posterior interval
  probabilities used for the coverage indicator. Defaults to
  `c(0.025, 0.975)`.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
of class `mlsim_recovery` with one row per model parameter:

- `type`: `"fixed"`, `"random_sd"`, `"random_cor"`, or `"rescor"`.
  Random-effect rows are reported on the standard deviation /
  correlation scale used by
  [`brms::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html).

- `group`: grouping factor for random-effect rows, `NA` otherwise.

- `parameter`: the brms-style parameter label.

- `truth`: the true generating value from the simulator
  ([`sqrt()`](https://rdrr.io/r/base/MathFun.html) of covariance
  diagonals for `random_sd` rows, correlations for `random_cor` rows).

- `estimate`, `est_error`, `lower`, `upper`: posterior summary of the
  fitted parameter at `probs`.

- `bias`: `estimate - truth`.

- `covered`: whether `truth` lies inside `[lower, upper]`.

- `simulator_name`: the simulator-side parameter identifier the row was
  derived from, for provenance.

## Details

The truth table is constructed by
[`prep_sim_analysis()`](https://florale.github.io/multilevelcoda/reference/prep_sim_analysis.md)
and stored in its `$truth` element; `sim_recovery()` matches every truth
row to exactly one fitted parameter and vice versa. Any parameter that
cannot be matched in either direction is an error, never a silently
wrong row: truth values are attached purely by name, and duplicate or
ambiguous names (for example response names that collide after brms
removes `_` and `.` characters) abort with an explanatory message.

All
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
families are supported, not only Gaussian models. Rows for a modeled
distributional parameter are labeled with the family's brms parameter
name (`sigma_*` for `"gaussian"`, `shape_*` for `"negbin"` and
`"gamma"`, `phi_*` for `"beta"`), and their truth values are on the log
link scale that both
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
and brms use for these parameters; location rows are on the family's
location link (identity, log, or logit), again matching on both sides.
`"poisson"` and `"binomial"` families have no distributional parameter,
so their tables contain no scale rows, and a binomial `trials()` term
contributes no parameter. Because `ar1()`, multivariate `mvbind()`
responses, and compositional outcomes are Gaussian-only simulator
features, recovery tables for the other families never contain `rescor`
rows.

Two situations yield no recovery row at all: cross-block random-effect
correlations when the analysis was prepared with
[`prep_sim_analysis()`](https://florale.github.io/multilevelcoda/reference/prep_sim_analysis.md)`(link_random = FALSE)`
(the fitted model never samples them), and analysis models whose
structure cannot be aligned with the simulator's single joint covariance
(recorded in `metadata$truth_unavailable_reason`).

Interpreting the `bias` and `covered` columns requires knowing what the
fitted model estimates.
[`prep_sim_analysis()`](https://florale.github.io/multilevelcoda/reference/prep_sim_analysis.md)
deliberately builds the pragmatic observed-data model that applied
analysts commonly fit, not the matched model for the generating process,
so some parameters estimate a related quantity rather than the one the
simulator recorded – most notably whenever the formula contains `ar1()`,
and for `between()`/[`within()`](https://rdrr.io/r/base/with.html) terms
under the default `centering = "manifest"`. Systematic bias in those
parameters is an expected property of the estimator, not an error in the
simulator. See the "Pragmatic default estimator" section of
[`prep_sim_analysis()`](https://florale.github.io/multilevelcoda/reference/prep_sim_analysis.md)
before drawing conclusions from a recovery table.

## See also

[`prep_sim_analysis()`](https://florale.github.io/multilevelcoda/reference/prep_sim_analysis.md)
for the analysis object and its `$truth` table,
[`gen_template()`](https://florale.github.io/multilevelcoda/reference/gen_template.md)
for authoring the truth values.

## Examples

``` r
# \donttest{
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  params <- list(
    location = list(beta = matrix(
      c(0.5, 0.2),
      nrow = 2,
      dimnames = list(c("(Intercept)", "x"), "y")
    )),
    scale = list(beta = matrix(
      log(0.5),
      nrow = 1,
      dimnames = list("(Intercept)", "y")
    ))
  )
  sim <- simulate_data(
    n = 200,
    seed = 1,
    generators = list(
      x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
      y = gen_outcome(y ~ x, scale = sigma ~ 1, params = params)
    )
  )
  analysis <- prep_sim_analysis(sim)
  analysis$truth

  fit <- brms::brm(
    analysis$formula,
    data = analysis$data,
    backend = "cmdstanr",
    refresh = 0
  )
  sim_recovery(fit, analysis)
}
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
# }
```
