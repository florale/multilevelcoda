#' Prepare Simulated Outcome Data for Analysis
#'
#' Convert an [simulate_data()] result with a [gen_outcome()] generator into an
#' analysis-ready data set and inferred `brms` formula. The helper recomputes
#' observed between- and within-person predictors from the raw simulated
#' columns, creates lagged outcome predictor columns for `ar1()` terms
#' (person-mean-centered by default; see `lag_center`), and
#' creates a [complr()] object when the outcome is compositional or when
#' `between()`/`within()` terms reference ILR coordinates of a compositional
#' predictor generator.
#'
#' The inferred analysis model is a deliberately *pragmatic default estimator*
#' built from observed data, not the matched model for the simulated
#' data-generating process. See the "Pragmatic default estimator" section
#' before interpreting parameter recovery results.
#'
#' @param sim An `mlsim_data` object returned by [simulate_data()].
#' @param outcome Optional character scalar naming the outcome generator. When
#'   `NULL`, the helper uses the only generator with
#'   `distribution = "outcome"`.
#' @param drop_lag_na Logical scalar. When `FALSE`, the default, first rows in
#'   each series are retained and lag-derived columns are left as `NA`. When
#'   `TRUE`, rows with missing lag-derived predictors are removed. If that
#'   removes every row of a group -- for example because the group's time
#'   spacing does not match the (inferred) `time_step` -- a warning lists the
#'   removed groups, which are also recorded in `metadata$dropped_groups`.
#' @param time_step Optional positive scalar giving the spacing between
#'   consecutive time points of the simulation `time_id`, used to build
#'   `ar1()` lag columns with [lag_by_time()]. When `NULL`, the default, the
#'   step is inferred as the smallest positive within-group time difference.
#'   If groups are spaced differently, the smallest step wins: every lag of a
#'   more widely spaced group is `NA` because that group has no observation
#'   one inferred step earlier. A warning then lists the affected groups;
#'   supply `time_step` explicitly to silence it when the chosen step is
#'   intended. For `Date` time the unit is days; for `POSIXct` time the unit
#'   is seconds. Ignored when the outcome formula has no `ar1()` term.
#' @param lag_center Character scalar controlling *whether* `ar1()` lag columns
#'   are centered; `centering` controls *which values* that centering uses.
#'   With `"within"`, the default, the lag is centered: under the default
#'   `centering = "manifest"` at the person's observed mean of the response,
#'   producing `lag_<response>_within` columns (cluster-mean centering, CMC),
#'   and under `centering = "latent"` at the latent conditional mean, which is
#'   the recorded latent residual state, producing `lag_<response>_latent`
#'   columns. With `"none"`, the raw lagged response is used as-is under either
#'   `centering`, producing `lag_<response>` columns (no centering, NC) -- the
#'   non-centered parametrization of Hamaker & Grasman (2014). Under `"none"`
#'   the model intercept estimates `(1 - phi_i) * mu_i` (multivariate:
#'   `(I - Phi_i) mu_i`) rather than the person mean `mu_i`, so location
#'   parameters lose their direct correspondence with the simulation truth;
#'   see the "Pragmatic default estimator" section. Ignored when the outcome
#'   formula has no `ar1()` term.
#' @param centering Character scalar controlling which values the analysis
#'   centers with, for `between()` and `within()` columns and for a centered
#'   `ar1()` lag alike. With `"manifest"`, the default, `between()`/`within()`
#'   columns are recomputed from the observed data as applied analysts would --
#'   person means for scalars and the closed arithmetic-mean composition for
#'   compositional predictors -- which are observed proxies for the generating
#'   components rather than the components themselves. With `"latent"`, they
#'   are taken from the generating components the simulator recorded, so they
#'   match the simulated quantity exactly. `"latent"` requires a grouped design
#'   and generators that emitted latent components, and is available only in
#'   simulated data; see the "Latent centering" section. When the formula
#'   contains `ar1()` and the lag is centered (`lag_center = "within"`, the
#'   default), `"latent"` also centers the lag with latent values, which makes
#'   it the recorded latent residual state; see the "Latent centering and
#'   `ar1()`" section for why even that oracle leaves the autoregressive
#'   parameters biased in models with group-level effects. With
#'   `lag_center = "none"` there is no lag centering
#'   for `centering` to act on, so the two settings compose freely. Ignored
#'   when the outcome formula has no `between()`, `within()`, or `ar1()` terms.
#' @param link_random Logical scalar. When `TRUE`, the default, grouping
#'   factors that appear in both the mean and the scale formulas emit brms
#'   ID-linked random effects (for example `(1 | p1 | ID)`), so the analysis
#'   model can estimate the cross-parameter random-effect correlations that
#'   the simulator generates from one joint covariance. Set to `FALSE` to
#'   emit separate, uncorrelated random-effect blocks instead; this drops
#'   those correlations from the analysis model but keeps large models
#'   tractable for approximate algorithms (variational inference and
#'   pathfinder often fail on large linked random-effect blocks).
#'
#' @return An `mlsim_analysis` object, a list with:
#'   \item{`data`}{Analysis data with derived between/within and lag columns.}
#'   \item{`formula`}{An inferred `brms` formula.}
#'   \item{`complr`}{A `complr` object covering every composition used by the
#'   analysis model (compositional predictors referenced through
#'   `between()`/`within()` and, when applicable, the compositional outcome),
#'   otherwise `NULL`. This object can be passed directly to [brmcoda()] and
#'   used with [substitution()].}
#'   \item{`truth`}{A [data.table::data.table] of the true generating
#'   parameter values labeled with the parameter names of the inferred
#'   analysis model, ready to be compared against the fitted model with
#'   [sim_recovery()]. `NULL` when the analysis model cannot be aligned with
#'   the simulator parameters (the reason is recorded in
#'   `metadata$truth_unavailable_reason`).}
#'   \item{`metadata`}{Preparation metadata, including derived column names and
#'   formula mappings.}
#'
#' @details
#' The analysis formula is inferred from the stored `gen_outcome()` formula.
#' Simulation terms `between(x)` and `within(x)` become observed-data columns
#' named `x_between` and `x_within`, computed from `x` by the simulation group
#' identifier. These columns are recomputed even when columns with the same
#' names already exist in `sim$data`.
#'
#' When `between()`/`within()` reference ILR coordinates of a compositional
#' [gen_mvn()] predictor, the helper instead builds a [complr()] object for
#' that predictor composition (using the simulator's parts, SBP, and total)
#' and maps the terms to the complr between/within coordinates (`bz<k>_<j>`
#' and `wz<k>_<j>`, where `k` is the coordinate and `j` the composition
#' index), exactly as in the standard `multilevelcoda` workflow. When the
#' outcome is also compositional, one multi-composition `complr` object
#' covers both the predictor and outcome compositions.
#'
#' The returned `brms` formula carries the family recorded by the outcome
#' generator (`gaussian()`, `poisson()`, `binomial()`,
#' [brms::negbinomial()], `Gamma(link = "log")`, or [brms::Beta()]). When the
#' generator has a scale model, it becomes the matching distributional
#' formula (`sigma`, `shape`, or `phi`). Binomial outcomes use
#' `y | trials(<outcome>_trials)` with the trials column generated by
#' [gen_outcome()]. When the same grouping factor appears in both the mean
#' and the scale formulas, the random effects are emitted with brms ID-linked
#' syntax (for example `(1 | p1 | ID)`) by default, because the simulator
#' draws these group-level effects from one joint covariance and the linked
#' syntax lets the analysis model estimate their correlation; see
#' `link_random` to opt out for large models fitted with approximate
#' algorithms.
#'
#' For dynamic formulas, `ar1()` is translated to lagged response columns whose
#' construction and name follow the two centering arguments:
#' \itemize{
#'   \item `lag_center = "within"`, `centering = "manifest"` (both defaults):
#'     the lagged observed response centered at each person's observed mean of
#'     all their values (not the mean of the lagged values),
#'     `lag_<response>_within`.
#'   \item `lag_center = "within"`, `centering = "latent"`: the lagged response
#'     centered at the latent conditional mean, which is the recorded latent
#'     residual state, `lag_<response>_latent`.
#'   \item `lag_center = "none"`, either `centering`: the raw lagged observed
#'     response, `lag_<response>`.
#' }
#' For compositional outcomes, the helper
#' rebuilds the ILR coordinates through [complr()] using the simulator's parts
#' and SBP metadata, then lags the generated `z` coordinates used by
#' [brmcoda()].
#'
#' The lags are built by [lag_by_time()] as *time-based* lags on the
#' simulation `time_id`: a row's lag comes from the observation at
#' `time - time_step`, not from the previous row. When a time point is
#' missing -- as can happen when analysing real data with skipped
#' observations -- the following row's lag is `NA` rather than the value from
#' before the gap. The autoregressive coefficient therefore refers to one
#' `time_step`. For data simulated with `ar1()` the time grid is complete and
#' equally spaced, so the result is identical to a positional shift. The
#' number of gap-affected rows is recorded in `metadata$lag_gaps`.
#'
#' Because a single `time_step` is applied to all groups, groups whose own
#' spacing is wider than the inferred step end up with only `NA` lags, and
#' `drop_lag_na = TRUE` would then remove them from the analysis data
#' entirely. `prep_sim_analysis()` warns in both situations and records the
#' diagnostics in the metadata: `metadata$time_step_by_group` (each group's
#' smallest positive spacing, `NA` for single-observation groups),
#' `metadata$time_step_heterogeneous` (whether those spacings differ),
#' `metadata$lag_na_groups` (groups without a single usable lagged row), and
#' `metadata$dropped_groups` (groups removed by `drop_lag_na`).
#'
#' @section Pragmatic default estimator:
#' [gen_outcome()] simulates latent residual AR/VAR dynamics around the
#' model-implied mean and resolves `between()`/`within()` from latent
#' generating components supplied by upstream predictor generators.
#' By default (`centering = "manifest"`) `prep_sim_analysis()` instead
#' constructs the model applied analysts commonly fit to observed data:
#' \itemize{
#'   \item `between(x)` and `within(x)` become `x_between` and `x_within`,
#'     recomputed from realised person means of the observed `x` (manifest
#'     centering), not from the latent generating components.
#'   \item `ar1()` becomes, by default, person-mean-centered lagged
#'     \emph{observed} response predictors (`lag_<response>_within`), not the
#'     latent residual state.
#'   \item For ILR coordinates of compositional predictors, `between(ilr)` and
#'     `within(ilr)` become the [complr()] coordinates `bz*`/`wz*`, where the
#'     between composition is the ILR of each person's closed
#'     \emph{arithmetic-mean} composition of the observed parts. The
#'     simulator's `between(ilr)`, in contrast, is the latent group-level ILR
#'     mean, whose back-transform is the closed \emph{geometric} (Aitchison)
#'     center. These are different estimands, and the gap has two parts: a
#'     structural Jensen gap between the arithmetic and geometric centers,
#'     which does \emph{not} shrink as the number of observations per group
#'     grows, plus ordinary sampling error, which does. Manifest centering is
#'     therefore biased for latent between-person effects even in long series,
#'     unlike the scalar case where only sampling error is involved.
#' }
#' These observed-data constructions target different estimands from the
#' simulation truth. Manifest person-mean centering and observed-score lagged
#' regression are known to yield biased estimates of between-person effects
#' and of inertia and cross-lag parameters relative to the latent generating
#' values, especially with short series (Ludtke et al. 2008; Hamaker &
#' Grasman 2014). This mismatch is intentional: it lets simulation studies
#' quantify the bias of the pragmatic estimator. Do not interpret systematic
#' discrepancies between estimates from this default analysis model and the
#' `gen_outcome()` truth parameters as errors in the simulator.
#'
#' @section Latent centering:
#' `centering = "latent"` replaces the manifest construction with the
#' generating components themselves, so for the affected `between()`/`within()`
#' terms the fitted model targets exactly the quantity the simulator recorded:
#' \itemize{
#'   \item Scalar predictors take `x_between` from the column the generator
#'     labelled with component `"between"` -- for `level = "level2"`
#'     predictors that is the variable itself -- and set
#'     `x_within <- x - x_between`.
#'   \item Compositional predictors keep the returned [complr()] object, but
#'     its between and within blocks (`bX`, `wX`, `bZ`, `wZ`, and the matching
#'     `b<part>`/`w<part>`/`bz*`/`wz*` columns) are rebuilt from the latent
#'     between composition emitted by the generator. Column names are
#'     unchanged, so [brmcoda()] and [substitution()] work as usual.
#' }
#' This is an oracle: the latent components exist only in simulated data, so
#' `"latent"` has no counterpart in a real analysis. Use it to separate
#' estimator bias from simulator error -- if a parameter is recovered under
#' `"latent"` but not under `"manifest"`, the gap is the estimator's, not the
#' simulator's. It requires a grouped design and a generator that emitted
#' latent components; `gen_custom()` generators that label roles without
#' emitting them can only be analysed with `"manifest"`.
#'
#' @section Latent centering and `ar1()`:
#' When the outcome formula contains `ar1()` and the lag is centered
#' (`lag_center = "within"`, the default), `centering = "latent"` centers the
#' lag with latent values, which replaces the lag predictor entirely: instead
#' of a lagged \emph{observed} response centered at the person's observed mean,
#' the analysis uses the latent residual state the simulator recorded, in
#' columns named `lag_<response>_latent`. The two constructions are the same
#' operation applied with different values, because the recorded residual is
#' `y_t - mu_t`, the response centered at its latent conditional mean. The
#' latent state needs no further centering: person-mean centering it on top
#' would reintroduce exactly the attenuation this construction avoids.
#'
#' With `lag_center = "none"` there is no lag centering for `centering` to act
#' on, so the raw lagged observed response is used under `"latent"` just as it
#' is under `"manifest"`. That combination is the point of keeping the two
#' arguments separate: a model can carry `between()`/`within()` predictors that
#' have nothing to do with the lag, and those can be built from the latent
#' generating components while the AR part keeps the non-centered
#' parametrization discussed below.
#'
#' A latent lag makes the AR, scale and `rescor` parameters estimate exactly
#' what the simulator recorded \strong{only when the analysis model has no
#' group-level effects}. The reason is not the proxy
#' problem that `lag_center` addresses. Writing the generating model as
#' `z_t = X_t beta + Phi e_{t-1} + eps_t`, the regressor `e_{t-1}` is
#' \emph{predetermined} but not \emph{strictly exogenous}: `Cov(eps_t, e_t)`
#' equals the innovation variance, and `e_t` is the next row's regressor. Any
#' estimator that mixes rows within a group -- which is what the random-effects
#' GLS transformation implied by a `(1 | group)` term does -- picks up that
#' cross-row covariance and incurs an O(1/T) dynamic-panel (Nickell) bias.
#' Substituting the true latent state removes the measurement problem but not
#' this one. So with any group-level effect, no parameter of a dynamic model is
#' guaranteed to keep its generating estimand, under either centering.
#'
#' For a matched-model recovery of a \emph{multilevel} dynamic model, the
#' appropriate analysis is not lag regression at all but a residual
#' autocorrelation structure, for example
#' `brms::bf(y ~ 1 + (1 | ID)) + brms::ar(time = day, gr = ID, p = 1, cov = TRUE)`.
#' That form is outside what this helper builds -- brms `ar()` supports only a
#' scalar unmoderated AR per response and parametrizes `sigma` as the marginal
#' rather than the innovation standard deviation -- so construct it by hand.
#'
#' Note that under `"manifest"` the latent columns of the same name are
#' overwritten in `analysis$data`; `sim$data` always keeps the originals.
#'
#' Whether the `ar1()` lag columns are centered at all is controlled by
#' `lag_center`, independently of `centering`.
#' Hamaker & Grasman (2014) show that in multilevel autoregressive models,
#' cluster-mean centering the lagged outcome (the `"within"` default)
#' attenuates the average autoregressive coefficient, whereas the
#' non-centered lagged outcome (`"none"`) recovers it nearly unbiased. The
#' trade-off is that under `"none"` the whole location mean structure is
#' reparametrized: the intercept estimates `(I - Phi_i) mu_i` rather than the
#' person mean, and any location covariate effects absorb omitted
#' `-Phi x_{t-1} beta` terms. Centering does not fully resolve this either:
#' the centered observed lag still carries the lagged mean structure, so
#' when a predictor varies within series its current coefficient can absorb
#' omitted lagged-covariate terms under `"within"` as well. Whenever the
#' formula contains `ar1()` and the lag is an observed response, then, no
#' parameter is guaranteed to keep its generating estimand, under either
#' `lag_center` setting. Use `"none"` when the average
#' autoregressive (inertia and cross-lag) coefficients are the estimands of
#' interest, and the default `"within"` when person means and their
#' predictors must remain interpretable; the data-simulation vignette
#' compares parameter recovery under both.
#'
#' @references
#' Ludtke, O., et al. (2008). The multilevel latent covariate model.
#' \emph{Psychological Methods}, 13(3), 203-229.
#'
#' Hamaker, E. L., & Grasman, R. P. P. P. (2014). To center or not to
#' center? Investigating inertia with a multilevel autoregressive model.
#' \emph{Frontiers in Psychology}, 5, 1492.
#'
#' @examples
#' params <- list(
#'   location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
#'   scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y")))
#' )
#' sim <- simulate_data(
#'   n = 5,
#'   seed = 1,
#'   generators = list(
#'     outcome = gen_outcome(
#'       y ~ 1,
#'       scale = sigma ~ 1,
#'       params = params
#'     )
#'   )
#' )
#' analysis <- prep_sim_analysis(sim)
#' analysis$formula
#'
#' @export
prep_sim_analysis <- function(sim, outcome = NULL, drop_lag_na = FALSE,
                              time_step = NULL,
                              lag_center = c("within", "none"),
                              centering = c("manifest", "latent"),
                              link_random = TRUE) {
  if (!inherits(sim, "mlsim_data")) {
    .mlsim_stop("`sim` must be an `mlsim_data` object returned by `simulate_data()`.")
  }
  drop_lag_na <- .mlsim_check_scalar_logical(drop_lag_na, "drop_lag_na")
  if (!is.null(time_step)) {
    if (inherits(time_step, "difftime")) {
      time_step <- as.numeric(time_step)
    }
    time_step <- .mlsim_check_scalar_number(time_step, "time_step")
    if (time_step <= 0) {
      .mlsim_stop("`time_step` must be a positive number.")
    }
  }
  lag_center <- match.arg(lag_center)
  centering <- match.arg(centering)
  link_random <- .mlsim_check_scalar_logical(link_random, "link_random")
  outcome <- .mlsim_analysis_resolve_outcome(sim, outcome)
  outcome_metadata <- sim$generator_metadata[[outcome]]
  parsed <- outcome_metadata$parsed %||% list()
  outcomes <- parsed$outcomes
  if (length(outcomes) == 0L) {
    .mlsim_stop("Outcome generator `%s` does not contain parsed outcome metadata.", outcome)
  }

  family <- parsed$family %||% "gaussian"
  dpar <- parsed$dpar %||% (if (identical(family, "gaussian")) "sigma" else NULL)
  trials_column <- parsed$trials_column

  data <- data.table::copy(sim$data)
  group_id <- sim$metadata$group_id
  grouped <- !is.null(sim$metadata$n_groups)
  time_id <- sim$metadata$time_id
  has_ar <- length(parsed$ar_terms %||% character()) > 0L
  # the two arguments have separate jobs: `lag_center` says whether the ar1()
  # lag is centered at all, `centering` says which values any centering uses.
  # Centering the lag with latent values means the recorded latent residual
  # state, y[t] - mu[t], the mirror of centering with the observed person mean
  lag_latent <- identical(centering, "latent") &&
    identical(lag_center, "within") &&
    isTRUE(has_ar)

  if (!is.null(trials_column) && !trials_column %in% names(data)) {
    .mlsim_stop(
      "Binomial analysis preparation requires the trials column `%s` generated by `gen_outcome()`.",
      trials_column
    )
  }

  if (isTRUE(grouped) && !group_id %in% names(data)) {
    .mlsim_stop("The active grouping column `%s` is missing from `sim$data`.", group_id)
  }
  latent_lag_sources <- NULL
  if (isTRUE(has_ar)) {
    if (is.null(time_id) || !time_id %in% names(data)) {
      .mlsim_stop("Preparing `ar1()` analysis terms requires the simulation `time_id` column.")
    }
    if (isTRUE(lag_latent)) {
      # attach before the reorder below: the recorded residual is in generation
      # row order, so setorderv() must carry each value with its own row
      latent_lag_sources <- .mlsim_analysis_attach_latent_residual(
        data = data,
        outcome_metadata = outcome_metadata,
        outcomes = outcomes,
        outcome = outcome
      )
    }
    # ordering keeps the output rows stable; lag correctness comes from the
    # time-based matching in lag_by_time(), not from row adjacency
    order_cols <- if (isTRUE(grouped)) c(group_id, time_id) else time_id
    data.table::setorderv(data, order_cols)
  }

  ilr_roles <- .mlsim_analysis_ilr_roles(sim, outcome_metadata)
  if (nrow(ilr_roles) > 0L && !isTRUE(grouped)) {
    .mlsim_stop("Preparing `between()` or `within()` analysis terms requires a grouped simulation design.")
  }
  derived_roles <- .mlsim_analysis_add_role_columns(
    data,
    outcome_metadata,
    grouped,
    group_id,
    exclude_variables = unique(ilr_roles$variable),
    centering = centering
  )
  compositional <- isTRUE(outcome_metadata$compositional)
  if (isTRUE(compositional) &&
      (is.null(outcome_metadata$parts) || is.null(outcome_metadata$sbp))) {
    .mlsim_stop("Compositional outcome `%s` is missing parts or SBP metadata.", outcome)
  }
  compositions <- .mlsim_analysis_compositions(sim, outcome_metadata, ilr_roles, compositional)
  response_map <- stats::setNames(outcomes, outcomes)
  special_term_map <- character()
  complr_object <- NULL

  if (length(compositions$parts) > 0L) {
    complr_object <- complr(
      data = data,
      parts = compositions$parts,
      sbp = compositions$sbp,
      total = compositions$total,
      idvar = if (isTRUE(grouped)) group_id else NULL
    )
    if (identical(centering, "latent")) {
      complr_object <- .mlsim_analysis_recenter_complr(
        complr_object = complr_object,
        compositions = compositions
      )
    }
    data <- data.table::copy(complr_object$dataout)
    if (!is.null(compositions$outcome_index)) {
      response_names <- colnames(complr_object$output[[compositions$outcome_index]]$Z)
      response_map <- stats::setNames(response_names, outcomes)
    }
    if (nrow(ilr_roles) > 0L) {
      special_term_map <- vapply(seq_len(nrow(ilr_roles)), function(i) {
        j <- match(ilr_roles$generator[[i]], compositions$generators)
        block <- if (identical(ilr_roles$component[[i]], "between")) "bZ" else "wZ"
        colnames(complr_object$output[[j]][[block]])[[ilr_roles$coordinate[[i]]]]
      }, character(1))
      names(special_term_map) <- sprintf("%s(%s)", ilr_roles$component, ilr_roles$variable)
    }
  }

  lag_info <- .mlsim_analysis_add_lag_columns(
    data = data,
    response_names = unname(response_map),
    has_ar = has_ar,
    grouped = grouped,
    group_id = group_id,
    time_id = time_id,
    time_step = time_step,
    drop_lag_na = drop_lag_na,
    lag_center = lag_center,
    # reordered to match `unname(response_map)` positionally: the residual is
    # keyed by simulator outcome name, the responses by analysis column name
    lag_source_names = if (isTRUE(lag_latent)) {
      unname(latent_lag_sources[names(response_map)])
    } else {
      NULL
    }
  )
  data <- lag_info$data

  if (!is.null(complr_object)) {
    complr_object$dataout <- data.table::copy(data)
    if (isTRUE(lag_info$dropped_rows > 0L)) {
      keep_rows <- lag_info$keep_rows
      complr_object$datain <- complr_object$datain[keep_rows]
      complr_object$output <- lapply(complr_object$output, function(out) {
        for (element in c("X", "bX", "wX", "Z", "bZ", "wZ", "dataout")) {
          if (!is.null(out[[element]])) {
            out[[element]] <- out[[element]][keep_rows, , drop = FALSE]
          }
        }
        out
      })
    }
  }

  formula <- .mlsim_analysis_build_brms_formula(
    formula = outcome_metadata$formula,
    scale = outcome_metadata$scale,
    response_map = response_map,
    lag_columns = lag_info$lag_columns,
    family = family,
    dpar = dpar,
    trials_column = trials_column,
    link_random = link_random,
    special_map = special_term_map
  )

  metadata <- list(
    outcome_generator = outcome,
    family = family,
    trials_column = trials_column,
    link_random = link_random,
    compositional = compositional,
    group_id = if (isTRUE(grouped)) group_id else NULL,
    time_id = time_id,
    drop_lag_na = drop_lag_na,
    original_formula = outcome_metadata$formula,
    original_scale = outcome_metadata$scale,
    response_map = response_map,
    derived_roles = derived_roles,
    composition_generators = compositions$generators,
    outcome_composition_index = compositions$outcome_index,
    special_term_map = special_term_map,
    lag_columns = lag_info$lag_columns,
    lag_center = lag_center,
    lag_latent = lag_latent,
    centering = centering,
    dropped_rows = lag_info$dropped_rows,
    time_step = lag_info$time_step,
    time_step_by_group = lag_info$time_step_by_group,
    time_step_heterogeneous = lag_info$time_step_heterogeneous,
    lag_na_groups = lag_info$lag_na_groups,
    dropped_groups = lag_info$dropped_groups,
    lag_gaps = lag_info$lag_gaps,
    lag_irregular = lag_info$lag_irregular,
    formula = formula
  )

  truth_info <- .mlsim_recovery_truth(outcome_metadata, metadata)
  metadata$truth_unavailable_reason <- truth_info$reason

  structure(
    list(
      data = data,
      formula = formula,
      complr = complr_object,
      truth = truth_info$truth,
      metadata = metadata
    ),
    class = "mlsim_analysis"
  )
}

#' @export
print.mlsim_analysis <- function(x, ...) {
  cat("<mlsim_analysis>\n")
  cat("  outcome: ", x$metadata$outcome_generator, "\n", sep = "")
  cat("  family: ", x$metadata$family %||% "gaussian", "\n", sep = "")
  cat("  rows: ", nrow(x$data), "\n", sep = "")
  cat("  compositional: ", if (isTRUE(x$metadata$compositional)) "yes" else "no", "\n", sep = "")
  if (length(x$metadata$derived_roles$column %||% character()) > 0L) {
    cat("  derived predictors: ", paste(unique(x$metadata$derived_roles$column), collapse = ", "), "\n", sep = "")
  }
  cat("  centering: ", x$metadata$centering %||% "manifest", "\n", sep = "")
  if (length(x$metadata$lag_columns) > 0L) {
    cat("  lag predictors: ", paste(x$metadata$lag_columns, collapse = ", "), "\n", sep = "")
    cat("  lag centering: ", x$metadata$lag_center %||% "within",
        if (isTRUE(x$metadata$lag_latent)) " (latent residual state)" else "",
        "\n", sep = "")
    if (isTRUE(x$metadata$time_step_heterogeneous)) {
      cat("  time step: heterogeneous across groups (see `$metadata$time_step_by_group`)\n")
    }
    if (length(x$metadata$dropped_groups %||% character()) > 0L) {
      cat("  groups removed by drop_lag_na: ",
          .mlsim_format_ids(x$metadata$dropped_groups), "\n", sep = "")
    } else if (length(x$metadata$lag_na_groups %||% character()) > 0L) {
      cat("  groups with all-NA lags: ",
          .mlsim_format_ids(x$metadata$lag_na_groups), "\n", sep = "")
    }
  }
  if (!is.null(x$truth)) {
    cat("  truth parameters: ", nrow(x$truth),
        " (see `$truth` and `sim_recovery()`)\n", sep = "")
  } else if (!is.null(x$metadata$truth_unavailable_reason)) {
    cat("  truth parameters: unavailable (", x$metadata$truth_unavailable_reason, ")\n", sep = "")
  }
  invisible(x)
}

.mlsim_analysis_resolve_outcome <- function(sim, outcome) {
  metadata <- sim$generator_metadata
  if (length(metadata) == 0L) {
    .mlsim_stop("Could not find an outcome generator in `sim$generator_metadata`.")
  }
  outcome_generators <- names(metadata)[vapply(metadata, function(x) {
    identical(x$distribution, "outcome")
  }, logical(1))]
  if (is.null(outcome)) {
    if (length(outcome_generators) != 1L) {
      .mlsim_stop(
        "Could not infer the outcome generator; found %d outcome generators. Supply `outcome`.",
        length(outcome_generators)
      )
    }
    return(outcome_generators[[1L]])
  }
  .mlsim_check_vars(outcome, 1L, "outcome")
  if (!outcome %in% names(metadata)) {
    .mlsim_stop("`outcome` must name a generator in `sim$generator_metadata`.")
  }
  if (!identical(metadata[[outcome]]$distribution, "outcome")) {
    .mlsim_stop("Generator `%s` is not a `gen_outcome()` outcome generator.", outcome)
  }
  outcome
}

.mlsim_analysis_ilr_roles <- function(sim, outcome_metadata) {
  empty <- data.table::data.table(
    variable = character(),
    component = character(),
    generator = character(),
    coordinate = integer()
  )
  roles <- outcome_metadata$selected_column_roles
  if (is.null(roles) || nrow(roles) == 0L || !"generator" %in% names(roles)) {
    return(empty)
  }
  roles <- unique(data.table::as.data.table(roles)[
    , c("variable", "component", "generator"), with = FALSE
  ])
  roles <- roles[roles$component %in% c("between", "within")]
  if (nrow(roles) == 0L) {
    return(empty)
  }
  coordinate <- vapply(seq_len(nrow(roles)), function(i) {
    gen_md <- sim$generator_metadata[[roles$generator[[i]]]]
    if (!isTRUE(gen_md$compositional)) {
      return(NA_integer_)
    }
    match(roles$variable[[i]], gen_md$ilr_vars %||% character())
  }, integer(1))
  roles$coordinate <- coordinate
  roles <- roles[!is.na(roles$coordinate)]
  if (nrow(roles) == 0L) {
    return(empty)
  }
  roles[]
}

.mlsim_analysis_compositions <- function(sim, outcome_metadata, ilr_roles, compositional_outcome) {
  generators <- intersect(sim$metadata$generator_order, unique(ilr_roles$generator))
  parts <- list()
  sbp <- list()
  total <- list()
  for (generator in generators) {
    gen_md <- sim$generator_metadata[[generator]]
    if (is.null(gen_md$parts) || is.null(gen_md$sbp)) {
      .mlsim_stop("Compositional predictor generator `%s` is missing parts or SBP metadata.", generator)
    }
    parts <- c(parts, list(gen_md$parts))
    sbp <- c(sbp, list(gen_md$sbp))
    total <- c(total, list(gen_md$total %||% 1))
  }
  outcome_index <- NULL
  if (isTRUE(compositional_outcome)) {
    parts <- c(parts, list(outcome_metadata$parts))
    sbp <- c(sbp, list(outcome_metadata$sbp))
    total <- c(total, list(outcome_metadata$total %||% 1))
    outcome_index <- length(parts)
  }
  list(
    parts = parts,
    sbp = sbp,
    total = total,
    generators = generators,
    outcome_index = outcome_index
  )
}

.mlsim_analysis_latent_residual_names <- function(outcomes) {
  paste0(".mlsim_", .mlsim_internal_name_piece(outcomes), "_ar_residual")
}

.mlsim_analysis_attach_latent_residual <- function(data, outcome_metadata, outcomes, outcome) {
  residual <- outcome_metadata$residual
  advice <- paste0(
    "Only Gaussian `gen_outcome()` models with `ar1()` record it; use ",
    "`centering = \"manifest\"` otherwise."
  )
  if (is.null(residual) || !is.matrix(residual)) {
    .mlsim_stop(
      paste0(
        "`centering = \"latent\"` requires the latent AR residual state recorded by ",
        "`gen_outcome()`, but `sim$generator_metadata[[\"%s\"]]$residual` is missing. %s"
      ),
      outcome, advice
    )
  }
  if (nrow(residual) != nrow(data)) {
    .mlsim_stop(
      paste0(
        "The latent AR residual state recorded for outcome generator `%s` has %d rows but ",
        "`sim$data` has %d. Refusing to build a misaligned latent lag."
      ),
      outcome, nrow(residual), nrow(data)
    )
  }
  missing_outcomes <- setdiff(outcomes, colnames(residual))
  if (length(missing_outcomes) > 0L) {
    .mlsim_stop(
      paste0(
        "The latent AR residual state recorded for outcome generator `%s` is missing ",
        "column(s) %s. %s"
      ),
      outcome,
      paste(sprintf("`%s`", missing_outcomes), collapse = ", "),
      advice
    )
  }
  source_names <- stats::setNames(
    .mlsim_analysis_latent_residual_names(outcomes),
    outcomes
  )
  existing <- intersect(source_names, names(data))
  if (length(existing) > 0L) {
    .mlsim_stop(
      "Cannot attach the latent AR residual state because column(s) %s already exist.",
      paste(sprintf("`%s`", existing), collapse = ", ")
    )
  }
  for (nm in outcomes) {
    data.table::set(data, j = source_names[[nm]], value = residual[, nm])
  }
  source_names
}

.mlsim_analysis_recenter_complr <- function(complr_object, compositions) {
  # predictor compositions only: a compositional outcome's own blocks are not
  # between/within predictors, so they keep the published arithmetic centering
  predictor_indices <- seq_along(compositions$generators)
  if (length(predictor_indices) == 0L) {
    return(complr_object)
  }
  if (!identical(complr_object$transform, "ilr")) {
    .mlsim_stop(
      "`centering = \"latent\"` is only supported for the `ilr` transform, not `%s`.",
      complr_object$transform
    )
  }
  # datain is the frame complr() was given, so its rows already match the
  # complr blocks even when ar1() reordered the analysis data
  datain <- complr_object$datain
  for (j in predictor_indices) {
    out <- complr_object$output[[j]]
    parts <- out$parts
    between_cols <- paste0(parts, "_between")
    missing_cols <- setdiff(between_cols, names(datain))
    if (length(missing_cols) > 0L) {
      .mlsim_stop(
        paste0(
          "`centering = \"latent\"` requires the latent between composition of generator ",
          "`%s`, but column(s) %s are missing. Only multilevel compositional generators ",
          "emit them; use `centering = \"manifest\"` otherwise."
        ),
        compositions$generators[[j]],
        paste(sprintf("`%s`", missing_cols), collapse = ", ")
      )
    }
    # mirrors complr(): close the between composition, take the perturbation
    # difference for within, then map both through the stored ilr basis
    bX <- acomp(datain[, between_cols, with = FALSE], total = out$total)
    colnames(bX) <- colnames(out$bX)
    wX <- out$X - bX
    colnames(wX) <- colnames(out$wX)
    bZ <- ilr(bX, V = out$psi)
    wZ <- ilr(wX, V = out$psi)
    colnames(bZ) <- colnames(out$bZ)
    colnames(wZ) <- colnames(out$wZ)

    out$bX <- bX
    out$wX <- wX
    out$bZ <- bZ
    out$wZ <- wZ
    out$dataout <- cbind(out$X, bX, wX, out$Z, bZ, wZ)
    complr_object$output[[j]] <- out

    recentered <- out$dataout
    for (nm in c(colnames(bX), colnames(wX), colnames(bZ), colnames(wZ))) {
      data.table::set(complr_object$dataout, j = nm, value = recentered[, nm])
    }
  }
  complr_object
}

.mlsim_analysis_add_role_columns <- function(data, outcome_metadata, grouped, group_id,
                                             exclude_variables = character(),
                                             centering = c("manifest", "latent")) {
  centering <- match.arg(centering)
  roles <- outcome_metadata$selected_column_roles
  if (is.null(roles) || nrow(roles) == 0L) {
    return(data.table::data.table(
      variable = character(),
      component = character(),
      column = character(),
      overwritten = logical()
    ))
  }
  role_cols <- intersect(c("variable", "component", "column", "generator"), names(roles))
  roles <- unique(roles[, role_cols, with = FALSE])
  roles <- roles[roles$component %in% c("between", "within")]
  roles <- roles[!roles$variable %in% exclude_variables]
  if (nrow(roles) == 0L) {
    return(data.table::data.table(
      variable = character(),
      component = character(),
      column = character(),
      overwritten = logical()
    ))
  }
  if (!isTRUE(grouped)) {
    .mlsim_stop("Preparing `between()` or `within()` analysis terms requires a grouped simulation design.")
  }
  if (!group_id %in% names(data)) {
    .mlsim_stop("The active grouping column `%s` is missing from `sim$data`.", group_id)
  }

  out <- data.table::data.table(
    variable = roles$variable,
    component = roles$component,
    column = paste(roles$variable, roles$component, sep = "_"),
    overwritten = paste(roles$variable, roles$component, sep = "_") %in% names(data)
  )
  for (variable in unique(roles$variable)) {
    if (!variable %in% names(data)) {
      .mlsim_stop("Cannot prepare analysis terms because observed predictor `%s` is missing.", variable)
    }
    if (!is.numeric(data[[variable]])) {
      .mlsim_stop("Observed predictor `%s` must be numeric to compute between/within analysis terms.", variable)
    }
    needed <- roles$component[roles$variable == variable]
    mean_col <- paste(variable, "between", sep = "_")
    within_col <- paste(variable, "within", sep = "_")
    mean_tmp <- paste0(".__mlsim_analysis_mean_", make.names(variable), "__")
    if (identical(centering, "latent")) {
      # the latent between column is the one the generator labelled with
      # component "between"; for `level = "level2"` predictors that is the
      # variable itself, so pasting "_between" would not find it
      source_col <- .mlsim_analysis_latent_between_column(roles, variable, data)
      data[, (mean_tmp) := get(source_col)]
    } else {
      data[, (mean_tmp) := mean(get(variable), na.rm = TRUE), by = group_id]
    }
    if ("between" %in% needed) {
      data[, (mean_col) := get(mean_tmp)]
    }
    if ("within" %in% needed) {
      data[, (within_col) := get(variable) - get(mean_tmp)]
    }
    data[, (mean_tmp) := NULL]
  }
  out
}

.mlsim_analysis_latent_between_column <- function(roles, variable, data) {
  is_between <- roles$variable == variable & roles$component == "between"
  column <- if ("column" %in% names(roles)) roles$column[is_between] else character()
  if (length(column) != 1L) {
    # only within() was referenced, so no between role row was recorded
    column <- paste(variable, "between", sep = "_")
  }
  if (!column %in% names(data)) {
    generator <- if ("generator" %in% names(roles)) {
      unique(roles$generator[roles$variable == variable])
    } else {
      character()
    }
    .mlsim_stop(
      paste0(
        "`centering = \"latent\"` requires the latent between component of `%s`, ",
        "but column `%s` is missing from `sim$data`%s. Generators that label ",
        "between/within roles without emitting the matching latent columns can only ",
        "be analysed with `centering = \"manifest\"`."
      ),
      variable,
      column,
      if (length(generator) == 1L) sprintf(" (generator `%s`)", generator) else ""
    )
  }
  column
}

.mlsim_analysis_group_time_steps <- function(data, group_id, time_id) {
  # per-group smallest positive spacing, mirroring the global inference in
  # .mlc_lag_infer_step() so the two are directly comparable
  data[, list(time_step = {
    d <- diff(sort(as.numeric(get(time_id))))
    d <- d[d > 0]
    if (length(d) == 0L) NA_real_ else min(d)
  }), by = group_id]
}

.mlsim_analysis_add_lag_columns <- function(data, response_names, has_ar, grouped,
                                            group_id, time_id, time_step, drop_lag_na,
                                            lag_center = "within",
                                            lag_source_names = NULL) {
  latent <- !is.null(lag_source_names)
  if (!isTRUE(has_ar)) {
    return(list(
      data = data,
      lag_columns = character(),
      keep_rows = rep(TRUE, nrow(data)),
      dropped_rows = 0L,
      time_step = NULL,
      time_step_by_group = NULL,
      time_step_heterogeneous = FALSE,
      lag_na_groups = character(),
      dropped_groups = character(),
      lag_gaps = 0L,
      lag_irregular = 0L
    ))
  }
  missing_responses <- response_names[!response_names %in% names(data)]
  if (length(missing_responses) > 0L) {
    .mlsim_stop(
      "Cannot prepare `ar1()` analysis terms because response columns are missing: %s.",
      paste(missing_responses, collapse = ", ")
    )
  }
  if (isTRUE(latent)) {
    if (length(lag_source_names) != length(response_names)) {
      .mlsim_stop(
        "Internal error: %d latent lag source column(s) for %d response(s).",
        length(lag_source_names), length(response_names)
      )
    }
    missing_sources <- lag_source_names[!lag_source_names %in% names(data)]
    if (length(missing_sources) > 0L) {
      .mlsim_stop(
        "Cannot prepare latent `ar1()` analysis terms because residual columns are missing: %s.",
        paste(missing_sources, collapse = ", ")
      )
    }
  }
  lag_columns <- if (isTRUE(latent)) {
    paste0("lag_", response_names, "_latent")
  } else if (identical(lag_center, "within")) {
    paste0("lag_", response_names, "_within")
  } else {
    paste0("lag_", response_names)
  }
  colliding <- lag_columns[lag_columns %in% names(data)]
  if (length(colliding) > 0L) {
    .mlsim_stop(
      "Preparing `ar1()` analysis terms would overwrite existing columns: %s. Rename or remove these columns from the simulated data first.",
      paste(sprintf("`%s`", colliding), collapse = ", ")
    )
  }
  group_steps <- NULL
  step_heterogeneous <- FALSE
  if (isTRUE(grouped)) {
    step_dt <- .mlsim_analysis_group_time_steps(data, group_id, time_id)
    group_steps <- stats::setNames(
      step_dt$time_step, as.character(step_dt[[group_id]])
    )
    known <- !is.na(group_steps)
    if (any(known)) {
      reference_step <- min(group_steps[known])
      tolerance <- 1e-8 * max(1, reference_step)
      wider <- known & (group_steps - reference_step) > tolerance
      step_heterogeneous <- any(wider)
      if (step_heterogeneous && is.null(time_step)) {
        .mlsim_warn(
          "Inferred `time_step` %s does not match every `%s` group: %d group(s) are spaced more widely (%s; steps %s), so all of their `ar1()` lag values will be `NA`. Supply `time_step` explicitly if a different step is intended.",
          format(reference_step), group_id, sum(wider),
          .mlsim_format_ids(names(group_steps)[wider]),
          .mlsim_format_ids(vapply(group_steps[wider], format, character(1L)))
        )
      }
    }
  }
  raw_suffix <- ".__mlsim_raw_lag__"
  # lag the latent residual state rather than the observed response under latent
  # centering; lag_by_time() applies the same time-based matching either way
  lag_inputs <- if (isTRUE(latent)) lag_source_names else response_names
  data <- lag_by_time(
    data,
    cols = lag_inputs,
    time = time_id,
    group = if (isTRUE(grouped)) group_id else NULL,
    time_step = time_step,
    suffix = raw_suffix
  )
  lag_meta <- attr(data, "lag_by_time")
  data.table::setattr(data, "lag_by_time", NULL)
  for (i in seq_along(response_names)) {
    response <- response_names[[i]]
    lag_col <- lag_columns[[i]]
    raw_col <- paste0(lag_inputs[[i]], raw_suffix)
    if (isTRUE(latent)) {
      # the latent residual is already zero-mean stationary; centering it would
      # reintroduce the very attenuation this construction avoids
      data[, (lag_col) := get(raw_col)]
    } else if (identical(lag_center, "none")) {
      data[, (lag_col) := get(raw_col)]
    } else if (isTRUE(grouped)) {
      data[, (lag_col) := get(raw_col) - mean(get(response), na.rm = TRUE), by = group_id]
    } else {
      data[, (lag_col) := get(raw_col) - mean(get(response), na.rm = TRUE)]
    }
    data[, (raw_col) := NULL]
  }
  if (isTRUE(latent)) {
    # drop the contemporaneous residual columns before the NA/keep_rows logic so
    # complete.cases() sees only lag columns, and so `y - e == mu` does not ship
    # a fully-determining oracle inside `analysis$data`
    data[, (lag_source_names) := NULL]
  }
  lag_na_groups <- character()
  if (isTRUE(grouped)) {
    usable <- data[, list(n_usable = sum(stats::complete.cases(.SD))),
                   by = group_id, .SDcols = lag_columns]
    lag_na_groups <- as.character(usable[[group_id]][usable$n_usable == 0L])
  }
  dropped_groups <- character()
  keep_rows <- rep(TRUE, nrow(data))
  if (isTRUE(drop_lag_na)) {
    keep_rows <- stats::complete.cases(data[, lag_columns, with = FALSE])
    # a group loses every row under drop_lag_na exactly when it has no
    # complete case over the lag columns, so the removed groups are the
    # all-NA-lag groups
    dropped_groups <- lag_na_groups
    if (length(dropped_groups) > 0L) {
      .mlsim_warn(
        "`drop_lag_na = TRUE` removed every row of %d `%s` group(s) because all of their `ar1()` lag values are `NA`: %s.",
        length(dropped_groups), group_id, .mlsim_format_ids(dropped_groups)
      )
    }
    data <- data[keep_rows]
  } else if (length(lag_na_groups) > 0L) {
    .mlsim_warn(
      "All `ar1()` lag values are `NA` for %d `%s` group(s): %s. These groups contribute no usable lagged rows; check their time spacing or supply `time_step`.",
      length(lag_na_groups), group_id, .mlsim_format_ids(lag_na_groups)
    )
  }
  list(
    data = data,
    lag_columns = lag_columns,
    keep_rows = keep_rows,
    dropped_rows = sum(!keep_rows),
    time_step = lag_meta$time_step,
    time_step_by_group = group_steps,
    time_step_heterogeneous = step_heterogeneous,
    lag_na_groups = lag_na_groups,
    dropped_groups = dropped_groups,
    lag_gaps = lag_meta$n_gaps,
    lag_irregular = lag_meta$n_irregular
  )
}

.mlsim_analysis_brms_family <- function(family) {
  switch(
    family,
    gaussian = stats::gaussian(),
    poisson = stats::poisson(),
    binomial = stats::binomial(),
    negbin = brms::negbinomial(),
    gamma = stats::Gamma(link = "log"),
    beta = brms::Beta(),
    .mlsim_stop("Unsupported outcome family `%s`.", family)
  )
}

.mlsim_analysis_build_brms_formula <- function(formula, scale, response_map, lag_columns,
                                               family = "gaussian", dpar = "sigma",
                                               trials_column = NULL,
                                               link_random = TRUE,
                                               special_map = character()) {
  rhs <- .mlsim_analysis_transform_expr(formula[[3L]], lag_columns, special_map)
  scale_rhs <- if (is.null(scale)) {
    NULL
  } else {
    .mlsim_analysis_transform_expr(scale[[3L]], character(), special_map)
  }

  # The simulator draws location and scale group-level effects from one joint
  # covariance, so grouping factors shared by both formulas use brms ID-linked
  # syntax by default so the analysis model can estimate their correlation.
  # `link_random = FALSE` opts out for large models where linked blocks are
  # intractable for approximate algorithms.
  shared_groups <- if (isTRUE(link_random)) {
    intersect(
      .mlsim_analysis_random_groups(rhs),
      if (is.null(scale_rhs)) character() else .mlsim_analysis_random_groups(scale_rhs)
    )
  } else {
    character()
  }
  if (length(shared_groups) > 0L) {
    labels <- stats::setNames(paste0("p", seq_along(shared_groups)), shared_groups)
    rhs <- .mlsim_analysis_link_random_terms(rhs, labels)
    scale_rhs <- .mlsim_analysis_link_random_terms(scale_rhs, labels)
  }

  brms_family <- .mlsim_analysis_brms_family(family)
  response_names <- unname(response_map)

  formulas <- lapply(response_names, function(response) {
    lhs <- if (is.null(trials_column)) {
      as.name(response)
    } else {
      call("|", as.name(response), call("trials", as.name(trials_column)))
    }
    if (is.null(scale_rhs)) {
      brms::bf(
        .mlsim_analysis_formula(lhs, rhs),
        family = brms_family
      )
    } else {
      brms::bf(
        .mlsim_analysis_formula(lhs, rhs),
        .mlsim_analysis_formula(as.name(dpar), scale_rhs),
        family = brms_family
      )
    }
  })
  if (length(formulas) == 1L) {
    return(formulas[[1L]])
  }
  Reduce(`+`, formulas) + brms::set_rescor(TRUE)
}

.mlsim_analysis_random_groups <- function(expr) {
  if (!is.call(expr)) {
    return(character())
  }
  if (identical(as.character(expr[[1L]]), "|") && is.symbol(expr[[3L]])) {
    return(unique(c(
      as.character(expr[[3L]]),
      .mlsim_analysis_random_groups(expr[[2L]])
    )))
  }
  groups <- unique(unlist(lapply(as.list(expr[-1L]), .mlsim_analysis_random_groups)))
  groups %||% character()
}

.mlsim_analysis_link_random_terms <- function(expr, labels) {
  if (!is.call(expr)) {
    return(expr)
  }
  if (identical(as.character(expr[[1L]]), "|") && is.symbol(expr[[3L]])) {
    group <- as.character(expr[[3L]])
    lhs <- .mlsim_analysis_link_random_terms(expr[[2L]], labels)
    if (group %in% names(labels)) {
      return(call("|", call("|", lhs, as.name(labels[[group]])), expr[[3L]]))
    }
    return(call("|", lhs, expr[[3L]]))
  }
  args <- lapply(as.list(expr[-1L]), .mlsim_analysis_link_random_terms, labels = labels)
  as.call(c(list(expr[[1L]]), args))
}

.mlsim_analysis_formula <- function(lhs, rhs) {
  out <- eval(call("~", lhs, rhs), envir = parent.frame())
  environment(out) <- parent.frame()
  out
}

.mlsim_analysis_transform_expr <- function(expr, lag_columns, special_map = character()) {
  if (is.symbol(expr) || is.atomic(expr)) {
    return(expr)
  }
  if (!is.call(expr)) {
    return(expr)
  }
  op <- as.character(expr[[1L]])
  if (identical(op, "(")) {
    return(call("(", .mlsim_analysis_transform_expr(expr[[2L]], lag_columns, special_map)))
  }
  if (identical(op, "between") || identical(op, "within")) {
    if (length(expr) != 2L || !is.symbol(expr[[2L]])) {
      .mlsim_stop("`%s()` terms must contain one variable name.", op)
    }
    key <- sprintf("%s(%s)", op, as.character(expr[[2L]]))
    if (key %in% names(special_map)) {
      return(as.name(special_map[[key]]))
    }
    return(as.name(paste(as.character(expr[[2L]]), op, sep = "_")))
  }
  if (identical(op, "ar1")) {
    if (length(expr) != 1L) {
      .mlsim_stop("`ar1()` does not take arguments.")
    }
    if (length(lag_columns) == 0L) {
      .mlsim_stop("Cannot translate `ar1()` because no lag columns were prepared.")
    }
    return(.mlsim_analysis_sum_expr(lapply(lag_columns, as.name)))
  }
  if (identical(op, "+")) {
    return(.mlsim_analysis_sum_expr(list(
      .mlsim_analysis_transform_expr(expr[[2L]], lag_columns, special_map),
      .mlsim_analysis_transform_expr(expr[[3L]], lag_columns, special_map)
    )))
  }
  if (identical(op, ":")) {
    left <- .mlsim_analysis_split_sum(.mlsim_analysis_transform_expr(expr[[2L]], lag_columns, special_map))
    right <- .mlsim_analysis_split_sum(.mlsim_analysis_transform_expr(expr[[3L]], lag_columns, special_map))
    expanded <- unlist(lapply(left, function(l) {
      lapply(right, function(r) call(":", l, r))
    }), recursive = FALSE)
    return(.mlsim_analysis_sum_expr(expanded))
  }
  if (identical(op, "|")) {
    return(call("|", .mlsim_analysis_transform_expr(expr[[2L]], lag_columns, special_map), expr[[3L]]))
  }
  args <- lapply(as.list(expr[-1L]), .mlsim_analysis_transform_expr,
                 lag_columns = lag_columns, special_map = special_map)
  as.call(c(list(expr[[1L]]), args))
}

.mlsim_analysis_split_sum <- function(expr) {
  if (is.call(expr) && identical(as.character(expr[[1L]]), "+")) {
    return(c(.mlsim_analysis_split_sum(expr[[2L]]), .mlsim_analysis_split_sum(expr[[3L]])))
  }
  list(expr)
}

.mlsim_analysis_sum_expr <- function(exprs) {
  exprs <- unlist(lapply(exprs, .mlsim_analysis_split_sum), recursive = FALSE)
  if (length(exprs) == 0L) {
    return(1)
  }
  Reduce(function(x, y) call("+", x, y), exprs)
}
