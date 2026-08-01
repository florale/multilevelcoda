#' Align Fitted Estimates with Simulation Truth
#'
#' `sim_recovery()` joins the posterior summaries of a model fitted to
#' simulated data with the true generating parameter values recorded by
#' [simulate_data()], using the brms parameter names on both sides.
#' It replaces the error-prone manual name matching otherwise needed to
#' compare [gen_outcome()] truth with [brms::fixef()] and
#' [brms::VarCorr()] output.
#'
#' The truth table is constructed by [prep_sim_analysis()] and stored in
#' its `$truth` element; `sim_recovery()` matches every truth row to
#' exactly one fitted parameter and vice versa. Any parameter that cannot
#' be matched in either direction is an error, never a silently wrong
#' row: truth values are attached purely by name, and duplicate or
#' ambiguous names (for example response names that collide after brms
#' removes `_` and `.` characters) abort with an explanatory message.
#'
#' All [gen_outcome()] families are supported, not only Gaussian models.
#' Rows for a modeled distributional parameter are labeled with the
#' family's brms parameter name (`sigma_*` for `"gaussian"`, `shape_*`
#' for `"negbin"` and `"gamma"`, `phi_*` for `"beta"`), and their truth
#' values are on the log link scale that both [gen_outcome()] and brms
#' use for these parameters; location rows are on the family's location
#' link (identity, log, or logit), again matching on both sides.
#' `"poisson"` and `"binomial"` families have no distributional
#' parameter, so their tables contain no scale rows, and a binomial
#' `trials()` term contributes no parameter. Because `ar1()`,
#' multivariate `mvbind()` responses, and compositional outcomes are
#' Gaussian-only simulator features, recovery tables for the other
#' families never contain `rescor` rows.
#'
#' Two situations yield no recovery row at all: cross-block random-effect
#' correlations when the analysis was prepared with
#' [prep_sim_analysis()]`(link_random = FALSE)` (the fitted model never
#' samples them), and analysis models whose structure cannot be aligned
#' with the simulator's single joint covariance (recorded in
#' `metadata$truth_unavailable_reason`).
#'
#' Interpreting the `bias` and `covered` columns requires knowing what the
#' fitted model estimates. [prep_sim_analysis()] deliberately builds the
#' pragmatic observed-data model that applied analysts commonly fit, not
#' the matched model for the generating process, so some parameters
#' estimate a related quantity rather than the one the simulator recorded
#' -- most notably whenever the formula contains `ar1()`, and for
#' `between()`/`within()` terms under the default
#' `centering = "manifest"`. Systematic bias in those parameters is an
#' expected property of the estimator, not an error in the simulator. See
#' the "Pragmatic default estimator" section of [prep_sim_analysis()]
#' before drawing conclusions from a recovery table.
#'
#' @param fit A [brmcoda()] model object or [brms::brmsfit] fitted to
#'   `analysis$data` (or `analysis$complr`) using `analysis$formula`.
#' @param analysis An `mlsim_analysis` object from [prep_sim_analysis()]
#'   with a non-`NULL` `$truth` element.
#' @param probs A length-two numeric vector of lower and upper posterior
#'   interval probabilities used for the coverage indicator.
#'   Defaults to `c(0.025, 0.975)`.
#'
#' @return A [data.table::data.table] of class `mlsim_recovery` with one
#'   row per model parameter:
#' \itemize{
#'   \item `type`: `"fixed"`, `"random_sd"`, `"random_cor"`, or
#'     `"rescor"`. Random-effect rows are reported on the standard
#'     deviation / correlation scale used by [brms::VarCorr()].
#'   \item `group`: grouping factor for random-effect rows, `NA`
#'     otherwise.
#'   \item `parameter`: the brms-style parameter label.
#'   \item `truth`: the true generating value from the simulator
#'     (`sqrt()` of covariance diagonals for `random_sd` rows,
#'     correlations for `random_cor` rows).
#'   \item `estimate`, `est_error`, `lower`, `upper`: posterior summary
#'     of the fitted parameter at `probs`.
#'   \item `bias`: `estimate - truth`.
#'   \item `covered`: whether `truth` lies inside `[lower, upper]`.
#'   \item `simulator_name`: the simulator-side parameter identifier the
#'     row was derived from, for provenance.
#' }
#'
#' @seealso [prep_sim_analysis()] for the analysis object and its
#'   `$truth` table, [gen_template()] for authoring the truth values.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   params <- list(
#'     location = list(beta = matrix(
#'       c(0.5, 0.2),
#'       nrow = 2,
#'       dimnames = list(c("(Intercept)", "x"), "y")
#'     )),
#'     scale = list(beta = matrix(
#'       log(0.5),
#'       nrow = 1,
#'       dimnames = list("(Intercept)", "y")
#'     ))
#'   )
#'   sim <- simulate_data(
#'     n = 200,
#'     seed = 1,
#'     generators = list(
#'       x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
#'       y = gen_outcome(y ~ x, scale = sigma ~ 1, params = params)
#'     )
#'   )
#'   analysis <- prep_sim_analysis(sim)
#'   analysis$truth
#'
#'   fit <- brms::brm(
#'     analysis$formula,
#'     data = analysis$data,
#'     backend = "cmdstanr",
#'     refresh = 0
#'   )
#'   sim_recovery(fit, analysis)
#' }
#' }
#'
#' @export
sim_recovery <- function(fit, analysis, probs = c(0.025, 0.975)) {
  if (!inherits(analysis, "mlsim_analysis")) {
    .mlsim_stop("`analysis` must be an `mlsim_analysis` object returned by `prep_sim_analysis()`.")
  }
  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs) ||
      probs[[1L]] <= 0 || probs[[2L]] >= 1 || probs[[1L]] >= probs[[2L]]) {
    .mlsim_stop("`probs` must be two increasing probabilities strictly between 0 and 1.")
  }
  truth <- analysis$truth
  if (is.null(truth)) {
    reason <- analysis$metadata$truth_unavailable_reason
    if (!is.null(reason)) {
      .mlsim_stop(
        "Simulation truth is not available for this analysis model: %s. Construct the recovery table by hand for this design.",
        reason
      )
    }
    .mlsim_stop(
      "`analysis` has no `$truth` table; re-run `prep_sim_analysis()` with this version of multilevelcoda."
    )
  }
  if (inherits(fit, "brmcoda")) {
    fit <- fit$model
  }
  if (!inherits(fit, "brmsfit")) {
    .mlsim_stop("`fit` must be a `brmcoda` or `brmsfit` model object.")
  }

  fixef_summary <- brms::fixef(fit, probs = probs)
  needs_random <- any(truth$type %in% c("random_sd", "random_cor"))
  needs_rescor <- any(truth$type == "rescor")
  variables <- if (needs_random || needs_rescor) brms::variables(fit) else character()
  varcorr <- if (needs_random) brms::VarCorr(fit, probs = probs) else NULL
  rescor <- if (needs_rescor) {
    # the table has a column named `truth`, so subset with a precomputed
    # index to keep data.table from resolving `truth` to that column
    rescor_index <- which(truth$type == "rescor")
    .mlsim_recovery_rescor_summary(fit, truth[rescor_index, ], variables, probs)
  } else {
    NULL
  }

  estimates <- .mlsim_recovery_estimates(fixef_summary, varcorr, rescor, variables, probs)
  out <- .mlsim_recovery_join(truth, estimates)
  data.table::setattr(out, "probs", probs)
  data.table::setattr(out, "class", c("mlsim_recovery", class(data.table::data.table())))
  out
}

#' @export
print.mlsim_recovery <- function(x, digits = 3, ...) {
  probs <- attr(x, "probs") %||% c(0.025, 0.975)
  tbl <- data.table::copy(x)
  data.table::setattr(tbl, "class", class(data.table::data.table()))
  print(tbl, digits = digits, ...)
  cat(sprintf(
    "Coverage of the %.0f%% interval: %d/%d.\n",
    100 * (probs[[2L]] - probs[[1L]]),
    sum(x$covered, na.rm = TRUE), nrow(x)
  ))
  invisible(x)
}

# truth construction -----------------------------------------------------

# Returns list(truth = data.table or NULL, reason = character or NULL).
# `stop()`s only on ambiguity that would mislabel truth (name collisions);
# structurally unsupported analysis models yield truth = NULL with a reason
# so prep_sim_analysis() keeps working for them.
.mlsim_recovery_truth <- function(outcome_metadata, analysis_metadata) {
  reason <- .mlsim_recovery_unsupported_reason(analysis_metadata)
  if (!is.null(reason)) {
    return(list(truth = NULL, reason = reason))
  }
  params <- outcome_metadata$params
  if (is.null(params)) {
    return(list(truth = NULL, reason = "the outcome generator metadata does not record `params`"))
  }

  parsed <- outcome_metadata$parsed %||% list()
  family <- parsed$family %||% "gaussian"
  family_info <- .mlsim_outcome_families()[[family]]
  dpar <- parsed$dpar %||% family_info$dpar
  scale_label <- family_info$scale_label
  response_map <- analysis_metadata$response_map
  special_term_map <- analysis_metadata$special_term_map %||% character()
  prefixes <- .mlsim_recovery_prefixes(response_map)
  lag_columns <- analysis_metadata$lag_columns %||% character()
  if (length(lag_columns) == length(response_map)) {
    lag_columns <- stats::setNames(lag_columns, names(response_map))
  }
  fixed <- .mlsim_recovery_fixed_truth(
    params, prefixes, dpar, lag_columns, special_term_map
  )
  random <- .mlsim_recovery_random_truth(
    params, prefixes, dpar, scale_label, lag_columns, special_term_map,
    isTRUE(analysis_metadata$link_random)
  )
  rescor <- .mlsim_recovery_rescor_truth(params, prefixes, family)

  truth <- data.table::rbindlist(list(fixed, random, rescor), use.names = TRUE)
  keys <- .mlsim_recovery_keys(truth)
  if (anyDuplicated(keys)) {
    duplicated_names <- unique(truth$parameter[keys %in% keys[duplicated(keys)]])
    .mlsim_stop(
      "Simulation truth maps more than one true value to the same model parameter name: %s. Refusing to build an ambiguous truth table; rename the affected simulator outcomes or terms.",
      paste(sprintf("`%s`", duplicated_names), collapse = ", ")
    )
  }
  list(truth = truth, reason = NULL)
}

.mlsim_recovery_fixed_truth <- function(params, prefixes, dpar, lag_columns,
                                        special_term_map) {
  rows <- list()
  add <- function(parameter, truth, simulator_name) {
    rows[[length(rows) + 1L]] <<- data.table::data.table(
      type = "fixed",
      group = NA_character_,
      parameter = parameter,
      coef1 = NA_character_,
      coef2 = NA_character_,
      truth = truth,
      simulator_name = simulator_name
    )
  }

  location <- params$location$beta
  for (outcome in colnames(location)) {
    for (term in rownames(location)) {
      add(
        .mlsim_recovery_paste_name(
          prefixes[[outcome]],
          .mlsim_recovery_term(term, special_term_map)
        ),
        location[term, outcome],
        sprintf("location:%s@%s", term, outcome)
      )
    }
  }
  scale_beta <- params$scale$beta
  if (!is.null(scale_beta)) {
    for (outcome in colnames(scale_beta)) {
      for (term in rownames(scale_beta)) {
        add(
          .mlsim_recovery_paste_name(
            dpar,
            prefixes[[outcome]],
            .mlsim_recovery_term(term, special_term_map)
          ),
          scale_beta[term, outcome],
          sprintf("scale:%s@%s", term, outcome)
        )
      }
    }
  }
  ar_beta <- params$ar$beta
  if (!is.null(ar_beta)) {
    for (ar_term in dimnames(ar_beta)[[1L]]) {
      for (to in dimnames(ar_beta)[[2L]]) {
        for (from in dimnames(ar_beta)[[3L]]) {
          add(
            .mlsim_recovery_paste_name(
              .mlsim_recovery_prefix(prefixes, to),
              .mlsim_recovery_ar_term(ar_term, from, lag_columns, special_term_map)
            ),
            ar_beta[ar_term, to, from],
            sprintf("ar:%s@to=%s,from=%s", ar_term, to, from)
          )
        }
      }
    }
  }
  data.table::rbindlist(rows)
}

.mlsim_recovery_random_truth <- function(params, prefixes, dpar, scale_label,
                                         lag_columns, special_term_map,
                                         link_random) {
  if (is.null(params$random) || length(params$random) == 0L) {
    return(NULL)
  }
  group <- names(params$random)[[1L]]
  covariance <- params$random[[group]]$covariance
  simulator_names <- rownames(covariance)
  parsed <- lapply(simulator_names, .mlsim_recovery_parse_random_name,
                   scale_label = scale_label)
  coef_names <- vapply(parsed, function(p) {
    if (identical(p$component, "ar")) {
      .mlsim_recovery_paste_name(
        .mlsim_recovery_prefix(prefixes, p$to),
        .mlsim_recovery_ar_term(p$term, p$from, lag_columns, special_term_map)
      )
    } else if (identical(p$component, "location")) {
      .mlsim_recovery_paste_name(
        .mlsim_recovery_prefix(prefixes, p$outcome),
        .mlsim_recovery_term(p$term, special_term_map)
      )
    } else {
      .mlsim_recovery_paste_name(
        dpar,
        .mlsim_recovery_prefix(prefixes, p$outcome),
        .mlsim_recovery_term(p$term, special_term_map)
      )
    }
  }, character(1))
  # brms estimates correlations within one covariance block; without ID-linked
  # random terms the location/AR (mu) and scale (dpar) formulas form separate
  # blocks, so cross-block truth correlations have no fitted counterpart.
  block <- vapply(parsed, function(p) {
    if (p$component %in% c("location", "ar")) "mu" else "dpar"
  }, character(1))

  sd_rows <- data.table::data.table(
    type = "random_sd",
    group = group,
    parameter = coef_names,
    coef1 = NA_character_,
    coef2 = NA_character_,
    truth = sqrt(diag(covariance)),
    simulator_name = simulator_names
  )

  sds <- sqrt(diag(covariance))
  denominator <- outer(sds, sds)
  correlation <- ifelse(denominator > 0, covariance / denominator, NA_real_)
  cor_rows <- list()
  n <- length(coef_names)
  if (n > 1L) {
    for (j in seq_len(n - 1L)) {
      for (i in seq((j + 1L), n)) {
        if (!link_random && !identical(block[[i]], block[[j]])) {
          next
        }
        cor_rows[[length(cor_rows) + 1L]] <- data.table::data.table(
          type = "random_cor",
          group = group,
          parameter = sprintf("cor(%s,%s)", coef_names[[j]], coef_names[[i]]),
          coef1 = coef_names[[j]],
          coef2 = coef_names[[i]],
          truth = correlation[i, j],
          simulator_name = sprintf("cor:%s,%s", simulator_names[[j]], simulator_names[[i]])
        )
      }
    }
  }
  data.table::rbindlist(c(list(sd_rows), cor_rows))
}

.mlsim_recovery_rescor_truth <- function(params, prefixes, family) {
  correlation <- params$scale$correlation
  if (!identical(family, "gaussian") || is.null(correlation) ||
      length(prefixes) < 2L) {
    return(NULL)
  }
  outcomes <- rownames(correlation)
  rows <- list()
  for (j in seq_len(length(outcomes) - 1L)) {
    for (i in seq((j + 1L), length(outcomes))) {
      coef1 <- .mlsim_recovery_prefix(prefixes, outcomes[[j]])
      coef2 <- .mlsim_recovery_prefix(prefixes, outcomes[[i]])
      rows[[length(rows) + 1L]] <- data.table::data.table(
        type = "rescor",
        group = NA_character_,
        parameter = sprintf("rescor(%s,%s)", coef1, coef2),
        coef1 = coef1,
        coef2 = coef2,
        truth = correlation[outcomes[[i]], outcomes[[j]]],
        simulator_name = sprintf("rescor:%s,%s", outcomes[[j]], outcomes[[i]])
      )
    }
  }
  data.table::rbindlist(rows)
}

.mlsim_recovery_unsupported_reason <- function(analysis_metadata) {
  formula <- analysis_metadata$original_formula
  if (!is.null(formula) &&
      .mlsim_recovery_count_random_terms(formula[[length(formula)]]) > 1L) {
    return("the location formula has more than one random-effect term, which brms fits as separate covariance blocks that cannot be aligned with the simulator's single joint covariance")
  }
  scale <- analysis_metadata$original_scale
  if (!is.null(scale) &&
      .mlsim_recovery_count_random_terms(scale[[length(scale)]]) > 1L) {
    return("the scale formula has more than one random-effect term, which brms fits as separate covariance blocks that cannot be aligned with the simulator's single joint covariance")
  }
  NULL
}

.mlsim_recovery_count_random_terms <- function(expr) {
  if (!is.call(expr)) {
    return(0L)
  }
  if (identical(as.character(expr[[1L]]), "|")) {
    return(1L)
  }
  sum(vapply(as.list(expr[-1L]), .mlsim_recovery_count_random_terms, integer(1)))
}

# name construction ------------------------------------------------------

# Multivariate brms models label parameters with sanitized response names
# (brms removes `_` and `.` after make.names()); univariate models use no
# response prefix. Sanitization can collide, which would silently mislabel
# truth, so collisions are an error rather than a rename.
.mlsim_recovery_prefixes <- function(response_map) {
  if (length(response_map) <= 1L) {
    return(stats::setNames(rep("", length(response_map)), names(response_map)))
  }
  sanitized <- gsub("\\.|_", "", make.names(unname(response_map)))
  if (anyDuplicated(sanitized)) {
    clash <- unique(sanitized[duplicated(sanitized)])
    detail <- vapply(clash, function(name) {
      sprintf(
        "%s -> `%s`",
        paste(sprintf("`%s`", unname(response_map)[sanitized == name]), collapse = ", "),
        name
      )
    }, character(1))
    .mlsim_stop(
      "Response names are ambiguous after brms removes `_` and `.` characters: %s. Rename the simulator outcomes or composition so the sanitized response names are unique.",
      paste(detail, collapse = "; ")
    )
  }
  stats::setNames(sanitized, names(response_map))
}

.mlsim_recovery_prefix <- function(prefixes, outcome) {
  if (!outcome %in% names(prefixes)) {
    .mlsim_stop(
      "Simulation truth references outcome `%s`, which is not one of the analysis responses (%s).",
      outcome,
      paste(sprintf("`%s`", names(prefixes)), collapse = ", ")
    )
  }
  prefixes[[outcome]]
}

.mlsim_recovery_paste_name <- function(...) {
  pieces <- c(...)
  paste(pieces[nzchar(pieces)], collapse = "_")
}

# Parses the simulator's structured random-effect names, e.g.
# "location|outcome=ilr1|term=(Intercept)" or
# "ar|term=ar1()|to=ilr1|from=ilr2". Values never contain `|`, and keys
# never contain `=`; splitting on the first `=` keeps future term text safe.
.mlsim_recovery_parse_random_name <- function(name, scale_label) {
  pieces <- strsplit(name, "|", fixed = TRUE)[[1L]]
  component <- pieces[[1L]]
  fields <- pieces[-1L]
  separator <- regexpr("=", fields, fixed = TRUE)
  expected_keys <- if (identical(component, "ar")) {
    c("term", "to", "from")
  } else {
    c("outcome", "term")
  }
  valid <- component %in% c("location", "ar", scale_label %||% character()) &&
    length(fields) == length(expected_keys) &&
    all(separator > 0L)
  if (valid) {
    keys <- substr(fields, 1L, separator - 1L)
    valid <- identical(keys, expected_keys)
  }
  if (!valid) {
    .mlsim_stop("Cannot parse the random-effect parameter name `%s`.", name)
  }
  values <- substring(fields, separator + 1L)
  c(list(component = component), stats::setNames(as.list(values), keys))
}

# Maps a simulator model-matrix term to the term name in the analysis
# model, mirroring .mlsim_analysis_transform_expr() at the name level.
.mlsim_recovery_term <- function(term, special_term_map) {
  if (identical(term, "(Intercept)")) {
    return("Intercept")
  }
  pieces <- strsplit(term, ":", fixed = TRUE)[[1L]]
  translated <- vapply(pieces, function(piece) {
    if (piece %in% names(special_term_map)) {
      return(unname(special_term_map[[piece]]))
    }
    for (role in c("between", "within")) {
      opening <- paste0(role, "(")
      if (startsWith(piece, opening) && endsWith(piece, ")")) {
        variable <- substr(piece, nchar(opening) + 1L, nchar(piece) - 1L)
        return(paste(variable, role, sep = "_"))
      }
    }
    piece
  }, character(1), USE.NAMES = FALSE)
  paste(translated, collapse = ":")
}

# An AR term like "treatmenttreatment:ar1()" becomes the analysis term with
# the ar1() factor replaced in place by the `from` response's lag column,
# e.g. "treatmenttreatment:lag_z1_1_within".
.mlsim_recovery_ar_term <- function(ar_term, from, lag_columns, special_term_map) {
  pieces <- strsplit(ar_term, ":", fixed = TRUE)[[1L]]
  is_ar <- pieces == "ar1()"
  if (sum(is_ar) != 1L) {
    .mlsim_stop(
      "Cannot align the AR parameter `%s`: expected exactly one `ar1()` factor.",
      ar_term
    )
  }
  lag_column <- lag_columns[[from]] %||% NA_character_
  if (is.na(lag_column)) {
    .mlsim_stop(
      "Cannot align the AR parameter `%s`: no lag column was prepared for outcome `%s`.",
      ar_term, from
    )
  }
  pieces <- vapply(seq_along(pieces), function(i) {
    if (is_ar[[i]]) lag_column else .mlsim_recovery_term(pieces[[i]], special_term_map)
  }, character(1))
  paste(pieces, collapse = ":")
}

# estimate extraction and joining ----------------------------------------

.mlsim_recovery_rescor_summary <- function(fit, rescor_truth, variables, probs) {
  candidates <- lapply(seq_len(nrow(rescor_truth)), function(i) {
    both <- c(
      sprintf("rescor__%s__%s", rescor_truth$coef1[[i]], rescor_truth$coef2[[i]]),
      sprintf("rescor__%s__%s", rescor_truth$coef2[[i]], rescor_truth$coef1[[i]])
    )
    found <- both[both %in% variables]
    if (length(found) == 0L) {
      return(NULL)
    }
    list(variable = found[[1L]], coef1 = rescor_truth$coef1[[i]],
         coef2 = rescor_truth$coef2[[i]])
  })
  candidates <- candidates[!vapply(candidates, is.null, logical(1))]
  if (length(candidates) == 0L) {
    return(NULL)
  }
  summary <- brms::posterior_summary(
    fit,
    variable = vapply(candidates, `[[`, character(1), "variable"),
    probs = probs
  )
  list(
    summary = summary,
    coef1 = vapply(candidates, `[[`, character(1), "coef1"),
    coef2 = vapply(candidates, `[[`, character(1), "coef2")
  )
}

.mlsim_recovery_summary_row <- function(summary_matrix, index, probs) {
  quantile_columns <- paste0("Q", probs * 100)
  missing_columns <- setdiff(c("Estimate", "Est.Error", quantile_columns),
                             colnames(summary_matrix))
  if (length(missing_columns) > 0L) {
    .mlsim_stop(
      "Posterior summary is missing the column(s) %s.",
      paste(sprintf("`%s`", missing_columns), collapse = ", ")
    )
  }
  list(
    estimate = summary_matrix[index, "Estimate"],
    est_error = summary_matrix[index, "Est.Error"],
    lower = summary_matrix[index, quantile_columns[[1L]]],
    upper = summary_matrix[index, quantile_columns[[2L]]]
  )
}

.mlsim_recovery_estimates <- function(fixef_summary, varcorr, rescor,
                                      variables, probs) {
  rows <- list()
  add <- function(type, group, parameter, coef1, coef2, stats) {
    rows[[length(rows) + 1L]] <<- data.table::data.table(
      type = type,
      group = group,
      parameter = parameter,
      coef1 = coef1,
      coef2 = coef2,
      estimate = stats$estimate,
      est_error = stats$est_error,
      lower = stats$lower,
      upper = stats$upper
    )
  }

  for (parameter in rownames(fixef_summary)) {
    add("fixed", NA_character_, parameter, NA_character_, NA_character_,
        .mlsim_recovery_summary_row(fixef_summary, parameter, probs))
  }

  for (group in setdiff(names(varcorr %||% list()), "residual__")) {
    sd_summary <- varcorr[[group]]$sd
    for (coefficient in rownames(sd_summary)) {
      add("random_sd", group, coefficient, NA_character_, NA_character_,
          .mlsim_recovery_summary_row(sd_summary, coefficient, probs))
    }
    cor_array <- varcorr[[group]]$cor
    if (!is.null(cor_array)) {
      coefficients <- dimnames(cor_array)[[1L]]
      quantile_columns <- paste0("Q", probs * 100)
      n <- length(coefficients)
      for (j in seq_len(max(n - 1L, 0L))) {
        for (i in seq((j + 1L), n)) {
          # VarCorr() zero-fills correlations the model never sampled (e.g.
          # across unlinked mu/dpar blocks); only pairs with a posterior
          # variable are real estimates.
          sampled <- any(sprintf(
            "cor_%s__%s__%s",
            group,
            c(coefficients[[j]], coefficients[[i]]),
            c(coefficients[[i]], coefficients[[j]])
          ) %in% variables)
          if (!sampled) {
            next
          }
          pair_summary <- cor_array[i, , j]
          add(
            "random_cor", group,
            sprintf("cor(%s,%s)", coefficients[[j]], coefficients[[i]]),
            coefficients[[j]], coefficients[[i]],
            list(
              estimate = pair_summary[["Estimate"]],
              est_error = pair_summary[["Est.Error"]],
              lower = pair_summary[[quantile_columns[[1L]]]],
              upper = pair_summary[[quantile_columns[[2L]]]]
            )
          )
        }
      }
    }
  }

  if (!is.null(rescor)) {
    for (i in seq_along(rescor$coef1)) {
      add("rescor", NA_character_,
          sprintf("rescor(%s,%s)", rescor$coef1[[i]], rescor$coef2[[i]]),
          rescor$coef1[[i]], rescor$coef2[[i]],
          .mlsim_recovery_summary_row(rescor$summary, i, probs))
    }
  }

  data.table::rbindlist(rows)
}

# The join key treats correlation pairs as unordered so truth and estimate
# rows match regardless of which coefficient brms lists first.
.mlsim_recovery_keys <- function(table) {
  if (is.null(table) || nrow(table) == 0L) {
    return(character())
  }
  pair <- !is.na(table$coef1) & !is.na(table$coef2)
  name <- table$parameter
  name[pair] <- paste0(
    pmin(table$coef1[pair], table$coef2[pair]),
    "::",
    pmax(table$coef1[pair], table$coef2[pair])
  )
  paste(table$type, ifelse(is.na(table$group), "", table$group), name, sep = "\r")
}

.mlsim_recovery_join <- function(truth, estimates) {
  truth_keys <- .mlsim_recovery_keys(truth)
  estimate_keys <- .mlsim_recovery_keys(estimates)
  missing_estimates <- truth$parameter[!truth_keys %in% estimate_keys]
  if (length(missing_estimates) > 0L) {
    .mlsim_stop(
      "Simulation truth parameters were not found in the fitted model: %s.",
      paste(sprintf("`%s`", missing_estimates), collapse = ", ")
    )
  }
  extra_estimates <- estimates$parameter[!estimate_keys %in% truth_keys]
  if (length(extra_estimates) > 0L) {
    .mlsim_stop(
      "Fitted parameters have no simulation truth: %s. The fitted model does not match `analysis$formula`.",
      paste(sprintf("`%s`", extra_estimates), collapse = ", ")
    )
  }
  matched <- estimates[match(truth_keys, estimate_keys), ]
  out <- data.table::data.table(
    type = truth$type,
    group = truth$group,
    parameter = truth$parameter,
    truth = truth$truth,
    estimate = matched$estimate,
    est_error = matched$est_error,
    lower = matched$lower,
    upper = matched$upper,
    bias = matched$estimate - truth$truth,
    covered = truth$truth >= matched$lower & truth$truth <= matched$upper,
    simulator_name = truth$simulator_name
  )
  out
}
