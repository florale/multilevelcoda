#' Generate Dynamic Gaussian and Compositional Outcomes
#'
#' Create an outcome generator for [simulate_data()]. Outcomes are simulated on
#' a Gaussian scale, optionally as ILR coordinates that are back-transformed to
#' strictly positive compositions. If `ar1()` appears in the location formula,
#' it defines a residual VAR(1) process, not an observed lagged-outcome
#' predictor.
#'
#' @param formula Outcome location formula. The left-hand side may be a single
#'   outcome (`y`) or `mvbind(y1, y2, ...)`. The right-hand side may include
#'   ordinary model terms, `between(x)`, `within(x)`, `ar1()`, interactions, and
#'   one brms/lme4-style grouping term.
#' @param scale Required scale formula with `sigma` on the left-hand side, for
#'   example `sigma ~ 1` or `sigma ~ treatment + (1 | ID)`. The scale model is
#'   on the log conditional standard-deviation scale.
#' @param params List of true parameters. Required components are
#'   `params$location$beta`, `params$scale$beta`, and, for multivariate
#'   outcomes, `params$scale$correlation`. When AR terms are present,
#'   `params$ar$beta` is required. When grouping terms are present,
#'   `params$random[[group]]$covariance` is required.
#' @param burnin Fixed non-negative integer burn-in length used when AR is
#'   active. Ignored when no AR terms are present.
#' @param composition List controlling optional ILR back-transformation. Use
#'   `parts` or `sbp` to request compositional output, `total` for the closure
#'   total, and `keep_ilr` to keep ILR coordinates alongside parts.
#' @param ar_stability Handling for unstable AR matrices: `"resample"`,
#'   `"shrink"`, or `"error"`.
#' @param max_stability_attempts Positive integer maximum number of resampling
#'   or shrinkage attempts.
#' @param shrink_target_radius Target spectral radius used by
#'   `ar_stability = "shrink"`.
#'
#' @return An `mlsim_generator_spec` object for use in [simulate_data()].
#'
#' @examples
#' beta_location <- matrix(
#'   c(0, 0, 0.2, -0.1),
#'   nrow = 2,
#'   dimnames = list(c("(Intercept)", "treatmenttreatment"), c("ilr1", "ilr2"))
#' )
#' beta_scale <- matrix(
#'   log(c(0.4, 0.35)),
#'   nrow = 1,
#'   dimnames = list("(Intercept)", c("ilr1", "ilr2"))
#' )
#' beta_ar <- array(
#'   c(0.25, 0.02, -0.01, 0.2),
#'   dim = c(1, 2, 2),
#'   dimnames = list("ar1()", c("ilr1", "ilr2"), c("ilr1", "ilr2"))
#' )
#' corr <- diag(2)
#' dimnames(corr) <- list(c("ilr1", "ilr2"), c("ilr1", "ilr2"))
#'
#' sim <- simulate_data(
#'   n_groups = 4,
#'   n_per_group = 4,
#'   group_id = "ID",
#'   time_id = "day",
#'   seed = 1,
#'   generators = list(
#'     treatment = gen_categorical(
#'       "treatment",
#'       level = "level2",
#'       categories = c("control", "treatment"),
#'       fixed_intercept = stats::qlogis(0.5)
#'     ),
#'     outcome = gen_outcome(
#'       mvbind(ilr1, ilr2) ~ treatment + ar1(),
#'       scale = sigma ~ 1,
#'       params = list(
#'         location = list(beta = beta_location),
#'         scale = list(beta = beta_scale, correlation = corr),
#'         ar = list(beta = beta_ar)
#'       ),
#'       burnin = 10,
#'       composition = list(parts = c("sleep", "active", "sedentary"), total = 24)
#'     )
#'   )
#' )
#' head(sim$data)
#'
#' @family predictor generators
#' @export
gen_outcome <- function(formula, scale, params, burnin,
                        composition = list(
                          parts = NULL,
                          total = 24,
                          sbp = NULL,
                          keep_ilr = TRUE
                        ),
                        ar_stability = c("resample", "shrink", "error"),
                        max_stability_attempts = 1000,
                        shrink_target_radius = 0.98) {
  if (missing(formula)) {
    .mlsim_stop("`formula` is required.")
  }
  if (missing(scale)) {
    .mlsim_stop("`scale` is required. Use `scale = sigma ~ 1` for constant conditional SDs.")
  }
  if (missing(params)) {
    .mlsim_stop("`params` is required.")
  }
  if (missing(burnin)) {
    .mlsim_stop("`burnin` is required and must be a fixed non-negative integer.")
  }
  ar_stability <- match.arg(ar_stability, c("resample", "shrink", "error"))
  controls <- .mlsim_prepare_outcome_controls(
    formula = formula,
    scale = scale,
    burnin = burnin,
    composition = composition,
    ar_stability = ar_stability,
    max_stability_attempts = max_stability_attempts,
    shrink_target_radius = shrink_target_radius
  )
  formula <- controls$formula
  scale <- controls$scale
  burnin <- controls$burnin
  composition <- controls$composition
  ar_stability <- controls$ar_stability
  max_stability_attempts <- controls$max_stability_attempts
  shrink_target_radius <- controls$shrink_target_radius

  simulate <- function(context) {
    spec <- .mlsim_build_outcome_spec(context, formula, scale, composition)
    .mlsim_simulate_outcome(
      context = context,
      spec = spec,
      params = params,
      burnin = burnin,
      ar_stability = ar_stability,
      max_stability_attempts = max_stability_attempts,
      shrink_target_radius = shrink_target_radius
    )
  }

  .mlsim_spec(
    "outcome",
    character(),
    "single",
    simulate,
    formula = formula,
    scale = scale,
    burnin = burnin,
    composition = composition,
    ar_stability = ar_stability,
    max_stability_attempts = max_stability_attempts,
    shrink_target_radius = shrink_target_radius
  )
}

#' Create a Parameter Template for `gen_outcome()`
#'
#' Build a generator that parses a [gen_outcome()] specification and records a
#' complete parameter template without simulating outcome values. Use this
#' inside [simulate_data()] after any generators that define predictors used by
#' the outcome formula, especially predictors referenced by `between()` or
#' `within()`.
#'
#' @inheritParams gen_outcome
#'
#' @details
#' The generated template is stored at
#' `sim$generator_metadata[[name]]$params`, where `name` is the name given to
#' this generator in the `generators` list. The template must be created with
#' the same simulation design, previous generators, factor levels, formula,
#' scale formula, and composition settings that will be used for the later
#' [gen_outcome()] call.
#'
#' Template values are initialized to zero for regression, AR, and group-level
#' covariance parameters, and to identity for residual ILR correlations. The
#' object is ready to edit and pass as `params` to [gen_outcome()].
#'
#' @return An `mlsim_generator_spec` object for use in [simulate_data()]. It
#'   emits no data columns and stores the parameter template in generator
#'   metadata.
#'
#' @examples
#' template_sim <- simulate_data(
#'   n_groups = 3,
#'   n_per_group = 2,
#'   group_id = "ID",
#'   seed = 1,
#'   generators = list(
#'     treatment = gen_categorical(
#'       "treatment",
#'       level = "level2",
#'       categories = c("control", "treatment"),
#'       fixed_intercept = stats::qlogis(0.5)
#'     ),
#'     outcome_template = gen_template(
#'       mvbind(ilr1, ilr2) ~ treatment,
#'       scale = sigma ~ 1,
#'       burnin = 0,
#'       composition = list(parts = c("sleep", "active", "sedentary"), total = 24)
#'     )
#'   )
#' )
#'
#' params <- template_sim$generator_metadata$outcome_template$params
#' params$location$beta["treatmenttreatment", "ilr1"] <- 0.2
#'
#' @family predictor generators
#' @export
gen_template <- function(formula, scale, burnin,
                         composition = list(
                           parts = NULL,
                           total = 24,
                           sbp = NULL,
                           keep_ilr = TRUE
                         ),
                         ar_stability = c("resample", "shrink", "error"),
                         max_stability_attempts = 1000,
                         shrink_target_radius = 0.98) {
  if (missing(formula)) {
    .mlsim_stop("`formula` is required.")
  }
  if (missing(scale)) {
    .mlsim_stop("`scale` is required. Use `scale = sigma ~ 1` for constant conditional SDs.")
  }
  if (missing(burnin)) {
    .mlsim_stop("`burnin` is required and must be a fixed non-negative integer.")
  }
  ar_stability <- match.arg(ar_stability, c("resample", "shrink", "error"))
  controls <- .mlsim_prepare_outcome_controls(
    formula = formula,
    scale = scale,
    burnin = burnin,
    composition = composition,
    ar_stability = ar_stability,
    max_stability_attempts = max_stability_attempts,
    shrink_target_radius = shrink_target_radius
  )
  formula <- controls$formula
  scale <- controls$scale
  burnin <- controls$burnin
  composition <- controls$composition
  ar_stability <- controls$ar_stability
  max_stability_attempts <- controls$max_stability_attempts
  shrink_target_radius <- controls$shrink_target_radius

  simulate <- function(context) {
    spec <- .mlsim_build_outcome_spec(context, formula, scale, composition)
    params <- .mlsim_outcome_template_params(spec)
    metadata <- list(
      distribution = "outcome_template",
      level = if (is.null(context$n_groups)) "single" else "multilevel",
      vars = character(),
      formula = spec$formula,
      scale = spec$scale,
      params = params,
      parsed = .mlsim_outcome_parsed_metadata(spec),
      term_map = spec$term_map,
      expected_parameter_names = spec$expected_parameter_names,
      selected_column_roles = spec$selected_column_roles,
      composition = spec$composition,
      ar_stability = ar_stability,
      max_stability_attempts = max_stability_attempts,
      shrink_target_radius = shrink_target_radius,
      burnin = list(
        burnin = burnin,
        burnin_type = "fixed",
        burnin_reference_regime = if (spec$has_ar) "first_observed_row" else NULL
      ),
      time_validation = list(required_complete_equally_spaced = spec$has_ar)
    )
    .mlsim_result(data.table::data.table(), character(), metadata)
  }

  .mlsim_spec(
    "template",
    character(),
    "single",
    simulate,
    formula = formula,
    scale = scale,
    burnin = burnin,
    composition = composition,
    ar_stability = ar_stability,
    max_stability_attempts = max_stability_attempts,
    shrink_target_radius = shrink_target_radius
  )
}

.mlsim_prepare_outcome_controls <- function(formula, scale, burnin, composition,
                                            ar_stability, max_stability_attempts,
                                            shrink_target_radius) {
  formula <- stats::as.formula(formula)
  scale <- stats::as.formula(scale)
  burnin <- .mlsim_check_burnin(burnin)
  ar_stability <- match.arg(ar_stability, c("resample", "shrink", "error"))
  max_stability_attempts <- .mlsim_recycle_integer(max_stability_attempts, 1L, "max_stability_attempts")
  if (max_stability_attempts < 1L) {
    .mlsim_stop("`max_stability_attempts` must be at least 1.")
  }
  shrink_target_radius <- .mlsim_check_open_prob(shrink_target_radius, "shrink_target_radius")
  if (length(shrink_target_radius) != 1L) {
    .mlsim_stop("`shrink_target_radius` must be a scalar in (0, 1).")
  }
  composition <- .mlsim_normalize_outcome_composition(composition)
  list(
    formula = formula,
    scale = scale,
    burnin = burnin,
    composition = composition,
    ar_stability = ar_stability,
    max_stability_attempts = max_stability_attempts,
    shrink_target_radius = shrink_target_radius
  )
}

.mlsim_check_burnin <- function(burnin) {
  burnin <- .mlsim_recycle_integer(burnin, 1L, "burnin")
  if (burnin < 0L) {
    .mlsim_stop("`burnin` must be non-negative.")
  }
  burnin
}

.mlsim_normalize_outcome_composition <- function(composition) {
  if (is.null(composition)) {
    composition <- list(parts = NULL, total = 24, sbp = NULL, keep_ilr = TRUE)
  }
  if (!is.list(composition)) {
    .mlsim_stop("`composition` must be a list.")
  }
  parts <- composition$parts
  sbp <- composition$sbp
  total <- composition$total %||% 24
  keep_ilr <- composition$keep_ilr %||% TRUE
  keep_ilr <- .mlsim_check_scalar_logical(keep_ilr, "composition$keep_ilr")
  total <- .mlsim_check_positive(total, "composition$total")
  if (length(total) != 1L) {
    .mlsim_stop("`composition$total` must be a scalar positive value.")
  }
  list(parts = parts, total = total, sbp = sbp, keep_ilr = keep_ilr)
}

.mlsim_build_outcome_spec <- function(context, formula, scale, composition) {
  data <- data.table::copy(context$data)
  roles <- .mlsim_collect_column_roles(context$generator_metadata, data)
  special_env <- .mlsim_special_term_env(data, roles, context)

  outcomes <- .mlsim_parse_outcome_lhs(formula[[2L]])
  .mlsim_check_vars(outcomes, NULL, "outcome names")
  d <- length(outcomes)

  location_parts <- .mlsim_split_formula_rhs(formula[[3L]], component = "location")
  location_expr <- .mlsim_transform_specials(
    location_parts$fixed %||% 1,
    special_env,
    allow_ar = TRUE,
    component = "location"
  )
  data <- special_env$data
  if (isTRUE(special_env$has_ar)) {
    data[[.mlsim_ar_placeholder()]] <- 1
    special_env$data[[.mlsim_ar_placeholder()]] <- 1
  }
  location_matrix <- .mlsim_build_model_matrix(
    location_expr,
    data,
    component = "location",
    term_map = special_env$term_map,
    allow_ar = TRUE
  )
  location <- .mlsim_split_ar_matrix(location_matrix)

  scale_parts <- .mlsim_split_scale_formula(scale)
  scale_expr <- .mlsim_transform_specials(
    scale_parts$fixed,
    special_env,
    allow_ar = FALSE,
    component = "scale"
  )
  data <- special_env$data
  scale_matrix <- .mlsim_build_model_matrix(
    scale_expr,
    data,
    component = "scale",
    term_map = special_env$term_map,
    allow_ar = FALSE
  )

  random <- .mlsim_build_outcome_random_spec(
    location_random = location_parts$random,
    scale_random = scale_parts$random,
    special_env = special_env,
    data = data,
    context = context,
    location_terms = colnames(location$X),
    ar_terms = colnames(location$Q),
    scale_terms = colnames(scale_matrix$matrix),
    outcomes = outcomes
  )

  composition <- .mlsim_prepare_outcome_composition(composition, outcomes)
  has_ar <- ncol(location$Q) > 0L
  if (has_ar) {
    .mlsim_validate_outcome_time(context, data)
  }

  expected <- list(
    location = colnames(location$X),
    scale = colnames(scale_matrix$matrix),
    ar = colnames(location$Q),
    random = random$names
  )

  list(
    formula = formula,
    scale = scale,
    outcomes = outcomes,
    X = location$X,
    Q = location$Q,
    W = scale_matrix$matrix,
    has_ar = has_ar,
    random = random,
    expected_parameter_names = expected,
    data = data,
    selected_column_roles = special_env$selected_roles,
    term_map = special_env$term_map,
    composition = composition,
    series = .mlsim_outcome_series(context, data, has_ar)
  )
}

.mlsim_parse_outcome_lhs <- function(lhs) {
  if (is.symbol(lhs)) {
    return(as.character(lhs))
  }
  if (is.call(lhs) && identical(as.character(lhs[[1L]]), "mvbind")) {
    args <- as.list(lhs[-1L])
    if (length(args) == 0L || !all(vapply(args, is.symbol, logical(1)))) {
      .mlsim_stop("`mvbind()` response terms must be variable names.")
    }
    return(vapply(args, as.character, character(1)))
  }
  .mlsim_stop("Outcome formula left-hand side must be `y` or `mvbind(y1, y2, ...)`.")
}

.mlsim_split_scale_formula <- function(scale) {
  if (length(scale) != 3L || !is.symbol(scale[[2L]]) || !identical(as.character(scale[[2L]]), "sigma")) {
    .mlsim_stop(
      "The scale formula must have `sigma` on the left-hand side, and ar1() is not allowed in the scale formula. The scale model defines row-wise conditional residual/innovation variability, not scale dynamics."
    )
  }
  .mlsim_split_formula_rhs(scale[[3L]], component = "scale")
}

.mlsim_split_formula_rhs <- function(expr, component) {
  .mlsim_reject_minus_terms(expr, component)
  random <- .mlsim_collect_random_terms(expr, component)
  fixed <- .mlsim_remove_random_terms(expr)
  if (is.null(fixed)) {
    fixed <- 1
  }
  list(fixed = fixed, random = random)
}

.mlsim_reject_minus_terms <- function(expr, component) {
  if (!is.call(expr)) {
    return(invisible(TRUE))
  }
  op <- as.character(expr[[1L]])
  if (identical(op, "-")) {
    .mlsim_stop("The %s formula must not remove terms or the intercept.", component)
  }
  lapply(as.list(expr[-1L]), .mlsim_reject_minus_terms, component = component)
  invisible(TRUE)
}

.mlsim_is_random_term <- function(expr) {
  if (is.call(expr) && identical(as.character(expr[[1L]]), "(")) {
    return(.mlsim_is_random_term(expr[[2L]]))
  }
  is.call(expr) && identical(as.character(expr[[1L]]), "|")
}

.mlsim_collect_random_terms <- function(expr, component) {
  if (.mlsim_is_random_term(expr)) {
    return(list(.mlsim_parse_random_term(expr, component)))
  }
  if (is.call(expr) && identical(as.character(expr[[1L]]), "+")) {
    return(c(
      .mlsim_collect_random_terms(expr[[2L]], component),
      .mlsim_collect_random_terms(expr[[3L]], component)
    ))
  }
  list()
}

.mlsim_remove_random_terms <- function(expr) {
  if (.mlsim_is_random_term(expr)) {
    return(NULL)
  }
  if (is.call(expr) && identical(as.character(expr[[1L]]), "+")) {
    left <- .mlsim_remove_random_terms(expr[[2L]])
    right <- .mlsim_remove_random_terms(expr[[3L]])
    if (is.null(left)) {
      return(right)
    }
    if (is.null(right)) {
      return(left)
    }
    return(call("+", left, right))
  }
  expr
}

.mlsim_parse_random_term <- function(expr, component) {
  if (is.call(expr) && identical(as.character(expr[[1L]]), "(")) {
    expr <- expr[[2L]]
  }
  if (!is.call(expr) || !identical(as.character(expr[[1L]]), "|")) {
    .mlsim_stop("Internal random-effect parser received a non-random term.")
  }
  if (is.call(expr[[2L]]) && identical(as.character(expr[[2L]][[1L]]), "|")) {
    .mlsim_stop("Three-part random-effect syntax is not supported; use two-part terms such as `(1 + x | ID)`.")
  }
  group <- expr[[3L]]
  if (!is.symbol(group)) {
    .mlsim_stop("Grouping factors in %s random effects must be simple variable names.", component)
  }
  list(component = component, expr = expr[[2L]], group = as.character(group))
}

.mlsim_special_term_env <- function(data, roles, context) {
  env <- new.env(parent = emptyenv())
  env$data <- data
  env$roles <- roles
  env$context <- context
  env$term_map <- character()
  env$selected_roles <- data.table::data.table(
    column = character(),
    variable = character(),
    component = character(),
    level = character(),
    generator = character()
  )
  env$has_ar <- FALSE
  env
}

.mlsim_ar_placeholder <- function() ".__mlsim_ar1__"

.mlsim_transform_specials <- function(expr, env, allow_ar, component) {
  if (is.symbol(expr) || is.atomic(expr)) {
    return(expr)
  }
  if (!is.call(expr)) {
    return(expr)
  }
  op <- as.character(expr[[1L]])
  .mlsim_reject_multiple_ar_calls(expr)
  if (identical(op, "(")) {
    return(call("(", .mlsim_transform_specials(expr[[2L]], env, allow_ar, component)))
  }
  if (op %in% c("+", ":", "*", "/")) {
    args <- lapply(as.list(expr[-1L]), .mlsim_transform_specials, env = env, allow_ar = allow_ar, component = component)
    return(as.call(c(list(as.name(op)), args)))
  }
  if (identical(op, "between") || identical(op, "within")) {
    if (length(expr) != 2L || !is.symbol(expr[[2L]])) {
      .mlsim_stop("`%s()` terms must contain one variable name.", op)
    }
    return(as.name(.mlsim_resolve_special_role(as.character(expr[[2L]]), op, env)))
  }
  if (identical(op, "ar1")) {
    if (!isTRUE(allow_ar)) {
      .mlsim_stop(
        "The scale formula must have `sigma` on the left-hand side, and ar1() is not allowed in the scale formula. The scale model defines row-wise conditional residual/innovation variability, not scale dynamics."
      )
    }
    if (length(expr) != 1L) {
      .mlsim_stop("`ar1()` does not take arguments.")
    }
    env$has_ar <- TRUE
    placeholder <- .mlsim_ar_placeholder()
    env$term_map[placeholder] <- "ar1()"
    return(as.name(placeholder))
  }
  if (.mlsim_expr_has_ar(expr)) {
    .mlsim_stop("`ar1()` may only appear as a main effect or inside interactions in the location formula.")
  }
  expr
}

.mlsim_expr_has_ar <- function(expr) {
  if (!is.call(expr)) {
    return(FALSE)
  }
  if (identical(as.character(expr[[1L]]), "ar1")) {
    return(TRUE)
  }
  any(vapply(as.list(expr[-1L]), .mlsim_expr_has_ar, logical(1)))
}

.mlsim_count_ar_calls <- function(expr) {
  if (!is.call(expr)) {
    return(0L)
  }
  count <- if (identical(as.character(expr[[1L]]), "ar1")) 1L else 0L
  count + sum(vapply(as.list(expr[-1L]), .mlsim_count_ar_calls, integer(1)))
}

.mlsim_reject_multiple_ar_calls <- function(expr) {
  if (!is.call(expr)) {
    return(invisible(TRUE))
  }
  op <- as.character(expr[[1L]])
  if (identical(op, "(")) {
    return(.mlsim_reject_multiple_ar_calls(expr[[2L]]))
  }
  if (identical(op, "+")) {
    lapply(as.list(expr[-1L]), .mlsim_reject_multiple_ar_calls)
    return(invisible(TRUE))
  }
  if (.mlsim_count_ar_calls(expr) > 1L) {
    .mlsim_stop("Terms containing more than one `ar1()` are not supported.")
  }
  invisible(TRUE)
}

.mlsim_resolve_special_role <- function(variable, component, env) {
  roles <- env$roles
  target_variable <- variable
  target_component <- component
  matches <- roles[roles[["variable"]] == target_variable & roles[["component"]] == target_component]
  if (nrow(matches) != 1L) {
    .mlsim_stop(
      "`%s(%s)` requires exactly one prior generator column role with variable `%s` and component `%s`; found %d.",
      component,
      variable,
      variable,
      component,
      nrow(matches)
    )
  }
  source_col <- matches$column[[1L]]
  safe_col <- sprintf(".__mlsim_%s_%s__", component, make.names(variable))
  env$data[[safe_col]] <- env$data[[source_col]]
  env$term_map[safe_col] <- sprintf("%s(%s)", component, variable)
  env$selected_roles <- unique(data.table::rbindlist(list(env$selected_roles, matches), fill = TRUE))
  safe_col
}

.mlsim_build_model_matrix <- function(expr, data, component, term_map, allow_ar) {
  form <- stats::as.formula(call("~", expr))
  trm <- stats::terms(form, data = data)
  if (identical(attr(trm, "intercept"), 0L)) {
    .mlsim_stop("The %s formula must include an intercept.", component)
  }
  mm <- tryCatch(
    stats::model.matrix(form, data = data),
    error = function(e) .mlsim_stop("Could not build the %s model matrix: %s", component, conditionMessage(e))
  )
  term_labels <- attr(trm, "term.labels")
  assign <- attr(mm, "assign")
  column_terms <- rep("(Intercept)", length(assign))
  non_intercept <- assign > 0L
  column_terms[non_intercept] <- term_labels[assign[non_intercept]]
  readable_columns <- .mlsim_restore_term_names(colnames(mm), term_map)
  readable_terms <- .mlsim_restore_term_names(column_terms, term_map)
  colnames(mm) <- readable_columns
  list(matrix = mm, column_terms = readable_terms)
}

.mlsim_restore_term_names <- function(x, term_map) {
  if (length(x) == 0L || length(term_map) == 0L) {
    return(x)
  }
  out <- x
  keys <- names(term_map)
  keys <- keys[order(nchar(keys), decreasing = TRUE)]
  for (key in keys) {
    out <- gsub(key, term_map[[key]], out, fixed = TRUE)
  }
  out
}

.mlsim_split_ar_matrix <- function(model) {
  ar_cols <- grepl("ar1\\(\\)", model$column_terms)
  X <- model$matrix[, !ar_cols, drop = FALSE]
  Q <- model$matrix[, ar_cols, drop = FALSE]
  colnames(Q) <- colnames(model$matrix)[ar_cols]
  if (ncol(X) == 0L) {
    X <- matrix(1, nrow(model$matrix), 1L, dimnames = list(NULL, "(Intercept)"))
  }
  list(X = X, Q = Q)
}

.mlsim_build_outcome_random_spec <- function(location_random, scale_random,
                                             special_env, data, context,
                                             location_terms, ar_terms,
                                             scale_terms, outcomes) {
  all_random <- c(location_random, scale_random)
  if (length(all_random) == 0L) {
    return(list(group = NULL, location_Z = NULL, ar_Z = NULL, scale_Z = NULL,
                location_terms = character(), ar_terms = character(),
                scale_terms = character(), names = character()))
  }
  groups <- unique(vapply(all_random, `[[`, character(1), "group"))
  if (length(groups) > 1L) {
    .mlsim_stop("Outcome simulation supports at most one grouping factor.")
  }
  group <- groups[[1L]]
  if (is.null(context$n_groups)) {
    .mlsim_stop("Outcome group-level effects require `simulate_data()` to use `n_groups`.")
  }
  if (!identical(group, context$group_id)) {
    .mlsim_stop("Outcome group-level effects must use the active grouping factor `%s`.", context$group_id)
  }

  loc_Z <- NULL
  ar_Z <- NULL
  sc_Z <- NULL
  for (term in location_random) {
    transformed <- .mlsim_transform_specials(term$expr, special_env, allow_ar = TRUE, component = "location")
    data <- special_env$data
    mat <- .mlsim_build_model_matrix(transformed, data, "location random-effect", special_env$term_map, allow_ar = TRUE)
    split <- .mlsim_split_ar_matrix(mat)
    loc_Z <- .mlsim_cbind_unique(loc_Z, split$X)
    ar_Z <- .mlsim_cbind_unique(ar_Z, split$Q)
  }
  for (term in scale_random) {
    transformed <- .mlsim_transform_specials(term$expr, special_env, allow_ar = FALSE, component = "scale")
    data <- special_env$data
    mat <- .mlsim_build_model_matrix(transformed, data, "scale random-effect", special_env$term_map, allow_ar = FALSE)
    sc_Z <- .mlsim_cbind_unique(sc_Z, mat$matrix)
  }

  loc_terms <- colnames(loc_Z) %||% character()
  ar_random_terms <- colnames(ar_Z) %||% character()
  sc_terms <- colnames(sc_Z) %||% character()
  .mlsim_check_random_terms_match(loc_terms, location_terms, "location")
  .mlsim_check_random_terms_match(ar_random_terms, ar_terms, "AR")
  .mlsim_check_random_terms_match(sc_terms, scale_terms, "scale")

  names <- c(
    .mlsim_random_effect_names("location", loc_terms, outcomes),
    .mlsim_random_effect_names("ar", ar_random_terms, outcomes),
    .mlsim_random_effect_names("scale", sc_terms, outcomes)
  )
  list(
    group = group,
    location_Z = loc_Z,
    ar_Z = ar_Z,
    scale_Z = sc_Z,
    location_terms = loc_terms,
    ar_terms = ar_random_terms,
    scale_terms = sc_terms,
    names = names
  )
}

.mlsim_cbind_unique <- function(x, y) {
  if (is.null(y) || ncol(y) == 0L) {
    return(x)
  }
  if (is.null(x)) {
    return(y)
  }
  new_cols <- colnames(y)[!colnames(y) %in% colnames(x)]
  if (length(new_cols) == 0L) {
    return(x)
  }
  cbind(x, y[, new_cols, drop = FALSE])
}

.mlsim_check_random_terms_match <- function(random_terms, population_terms, component) {
  missing <- random_terms[!random_terms %in% population_terms]
  if (length(missing) > 0L) {
    .mlsim_stop(
      "%s group-level terms must have matching population-level terms in the same component; missing: %s.",
      component,
      paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}

.mlsim_random_effect_names <- function(component, terms, outcomes) {
  if (length(terms) == 0L) {
    return(character())
  }
  if (identical(component, "ar")) {
    return(unlist(lapply(terms, function(term) {
      unlist(lapply(outcomes, function(to) {
        sprintf("ar|term=%s|to=%s|from=%s", term, to, outcomes)
      }), use.names = FALSE)
    }), use.names = FALSE))
  }
  as.vector(vapply(terms, function(term) {
    sprintf("%s|outcome=%s|term=%s", component, outcomes, term)
  }, character(length(outcomes))))
}

.mlsim_prepare_outcome_composition <- function(composition, outcomes) {
  active <- !is.null(composition$parts) || !is.null(composition$sbp)
  if (!isTRUE(active)) {
    return(c(composition, list(active = FALSE, resolved = NULL)))
  }
  resolved <- .mlsim_resolve_composition(composition$parts, composition$sbp, length(outcomes))
  .mlsim_check_composition_output_names(outcomes, resolved$parts, composition$keep_ilr)
  c(composition, list(active = TRUE, resolved = resolved))
}

.mlsim_validate_outcome_time <- function(context, data) {
  if (is.null(context$time_id)) {
    .mlsim_stop("`ar1()` outcome simulation requires `simulate_data()` to supply `time_id`.")
  }
  time <- data[[context$time_id]]
  if (!(is.numeric(time) || inherits(time, "Date") || inherits(time, "POSIXt"))) {
    .mlsim_stop("`ar1()` outcome simulation requires numeric, Date, or POSIXct time values.")
  }
  groups <- if (is.null(context$group_index)) rep(1L, nrow(data)) else context$group_index
  for (g in unique(groups)) {
    idx <- which(groups == g)
    ordered_time <- sort(time[idx])
    if (anyDuplicated(ordered_time)) {
      .mlsim_stop("`ar1()` outcome simulation requires unique time values within each group.")
    }
    if (length(ordered_time) > 2L) {
      diffs <- diff(as.numeric(ordered_time))
      if (any(diffs <= 0) || max(abs(diffs - diffs[[1L]])) > 1e-8) {
        .mlsim_stop("`ar1()` outcome simulation requires complete equally spaced time within each group.")
      }
    }
  }
  invisible(TRUE)
}

.mlsim_outcome_series <- function(context, data, has_ar) {
  groups <- if (is.null(context$group_index)) rep(1L, nrow(data)) else context$group_index
  if (isTRUE(has_ar)) {
    order_index <- if (is.null(context$group_index)) {
      order(data[[context$time_id]])
    } else {
      order(groups, data[[context$time_id]])
    }
  } else {
    order_index <- seq_len(nrow(data))
  }
  ordered_groups <- groups[order_index]
  rows <- split(order_index, ordered_groups)
  list(group_index = groups, order = order_index, rows = rows)
}

.mlsim_simulate_outcome <- function(context, spec, params, burnin, ar_stability,
                                    max_stability_attempts, shrink_target_radius) {
  checked <- .mlsim_validate_outcome_params(spec, params)
  random <- .mlsim_draw_outcome_random(
    spec,
    checked,
    ar_stability,
    max_stability_attempts,
    shrink_target_radius
  )
  mu <- .mlsim_outcome_location(spec, checked, random$draws)
  scale_eta <- .mlsim_outcome_scale(spec, checked, random$draws)
  sigma <- exp(scale_eta)
  phi <- .mlsim_outcome_phi(spec, checked, random$draws)
  sim <- .mlsim_draw_outcome_values(spec, mu, sigma, checked$scale$correlation, phi, burnin)
  values <- sim$z
  colnames(values) <- spec$outcomes
  output <- values
  metadata_composition <- list(compositional = FALSE)
  if (isTRUE(spec$composition$active)) {
    comp <- .mlsim_ilr_inverse(
      values,
      parts = spec$composition$resolved$parts,
      total = spec$composition$total,
      coordinate_names = spec$outcomes,
      sbp = spec$composition$resolved$sbp
    )
    output <- if (isTRUE(spec$composition$keep_ilr)) cbind(values, comp$values) else comp$values
    metadata_composition <- list(
      compositional = TRUE,
      parts = colnames(comp$values),
      sbp = comp$sbp,
      basis = comp$basis,
      basis_matrix = comp$basis,
      basis_source = "sbp",
      ilr_coordinate_map = comp$coordinate_map,
      total = comp$total,
      keep_ilr = spec$composition$keep_ilr,
      ilr_vars = spec$outcomes
    )
  }
  output <- as.matrix(output)

  metadata <- c(
    list(
      distribution = "outcome",
      level = if (is.null(context$n_groups)) "single" else "multilevel",
      vars = colnames(output),
      formula = spec$formula,
      scale = spec$scale,
      params = checked$raw,
      parsed = .mlsim_outcome_parsed_metadata(spec),
      term_map = spec$term_map,
      expected_parameter_names = spec$expected_parameter_names,
      selected_column_roles = spec$selected_column_roles,
      group_level_effects = random$draws,
      proposal_group_level_effect_covariance = random$proposal_covariance,
      realised_group_level_effect_covariance = random$realised_covariance,
      ar = list(
        beta = checked$ar$beta,
        phi_by_row = phi$phi,
        fixed_phi_by_row = phi$fixed_phi,
        stability = random$stability
      ),
      mu = mu,
      residual = sim$residual,
      innovation = sim$innovation,
      sigma = sigma,
      ilr_conditional_correlation = checked$scale$correlation,
      Sigma_epsilon = .mlsim_mvn_row_covariance_array(sigma, checked$scale$correlation),
      burnin = list(
        burnin = burnin,
        burnin_type = "fixed",
        burnin_reference_regime = if (spec$has_ar) "first_observed_row" else NULL
      ),
      time_validation = list(required_complete_equally_spaced = spec$has_ar)
    ),
    metadata_composition
  )

  .mlsim_result(output, colnames(output), metadata)
}

.mlsim_outcome_template_params <- function(spec) {
  location_beta <- matrix(
    0,
    nrow = length(spec$expected_parameter_names$location),
    ncol = length(spec$outcomes),
    dimnames = list(spec$expected_parameter_names$location, spec$outcomes)
  )
  scale_beta <- matrix(
    0,
    nrow = length(spec$expected_parameter_names$scale),
    ncol = length(spec$outcomes),
    dimnames = list(spec$expected_parameter_names$scale, spec$outcomes)
  )
  correlation <- diag(length(spec$outcomes))
  dimnames(correlation) <- list(spec$outcomes, spec$outcomes)

  params <- list(
    location = list(beta = location_beta),
    scale = list(beta = scale_beta, correlation = correlation)
  )

  if (isTRUE(spec$has_ar)) {
    params$ar <- list(
      beta = array(
        0,
        dim = c(
          length(spec$expected_parameter_names$ar),
          length(spec$outcomes),
          length(spec$outcomes)
        ),
        dimnames = list(
          spec$expected_parameter_names$ar,
          spec$outcomes,
          spec$outcomes
        )
      )
    )
  }

  if (length(spec$random$names) > 0L) {
    covariance <- matrix(
      0,
      nrow = length(spec$random$names),
      ncol = length(spec$random$names),
      dimnames = list(spec$random$names, spec$random$names)
    )
    params$random <- list()
    params$random[[spec$random$group]] <- list(covariance = covariance)
  }

  params
}

.mlsim_outcome_parsed_metadata <- function(spec) {
  list(
    outcomes = spec$outcomes,
    location_terms = colnames(spec$X),
    scale_terms = colnames(spec$W),
    ar_terms = colnames(spec$Q),
    random_effect_names = spec$random$names
  )
}

.mlsim_validate_outcome_params <- function(spec, params) {
  if (!is.list(params)) {
    .mlsim_stop("`params` must be a list.")
  }
  location_beta <- .mlsim_check_parameter_matrix(
    params$location$beta,
    spec$expected_parameter_names$location,
    spec$outcomes,
    "params$location$beta"
  )
  scale_beta <- .mlsim_check_parameter_matrix(
    params$scale$beta,
    spec$expected_parameter_names$scale,
    spec$outcomes,
    "params$scale$beta"
  )
  correlation <- .mlsim_check_outcome_correlation(
    params$scale$correlation,
    spec$outcomes,
    "params$scale$correlation"
  )
  ar_beta <- NULL
  if (isTRUE(spec$has_ar)) {
    ar_beta <- .mlsim_check_ar_array(params$ar$beta, spec$expected_parameter_names$ar, spec$outcomes)
  }
  random_cov <- NULL
  if (length(spec$random$names) > 0L) {
    random_cov <- params$random[[spec$random$group]]$covariance
    random_cov <- .mlsim_check_exact_covariance(
      random_cov,
      spec$random$names,
      sprintf("params$random$%s$covariance", spec$random$group)
    )
  }
  list(
    raw = params,
    location = list(beta = location_beta),
    scale = list(beta = scale_beta, correlation = correlation),
    ar = list(beta = ar_beta),
    random = list(covariance = random_cov)
  )
}

.mlsim_check_parameter_matrix <- function(x, row_names, col_names, arg) {
  if (is.null(x)) {
    .mlsim_stop("`%s` is required.", arg)
  }
  x <- as.matrix(x)
  if (!identical(dim(x), c(length(row_names), length(col_names)))) {
    .mlsim_stop("`%s` must be a %d by %d matrix.", arg, length(row_names), length(col_names))
  }
  .mlsim_require_exact_dimnames(x, row_names, col_names, arg)
  .mlsim_check_finite_numeric(x, arg)
}

.mlsim_require_exact_dimnames <- function(x, row_names, col_names, arg) {
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    .mlsim_stop("`%s` must have row and column names in the required order.", arg)
  }
  if (!identical(rownames(x), row_names) || !identical(colnames(x), col_names)) {
    .mlsim_stop(
      "`%s` names must match exactly in this order. Expected rows: %s. Observed rows: %s. Expected columns: %s. Observed columns: %s.",
      arg,
      paste(row_names, collapse = ", "),
      paste(rownames(x), collapse = ", "),
      paste(col_names, collapse = ", "),
      paste(colnames(x), collapse = ", ")
    )
  }
  invisible(TRUE)
}

.mlsim_check_outcome_correlation <- function(x, outcomes, arg) {
  if (is.null(x) && length(outcomes) == 1L) {
    x <- matrix(1, 1, 1, dimnames = list(outcomes, outcomes))
  }
  if (length(outcomes) == 1L && length(x) == 1L && !is.matrix(x)) {
    x <- matrix(x, 1, 1, dimnames = list(outcomes, outcomes))
  }
  x <- .mlsim_check_exact_covariance(x, outcomes, arg)
  if (any(abs(diag(x) - 1) > 1e-8)) {
    .mlsim_stop("`%s` must have ones on the diagonal.", arg)
  }
  if (any(x < -1 - 1e-8 | x > 1 + 1e-8)) {
    .mlsim_stop("`%s` values must be correlations in [-1, 1].", arg)
  }
  x
}

.mlsim_check_exact_covariance <- function(x, names, arg) {
  if (is.null(x)) {
    .mlsim_stop("`%s` is required.", arg)
  }
  x <- as.matrix(x)
  if (!identical(dim(x), c(length(names), length(names)))) {
    .mlsim_stop("`%s` must be a %d by %d matrix.", arg, length(names), length(names))
  }
  .mlsim_require_exact_dimnames(x, names, names, arg)
  .mlsim_check_cov(x, length(names), arg, allow_zero = TRUE)
}

.mlsim_check_ar_array <- function(x, ar_terms, outcomes) {
  if (is.null(x)) {
    .mlsim_stop("`params$ar$beta` is required when `ar1()` appears in the outcome formula.")
  }
  if (!is.array(x) || length(dim(x)) != 3L) {
    .mlsim_stop("`params$ar$beta` must be a three-dimensional array: ar_term x to_outcome x from_outcome.")
  }
  expected_dim <- c(length(ar_terms), length(outcomes), length(outcomes))
  if (!identical(dim(x), expected_dim)) {
    .mlsim_stop("`params$ar$beta` must have dimensions %s.", paste(expected_dim, collapse = " x "))
  }
  dn <- dimnames(x)
  if (is.null(dn) || length(dn) != 3L ||
      !identical(dn[[1L]], ar_terms) ||
      !identical(dn[[2L]], outcomes) ||
      !identical(dn[[3L]], outcomes)) {
    .mlsim_stop("`params$ar$beta` dimnames must match AR terms, target outcomes, and source outcomes in exact order.")
  }
  .mlsim_check_finite_numeric(x, "params$ar$beta")
}

.mlsim_draw_outcome_random <- function(spec, checked, ar_stability,
                                       max_stability_attempts, shrink_target_radius) {
  names <- spec$random$names
  n_random <- length(names)
  if (n_random == 0L) {
    fixed <- .mlsim_outcome_phi(spec, checked, NULL)
    .mlsim_check_fixed_phi_stability(spec, fixed$fixed_phi)
    stability <- .mlsim_make_stability_metadata(spec, fixed$phi, NULL, ar_stability)
    stability <- .mlsim_finalize_stability_metadata(spec, stability)
    return(list(
      draws = NULL,
      proposal_covariance = NULL,
      realised_covariance = NULL,
      stability = stability
    ))
  }
  covariance <- checked$random$covariance
  n_groups <- max(spec$series$group_index)
  draws <- matrix(0, n_groups, n_random, dimnames = list(as.character(seq_len(n_groups)), names))
  attempts <- rep(1L, n_groups)
  shrink_factor <- rep(1, n_groups)
  before <- after <- rep(NA_real_, n_groups)

  fixed_phi <- .mlsim_outcome_phi(spec, checked, NULL)$fixed_phi
  .mlsim_check_fixed_phi_stability(spec, fixed_phi)
  for (g in seq_len(n_groups)) {
    if (identical(ar_stability, "resample")) {
      accepted <- FALSE
      for (attempt in seq_len(max_stability_attempts)) {
        candidate <- .mlsim_rmvnorm(1L, rep(0, n_random), covariance)
        colnames(candidate) <- names
        candidate_draws <- draws
        candidate_draws[g, ] <- candidate[1L, ]
        phi <- .mlsim_outcome_phi(spec, checked, candidate_draws)$phi
        radius <- .mlsim_group_max_radius(spec, phi, g)
        if (radius < 1 - 1e-10) {
          draws[g, ] <- candidate[1L, ]
          attempts[[g]] <- attempt
          after[[g]] <- radius
          accepted <- TRUE
          break
        }
      }
      if (!accepted) {
        .mlsim_stop("Could not draw stable AR group-level effects for group %d after %d attempts.", g, max_stability_attempts)
      }
    } else {
      candidate <- .mlsim_rmvnorm(1L, rep(0, n_random), covariance)
      colnames(candidate) <- names
      draws[g, ] <- candidate[1L, ]
      phi <- .mlsim_outcome_phi(spec, checked, draws)$phi
      radius <- .mlsim_group_max_radius(spec, phi, g)
      before[[g]] <- radius
      if (radius >= 1 - 1e-10 && identical(ar_stability, "error")) {
        .mlsim_stop("Unstable AR matrix for group %d; spectral radius %.3f.", g, radius)
      }
      if (radius >= shrink_target_radius && identical(ar_stability, "shrink")) {
        lambda <- 1
        for (attempt in seq_len(max_stability_attempts)) {
          candidate_draws <- .mlsim_shrink_ar_draws(spec, draws, g, lambda)
          phi <- .mlsim_outcome_phi(spec, checked, candidate_draws)$phi
          radius <- .mlsim_group_max_radius(spec, phi, g)
          if (radius < shrink_target_radius) {
            draws <- candidate_draws
            shrink_factor[[g]] <- lambda
            attempts[[g]] <- attempt
            after[[g]] <- radius
            break
          }
          lambda <- lambda * 0.8
        }
        if (radius >= shrink_target_radius) {
          .mlsim_stop("Could not shrink AR group-level effects for group %d below spectral radius %.3f.", g, shrink_target_radius)
        }
      } else {
        after[[g]] <- radius
      }
    }
  }
  realised_cov <- if (nrow(draws) > 1L) stats::cov(draws) else matrix(0, n_random, n_random, dimnames = list(names, names))
  final_phi <- .mlsim_outcome_phi(spec, checked, draws)$phi
  stability <- .mlsim_make_stability_metadata(spec, final_phi, attempts, ar_stability)
  stability$shrink_factor_by_group_level <- if (identical(ar_stability, "shrink")) shrink_factor else NULL
  stability$spectral_radius_before_shrinkage <- if (identical(ar_stability, "shrink")) before else NULL
  stability$spectral_radius_after_shrinkage <- if (identical(ar_stability, "shrink")) after else NULL
  stability <- .mlsim_finalize_stability_metadata(spec, stability)
  list(
    draws = draws,
    proposal_covariance = covariance,
    realised_covariance = realised_cov,
    stability = stability
  )
}

.mlsim_finalize_stability_metadata <- function(spec, stability) {
  if (isTRUE(spec$has_ar) &&
      max(stability$max_spectral_radius_by_group_level, na.rm = TRUE) > 0.95) {
    warning("Accepted AR matrices have spectral radius above 0.95; consider a longer burn-in or a sensitivity check.", call. = FALSE)
    stability$high_persistence_warning <- TRUE
  } else {
    stability$high_persistence_warning <- FALSE
  }
  stability
}

.mlsim_check_fixed_phi_stability <- function(spec, fixed_phi) {
  if (!isTRUE(spec$has_ar)) {
    return(invisible(TRUE))
  }
  radii <- .mlsim_phi_radii(fixed_phi)
  if (any(radii >= 1 - 1e-10)) {
    idx <- which.max(radii)
    .mlsim_stop("Population-level AR matrix is unstable at row %d; spectral radius %.3f.", idx, radii[[idx]])
  }
  invisible(TRUE)
}

.mlsim_shrink_ar_draws <- function(spec, draws, group, lambda) {
  out <- draws
  if (length(spec$random$ar_terms) == 0L) {
    return(out)
  }
  ar_names <- .mlsim_random_effect_names("ar", spec$random$ar_terms, spec$outcomes)
  out[group, ar_names] <- lambda * out[group, ar_names]
  out
}

.mlsim_make_stability_metadata <- function(spec, phi, attempts, rule) {
  if (!isTRUE(spec$has_ar)) {
    return(list(stability_rule = rule, stability_attempts_by_group_level = NULL))
  }
  radii <- .mlsim_phi_radii(phi)
  radius_table <- data.table::data.table(
    group = spec$series$group_index,
    row = seq_along(spec$series$group_index),
    spectral_radius = radii
  )
  max_by_group <- as.numeric(tapply(radius_table$spectral_radius, radius_table$group, max))
  list(
    stability_rule = rule,
    stability_attempts_by_group_level = attempts,
    stability_acceptance_rate = if (is.null(attempts)) NULL else mean(1 / attempts),
    spectral_radius_by_group_level_and_row = radius_table,
    max_spectral_radius_by_group_level = max_by_group,
    max_spectral_radius_overall = max(radii, na.rm = TRUE)
  )
}

.mlsim_group_max_radius <- function(spec, phi, group) {
  idx <- which(spec$series$group_index == group)
  max(.mlsim_phi_radii(phi[idx, , , drop = FALSE]), na.rm = TRUE)
}

.mlsim_phi_radii <- function(phi) {
  if (is.null(phi) || length(phi) == 0L) {
    return(numeric())
  }
  vapply(seq_len(dim(phi)[[1L]]), function(i) {
    max(Mod(eigen(phi[i, , ], only.values = TRUE)$values))
  }, numeric(1))
}

.mlsim_outcome_location <- function(spec, checked, draws) {
  mu <- spec$X %*% checked$location$beta
  colnames(mu) <- spec$outcomes
  if (!is.null(draws) && length(spec$random$location_terms) > 0L) {
    for (term in spec$random$location_terms) {
      z <- spec$random$location_Z[, term]
      for (outcome in spec$outcomes) {
        name <- sprintf("location|outcome=%s|term=%s", outcome, term)
        mu[, outcome] <- mu[, outcome] + z * draws[spec$series$group_index, name]
      }
    }
  }
  mu
}

.mlsim_outcome_scale <- function(spec, checked, draws) {
  eta <- spec$W %*% checked$scale$beta
  colnames(eta) <- spec$outcomes
  if (!is.null(draws) && length(spec$random$scale_terms) > 0L) {
    for (term in spec$random$scale_terms) {
      z <- spec$random$scale_Z[, term]
      for (outcome in spec$outcomes) {
        name <- sprintf("scale|outcome=%s|term=%s", outcome, term)
        eta[, outcome] <- eta[, outcome] + z * draws[spec$series$group_index, name]
      }
    }
  }
  eta
}

.mlsim_outcome_phi <- function(spec, checked, draws) {
  n <- nrow(spec$X)
  k <- length(spec$outcomes)
  fixed_phi <- array(0, dim = c(n, k, k), dimnames = list(NULL, spec$outcomes, spec$outcomes))
  phi <- fixed_phi
  if (isTRUE(spec$has_ar)) {
    for (term in colnames(spec$Q)) {
      beta_term <- matrix(
        checked$ar$beta[term, , ],
        nrow = k,
        ncol = k,
        dimnames = list(spec$outcomes, spec$outcomes)
      )
      for (i in seq_len(n)) {
        fixed_phi[i, , ] <- fixed_phi[i, , ] + spec$Q[i, term] * beta_term
      }
    }
    phi <- fixed_phi
    if (!is.null(draws) && length(spec$random$ar_terms) > 0L) {
      for (term in spec$random$ar_terms) {
        z <- spec$random$ar_Z[, term]
        for (to in spec$outcomes) {
          for (from in spec$outcomes) {
            name <- sprintf("ar|term=%s|to=%s|from=%s", term, to, from)
            phi[, to, from] <- phi[, to, from] + z * draws[spec$series$group_index, name]
          }
        }
      }
    }
  }
  list(phi = phi, fixed_phi = fixed_phi)
}

.mlsim_draw_outcome_values <- function(spec, mu, sigma, correlation, phi, burnin) {
  n <- nrow(mu)
  k <- ncol(mu)
  z <- residual <- innovation <- matrix(NA_real_, n, k, dimnames = list(NULL, colnames(mu)))
  sigma_array <- .mlsim_mvn_row_covariance_array(sigma, correlation)
  for (rows in spec$series$rows) {
    e_prev <- rep(0, k)
    if (isTRUE(spec$has_ar) && burnin > 0L) {
      first <- rows[[1L]]
      for (b in seq_len(burnin)) {
        eps <- as.vector(.mlsim_rmvnorm(1L, rep(0, k), sigma_array[first, , ]))
        e_prev <- as.vector(phi$phi[first, , ] %*% e_prev + eps)
      }
    }
    for (r in rows) {
      eps <- as.vector(.mlsim_rmvnorm(1L, rep(0, k), sigma_array[r, , ]))
      e <- if (isTRUE(spec$has_ar)) {
        as.vector(phi$phi[r, , ] %*% e_prev + eps)
      } else {
        eps
      }
      innovation[r, ] <- eps
      residual[r, ] <- e
      z[r, ] <- mu[r, ] + e
      e_prev <- e
    }
  }
  list(z = z, residual = residual, innovation = innovation)
}
