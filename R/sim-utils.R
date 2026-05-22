#' Internal Validation, Recycling, and Random-Parameter Helpers
#'
#' Validate user inputs, recycle scalar parameters, draw covariance-structured
#' random values, preserve random seeds, and assemble metadata shared by
#' generator implementations.
#'
#' @param level Simulation level.
#' @param ... Values passed to formatting or error helpers.
#' @param call. Logical flag passed to [stop()].
#' @param args,used,ignored,supplied,use Character vectors or named logicals
#'   used to build diagnostics.
#' @param generator Generator label used in diagnostic messages.
#' @param vars,n,names,arg,arg_name Variable names, target lengths, or argument
#'   labels.
#' @param allow_empty Logical flag allowing an empty name vector.
#' @param new,existing Candidate and existing generated names.
#' @param x,y Generic vectors, matrices, formulas, or fallback values.
#' @param categories,reference,labels Category values and reference labels.
#' @param ordered,output Categorical output controls.
#' @param logits,nonreference Baseline-category logit parameters.
#' @param random_cov,cov,residual_cor,residual_sd Covariance, correlation, or
#'   standard-deviation inputs.
#' @param mean Mean vector for multivariate normal draws.
#' @param prob,eta,reference_index,codes Category probabilities, linear
#'   predictors, and sampled category codes.
#' @param size Count-distribution size specification.
#' @param context Simulation context supplied by [simulate_data()].
#' @param values,result Generated values or raw custom-generator result.
#' @param scale Logical flag for joint location/scale random effects.
#' @param fixed_intercept,scale_fixed_intercept Location and scale intercepts.
#' @param expected_names Optional names used to align named vector and matrix
#'   inputs.
#' @param random Random-effect draw metadata.
#' @param scale_linear_predictor Row-level scale linear predictor.
#' @param metadata,wrapper_metadata Metadata lists.
#' @param allow_zero Logical flag allowing positive semi-definite covariance
#'   matrices.
#' @param seed Optional random seed.
#' @param expr Expression evaluated under a temporary seed.
#' @param spec Count-distribution specification list.
#' @param type Distribution type.
#'
#' @return Internal helper outputs used by generator implementations.
#'
#' @examples
#' multilevelcoda:::.mlsim_recycle(1, 3, "x")
#' multilevelcoda:::.mlsim_check_cov(diag(2), 2, "cov")
#' multilevelcoda:::.mlsim_mean_from_eta("poisson", log(2))
#'
#' @keywords internal
#' @importFrom extraoperators %a==% %ag% %age% %agele% %agl% %ale%
#' @importFrom extraoperators %any!in% %any==% %anyg% %anyin% %sin%
#' @name multilevelcoda-internal-utils
NULL

#' @rdname multilevelcoda-internal-utils
.mlsim_match_level <- function(level) {
  match.arg(level, c("single", "level2", "multilevel"))
}

#' @rdname multilevelcoda-internal-utils
.mlsim_stop <- function(..., call. = FALSE) {
  stop(sprintf(...), call. = call.)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_format_args <- function(args) {
  paste(sprintf("`%s`", args), collapse = ", ")
}

#' @rdname multilevelcoda-internal-utils
.mlsim_warn_ignored_args <- function(generator, used, ignored) {
  ignored <- ignored[!is.na(ignored) & ignored != ""]
  used <- used[!is.na(used) & used != ""]
  if (length(ignored) == 0L || length(used) == 0L) {
    return(invisible(NULL))
  }

  warning(
    sprintf(
      "`%s` generator will use %s and ignore %s.",
      generator,
      .mlsim_format_args(used),
      .mlsim_format_args(ignored)
    ),
    call. = FALSE
  )
  invisible(NULL)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_multilevel_unused <- function(level, supplied, generator, use) {
  if (!identical(level, "multilevel")) {
    return(invisible(NULL))
  }
  unused <- names(supplied)[as.logical(supplied)]
  if (length(unused) > 0L) {
    .mlsim_stop(
      "Multilevel `%s` generator does not accept %s; use %s.",
      generator,
      .mlsim_format_args(unused),
      use
    )
  }
  invisible(NULL)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_multilevel_required <- function(level, supplied, generator) {
  if (!identical(level, "multilevel")) {
    return(invisible(NULL))
  }
  missing <- names(supplied)[vapply(supplied, is.null, logical(1))]
  if (length(missing) > 0L) {
    .mlsim_stop(
      "Multilevel `%s` generator requires %s.",
      generator,
      .mlsim_format_args(missing)
    )
  }
  invisible(NULL)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_vars <- function(vars, n = NULL, arg = "vars", allow_empty = FALSE) {
  if (!is.character(vars) || anyNA(vars) || vars %any==% "") {
    .mlsim_stop("`%s` must be a character vector of non-empty names.", arg)
  }
  if (length(vars) == 0L && !isTRUE(allow_empty)) {
    .mlsim_stop("`%s` must contain at least one name.", arg)
  }
  if (anyDuplicated(vars)) {
    .mlsim_stop("`%s` must not contain duplicate names.", arg)
  }
  syntactic <- make.names(vars)
  if (anyDuplicated(syntactic)) {
    collision <- vars[duplicated(syntactic) | duplicated(syntactic, fromLast = TRUE)]
    .mlsim_stop(
      "`%s` must not contain names that collide after `make.names()`: %s.",
      arg,
      paste(unique(collision), collapse = ", ")
    )
  }
  if (!is.null(n) && length(vars) != n) {
    .mlsim_stop("`%s` must have length %d.", arg, n)
  }
  vars
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_new_names_compatible <- function(new, existing, arg = "generated names") {
  if (length(new) == 0L || length(existing) == 0L) {
    return(invisible(TRUE))
  }
  existing_syntactic <- make.names(existing)
  new_syntactic <- make.names(new)
  collision <- new[new_syntactic %in% existing_syntactic]
  if (length(collision) > 0L) {
    .mlsim_stop(
      "`%s` collide with existing names after `make.names()`: %s.",
      arg,
      paste(unique(collision), collapse = ", ")
    )
  }
  invisible(TRUE)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_finite_numeric <- function(x, arg) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    .mlsim_stop("`%s` must contain finite numeric values.", arg)
  }
  x
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_prob <- function(x, arg = "prob") {
  if (!is.numeric(x) || anyNA(x) || !(x %agele% c(0, 1))) {
    .mlsim_stop("`%s` must be numeric probabilities in [0, 1].", arg)
  }
  x
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_open_prob <- function(x, arg = "prob") {
  if (!is.numeric(x) || anyNA(x) || !(x %agl% c(0, 1))) {
    .mlsim_stop("`%s` must be numeric probabilities in (0, 1).", arg)
  }
  x
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_categories <- function(categories) {
  if (is.null(categories)) {
    categories <- c(0L, 1L)
  }
  if (is.matrix(categories) || !is.atomic(categories) || length(categories) < 2L || anyNA(categories)) {
    .mlsim_stop("`categories` must be a vector of at least two non-missing values.")
  }
  labels <- as.character(categories)
  if (anyNA(labels) || labels %any==% "") {
    .mlsim_stop("`categories` must have non-empty labels.")
  }
  if (anyDuplicated(labels)) {
    .mlsim_stop("`categories` must not contain duplicate labels.")
  }
  list(values = categories, labels = labels)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_reference_index <- function(reference, labels) {
  if (is.null(reference)) {
    return(1L)
  }
  reference <- as.character(reference)
  if (length(reference) != 1L || is.na(reference) || reference == "") {
    .mlsim_stop("`reference` must be one category value.")
  }
  index <- match(reference, labels)
  if (is.na(index)) {
    .mlsim_stop("`reference` must match one of `categories`.")
  }
  index
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_categorical_ordered <- function(ordered, output) {
  ordered <- .mlsim_check_scalar_logical(ordered, "ordered")
  if (isTRUE(ordered) && !identical(output, "factor")) {
    .mlsim_stop("`ordered = TRUE` requires `output = \"factor\"`.")
  }
  ordered
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_categorical_logits <- function(logits, nonreference, arg) {
  if (is.null(logits)) {
    return(NULL)
  }
  d <- length(nonreference)
  if (!is.numeric(logits) || anyNA(logits) || any(!is.finite(logits))) {
    .mlsim_stop("`%s` must contain finite numeric baseline-category logits.", arg)
  }
  if (!is.null(names(logits))) {
    if (names(logits) %any==% "" || !setequal(names(logits), nonreference)) {
      .mlsim_stop("`%s` names must match the non-reference category labels.", arg)
    }
    logits <- logits[nonreference]
  }
  logits <- .mlsim_recycle(logits, d, arg)
  names(logits) <- nonreference
  logits
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_categorical_random_cov <- function(random_cov, nonreference) {
  d <- length(nonreference)
  if (d > 1L && length(random_cov) == 1L && !is.matrix(random_cov) && !is.null(random_cov)) {
    .mlsim_stop("`random_cov` must be a %d by %d covariance matrix.", d, d)
  }
  .mlsim_check_named_square_cov(
    random_cov,
    nonreference,
    "random_cov",
    default = "zero",
    allow_scalar = FALSE
  )
}

#' @rdname multilevelcoda-internal-utils
.mlsim_uniform_category_prob <- function(n, labels) {
  prob <- matrix(1 / length(labels), nrow = n, ncol = length(labels))
  colnames(prob) <- labels
  prob
}

#' @rdname multilevelcoda-internal-utils
.mlsim_resolve_category_prob <- function(prob, n, labels, arg = "prob") {
  k <- length(labels)
  if (is.matrix(prob) || is.data.frame(prob)) {
    prob <- as.matrix(prob)
    if (ncol(prob) != k) {
      .mlsim_stop("`%s` must have one column per category.", arg)
    }
    prob_names <- colnames(prob)
    if (!is.null(prob_names)) {
      if (prob_names %any==% "" || anyDuplicated(prob_names) || !setequal(prob_names, labels)) {
        .mlsim_stop("`%s` column names must match the category labels.", arg)
      }
      prob <- prob[, labels, drop = FALSE]
    }
    if (nrow(prob) == 1L) {
      prob <- prob[rep(1L, n), , drop = FALSE]
    } else if (nrow(prob) != n) {
      .mlsim_stop("`%s` must have 1 or %d rows.", arg, n)
    }
    colnames(prob) <- labels
    return(.mlsim_check_category_prob_matrix(prob, arg))
  }

  if (!is.numeric(prob) || anyNA(prob) || any(!is.finite(prob))) {
    .mlsim_stop("`%s` must contain finite numeric probabilities.", arg)
  }

  if (!is.null(names(prob)) && length(prob) != 1L) {
    if (length(prob) != k || names(prob) %any==% "" || !setequal(names(prob), labels)) {
      .mlsim_stop("`%s` names must match the category labels.", arg)
    }
    prob <- prob[labels]
  }

  if (k == 2L && length(prob) == 1L) {
    success_prob <- .mlsim_check_prob(prob, arg)
    out <- cbind(1 - success_prob, success_prob)
    colnames(out) <- labels
    return(out[rep(1L, n), , drop = FALSE])
  }

  if (length(prob) != k) {
    .mlsim_stop(
      "`%s` must be a scalar binary success probability, a vector with one value per category, or a probability matrix.",
      arg
    )
  }
  out <- matrix(rep(prob, n), nrow = n, byrow = TRUE)
  colnames(out) <- labels
  .mlsim_check_category_prob_matrix(out, arg)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_category_prob_matrix <- function(prob, arg = "prob") {
  if (!is.numeric(prob) || anyNA(prob) || any(!is.finite(prob)) || any(prob < 0 | prob > 1)) {
    .mlsim_stop("`%s` must be numeric probabilities in [0, 1].", arg)
  }
  row_sums <- rowSums(prob)
  if (abs(row_sums - 1) %anyg% 1e-8) {
    .mlsim_stop("Rows of `%s` must sum to 1.", arg)
  }
  prob
}

#' @rdname multilevelcoda-internal-utils
.mlsim_baseline_category_prob <- function(eta, reference_index, labels) {
  eta <- as.matrix(eta)
  full_eta <- matrix(0, nrow = nrow(eta), ncol = length(labels))
  full_eta[, -reference_index] <- eta
  shifted <- sweep(full_eta, 1L, apply(full_eta, 1L, max), "-")
  exp_eta <- exp(shifted)
  prob <- exp_eta / rowSums(exp_eta)
  colnames(prob) <- labels
  prob
}

#' @rdname multilevelcoda-internal-utils
.mlsim_sample_categories <- function(prob) {
  if (nrow(prob) == 0L) {
    return(integer())
  }
  cumulative <- t(apply(prob, 1L, cumsum))
  rowSums(sweep(cumulative, 1L, stats::runif(nrow(prob)), "<")) + 1L
}

#' @rdname multilevelcoda-internal-utils
.mlsim_encode_categories <- function(codes, labels, output, ordered) {
  switch(
    output,
    factor = factor(labels[codes], levels = labels, ordered = ordered),
    character = labels[codes],
    integer = as.integer(codes - 1L),
    .mlsim_stop("Unknown categorical output `%s`.", output)
  )
}

#' @rdname multilevelcoda-internal-utils
.mlsim_category_code_map <- function(labels) {
  data.frame(code = seq_along(labels) - 1L, category = labels, stringsAsFactors = FALSE)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_positive <- function(x, arg) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || !(x %ag% 0)) {
    .mlsim_stop("`%s` must contain positive finite numbers.", arg)
  }
  x
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_nonnegative <- function(x, arg) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || !(x %age% 0)) {
    .mlsim_stop("`%s` must contain non-negative finite numbers.", arg)
  }
  x
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_scalar_number <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    .mlsim_stop("`%s` must be a scalar finite number.", arg)
  }
  x
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_scalar_logical <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .mlsim_stop("`%s` must be TRUE or FALSE.", arg)
  }
  x
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_whole_nonnegative <- function(x, arg) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || !(x %age% 0) ||
      !(abs(x - round(x)) %ale% sqrt(.Machine$double.eps))) {
    .mlsim_stop("`%s` must contain non-negative whole numbers.", arg)
  }
  x
}

#' @rdname multilevelcoda-internal-utils
.mlsim_recycle <- function(x, n, arg) {
  if (length(x) == 1L) {
    rep(x, n)
  } else if (length(x) == n) {
    x
  } else {
    .mlsim_stop("`%s` must have length 1 or %d.", arg, n)
  }
}

#' @rdname multilevelcoda-internal-utils
.mlsim_align_named_vector <- function(x, expected_names, arg) {
  supplied_names <- names(x)
  if (is.null(supplied_names)) {
    out <- .mlsim_recycle(x, length(expected_names), arg)
    names(out) <- expected_names
    return(out)
  }

  if (length(x) != length(expected_names) || anyNA(supplied_names) ||
      supplied_names %any==% "" || anyDuplicated(supplied_names) ||
      !setequal(supplied_names, expected_names)) {
    .mlsim_stop("`%s` names must match the expected names.", arg)
  }

  out <- x[expected_names]
  names(out) <- expected_names
  out
}

#' @rdname multilevelcoda-internal-utils
.mlsim_recycle_integer <- function(x, n, arg) {
  x <- .mlsim_recycle(x, n, arg)
  x <- .mlsim_check_whole_nonnegative(x, arg)
  as.integer(x)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_as_matrix <- function(x, n, arg) {
  if (is.null(x)) {
    return(diag(n))
  }
  if (length(x) == 1L && is.numeric(x) && !is.matrix(x)) {
    x <- diag(as.numeric(x), n)
  } else {
    x <- as.matrix(x)
  }
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || !(dim(x) %a==% c(n, n))) {
    .mlsim_stop("`%s` must be a finite numeric %d by %d matrix.", arg, n, n)
  }
  x
}

#' @rdname multilevelcoda-internal-utils
.mlsim_align_square_matrix_dimnames <- function(x, expected_names, arg) {
  row_names <- rownames(x)
  col_names <- colnames(x)

  if (is.null(row_names) && is.null(col_names)) {
    dimnames(x) <- list(expected_names, expected_names)
    return(x)
  }

  if (is.null(row_names) || is.null(col_names)) {
    .mlsim_stop("`%s` row and column names must either both be supplied or both be omitted.", arg)
  }
  if (anyNA(row_names) || anyNA(col_names) ||
      row_names %any==% "" || col_names %any==% "" ||
      anyDuplicated(row_names) || anyDuplicated(col_names) ||
      !setequal(row_names, expected_names) || !setequal(col_names, expected_names)) {
    .mlsim_stop("`%s` row and column names must match the expected names.", arg)
  }

  x <- x[expected_names, expected_names, drop = FALSE]
  dimnames(x) <- list(expected_names, expected_names)
  x
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_named_square_cov <- function(cov, expected_names, arg,
                                          default = c("identity", "zero"),
                                          allow_scalar = TRUE,
                                          allow_zero = TRUE) {
  default <- match.arg(default)
  d <- length(expected_names)
  if (is.null(cov)) {
    cov <- if (identical(default, "zero")) matrix(0, d, d) else diag(d)
  } else if (!isTRUE(allow_scalar) && d > 1L && length(cov) == 1L && !is.matrix(cov)) {
    .mlsim_stop("`%s` must be a %d by %d covariance matrix.", arg, d, d)
  }

  cov <- .mlsim_as_matrix(cov, d, arg)
  cov <- .mlsim_align_square_matrix_dimnames(cov, expected_names, arg)
  .mlsim_check_cov(cov, d, arg, allow_zero = allow_zero)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_cov <- function(x, n, arg, allow_zero = TRUE) {
  x <- .mlsim_as_matrix(x, n, arg)
  if (!isTRUE(all.equal(x, t(x), tolerance = sqrt(.Machine$double.eps)))) {
    .mlsim_stop("`%s` must be symmetric.", arg)
  }
  eig <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  threshold <- if (allow_zero) -1e-8 else 1e-8
  if (!(eig %age% threshold)) {
    .mlsim_stop("`%s` must be positive %s.", arg, if (allow_zero) "semi-definite" else "definite")
  }
  x
}

#' @rdname multilevelcoda-internal-utils
.mlsim_rmvnorm <- function(n, mean, cov) {
  d <- length(mean)
  if (n == 0L) {
    return(matrix(numeric(), nrow = 0L, ncol = d))
  }
  cov <- .mlsim_check_cov(cov, d, "cov")
  eig <- eigen(cov, symmetric = TRUE)
  transform <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = d)
  z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  sweep(z %*% t(transform), 2L, mean, "+")
}

#' @rdname multilevelcoda-internal-utils
.mlsim_with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(force(expr))
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(expr)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_n_for_level <- function(level, context) {
  if (identical(level, "single")) context$n_rows else context$n_groups
}

#' @rdname multilevelcoda-internal-utils
.mlsim_expand_level2 <- function(values, level, context) {
  if (identical(level, "level2")) {
    values[context$group_index, , drop = FALSE]
  } else {
    values
  }
}

#' @rdname multilevelcoda-internal-utils
.mlsim_expand_level2_vector <- function(values, level, context) {
  if (identical(level, "level2")) values[context$group_index] else values
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_multilevel_only_args <- function(level, supplied) {
  if (identical(level, "multilevel")) {
    return(invisible(NULL))
  }
  supplied <- names(supplied)[as.logical(supplied)]
  if (length(supplied) > 0L) {
    .mlsim_stop("%s require `level = \"multilevel\"`.", .mlsim_format_args(supplied))
  }
  invisible(NULL)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_random_intercept_names <- function(vars, scale = FALSE) {
  location <- paste0(vars, ":(Intercept)")
  if (!isTRUE(scale)) {
    return(location)
  }
  c(location, paste0("scale_", vars, ":(Intercept)"))
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_position_cov <- function(cov, names, arg = "random_cov") {
  d <- length(names)
  if (d > 1L && length(cov) == 1L && !is.matrix(cov) && !is.null(cov)) {
    .mlsim_stop("`%s` must be a %d by %d covariance matrix.", arg, d, d)
  }
  .mlsim_check_named_square_cov(
    cov,
    names,
    arg,
    default = "zero",
    allow_scalar = FALSE
  )
}

#' @rdname multilevelcoda-internal-utils
.mlsim_joint_random_intercepts <- function(vars, random_cov, context, scale = FALSE) {
  names <- .mlsim_random_intercept_names(vars, scale)
  random_cov <- .mlsim_check_position_cov(random_cov, names, "random_cov")
  draws <- .mlsim_rmvnorm(context$n_groups, rep(0, length(names)), random_cov)
  colnames(draws) <- names

  n_vars <- length(vars)
  random_intercepts <- draws[, seq_len(n_vars), drop = FALSE]
  colnames(random_intercepts) <- vars

  scale_random_intercepts <- NULL
  if (isTRUE(scale)) {
    scale_random_intercepts <- draws[, n_vars + seq_len(n_vars), drop = FALSE]
    colnames(scale_random_intercepts) <- vars
  }

  list(
    random_cov = random_cov,
    random_intercepts = random_intercepts,
    scale_random_intercepts = scale_random_intercepts,
    joint_random_intercepts = draws
  )
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_scale_fixed_intercept <- function(scale_fixed_intercept, n = 1L, expected_names = NULL) {
  if (is.null(scale_fixed_intercept)) {
    return(NULL)
  }
  scale_fixed_intercept <- if (is.null(expected_names)) {
    .mlsim_recycle(scale_fixed_intercept, n, "scale_fixed_intercept")
  } else {
    .mlsim_align_named_vector(scale_fixed_intercept, expected_names, "scale_fixed_intercept")
  }
  if (!is.numeric(scale_fixed_intercept) || anyNA(scale_fixed_intercept) || any(!is.finite(scale_fixed_intercept))) {
    .mlsim_stop("`scale_fixed_intercept` must contain finite numbers.")
  }
  scale_fixed_intercept
}

#' @rdname multilevelcoda-internal-utils
.mlsim_random_eta <- function(type, fixed_intercept, random_cov, context, vars,
                              scale_fixed_intercept = NULL) {
  if (is.null(fixed_intercept)) {
    .mlsim_stop("Multilevel `%s` generator requires `fixed_intercept`.", type)
  }
  fixed_intercept <- .mlsim_check_scalar_number(fixed_intercept, "fixed_intercept")
  scale_fixed_intercept <- .mlsim_check_scale_fixed_intercept(scale_fixed_intercept, 1L)

  random <- .mlsim_joint_random_intercepts(
    vars,
    random_cov,
    context,
    scale = !is.null(scale_fixed_intercept)
  )

  eta <- fixed_intercept + random$random_intercepts[context$group_index, 1L]
  metadata <- list(
    fixed_parameters = list(intercept = fixed_intercept),
    group_parameters = list(
      random_cov = random$random_cov,
      random_intercepts = random$random_intercepts
    ),
    fixed_intercept = fixed_intercept,
    random_cov = random$random_cov,
    random_intercepts = random$random_intercepts
  )

  scale_eta <- NULL
  if (!is.null(scale_fixed_intercept)) {
    scale_eta <- scale_fixed_intercept + random$scale_random_intercepts[context$group_index, 1L]
    metadata$fixed_parameters$scale_intercept <- scale_fixed_intercept
    metadata$group_parameters$scale_random_intercepts <- random$scale_random_intercepts
    metadata$scale_fixed_intercept <- scale_fixed_intercept
    metadata$scale_random_intercepts <- random$scale_random_intercepts
    metadata$scale_linear_predictor <- scale_eta
  }

  list(
    eta = eta,
    scale_eta = scale_eta,
    random = random,
    metadata = metadata
  )
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_simple_residual_cor <- function(residual_cor, vars) {
  d <- length(vars)
  if (d > 1L && length(residual_cor) == 1L && !is.matrix(residual_cor) && !is.null(residual_cor)) {
    .mlsim_stop("`residual_cor` must be a %d by %d correlation matrix.", d, d)
  }
  residual_cor <- .mlsim_check_named_square_cov(
    residual_cor,
    vars,
    "residual_cor",
    default = "identity",
    allow_scalar = FALSE
  )
  if (any(abs(diag(residual_cor) - 1) > 1e-8)) {
    .mlsim_stop("`residual_cor` must have ones on the diagonal.")
  }
  if (any(residual_cor < -1 - 1e-8 | residual_cor > 1 + 1e-8)) {
    .mlsim_stop("`residual_cor` values must be correlations in [-1, 1].")
  }
  residual_cor
}

#' @rdname multilevelcoda-internal-utils
.mlsim_mvn_row_residual <- function(residual_sd, residual_cor) {
  residual_sd <- as.matrix(residual_sd)
  base <- .mlsim_rmvnorm(nrow(residual_sd), rep(0, ncol(residual_sd)), residual_cor)
  base * residual_sd
}

#' @rdname multilevelcoda-internal-utils
.mlsim_mvn_row_covariance_array <- function(residual_sd, residual_cor) {
  d <- ncol(residual_sd)
  out <- array(
    NA_real_,
    dim = c(nrow(residual_sd), d, d),
    dimnames = list(NULL, colnames(residual_sd), colnames(residual_sd))
  )
  for (i in seq_len(nrow(residual_sd))) {
    out[i, , ] <- diag(residual_sd[i, ], d) %*% residual_cor %*% diag(residual_sd[i, ], d)
  }
  out
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_scale_or_arg <- function(arg, arg_name, scale_fixed_intercept, generator) {
  if (!is.null(arg) && !is.null(scale_fixed_intercept)) {
    .mlsim_stop("`%s` cannot be supplied with `scale_fixed_intercept`.", arg_name)
  }
  if (is.null(arg) && is.null(scale_fixed_intercept)) {
    .mlsim_stop("Multilevel `%s` generator requires `%s` or `scale_fixed_intercept`.", generator, arg_name)
  }
  invisible(NULL)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_check_positive_or_null <- function(x, arg) {
  if (is.null(x)) {
    return(NULL)
  }
  .mlsim_check_positive(x, arg)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_simple_random_metadata <- function(fixed_intercept, random,
                                          scale_fixed_intercept = NULL,
                                          scale_linear_predictor = NULL) {
  metadata <- list(
    fixed_parameters = list(intercept = fixed_intercept),
    group_parameters = list(
      random_cov = random$random_cov,
      random_intercepts = random$random_intercepts
    ),
    fixed_intercept = fixed_intercept,
    random_cov = random$random_cov,
    random_intercepts = random$random_intercepts
  )
  if (!is.null(scale_fixed_intercept)) {
    metadata$fixed_parameters$scale_intercept <- scale_fixed_intercept
    metadata$group_parameters$scale_random_intercepts <- random$scale_random_intercepts
    metadata$scale_fixed_intercept <- scale_fixed_intercept
    metadata$scale_random_intercepts <- random$scale_random_intercepts
    metadata$scale_linear_predictor <- scale_linear_predictor
  }
  metadata
}

#' @rdname multilevelcoda-internal-utils
.mlsim_validate_scale_multilevel_args <- function(level, scale_fixed_intercept,
                                                  random_cov = NULL,
                                                  residual_cor = NULL) {
  .mlsim_check_multilevel_only_args(
    level,
    c(
      random_cov = !is.null(random_cov),
      scale_fixed_intercept = !is.null(scale_fixed_intercept),
      residual_cor = !is.null(residual_cor)
    )
  )
}

#' @rdname multilevelcoda-internal-utils
.mlsim_mean_from_eta <- function(type, eta) {
  switch(
    type,
    bernoulli = stats::plogis(eta),
    binomial = stats::plogis(eta),
    poisson = exp(eta),
    negbin = exp(eta),
    gamma = exp(eta),
    beta = stats::plogis(eta),
    .mlsim_stop("No mean function is defined for `%s`.", type)
  )
}

#' @rdname multilevelcoda-internal-utils
.mlsim_parameter_metadata <- function(level, context, ...) {
  parameter_level <- if (identical(level, "level2")) "group" else "row"
  parameter_count <- if (identical(parameter_level, "group")) context$n_groups else context$n_rows
  parameters_name <- paste0(parameter_level, "_parameters")
  parameters <- setNames(list(list(...)), parameters_name)
  c(parameters, list(parameter_level = parameter_level, parameter_count = parameter_count))
}

#' @rdname multilevelcoda-internal-utils
.mlsim_resolve_size <- function(size, n, context) {
  if (is.function(size)) {
    size <- size(n)
  } else if (is.list(size)) {
    size <- .mlsim_draw_count_distribution(size, n)
  }
  .mlsim_recycle_integer(size, n, context)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_result <- function(values, names, metadata) {
  data <- data.table::as.data.table(values)
  names <- .mlsim_check_vars(names, ncol(data), "generated names", allow_empty = TRUE)
  data.table::setnames(data, names)
  list(data = data, names = names, metadata = metadata)
}

#' @rdname multilevelcoda-internal-utils
.mlsim_custom_result <- function(result, vars, level, context, wrapper_metadata) {
  if (!is.list(result) || is.null(result$data) || is.null(result$names)) {
    .mlsim_stop("Custom generators must return a list with `data` and `names` elements.")
  }
  names <- .mlsim_check_vars(result$names, length(vars), "custom result names")
  if (!identical(names, vars)) {
    .mlsim_stop("Custom result `names` must exactly match `vars`.")
  }

  metadata <- result$metadata %||% list()
  if (!is.list(metadata)) {
    .mlsim_stop("Custom result `metadata` must be a list when supplied.")
  }

  data <- data.table::as.data.table(result$data)
  if (ncol(data) != length(vars)) {
    .mlsim_stop("Custom result `data` must have %d columns.", length(vars))
  }

  output_level <- "row"
  if (identical(level, "level2") && nrow(data) == context$n_groups) {
    group_keys <- data.table::copy(context$group_data)
    row_keys <- data.table::copy(context$design_data[, context$group_id, with = FALSE])
    value_names <- paste0("__mlsim_custom_", seq_along(vars))
    data.table::setnames(data, value_names)

    group_values <- cbind(group_keys, data)
    data.table::setkeyv(group_values, context$group_id)
    data.table::setkeyv(row_keys, context$group_id)
    data <- group_values[row_keys]
    data[, (context$group_id) := NULL]
    data.table::setnames(data, value_names, vars)
    output_level <- "group"
  } else if (nrow(data) != context$n_rows) {
    expected <- if (identical(level, "level2")) {
      sprintf("%d rows or %d rows", context$n_groups, context$n_rows)
    } else {
      sprintf("%d rows", context$n_rows)
    }
    .mlsim_stop("Custom result `data` must have %s for level `%s`.", expected, level)
  }

  .mlsim_result(
    data,
    vars,
    c(
      list(
        distribution = "custom",
        level = level,
        vars = vars,
        custom_output_level = output_level,
        wrapper_metadata = wrapper_metadata
      ),
      metadata
    )
  )
}

#' @rdname multilevelcoda-internal-utils
.mlsim_draw_count_distribution <- function(spec, n) {
  distribution <- spec$distribution %||% spec$dist
  if (is.null(distribution) || !is.character(distribution) || length(distribution) != 1L) {
    .mlsim_stop("Count distribution specs need a scalar `distribution` or `dist` name.")
  }

  distribution <- match.arg(distribution, c("constant", "uniform", "poisson", "negative_binomial"))
  values <- switch(
    distribution,
    constant = {
      value <- spec$value %||% spec$n
      if (is.null(value)) {
        .mlsim_stop("Constant count specs need `value`.")
      }
      .mlsim_recycle_integer(value, n, "value")
    },
    uniform = {
      min_value <- spec$min %||% spec$minimum
      max_value <- spec$max %||% spec$maximum
      if (is.null(min_value) || is.null(max_value)) {
        .mlsim_stop("Uniform count specs need `min`/`max` or `minimum`/`maximum`.")
      }
      min_value <- .mlsim_recycle_integer(min_value, 1L, "min")
      max_value <- .mlsim_recycle_integer(max_value, 1L, "max")
      if (min_value > max_value) {
        .mlsim_stop("Uniform count specs require `min` to be less than or equal to `max`.")
      }
      sample(seq.int(as.integer(min_value), as.integer(max_value)), n, replace = TRUE)
    },
    poisson = {
      lambda <- spec$lambda
      if (is.null(lambda)) {
        .mlsim_stop("Poisson count specs need `lambda`.")
      }
      lambda <- .mlsim_check_nonnegative(.mlsim_recycle(lambda, n, "lambda"), "lambda")
      stats::rpois(n, lambda)
    },
    negative_binomial = {
      mu <- spec$mu
      size <- spec$size
      if (is.null(mu) || is.null(size)) {
        .mlsim_stop("Negative-binomial count specs need `mu` and `size`.")
      }
      mu <- .mlsim_check_nonnegative(.mlsim_recycle(mu, n, "mu"), "mu")
      size <- .mlsim_check_positive(.mlsim_recycle(size, n, "size"), "size")
      stats::rnbinom(n, mu = mu, size = size)
    }
  )

  minimum <- spec$minimum %||% spec$min_count
  maximum <- spec$maximum %||% spec$max_count
  if (!is.null(minimum)) {
    minimum <- .mlsim_recycle_integer(minimum, 1L, "minimum")
    values <- pmax(values, minimum)
  }
  if (!is.null(maximum)) {
    maximum <- .mlsim_recycle_integer(maximum, 1L, "maximum")
    values <- pmin(values, maximum)
  }
  if (!is.null(minimum) && !is.null(maximum) && minimum > maximum) {
    .mlsim_stop("Count distribution bounds require `minimum` to be less than or equal to `maximum`.")
  }
  .mlsim_recycle_integer(values, n, "generated counts")
}

#' @rdname multilevelcoda-internal-utils
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
