# Outcome fitting-metadata and helper-column internals.

#' @rdname multilevelcoda-internal-outcome
.mlsim_outcome_template_metadata <- function(setup) {
  list(
    distribution = if (length(setup$response) == 1L) "gaussian" else "mvn",
    level = if (is.null(setup$context$n_groups)) "single" else "multilevel",
    vars = setup$output_names,
    response = setup$response,
    formula = setup$formula,
    brms_formula = setup$brms_formula,
    fit_formula = setup$fit_formula,
    fit_fixed_formula = setup$fit_fixed_formula,
    fit_term_map = setup$fit_term_map,
    fit_helpers = setup$fit_helpers,
    helper_columns = setup$helper_columns,
    coefficients = .mlsim_default_coefficients(setup),
    residual_cov = if (isTRUE(setup$scale$enabled)) NULL else .mlsim_default_residual_cov(setup),
    residual_cor = if (isTRUE(setup$scale$enabled)) .mlsim_default_residual_cor(setup) else NULL,
    scale_formula = setup$scale$formula,
    scale_fit_formula = setup$scale$fit_formula,
    scale_fit_term_map = setup$scale$fit_term_map,
    scale_helper_columns = setup$scale$helper_columns,
    scale_coefficients = .mlsim_default_scale_coefficients(setup),
    random_cov = .mlsim_default_random_cov(setup),
    ar_terms = setup$ar_info$terms,
    burnin = setup$ar_info$burnin,
    compositional = setup$composition$compositional,
    parts = setup$composition$parts,
    sbp = setup$composition$sbp,
    total = setup$composition$total,
    keep_ilr = setup$composition$keep_ilr
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_scale_fit_formula <- function(fixed_expr, random_terms, state, env) {
  fit_expr <- .mlsim_transform_for_fit(fixed_expr %||% 1, state)
  for (term in random_terms) {
    random_expr <- .mlsim_fit_random_term_expr(term, state)
    fit_expr <- .mlsim_plus_expr(fit_expr, random_expr)
  }
  stats::as.formula(
    as.call(list(as.name("~"), as.name("sigma"), fit_expr)),
    env = env
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_scale_fit_term_map <- function(fixed_names, fixed_expr, state, env) {
  fit_expr <- .mlsim_replace_symbols(fixed_expr %||% 1, state$fit_map)
  fit_formula <- .mlsim_rhs_formula(fit_expr, env)
  fit_data <- .mlsim_fit_data_from_state(state)
  fit_matrix <- .mlsim_model_matrix(fit_formula, data = fit_data)
  data.frame(
    component = "scale",
    sim_term = fixed_names,
    fit_term = colnames(fit_matrix),
    stringsAsFactors = FALSE
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_scale_helper_columns <- function(fit_formula, state) {
  helpers <- .mlsim_fit_helpers_dataframe(state)
  if (nrow(helpers) == 0L || is.null(fit_formula)) {
    return(character())
  }
  unique(all.vars(fit_formula) %sin% helpers$column)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_register_fit_helper <- function(state, type, source, internal, column,
                                       sim_term, response_derived) {
  state$fit_map[[internal]] <- column
  if (is.null(state$fit_helpers[[internal]])) {
    state$fit_helpers[[internal]] <- list(
      type = type,
      source = source,
      internal = internal,
      column = column,
      sim_term = sim_term,
      response_derived = isTRUE(response_derived)
    )
  }
  invisible(state)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_fit_column_name <- function(type, source) {
  switch(
    type,
    ar1 = paste0("lag_", source),
    lag1 = paste0("lag_", source),
    between = paste0(source, "_between"),
    within = paste0(source, "_within"),
    .mlsim_stop("No fit column naming rule is defined for `%s`.", type)
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_fit_helpers_dataframe <- function(state) {
  if (length(state$fit_helpers) == 0L) {
    return(data.frame(
      type = character(),
      source = character(),
      internal = character(),
      column = character(),
      sim_term = character(),
      response_derived = logical(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(
    rbind,
    lapply(state$fit_helpers, function(x) {
      data.frame(
        type = x$type,
        source = x$source,
        internal = x$internal,
        column = x$column,
        sim_term = x$sim_term,
        response_derived = x$response_derived,
        stringsAsFactors = FALSE
      )
    })
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_fit_data_from_state <- function(state) {
  data <- state$data
  helpers <- .mlsim_fit_helpers_dataframe(state)
  for (i in seq_len(nrow(helpers))) {
    data[[helpers$column[[i]]]] <- data[[helpers$internal[[i]]]]
  }
  data
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_fit_helper_values <- function(setup, values) {
  helpers <- setup$fit_helpers
  out <- data.frame(row.names = seq_len(nrow(setup$data)))
  if (nrow(helpers) == 0L) {
    return(out)
  }
  for (i in seq_len(nrow(helpers))) {
    helper <- helpers[i, ]
    out[[helper$column]] <- if (isTRUE(helper$response_derived)) {
      .mlsim_lag1(values[, helper$source], setup$context)
    } else {
      setup$data[[helper$internal]]
    }
  }
  out
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_brms_formula <- function(formula, fixed_expr, random_terms, state) {
  brms_expr <- .mlsim_transform_for_brms(fixed_expr %||% 1, state)
  for (term in random_terms) {
    random_expr <- .mlsim_random_term_expr(term, state)
    brms_expr <- .mlsim_plus_expr(brms_expr, random_expr)
  }
  stats::as.formula(
    as.call(list(as.name("~"), formula[[2L]], brms_expr)),
    env = environment(formula)
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_transform_for_brms <- function(expr, state) {
  transformed <- .mlsim_transform_outcome_expr(expr, state)$expr
  .mlsim_replace_symbols(transformed, state$brms_map)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_fit_formula <- function(formula, fixed_expr, random_terms, state) {
  fit_expr <- .mlsim_transform_for_fit(fixed_expr %||% 1, state)
  for (term in random_terms) {
    random_expr <- .mlsim_fit_random_term_expr(term, state)
    fit_expr <- .mlsim_plus_expr(fit_expr, random_expr)
  }
  stats::as.formula(
    as.call(list(as.name("~"), formula[[2L]], fit_expr)),
    env = environment(formula)
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_fit_fixed_formula <- function(formula, fixed_expr, state) {
  fit_expr <- .mlsim_transform_for_fit(fixed_expr %||% 1, state)
  stats::as.formula(
    as.call(list(as.name("~"), formula[[2L]], fit_expr)),
    env = environment(formula)
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_transform_for_fit <- function(expr, state) {
  transformed <- .mlsim_transform_outcome_expr(expr, state)$expr
  .mlsim_replace_symbols(transformed, state$fit_map)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_fit_random_term_expr <- function(term, state) {
  effects <- .mlsim_transform_for_fit(term$effects, state)
  if (is.null(term$id)) {
    as.call(list(as.name("("), as.call(list(as.name("|"), effects, as.name(term$group)))))
  } else {
    inner <- as.call(list(as.name("|"), effects, as.name(term$id)))
    as.call(list(as.name("("), as.call(list(as.name("|"), inner, as.name(term$group)))))
  }
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_fit_term_map <- function(fixed_names, fixed_expr, state, env) {
  fit_expr <- .mlsim_transform_for_fit(fixed_expr %||% 1, state)
  fit_formula <- .mlsim_rhs_formula(fit_expr, env)
  fit_data <- .mlsim_fit_data_from_state(state)
  fit_matrix <- .mlsim_model_matrix(fit_formula, data = fit_data)
  data.frame(
    component = "fixed",
    sim_term = fixed_names,
    fit_term = colnames(fit_matrix),
    stringsAsFactors = FALSE
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_random_term_expr <- function(term, state) {
  effects <- .mlsim_transform_for_brms(term$effects, state)
  if (is.null(term$id)) {
    as.call(list(as.name("("), as.call(list(as.name("|"), effects, as.name(term$group)))))
  } else {
    inner <- as.call(list(as.name("|"), effects, as.name(term$id)))
    as.call(list(as.name("("), as.call(list(as.name("|"), inner, as.name(term$group)))))
  }
}


