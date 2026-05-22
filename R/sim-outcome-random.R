# Outcome random-effect specification and contribution internals.

#' @rdname multilevelcoda-internal-outcome
.mlsim_extract_random_terms <- function(expr) {
  out <- .mlsim_extract_random_terms_rec(expr)
  out$random <- out$random %||% list()
  out
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_extract_random_terms_rec <- function(expr) {
  if (!is.call(expr)) {
    return(list(fixed_expr = expr, random = list()))
  }
  if (identical(.mlsim_call_name(expr), "(")) {
    inner <- expr[[2L]]
    if (.mlsim_is_random_call(inner)) {
      return(list(fixed_expr = NULL, random = list(.mlsim_parse_random_term(inner))))
    }
    out <- .mlsim_extract_random_terms_rec(inner)
    if (is.null(out$fixed_expr)) {
      return(out)
    }
    out$fixed_expr <- as.call(list(as.name("("), out$fixed_expr))
    return(out)
  }
  if (.mlsim_is_random_call(expr)) {
    return(list(fixed_expr = NULL, random = list(.mlsim_parse_random_term(expr))))
  }
  if (identical(.mlsim_call_name(expr), "+")) {
    left <- .mlsim_extract_random_terms_rec(expr[[2L]])
    right <- .mlsim_extract_random_terms_rec(expr[[3L]])
    return(list(
      fixed_expr = .mlsim_plus_expr(left$fixed_expr, right$fixed_expr),
      random = c(left$random, right$random)
    ))
  }
  list(fixed_expr = expr, random = list())
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_parse_random_term <- function(expr) {
  expr <- .mlsim_unwrap_parens(expr)
  if (!.mlsim_is_random_call(expr)) {
    .mlsim_stop("Invalid random-effect term.")
  }
  lhs <- .mlsim_unwrap_parens(expr[[2L]])
  group <- .mlsim_deparse_one(expr[[3L]])
  id <- NULL
  effects <- lhs
  if (.mlsim_is_random_call(lhs)) {
    effects <- .mlsim_unwrap_parens(lhs[[2L]])
    id <- .mlsim_deparse_one(lhs[[3L]])
  }
  list(effects = effects, group = group, id = id)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_random_specs <- function(random_terms, state, env, response) {
  if (length(random_terms) == 0L) {
    return(list())
  }
  specs <- list()
  for (term in random_terms) {
    if (term$group %any!in% names(state$data)) {
      .mlsim_stop("Random-effect grouping variable `%s` is not available.", term$group)
    }
    term_state <- state
    effects <- .mlsim_transform_outcome_expr(term$effects, term_state)$expr
    z_formula <- .mlsim_rhs_formula(effects, env)
    z <- .mlsim_model_matrix(z_formula, data = term_state$data)
    z_xlevels <- .mlsim_model_xlevels(z_formula, term_state$data)
    z_contrasts <- attr(z, "contrasts")
    z_names <- .mlsim_pretty_design_names(colnames(z), term_state)
    if (anyDuplicated(z_names)) {
      .mlsim_stop("Random-effect term names must be unique within a grouping block.")
    }
    group_values <- state$data[[term$group]]
    group_levels <- unique(group_values)

    if (is.null(term$id) && length(response) > 1L) {
      for (resp in response) {
        key <- paste(term$group, resp, sep = "::")
        specs[[length(specs) + 1L]] <- list(
          key = key, group = term$group, id = NULL, response = resp,
          responses = resp, z = z, z_formula = z_formula,
          z_xlevels = z_xlevels, z_contrasts = z_contrasts,
          z_names = z_names,
          group_values = group_values, group_levels = group_levels,
          coef_names = paste(resp, z_names, sep = ":")
        )
      }
    } else {
      key_id <- term$id %||% "default"
      key <- paste(term$group, key_id, sep = "::")
      specs[[length(specs) + 1L]] <- list(
        key = key, group = term$group, id = term$id, response = NULL,
        responses = response, z = z, z_formula = z_formula,
        z_xlevels = z_xlevels, z_contrasts = z_contrasts,
        z_names = z_names,
        group_values = group_values, group_levels = group_levels,
        coef_names = unlist(lapply(response, function(resp) paste(resp, z_names, sep = ":")), use.names = FALSE)
      )
    }
  }
  .mlsim_merge_random_specs(specs)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_merge_random_specs <- function(specs) {
  merged <- list()
  for (spec in specs) {
    if (is.null(merged[[spec$key]])) {
      merged[[spec$key]] <- spec
    } else {
      current <- merged[[spec$key]]
      if (!identical(current$responses, spec$responses) ||
          !identical(current$group_values, spec$group_values)) {
        .mlsim_stop("Random-effect terms sharing an ID must use the same grouping structure.")
      }
      current$z <- cbind(current$z, spec$z)
      current$z_formula <- .mlsim_plus_formula(current$z_formula, spec$z_formula)
      current$z_xlevels <- c(current$z_xlevels, spec$z_xlevels)
      current$z_contrasts <- c(current$z_contrasts, spec$z_contrasts)
      current$z_names <- c(current$z_names, spec$z_names)
      if (anyDuplicated(current$z_names)) {
        .mlsim_stop("Random-effect terms sharing an ID must have unique names.")
      }
      current$coef_names <- unlist(lapply(current$responses, function(resp) {
        paste(resp, current$z_names, sep = ":")
      }), use.names = FALSE)
      merged[[spec$key]] <- current
    }
  }
  merged
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_scale_random_specs <- function(random_terms, state, env, response) {
  if (length(random_terms) == 0L) {
    return(list())
  }
  specs <- list()
  for (term in random_terms) {
    if (term$group %any!in% names(state$data)) {
      .mlsim_stop("Scale random-effect grouping variable `%s` is not available.", term$group)
    }
    effects <- .mlsim_transform_outcome_expr(term$effects, state, allow_ar1 = FALSE)$expr
    z_formula <- .mlsim_rhs_formula(effects, env)
    z <- .mlsim_model_matrix(z_formula, data = state$data)
    z_xlevels <- .mlsim_model_xlevels(z_formula, state$data)
    z_contrasts <- attr(z, "contrasts")
    z_names <- .mlsim_pretty_design_names(colnames(z), state)
    if (!identical(z_names, "(Intercept)")) {
      .mlsim_stop("Scale random effects currently support intercept-only terms.")
    }
    group_values <- state$data[[term$group]]
    group_levels <- unique(group_values)

    if (is.null(term$id) && length(response) > 1L) {
      for (resp in response) {
        key <- paste(term$group, paste0("sigma_", resp), sep = "::")
        specs[[length(specs) + 1L]] <- list(
          key = key, group = term$group, id = NULL, response = resp,
          responses = resp, z = z, z_formula = z_formula,
          z_xlevels = z_xlevels, z_contrasts = z_contrasts,
          z_names = z_names,
          group_values = group_values, group_levels = group_levels,
          coef_names = paste0("sigma_", resp, ":(Intercept)")
        )
      }
    } else {
      key_id <- term$id %||% "default"
      key <- paste(term$group, key_id, sep = "::")
      specs[[length(specs) + 1L]] <- list(
        key = key, group = term$group, id = term$id, response = NULL,
        responses = response, z = z, z_formula = z_formula,
        z_xlevels = z_xlevels, z_contrasts = z_contrasts,
        z_names = z_names,
        group_values = group_values, group_levels = group_levels,
        coef_names = paste0("sigma_", response, ":(Intercept)")
      )
    }
  }
  .mlsim_merge_random_specs(specs)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_joint_random_specs <- function(location_specs, scale_specs) {
  keys <- unique(c(names(location_specs), names(scale_specs)))
  out <- list()
  for (key in keys) {
    location <- location_specs[[key]]
    scale <- scale_specs[[key]]
    spec <- location %||% scale
    if (!is.null(location) && !is.null(scale)) {
      if (!identical(location$group, scale$group) ||
          !identical(location$group_values, scale$group_values)) {
        .mlsim_stop("Location and scale random effects sharing an ID must use the same grouping structure.")
      }
    }
    out[[key]] <- list(
      key = key,
      group = spec$group,
      group_values = spec$group_values,
      group_levels = spec$group_levels,
      coef_names = c(location$coef_names %||% character(), scale$coef_names %||% character())
    )
  }
  out
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_random_contribution <- function(setup, random_cov) {
  n <- nrow(setup$data)
  d <- length(setup$response)
  contribution <- matrix(0, nrow = n, ncol = d, dimnames = list(NULL, setup$response))
  scale_contribution <- matrix(0, nrow = n, ncol = d, dimnames = list(NULL, setup$response))
  effects <- list()
  if (length(setup$joint_random_specs) == 0L) {
    return(list(contribution = contribution, scale_contribution = scale_contribution, effects = effects))
  }

  for (key in names(setup$joint_random_specs)) {
    spec <- setup$joint_random_specs[[key]]
    cov <- random_cov[[key]]
    draws <- .mlsim_rmvnorm(length(spec$group_levels), rep(0, nrow(cov)), cov)
    colnames(draws) <- spec$coef_names
    rownames(draws) <- as.character(spec$group_levels)
    effects[[key]] <- draws
  }
  contribution <- .mlsim_random_contribution_from_z(
    setup,
    effects,
    lapply(setup$random_specs, function(spec) spec$z)
  )
  scale_contribution <- .mlsim_scale_random_contribution(setup, effects)
  list(contribution = contribution, scale_contribution = scale_contribution, effects = effects)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_random_contribution_data <- function(setup, effects, data) {
  .mlsim_random_contribution_from_z(
    setup,
    effects,
    lapply(setup$random_specs, function(spec) .mlsim_random_model_matrix(spec, data))
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_random_contribution_from_z <- function(setup, effects, z_by_key) {
  n <- nrow(setup$data)
  d <- length(setup$response)
  contribution <- matrix(0, nrow = n, ncol = d, dimnames = list(NULL, setup$response))
  if (length(setup$random_specs) == 0L) {
    return(contribution)
  }

  for (key in names(setup$random_specs)) {
    spec <- setup$random_specs[[key]]
    draws <- effects[[key]]
    z <- z_by_key[[key]]
    group_match <- match(spec$group_values, spec$group_levels)
    draw_rows <- draws[group_match, spec$coef_names, drop = FALSE]

    if (length(spec$responses) == 1L) {
      resp_index <- match(spec$responses, setup$response)
      contribution[, resp_index] <- contribution[, resp_index] +
        rowSums(z * draw_rows)
    } else {
      q <- length(spec$z_names)
      for (resp_index in seq_along(setup$response)) {
        draw_cols <- (resp_index - 1L) * q + seq_len(q)
        contribution[, resp_index] <- contribution[, resp_index] +
          rowSums(z * draw_rows[, draw_cols, drop = FALSE])
      }
    }
  }
  contribution
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_random_contribution_row <- function(setup, effects, row, row_data) {
  contribution <- rep(0, length(setup$response))
  names(contribution) <- setup$response
  if (length(setup$random_specs) == 0L) {
    return(contribution)
  }

  for (key in names(setup$random_specs)) {
    spec <- setup$random_specs[[key]]
    draws <- effects[[key]]
    group <- as.character(spec$group_values[[row]])
    z <- .mlsim_random_model_matrix_row(spec, row_data)
    draw <- draws[group, spec$coef_names, drop = TRUE]
    if (length(spec$responses) == 1L) {
      resp_index <- match(spec$responses, setup$response)
      contribution[[resp_index]] <- contribution[[resp_index]] + sum(z * draw)
    } else {
      q <- length(spec$z_names)
      b <- matrix(draw, nrow = q, ncol = length(setup$response))
      contribution <- contribution + as.numeric(z %*% b)
    }
  }
  contribution
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_scale_random_contribution <- function(setup, effects) {
  n <- nrow(setup$data)
  d <- length(setup$response)
  contribution <- matrix(0, nrow = n, ncol = d, dimnames = list(NULL, setup$response))
  if (!isTRUE(setup$scale$enabled) || length(setup$scale$random_specs) == 0L) {
    return(contribution)
  }

  for (key in names(setup$scale$random_specs)) {
    spec <- setup$scale$random_specs[[key]]
    draws <- effects[[key]]
    z <- spec$z
    group_match <- match(spec$group_values, spec$group_levels)
    draw_rows <- draws[group_match, spec$coef_names, drop = FALSE]

    if (length(spec$responses) == 1L) {
      resp_index <- match(spec$responses, setup$response)
      contribution[, resp_index] <- contribution[, resp_index] +
        rowSums(z * draw_rows)
    } else {
      q <- length(spec$z_names)
      for (resp_index in seq_along(setup$response)) {
        draw_cols <- (resp_index - 1L) * q + seq_len(q)
        contribution[, resp_index] <- contribution[, resp_index] +
          rowSums(z * draw_rows[, draw_cols, drop = FALSE])
      }
    }
  }
  contribution
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_random_model_matrix_row <- function(spec, row_data) {
  .mlsim_random_model_matrix(spec, row_data)
}


