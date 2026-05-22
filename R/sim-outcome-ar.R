# Outcome autoregressive simulation internals.

#' @rdname multilevelcoda-internal-outcome
.mlsim_simulate_ar_outcome <- function(setup, coefficients, random, residual, residual_info, ar) {
  if (!isTRUE(ar$linear)) {
    return(.mlsim_simulate_ar_outcome_rowwise(setup, coefficients, random, residual, residual_info))
  }
  n <- nrow(setup$data)
  d <- length(setup$response)
  values <- matrix(NA_real_, nrow = n, ncol = d, dimnames = list(NULL, setup$response))
  groups <- if (is.null(setup$context$n_groups)) rep(1L, n) else setup$data[[setup$context$group_id]]
  group_indices <- split(seq_len(n), groups)

  for (idx in group_indices) {
    prior <- .mlsim_burnin_prior_ar(setup, ar, residual_info = residual_info, first_rows = idx)
    for (row in idx) {
      if (!isTRUE(setup$complete_rows[[row]])) {
        next
      }
      transition <- .mlsim_ar_transition_row(ar, row, d, setup$response)
      current <- as.numeric(ar$base[row, ]) + as.numeric(prior %*% transition) + residual[row, ]
      values[row, ] <- current
      prior <- current
    }
  }
  values
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_simulate_ar_outcome_rowwise <- function(setup, coefficients, random, residual, residual_info) {
  n <- nrow(setup$data)
  d <- length(setup$response)
  values <- matrix(NA_real_, nrow = n, ncol = d, dimnames = list(NULL, setup$response))
  groups <- if (is.null(setup$context$n_groups)) rep(1L, n) else setup$data[[setup$context$group_id]]
  group_indices <- split(seq_len(n), groups)

  for (idx in group_indices) {
    prior <- .mlsim_burnin_prior(setup, coefficients, random, residual_info = residual_info, first_rows = idx)
    for (row in idx) {
      if (!isTRUE(setup$complete_rows[[row]])) {
        next
      }
      row_data <- setup$data[row, , drop = FALSE]
      for (j in seq_along(setup$ar_internal)) {
        row_data[[setup$ar_internal[[j]]]] <- prior[[j]]
      }
      x <- .mlsim_model_matrix_row(setup, row_data)
      random_row <- .mlsim_random_contribution_row(setup, random$effects, row, row_data)
      current <- as.numeric(x %*% coefficients) + random_row + residual[row, ]
      values[row, ] <- current
      prior <- current
    }
  }
  values
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_burnin_prior_ar <- function(setup, ar, residual_info, first_rows) {
  d <- length(setup$response)
  burnin <- setup$ar_info$burnin
  prior <- rep(0, d)
  if (burnin == 0L) {
    return(prior)
  }
  first_complete <- first_rows[setup$complete_rows[first_rows]][1L]
  if (is.na(first_complete)) {
    return(prior)
  }
  transition <- .mlsim_ar_transition_row(ar, first_complete, d, setup$response)
  for (i in seq_len(burnin)) {
    prior <- as.numeric(ar$base[first_complete, ]) + as.numeric(prior %*% transition) +
      .mlsim_burnin_residual_draw(residual_info, first_complete)
  }
  prior
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_burnin_prior <- function(setup, coefficients, random, residual_info, first_rows) {
  d <- length(setup$response)
  burnin <- setup$ar_info$burnin
  prior <- rep(0, d)
  if (burnin == 0L) {
    return(prior)
  }
  first_complete <- first_rows[setup$complete_rows[first_rows]][1L]
  if (is.na(first_complete)) {
    return(prior)
  }
  for (i in seq_len(burnin)) {
    row_data <- setup$data[first_complete, , drop = FALSE]
    for (j in seq_along(setup$ar_internal)) {
      row_data[[setup$ar_internal[[j]]]] <- prior[[j]]
    }
    x <- .mlsim_model_matrix_row(setup, row_data)
    random_row <- .mlsim_random_contribution_row(setup, random$effects, first_complete, row_data)
    prior <- as.numeric(x %*% coefficients) + random_row +
      .mlsim_burnin_residual_draw(residual_info, first_complete)
  }
  prior
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_ar_info <- function(setup) {
  if (length(setup$ar_internal) == 0L) {
    return(list(terms = character(), burnin = 0L, max_rho = 0))
  }
  burnin <- setup$burnin
  if (identical(burnin, "auto")) {
    burnin <- 50L
  }
  list(terms = setup$ar_names, burnin = burnin, max_rho = NA_real_)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_ar_components <- function(setup, coefficients, random) {
  n <- nrow(setup$data)
  d <- length(setup$response)
  zero_data <- .mlsim_ar_data(setup, 0)
  base <- setup$fixed_matrix %*% coefficients + random$contribution
  transition <- array(
    0,
    dim = c(n, d, d),
    dimnames = list(NULL, setup$response, setup$response)
  )

  for (j in seq_along(setup$ar_internal)) {
    one_data <- zero_data
    one_data[[setup$ar_internal[[j]]]] <- 1
    delta <- .mlsim_ar_eta(setup, coefficients, random, one_data) - base
    transition[, j, ] <- delta
  }

  list(
    base = base,
    transition = transition,
    linear = .mlsim_ar_components_are_linear(setup, coefficients, random, zero_data, base, transition)
  )
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_ar_eta <- function(setup, coefficients, random, data) {
  .mlsim_model_matrix_data(setup, data) %*% coefficients +
    .mlsim_random_contribution_data(setup, random$effects, data)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_ar_data <- function(setup, value) {
  data <- setup$data
  for (internal in setup$ar_internal) {
    data[[internal]] <- value
  }
  data
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_ar_components_are_linear <- function(setup, coefficients, random, zero_data, base, transition) {
  n <- nrow(setup$data)
  d <- length(setup$response)
  complete_rows <- setup$complete_rows
  if (!any(complete_rows)) {
    return(TRUE)
  }

  for (j in seq_along(setup$ar_internal)) {
    two_data <- zero_data
    two_data[[setup$ar_internal[[j]]]] <- 2
    expected <- base + 2 * .mlsim_ar_transition_delta(transition, j, n, d)
    actual <- .mlsim_ar_eta(setup, coefficients, random, two_data)
    if (!.mlsim_same_matrix(actual, expected, complete_rows)) {
      return(FALSE)
    }
  }

  if (length(setup$ar_internal) > 1L) {
    for (j in seq_len(length(setup$ar_internal) - 1L)) {
      for (k in seq.int(j + 1L, length(setup$ar_internal))) {
        pair_data <- zero_data
        pair_data[[setup$ar_internal[[j]]]] <- 1
        pair_data[[setup$ar_internal[[k]]]] <- 1
        expected <- base +
          .mlsim_ar_transition_delta(transition, j, n, d) +
          .mlsim_ar_transition_delta(transition, k, n, d)
        actual <- .mlsim_ar_eta(setup, coefficients, random, pair_data)
        if (!.mlsim_same_matrix(actual, expected, complete_rows)) {
          return(FALSE)
        }
      }
    }
  }
  TRUE
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_ar_transition_delta <- function(transition, j, n, d) {
  matrix(transition[, j, ], nrow = n, ncol = d)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_same_matrix <- function(actual, expected, rows, tolerance = 1e-8) {
  rows <- which(rows)
  if (length(rows) == 0L) {
    return(TRUE)
  }
  delta <- actual[rows, , drop = FALSE] - expected[rows, , drop = FALSE]
  all(is.finite(delta)) && max(abs(delta)) <= tolerance
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_validate_ar_stability <- function(setup, ar = NULL) {
  if (is.null(ar) || length(setup$ar_internal) == 0L) {
    return(invisible(TRUE))
  }
  d <- length(setup$response)
  for (row in which(setup$complete_rows)) {
    ar_matrix <- .mlsim_ar_transition_row(ar, row, d, setup$response)
    rho <- max(Mod(eigen(ar_matrix, only.values = TRUE)$values))
    if (!is.finite(rho) || rho >= 1) {
      .mlsim_stop("Every observed AR matrix must have spectral radius less than 1.")
    }
  }
  invisible(TRUE)
}


#' @rdname multilevelcoda-internal-outcome
.mlsim_ar_transition_row <- function(ar, row, d, response) {
  matrix(
    ar$transition[row, , ],
    nrow = d,
    ncol = d,
    dimnames = list(response, response)
  )
}


