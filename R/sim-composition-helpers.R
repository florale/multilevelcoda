#' Internal Composition Helpers
#'
#' Validate sequential binary partitions, resolve composition part names, and
#' transform ILR coordinates back to closed composition parts.
#'
#' @param ilr Matrix or data frame of ILR coordinates.
#' @param parts Character vector naming composition parts.
#' @param total Positive scalar total used to close compositions.
#' @param coordinate_names Character vector naming ILR coordinates.
#' @param sbp Sequential binary partition matrix with part names in columns.
#' @param d Number of ILR coordinates.
#' @param keep_ilr Logical flag indicating whether ILR coordinates are emitted
#'   alongside composition parts.
#'
#' @return Internal composition metadata or transformed values.
#'
#' @examples
#' resolved <- multilevelcoda:::.mlsim_resolve_composition(
#'   parts = c("a", "b", "c"),
#'   d = 2
#' )
#' resolved$parts
#'
#' inv <- multilevelcoda:::.mlsim_ilr_inverse(
#'   matrix(c(0, 0, 0.1, -0.1), nrow = 2, byrow = TRUE),
#'   parts = c("a", "b", "c")
#' )
#' rowSums(inv$values)
#'
#' @keywords internal
#' @name multilevelcoda-internal-composition
NULL

#' @rdname multilevelcoda-internal-composition
.mlsim_ilr_inverse <- function(ilr, parts = NULL, total = 1, coordinate_names = NULL, sbp = NULL) {
  ilr <- as.matrix(ilr)
  composition <- .mlsim_resolve_composition(parts, sbp, ncol(ilr))
  parts <- composition$parts
  coordinate_names <- coordinate_names %||% colnames(ilr) %||% paste0("ilr_", seq_len(ncol(ilr)))
  coordinate_names <- .mlsim_check_vars(coordinate_names, ncol(ilr), "coordinate_names")
  total <- .mlsim_check_positive(total, "total")
  if (length(total) != 1L) {
    .mlsim_stop("`total` must be a scalar positive value.")
  }

  sbp <- composition$sbp
  basis <- compositions::gsi.buildilrBase(t(sbp))
  values <- as.matrix(compositions::clo(compositions::ilrInv(ilr, V = basis), total = total))
  expected_dim <- c(nrow(ilr), length(parts))
  if (!identical(dim(values), expected_dim)) {
    transposed <- t(values)
    if (identical(dim(transposed), expected_dim)) {
      values <- transposed
    } else {
      .mlsim_stop("Internal ILR inverse produced an unexpected composition shape.")
    }
  }
  colnames(values) <- parts

  list(
    values = values,
    sbp = sbp,
    basis = basis,
    coordinate_map = .mlsim_ilr_coordinate_map(sbp, coordinate_names),
    total = total
  )
}

#' @rdname multilevelcoda-internal-composition
.mlsim_resolve_composition <- function(parts = NULL, sbp = NULL, d) {
  if (is.null(parts) && is.null(sbp)) {
    .mlsim_stop("Compositional specs need `parts` or `sbp`.")
  }

  if (is.null(sbp)) {
    parts <- .mlsim_check_vars(parts, d + 1L, "parts")
    sbp <- build.sbp(parts)
  } else {
    sbp <- as.matrix(sbp)
    if (!is.numeric(sbp) || anyNA(sbp) || !(dim(sbp) %a==% c(d, d + 1L))) {
      .mlsim_stop("`sbp` must be a numeric %d by %d matrix.", d, d + 1L)
    }
    if (sbp %any!in% c(-1, 0, 1)) {
      .mlsim_stop("`sbp` must contain only -1, 0, and 1.")
    }
    if (is.null(colnames(sbp))) {
      if (is.null(parts)) {
        .mlsim_stop("Explicit `sbp` needs column names or `parts`.")
      }
      parts <- .mlsim_check_vars(parts, d + 1L, "parts")
      colnames(sbp) <- parts
    } else {
      sbp_parts <- .mlsim_check_vars(colnames(sbp), d + 1L, "colnames(sbp)")
      if (is.null(parts)) {
        parts <- sbp_parts
      } else {
        parts <- .mlsim_check_vars(parts, d + 1L, "parts")
        if (!identical(parts, sbp_parts)) {
          .mlsim_stop("`parts` must match `colnames(sbp)`.")
        }
      }
    }
  }

  rownames(sbp) <- rownames(sbp) %||% paste0("balance_", seq_len(d))
  .mlsim_validate_sbp(sbp, d)
  list(parts = parts, sbp = sbp)
}

#' @rdname multilevelcoda-internal-composition
.mlsim_validate_sbp <- function(sbp, d) {
  row_has_positive <- rowSums(sbp == 1) > 0L
  row_has_negative <- rowSums(sbp == -1) > 0L
  if (!all(row_has_positive & row_has_negative)) {
    .mlsim_stop("Every row of `sbp` must contain at least one 1 and one -1.")
  }

  column_participates <- colSums(sbp != 0) > 0L
  if (!all(column_participates)) {
    .mlsim_stop("Every part column in `sbp` must participate in at least one balance.")
  }

  basis <- compositions::gsi.buildilrBase(t(sbp))
  if (!is.numeric(basis) || anyNA(basis) || any(!is.finite(basis)) ||
      !(dim(basis) %a==% c(d + 1L, d))) {
    .mlsim_stop("`sbp` must define a finite ILR basis.")
  }
  if (max(abs(crossprod(basis) - diag(d))) > 1e-8) {
    .mlsim_stop(
      "`sbp` must be a nested sequential binary partition: the implied ILR basis is not orthonormal, so simulated ILR parameters would not be recovered by `ilr()`."
    )
  }
  invisible(TRUE)
}

#' @rdname multilevelcoda-internal-composition
.mlsim_check_composition_output_names <- function(coordinate_names, parts, keep_ilr) {
  if (!isTRUE(keep_ilr)) {
    return(invisible(TRUE))
  }
  overlap <- coordinate_names[coordinate_names %in% parts]
  if (length(overlap) > 0L) {
    .mlsim_stop(
      "Compositional output names must be unique; ILR coordinate names overlap with `parts`: %s.",
      paste(unique(overlap), collapse = ", ")
    )
  }
  invisible(TRUE)
}

#' @rdname multilevelcoda-internal-composition
.mlsim_ilr_coordinate_map <- function(sbp, coordinate_names) {
  data.table::data.table(
    ilr = coordinate_names,
    sbp_row = seq_len(nrow(sbp)),
    positive_parts = I(lapply(seq_len(nrow(sbp)), function(i) colnames(sbp)[sbp[i, ] == 1])),
    negative_parts = I(lapply(seq_len(nrow(sbp)), function(i) colnames(sbp)[sbp[i, ] == -1]))
  )
}
