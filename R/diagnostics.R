#' Generic diagnostics function for \pkg{multilevelcoda}
#' 
#' This is a generic function for diagnostics methods in \pkg{multilevelcoda}. See
#' \code{\link{diagnostics.complr}} for the method for \code{\link{complr}} objects.
#' 
#' @param x An object to calculate diagnostics, primarily around extreme values
#'   and assumed distributions.
#' @param ... Additional arguments passed to methods.
#' @return \code{diagnostics} class object
#' @export
#' @keywords multivariate
#' @rdname diagnostics
diagnostics <- function(x, ...) {
  UseMethod("diagnostics", x)
}

#' Checks if argument is a \code{diagnostics} object
#'
#' @param x An object of class \code{diagnostics}.
#'
#' @export
is.diagnostics <- function(x) {
  inherits(x, "diagnostics")
}

#' Coerce a list to a \code{diagnostics} object
#'
#' This is a constructor function for a \code{diagnostics} object.
#' It checks that the input list has the necessary structure and
#' components to be considered a \code{diagnostics} object,
#' and if so, it assigns the class "diagnostics" to it.
#' This allows for method dispatch on \code{diagnostics} objects.
#'
#' @param x A list with elements \code{x}, \code{distance}, \code{cutoff},
#' \code{ev.perc}, \code{extremevalues}, and \code{levels}.
#'
#' @return A \code{diagnostics} object.
#' @importFrom extraoperators %ain%
#'
#' @export
as.diagnostics <- function(x) {
  if (is.diagnostics(x)) {
    return(x)
  }

  stopifnot(is.list(x))

  required_names <- c("x", "distance", "cutoff", "ev.perc", "extremevalues", "levels")

  if (!required_names %ain% names(x) || !names(x) %ain% required_names) {
    stop("x must contain exactly the following elements: ",
      paste(required_names, collapse = ", "))
  }

  if (!is.matrix(x$x)) {
    stop(sprintf("x must be a matrix but has class %s.", class(x$x)))
  }

  if (!inherits(x$distance, "data.table")) {
    stop(sprintf("distance must be a data.table but has class %s.", class(x$distance)))
  }

  if (!is.numeric(x$cutoff) || length(x$cutoff) != 1L) {
    stop(sprintf(
      "cutoff must be a length 1 numeric value, but has class %s and length %d.",
      class(x$cutoff), length(x$cutoff)
    ))
  }

  if (!is.numeric(x$ev.perc) || length(x$ev.perc) != 1L) {
    stop(sprintf(
      "ev.perc must be a length 1 numeric value, but has class %s and length %d.",
      class(x$ev.perc), length(x$ev.perc)
    ))
  }

  if (!is.character(x$extremevalues) || length(x$extremevalues) != 1L) {
    stop(sprintf(
      "extremevalues must be a length 1 character string, but has class %s and length %d.",
      class(x$extremevalues), length(x$extremevalues)
    ))
  }

  if (!is.character(x$levels) || length(x$levels) != 1L) {
    stop(sprintf(
      "levels must be a length 1 character string, but has class %s and length %d.",
      class(x$levels), length(x$levels)
    ))
  }

  x <- x[required_names]  # ensure correct order of elements

  structure(x, class = "diagnostics")
}

#' Robust multivariate normal diagnostics for \code{complr} coordinates
#'
#' Uses robust minimum covariance determinant estimates to calculate squared
#' Mahalanobis distances for the total, between, or within log-ratio
#' coordinates stored in a \code{\link{complr}} object.
#' This can be used to identify extreme values (outliers).
#'
#' @param x A \code{\link{complr}} object.
#' @param level A character string indicating which coordinates to diagnose:
#' \code{"total"}, \code{"between"}, or \code{"within"}.
#' @param parts Optional character vector or integer specifying which set of
#' compositional parts should be used. Defaults to the first composition in the
#' \code{complr} object.
#' @param ev.perc Proportion in the upper tail used to flag extreme values.
#' Defaults to \code{.001}.
#' @param extremevalues Character string indicating whether extreme values
#' should be flagged using a \code{"theoretical"} chi-squared cutoff or an
#' \code{"empirical"} cutoff from the observed robust distances.
#' @param ... Additional arguments passed to \code{\link[robustbase:covMcd]{covMcd}}.
#'
#' @return A \code{diagnostics} object with the raw composition matrix in
#' \code{x}, a \code{data.table} of robust distances in \code{distance}, the
#' extreme-value cutoff in \code{cutoff}, the requested upper-tail proportion
#' in \code{ev.perc}, the cutoff type in \code{extremevalues}, and the
#' diagnosed coordinate level in \code{levels}. The \code{distance} element
#' includes an \code{is_extremevalue} column indicating whether each fitted
#' observation is flagged as an extreme value. For \code{level = "between"},
#' diagnostics are fitted using one observation per \code{idvar}.
#'
#' @importFrom data.table data.table
#' @importFrom robustbase covMcd
#' @importFrom extraoperators %gele%
#' @importFrom stats mahalanobis qchisq ppoints resid lm
#' @export
#'
#' @examples
#' data(mcompd)
#' data(sbp)
#'
#' ids <- unique(mcompd$ID)[1:20]
#' cilr <- complr(
#'   data = mcompd[ID %in% ids, .SD[1:3], by = ID],
#'   sbp = sbp,
#'   parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'   idvar = "ID",
#'   total = 1440
#' )
#'
#' # One diagnostic object per coordinate level
#' total_diag <- diagnostics(cilr, level = "total")
#' between_diag <- diagnostics(cilr, level = "between")
#' within_diag <- diagnostics(cilr, level = "within")
#'
#' # Use an empirical cutoff and a larger tail proportion
#' empirical_diag <- diagnostics(
#'   cilr,
#'   level = "between",
#'   ev.perc = .05,
#'   extremevalues = "empirical"
#' )
#'
#' is.diagnostics(between_diag)
#' head(between_diag$distance)
diagnostics.complr <- function(x,
                       level = c("total", "between", "within"),
                       parts = 1,
                       ev.perc = .001,
                       extremevalues = c("theoretical", "empirical"),
                       ...) {
  stopifnot(is.complr(x))

  level <- match.arg(level)
  extremevalues <- match.arg(extremevalues)

  if (length(ev.perc) != 1L || is.na(ev.perc) || !ev.perc %gele% c(0, 1)) {
    stop("ev.perc must be a single numeric value between 0 and 1.")
  }

  if (identical(level, "between") && is.null(x$idvar)) {
    stop("between-level diagnostics require a complr object with an idvar.")
  }

  coord_type <- switch(
    level,
    total = "Z",
    between = "bZ",
    within = "wZ"
  )

  raw_type <- switch(
    level,
    total = "X",
    between = "bX",
    within = "wX"
  )

  coords <- x$output[[parts]][[coord_type]]
  raw <- x$output[[parts]][[raw_type]]

  if (is.null(coords)) {
    stop(sprintf(
      "%s-level coordinates are not available in this complr object.",
      level
    ))
  }

  if (is.null(raw)) {
    stop(sprintf(
      "%s-level raw values are not available in this complr object.",
      level
    ))
  }

  coords <- as.matrix(coords)
  raw <- as.matrix(raw)

  stopifnot(all(!is.na(coords)))
  stopifnot(all(!is.na(raw)))

  fit_index <- rep(TRUE, nrow(coords))

  if (identical(level, "between")) {
    ids <- x$dataout[[x$idvar]]
    fit_index <- !duplicated(ids)
  }

  coords <- coords[fit_index, , drop = FALSE]
  raw <- raw[fit_index, , drop = FALSE]

  n_obs <- nrow(coords)
  n_vars <- ncol(coords)

  if (n_obs <= n_vars) {
    stop(sprintf(
      "At least %d complete observations are required for %d coordinates.",
      n_vars + 1L,
      n_vars
    ))
  }

  desc <- covMcd(coords, ...)

  ## if any variances near zero, set to at least 1e-10
  if (any(diag(desc$cov) < 1e-10)) {
    diag(desc$cov) <- pmax(diag(desc$cov), 1e-10)
  }
  ## if covariance matrix is singular, inflate diagonal (variances)
  ## 2.5 percent increase per iteration, up to 20 iterations
  i <- 1L
  while(isTRUE(tryCatch(solve(desc$cov), error = function(e) TRUE)) && i <= 20) {
    diag(desc$cov) <- diag(desc$cov) * 1.025
    i <- i + 1L
  }

  distance <- mahalanobis(coords, desc$center, desc$cov)

  expected_distance <- stats::qchisq(
    stats::ppoints(n_obs),
    df = n_vars
  )[order(order(distance))]

  deviates <- resid(lm(distance ~ 0 + offset(expected_distance)))

  cutoff <- switch(
    extremevalues,
    theoretical = stats::qchisq(1 - ev.perc, df = n_vars),
    empirical = stats::quantile(
      distance,
      probs = 1 - ev.perc,
      na.rm = TRUE,
      names = FALSE
    )
  )

  dist <- data.table(
    distance = distance,
    expected_distance = expected_distance,
    deviates = deviates,
    is_extremevalue = distance > cutoff
  )

  as.diagnostics(list(
    x = raw,
    distance = dist,
    cutoff = cutoff,
    ev.perc = ev.perc,
    extremevalues = extremevalues,
    levels = level))
}
