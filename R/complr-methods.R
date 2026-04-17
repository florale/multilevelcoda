#' Checks if argument is a \code{complr} object
#'
#' @param x An object of class \code{complr}.
#'
#' @export
is.complr <- function(x) {
  inherits(x, "complr")
}

#' Coerce a list to a \code{complr} object
#' 
#' This is a constructor function for a \code{complr} object.
#' It checks that the input list has the necessary structure and
#' components to be considered a \code{complr} object,
#' and if so, it assigns the class "complr" to it.
#' This allows for method dispatch on \code{complr} objects.
#'
#' @param x A list with elements \code{output}, \code{datain}, \code{dataout},
#' \code{transform}, and \code{idvar}.
#' 
#' @return A \code{complr} object.
#' @importFrom data.table is.data.table
#' @importFrom extraoperators %nin% %ain%
#'
#' @export
as.complr <- function(x) {
  if (is.complr(x)) {
    return(x)
  }

  stopifnot(is.list(x))

  required_names <- c("output", "datain", "dataout", "transform", "idvar")

  if (!required_names %ain% names(x) || !names(x) %ain% required_names) {
    stop("x must contain exactly the following elements: ",
      paste(required_names, collapse = ", "))
  }

  if (!is.list(x$output)) {
    stop(sprintf("output must be a list but has class %s.", class(x$output)))
  }

  if (!is.data.table(x$datain) || !is.data.table(x$dataout)) {
    stop(sprintf("datain and dataout must be data.tables but have classes %s and %s.",
      class(x$datain), class(x$dataout)))
  }

  if (!is.character(x$transform) || length(x$transform) != 1L) {
    stop(sprintf(
      "transform must be a length 1 character string, but has class %s and length %d.",
      class(x$transform), length(x$transform)
    ))
  }

  if (!is.null(x$idvar)) {
    if (!is.character(x$idvar) || length(x$idvar) != 1L) {
     stop(sprintf(
        "idvar must be a length 1 character string, but has class %s and length %d.",
        class(x$idvar), length(x$idvar)
      ))
    }
    if (x$idvar %nin% names(x$dataout)) {
      stop(sprintf(
        "idvar '%s' is not a column in dataout.",
        x$idvar
      ))
    }
  }

  x <- x[required_names]  # ensure correct order of elements

  structure(x, class = "complr")
}

#' Mean and Variance of compositions presented in a \code{complr} object.
#' 
#' Internal function only.
#'
#' @inheritParams mean.complr
.meanvar.complr <- function(x,
                       weight = c("equal", "proportional"),
                       parts = 1,
                       ...) {
  
  parts <- .get_parts(object = x, parts = parts)
  idx   <- which(vapply(lapply(x$output,
                               function(x) x$parts),
                        function(p) identical(sort(parts), sort(p)), logical(1)))
  X  <- x$output[[idx]]$X
  bX <- x$output[[idx]]$bX
  wX <- x$output[[idx]]$wX
  Z  <- x$output[[idx]]$Z
  bZ <- x$output[[idx]]$bZ
  wZ <- x$output[[idx]]$wZ
  
  if (identical(weight, "proportional") || is.null(x$idvar)) {
    out <- list(
      if(!is.null(X)) (summary(X, robust = TRUE)) else (NULL),
      if(!is.null(bX)) (summary(bX, robust = TRUE)) else (NULL),
      if(!is.null(wX)) (summary(wX, robust = TRUE)) else (NULL),
      if(!is.null(Z)) (data.frame(summary(Z))) else (NULL),
      if(!is.null(bZ)) (data.frame(summary(bZ))) else (NULL),
      if(!is.null(wZ)) (data.frame(summary(wZ))) else (NULL)
    )
  } else {
    data <- cbind(x$dataout[, x$idvar, with = FALSE], data.frame(X, bX, wX, Z, bZ, wZ))
    out  <- list(
      # average X by ID in data
      summary(acomp(data[, lapply(.SD, mean), by = get(x$idvar), .SDcols = colnames(X)][, colnames(X), with = FALSE]), robust = TRUE),
      summary(acomp(data[, lapply(.SD, mean), by = get(x$idvar), .SDcols = colnames(bX)][, colnames(bX), with = FALSE]), robust = TRUE),
      summary(acomp(data[, lapply(.SD, mean), by = get(x$idvar), .SDcols = colnames(wX)][, colnames(wX), with = FALSE]), robust = TRUE),
      
      # average Z by ID in data
      summary(rmult(data[, lapply(.SD, mean), by = get(x$idvar), .SDcols = colnames(Z)][, colnames(Z), with = FALSE])),
      summary(rmult(data[, lapply(.SD, mean), by = get(x$idvar), .SDcols = colnames(bZ)][, colnames(bZ), with = FALSE])),
      summary(rmult(data[, lapply(.SD, mean), by = get(x$idvar), .SDcols = colnames(wZ)][, colnames(wZ), with = FALSE]))
    )
  }
  names(out) <- c("X", "bX", "wX", "Z", "bZ", "wZ")

  out
}


#' Mean amounts and mean compositions presented in a \code{complr} object.
#'
#' @param x An object of class \code{complr}.
#' @param weight A character value specifying the weight to use in calculation of the reference composition.
#' If \code{"equal"}, give equal weight to units (e.g., individuals).
#' If \code{"proportional"}, weights in proportion to the frequencies of units being averaged
#' (e.g., observations across individuals)
#' Default is \code{equal}.
#' @param parts A optional character string specifying names of compositional parts that should be considered
#' in the substitution analysis. This should correspond to a single set of names of compositional parts specified
#' in the \code{complr} object. Default to the first composition in the \code{complr} object.
#' @param ... generic argument, not in use.
#'
#' @method mean complr
#'
#' @examples
#' x <- complr(data = mcompd, sbp = sbp,
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                 idvar = "ID")
#' mean(x)
#' @export
mean.complr <- function(x,
                        weight = c("equal", "proportional"),
                        parts = 1,
                        ...) {
  out <- .meanvar.complr(x = x, weight = weight, parts = parts, ...)

  out <- list(X  = out$X$mean,
              bX = if(!is.null(out$bX$mean)) out$bX$mean else NULL,
              wX = out$wX$mean,
              Z  = out$Z["Mean", ],
              bZ = out$bZ["Mean", ],
              wZ = out$wZ["Mean", ]
  )
  # remove null elments
  out <- out[!vapply(out, is.null, logical(1))]

  return(out)
}

#' Variance of compositions presented in a \code{complr} object.
#'
#' @inheritParams mean.complr
#'
#' @importFrom compositions var
#'
#' @method var complr
#' @examples
#' x <- complr(data = mcompd, sbp = sbp,
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                 idvar = "ID")
#' ## ensure dispatch to the s3 generic var function
#' ## defined in the compositions package
#' compositions::var(x)
#' @export
var.complr <- function(x,
                       weight = c("equal", "proportional"),
                       parts = 1,
                       ...) {
  
  out <- .meanvar.complr(x = x, weight = weight, parts = parts, ...)

  out <- list(X  = out$X$variation,
              bX = out$bX$variation,
              wX = out$wX$variation
  )

  return(out)
}

#' Extract amounts and compositions in conventional formats
#' as data.frames, matrices, or arrays.
#'
#' @inheritParams mean.complr
#' @param row.names,optional Unused and only added for consistency with
#' the \code{\link[base:as.data.frame]{as.data.frame}} generic.
#'
#' @export
as.data.frame.complr <- function(x,
                                 row.names = NULL, optional = FALSE,
                                 ...) {
  do.call(cbind, lapply(x$output, function(y) {
    data.frame(y$X, y$bX, y$wX, y$Z, y$bZ, y$wZ)
  }))
}

#' @rdname as.data.frame.complr
#' @export
as.matrix.complr <- function(x, ...) {
  as.matrix(as.data.frame(x, ...))
}
