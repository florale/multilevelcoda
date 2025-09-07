#' Checks if argument is a \code{complr} object
#'
#' @param x An object of class \code{complr}.
#'
#' @export
is.complr <- function(x) {
  inherits(x, "complr")
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
  
  parts <- get_parts(object = x, parts = parts)
  idx   <- which(vapply(lapply(x$output, 
                               function(x) x$parts), 
                        function(p) identical(sort(parts), sort(p)), logical(1)))
  X  <- x$output[[idx]]$X
  bX <- x$output[[idx]]$bX
  wX <- x$output[[idx]]$wX
  Z  <- x$output[[idx]]$Z
  bZ <- x$output[[idx]]$bZ
  wZ <- x$output[[idx]]$wZ
  
  if (identical(weight, "proportional") | is.null(x$idvar)) {
    out <- list(
      if(!is.null(X))  (summary(X, robust = TRUE)) else (NULL),
      if(!is.null(bX)) (summary(bX, robust = TRUE)) else (NULL),
      if(!is.null(wX)) (summary(wX, robust = TRUE)) else (NULL),
      if(!is.null(Z))  (data.frame(summary(Z))) else (NULL),
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
  
  out <- list(X  = out$X$mean, 
              bX = out$bX$mean,
              wX = out$wX$mean,
              Z  = out$Z["Mean", ],
              bZ = out$bZ["Mean", ],
              wZ = out$wZ["Mean", ]
  )
  print(lapply(out, function(x) {
    as.data.frame(t(x))
  }), row.names = FALSE)
  
  invisible(out)
}

#' Variance of compositions presented in a \code{complr} object.
#'
#' @inheritParams mean.complr
#'
#' @importFrom compositions var
#'
#' @method var complr
#' @export
var.complr <- function(x,
                       weight = c("equal", "proportional"),
                       parts = 1,
                       ...) {
  
  parts <- get_parts(object = x, parts = parts)
  idx   <- which(vapply(lapply(x$output, 
                               function(x) x$parts), 
                        function(p) identical(sort(parts), sort(p)), logical(1)))
  X  <- x$output[[idx]]$X
  bX <- x$output[[idx]]$bX
  wX <- x$output[[idx]]$wX
  Z  <- x$output[[idx]]$Z
  bZ <- x$output[[idx]]$bZ
  wZ <- x$output[[idx]]$wZ
  
  if (identical(weight, "proportional") | is.null(x$idvar)) {
    out <- list(
      if(!is.null(X)) (summary(X, robust = TRUE)) else (NULL),
      if(!is.null(bX)) (summary(bX, robust = TRUE)) else (NULL),
      if(!is.null(wX)) (summary(wX, robust = TRUE)) else (NULL),
      if(!is.null(Z)) (data.frame(summary(Z))) else (NULL),
      if(!is.null(bZ)) (data.frame(summary(bZ))) else (NULL),
      if(!is.null(wZ)) (data.frame(summary(wZ))) else (NULL)
    )
  } else {
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
  
  out <- list(X  = out$X$variation, 
              bX = out$bX$variation,
              wX = out$wX$variation
  )
  print(lapply(out, function(x) {
    as.data.frame(t(x))
  }), row.names = FALSE)
  
  invisible(out)
}

#' Extract Compositional Data from \code{complr} object.
#'
#' Extract amounts and compositions in conventional formats
#' as data.frames, matrices, or arrays.
#'
#' @inheritParams mean.complr
#' @param row.names,optional Unused and only added for consistency with
#' the \code{\link[base:as.data.frame]{as.data.frame}} generic.
#'
#' @export
as.data.frame.complr <- function(x,
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

#' Extract variable names from a \code{complr} object.
#' @param object A \code{complr} object
#' 
#' @method get_variables complr
#' 
#' @export
get_variables = function(object) UseMethod("get_variables")

#' @rdname get_variables
#' @export
get_variables.complr <- function(object) {
  
  out <- lapply(object$output, function(x) {
    list(names(x$X),
         names(x$bX),
         names(x$wX),
         names(x$Z),
         names(x$bZ),
         names(x$wZ))
  })
  
  out <- as.data.frame(do.call(cbind, out))
  rownames(out) <- c("X", "bX", "wX", "Z", "bZ", "wZ")
  colnames(out) <- c(paste0("composition_", seq_along(object$out)))

  out
}

#' Extract names parts  from a \code{complr} object.
#' 
#' Internal function used to validate and extract the names of compositional parts.
#' 
#' @param object A \code{complr} object
#' @param parts A optional character string specifying names of compositional parts that should be considered
#' in the substitution analysis. This should correspond to a single set of names of compositional parts specified
#' in the \code{complr} object. Default to the first composition in the \code{complr} object.
#' 
#' @keywords internal
#' @noRd
get_parts <- function(object, parts = 1) {
  
  if (isFALSE(inherits(object, "complr"))) {
    stop(sprintf(
      "Can't handle an object of class (%s)
  It should be a 'complr' object
  See ?complr for details.",
      class(object)))
  }
  
  if (is.numeric(parts)) {
    if (length(parts) > 1) {
      stop(" 'parts' should be a single numeric value indicating which set of compositional parts to use.")
    }
    if (parts < 1 || parts > length(object$output)) {
      stop(sprintf(
        " 'parts' should be a single numeric value between 1 and %s, corresponding to the number of sets of compositional parts in the 'complr' object.",
        length(object$output)))
    }
    parts <- object$output[[parts]]$parts
    
  } else {
    if (isFALSE(inherits(parts, "character"))) {
      stop(" 'parts' should be a character vector of compositional parts.")
    }
    ## parts should be identical with either one of the parts presented in output of complr
    if (isFALSE((any(vapply(lapply(object$output, function(x) x$parts), function(p) identical(sort(parts), sort(p)), logical(1)))))) {
      stop(sprintf(
        "The specified 'parts' (%s) are not found in the complr object.",
        "  It should corespond to one set of compositional parts, either one of the following:",
        "%s",
        paste(parts, collapse = ", "),
        invisible(lapply(object$output, function(x) cat(paste(x$parts, collapse = ", "), "\n"))),
        sep = "\n"))
    }
  }
  parts
}