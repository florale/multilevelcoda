#' Summary \code{\link{substitution}} 
#' 
#' This method allows for summarising an existing \code{\link{substitution}} object.
#' 
#' @param object A \code{\link{substitution}} object containing the output of substitution models.
#' @param delta A integer, numeric value or vector indicating the desired \code{delta} 
#' at which substitution results should be summarised.
#' Default to all \code{delta} available in the \code{\link{substitution}} object.
#' @param to and @param from A character value or vector specifying the names of the compositional parts
#' that were reallocated to/from in the model.
#' @param ref Either a character value or vector ((\code{grandmean} and/or \code{clustermean} or \code{users}),
#' Default to all \code{ref} available in the \code{\link{substitution}} object .
#' @param level A character string or vector (\code{between} and/or \code{within}).
#' Default to all \code{level} available in the \code{\link{substitution}} object .
#' @param digits A integer value used for number formatting. Default to 3.
#' @param ... generic argument, not in use.
#' 
#' @return A summary of \code{\link{substitution}} object with twelve elements.
#' \itemize{
#'   \item{\code{BetweenComp}}{ A vector of class \code{acomp} representing one closed between-person composition
#'   or a matrix of class \code{acomp} representing multiple closed between-person compositions each in one row.}
#'   \item{\code{WithinComp}}{ A vector of class \code{acomp} representing one closed within-person composition
#'   or a matrix of class \code{acomp} representing multiple closed within-person compositions each in one row.}
#'   \item{\code{TotalComp}}{ A vector of class \code{acomp} representing one closed total composition
#'   or a matrix of class \code{acomp} representing multiple closed total compositions each in one row.}
#'   \item{\code{BetweenILR}}{ Isometric log ratio transform of between-person composition.}
#'   \item{\code{WithinILR}}{ Isometric log ratio transform of within-person composition.}
#'   \item{\code{TotalILR}}{ Isometric log ratio transform of total composition.}
#'   \item{\code{data}}{ The user's dataset or imputed dataset if the input data contains zeros.}
#'   \item{\code{psi}}{ A ILR matrix associated with user-defined partition structure.}
#'   \item{\code{sbp}}{ The user-defined sequential binary partition matrix.}
#'   \item{\code{parts}}{ Names of compositional variables.}
#'   \item{\code{idvar}}{ Name of the variable containing IDs.}
#'   \item{\code{total}}{ Total amount to which the compositions is closed.}
#' }
#' 
#' @exportS3Method summary substitution
summary.substitution <- function(object, delta, to, from,
                                 ref, level,
                                 digits = 2, ...) { 
  
  if (isTRUE(missing(delta))) {
    delta <- object$delta
  } else {
    delta <- as.integer(delta)
  }
  
  if (isTRUE(missing(ref))) {
    ref <- object$ref
  } else {
    if (isFALSE(any(c("grandmean", "clustermean", "users") %in% ref))) {
      stop("'ref' should be grandmean and/or clustermean or users.")
    }
    ref <- as.character(ref)
  }
  if (isTRUE(missing(level))) {
    level <- object$level
  } else {
    if (isFALSE(any(c("between", "within") %in% level))) {
      stop("'level' should be between and/or within.")
    }
    level <- as.character(level)
  }
  
  if (isTRUE(missing(to))) {
    to <- object$parts
  } else {
    if (isFALSE(any(object$parts %in% to))) {
      stop("'to' should be names of one or more compositional parts present in the substitution model.")
    }
    to <- as.character(to)
  }
  if (isTRUE(missing(from))) {
    from <- object$parts
  } else {
    if (isFALSE(any(object$parts %in% from))) {
      stop("'from' should be names of one or more compositional parts present in the substitution model.")
    }
    from <- as.character(from)
  }
  
  out <- lapply(lapply(object[1:4], function(x) {
    lapply(x, function(y) {
      y <- y[Delta %in% delta & Level %in% level & Reference %in% ref]
      y[, 1:3] <- round(y[, 1:3], digits)
      y
    })
  }), rbindlist)
  
  out <- rbindlist(out, use.names = TRUE)
  out <- out[To %in% to & From %in% from]
  
  out
  
}
