#' Create a summary of a fitted \code{brmsfit} model in a \code{brmcoda} object
#' 
#' @param object An object of class \code{brmcoda}.
#' @param ... Other arguments passed to \code{summary.brmsfit}
#' 
#' @method summary brmcoda
#' @export
summary.brmcoda <- function(object, ...) {
  if (inherits(object$Model, "brmsfit")) {
    summary(object$Model, ...)
  }
}

#' Print a summary for a fitted \code{brmsfit} model in a \code{brmcoda} object
#' 
#' @param x A \code{brmcoda} object.
#' @param ... Additional arguments to be passed to to method \code{summary} of \code{brmsfit}.
#' 
#' @seealso \code{\link{summary.brmcoda}}
#' 
#' @export
print.brmcoda <- function(x, ...) {
  print(summary(x$Model, ...), ...)
  
}

#' Create a summary of a substitution model represented by a \code{\link{substitution}} object
#' 
#' @param object A \code{\link{substitution}} class object.
#' @param delta A integer, numeric value or vector indicating the desired \code{delta} 
#' at which substitution results should be summarised.
#' Default to all \code{delta} available in the \code{\link{substitution}} object.
#' @param to A character value or vector specifying the names of the compositional parts
#' that were reallocated to in the model.
#' @param from A character value or vector specifying the names of the compositional parts
#' that were reallocated from in the model.
#' @param ref Either a character value or vector ((\code{grandmean} and/or \code{clustermean} or \code{users}),
#' Default to all \code{ref} available in the \code{\link{substitution}} object .
#' @param level A character string or vector (\code{between} and/or \code{within}).
#' Default to all \code{level} available in the \code{\link{substitution}} object.
#' @param digits A integer value used for number formatting. Default to 3.
#' @param ... generic argument, not in use.
#' 
#' @return A summary of \code{\link{substitution}} object.
#' \itemize{
#'   \item{\code{Mean}}{ Posterior means.}
#'   \item{\code{CI_low} and \code{CI_high}}{ 95% credible intervals.}
#'   \item{\code{Delta}}{ Amount substituted across compositional parts.}
#'   \item{\code{From}}{ Compositional part that is substituted from.}
#'   \item{\code{To}}{ Compositional parts that is substituted to.}
#'   \item{\code{Level}}{ Level where changes in composition takes place. Either \code{between} or \code{within}.}
#'   \item{\code{Reference}}{ Either \code{grandmean}, \code{clustermean}, or \code{users}.}
#' }
#' 
#' @method summary substitution
#' @export
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
  
  out <- lapply(object[1:4], rbindlist)
  out <- rbindlist(out, use.names = TRUE)
  out <- out[Delta %in% delta & Level %in% level & Reference %in% ref & To %in% to & From %in% from]
  
  if(isTRUE(dim(out)[1] == 0)) {
    stop("An empty data.table returned. Are you sure your arguments are correct?")
  }
    
  if(isFALSE(digits == "asis")) {
    out[, 1:3] <- round(out[, 1:3], digits)
  }
  
  out
  
}

#' Print a summary for a \code{substitution} object
#' 
#' @param x A \code{substitution} object.
#' @param ... Additional arguments to be passed to to method \code{summary} of \code{substitution}.
#' 
#' @seealso \code{\link{summary.substitution}}
#' 
#' @export
print.substitution <- function(x, ...) {
  print(summary(object, ...), ...)
}