#' Create a Summary of a \code{complr} object
#' 
#' @param object An object of class \code{complr}.
#' @param weight A character value specifying the weight to use in calculation of the reference composition.
#' If \code{"equal"}, give equal weight to units (e.g., individuals).
#' If \code{"proportional"}, weights in proportion to the frequencies of units being averaged 
#' (e.g., observations across individuals).
#' Default is \code{equal}.
#' @param ... generic argument, not in use.
#' 
#' @importFrom compositions summary.acomp summary.rmult clo acomp rmult
#' 
#' @method summary complr
#' 
#' @examples
#' cilr <- complr(data = mcompd, sbp = sbp, 
#'                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), 
#'                idvar = "ID")
#' summary(cilr)
#' @export
summary.complr <- function(object,
                           weight = c("equal", "proportional"),
                           ...) {
  
  out <- list()
  
  out$mean_comp <- mean.complr(object, weight = weight)$mean_comp
  out$mean_lr   <- mean.complr(object, weight = weight)$mean_lr
  out$variation_comp <- var.complr(object, weight = weight)
  
  out$transform_method <- list(
    transform_type = object$transform,
    psi = object$psi,
    total = object$total,
    sequential_binary_partition = object$sbp
  )
  out$data <- list(composition_parts = object$parts,
                   logratios = paste0(object$transform, seq_len(length(object$parts) - 1)),
                   idvar = if(exists("object$idvar")) (object$idvar) else (NULL),
                   nobs = nrow(object$data),
                   ngrps = length(unique(object$data[[object$idvar]]))
  )
  
  out$geometry <- list(composition_geometry = class(object$comp),
                       logratio_class = class(object$logratio))
  
  class(out) <- "summary.complr"
  
  return(out)
}

#' Summary for a \code{complr} object
#' 
#' @param x An object of class \code{summary.complr}.
#' @param ... Other arguments passed to \code{\link{summary.complr}}.
#' 
#' @seealso \code{\link{summary.complr}}
#' 
#' @method print summary.complr
#' 
#' @examples
#' 
#' cilr <- complr(data = mcompd, sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), 
#'                 idvar = "ID")
#' print(cilr)
#' @export
print.summary.complr <- function(x, ...) {
  
  out <- as.data.frame(as.array(append(append(x$data, x$transform_method[c("transform_type", "total")]), x$geometry)))
  names(out) <- NULL
  
  print(out, ...)
}

#' Print a Summary for a \code{complr} object
#' 
#' @param x An object of class \code{complr}.
#' @param ... Other arguments passed to \code{\link{summary.complr}}.
#' 
#' @seealso \code{\link{summary.complr}}
#' 
#' @examples
#' 
#' cilr <- complr(data = mcompd, sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), 
#'                 idvar = "ID")
#' print(cilr)
#' @export
print.complr <- function(x, ...) {
  summary(x, ...)
}

#' Create a Summary of a fitted \code{brmsfit} model in a \code{brmcoda} object
#' 
#' @param object An object of class \code{brmcoda}.
#' @param ... Other arguments passed to \code{\link{summary.brmsfit}}.
#' 
#' @method summary brmcoda
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                  idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'   
#'   summary(m)
#' }}
#' @export
summary.brmcoda <- function(object, ...) {
  if (inherits(object$model, "brmsfit")) {
    summary(object$model, ...)
  }
}

#' Print a Summary for a fitted \code{brmsfit} model in a \code{brmcoda} object
#' 
#' @param x An object of class \code{brmcoda}.
#' @param ... Other arguments passed to \code{summary.brmcoda}.
#' 
#' @seealso \code{\link{summary.brmcoda}}
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'   
#'  print(m)
#' }}
#' @export
print.brmcoda <- function(x, ...) {
  print(summary(x$model, ...), ...)
  
}

#' Create a Summary of a Substitution Model represented by a \code{substitution} object
#' 
#' @param object A \code{substitution} class object.
#' @param delta A integer, numeric value or vector indicating the desired \code{delta} 
#' at which substitution results should be summarised.
#' Default to all \code{delta} available in the \code{substitution} object.
#' @param to A character value or vector specifying the names of the compositional parts
#' that were reallocated to in the model.
#' @param from A character value or vector specifying the names of the compositional parts
#' that were reallocated from in the model.
#' @param ref Either a character value or vector ((\code{"grandmean"} and/or \code{"clustermean"} or \code{"users"}),
#' Default to all \code{ref} available in the \code{substitution} object.
#' @param level A character string or vector (\code{"between"} and/or \code{"within"}).
#' Default to all \code{level} available in the \code{substitution} object.
#' @param digits A integer value used for number formatting. Default is \code{2}.
#' @param ... generic argument, not in use.
#' 
#' @return A summary of \code{substitution} object.
#'   \item{\code{Mean}}{ Posterior means.}
#'   \item{\code{CI_low} and \code{CI_high}}{ 95% credible intervals.}
#'   \item{\code{Delta}}{ Amount substituted across compositional parts.}
#'   \item{\code{From}}{ Compositional part that is substituted from.}
#'   \item{\code{To}}{ Compositional parts that is substituted to.}
#'   \item{\code{Level}}{ Level where changes in composition takes place. Either \code{between} or \code{within}.}
#'   \item{\code{Reference}}{ Either \code{grandmean}, \code{clustermean}, or \code{users}.}
#' 
#' @method summary substitution
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   ## fit a model with compositional predictor at between and between-person levels
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'   
#'   subm <- substitution(object = m, delta = 5)
#'   summary(subm)
#' }}
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
      stop("'ref' should be either one of the following: \"grandmean\", \"clustermean\", or \"users\".")
    }
    ref <- as.character(ref)
  }
  if (isTRUE(missing(level))) {
    level <- object$level
  } else {
    if (isFALSE(any(c("between", "within", "aggregate") %in% level))) {
      stop("'level' should be either one of the following: \"between\", \"within\", \"aggregate\".")
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
  
  out <- lapply(object[c("between_simple_sub", "within_simple_sub", "simple_sub",
                         "between_avg_sub", "within_avg_sub", "avg_sub")],
                rbindlist)
  out <- rbindlist(out, use.names = TRUE)
  out <- out[Delta %in% delta & Level %in% level & Reference %in% ref & To %in% to & From %in% from]
  
  if(isTRUE(dim(out)[1] == 0)) {
    stop("An empty data.table returned. Please check that the arguments match with your substitution object.")
  }
  
  if(isFALSE(digits == "asis")) {
    # out[, 1:3] <- round(out[, 1:3], digits)
    out[] <- lapply(out, function(X) if(is.numeric(X)) round(X, digits) else X)
  }
  
  out
}

#' Print a Summary for a \code{substitution} object
#' 
#' @param x A \code{substitution} object.
#' @param ... Additional arguments to be passed to to method \code{summary} of \code{substitution}.
#' 
#' @seealso \code{\link{summary.substitution}}
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   ## fit a model with compositional predictor at between and between-person levels
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'   
#'   subm <- substitution(object = m, delta = 5)
#'   print(subm)
#' }}
#' @export
print.substitution <- function(x, ...) {
  print(summary(x, ...), ...)
}