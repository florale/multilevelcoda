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
#' @param ... generic argument, not in use.
#' 
#' @method mean complr
#' 
#' @examples
#' cilr <- complr(data = mcompd, sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), 
#'                 idvar = "ID")
#' mean(cilr)
#' @export
mean.complr <- function(x,
                        weight = c("equal", "proportional"),
                        ...) {
  
  stats <- .get.complr(object = x,
                       weight = weight)
  
  mean_comp <- do.call(rbind, lapply(stats[c("comp", "between_comp", "within_comp")], "[[", "mean"))
  mean_lr   <- do.call(rbind, lapply(stats[c("logratio", "between_logratio", "within_logratio")], function(x) x["Mean",]))
  
  out <- list(mean_comp = mean_comp,
              mean_lr = mean_comp)
  out
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
                       ...) {

  stats <- .get.complr(object = x,
                       weight = weight)
  
  out <- lapply(stats[c("comp", "between_comp", "within_comp")], "[[", "variation")
  out
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
                                 row.names = NULL, optional = TRUE,
                                 ...) {
  as.data.frame(do.call(cbind,
                        x[c("comp", "between_comp", "within_comp",
                            "logratio", "between_logratio", "within_logratio")]))
}

#' @rdname as.data.frame.complr
#' @export
as.matrix.complr <- function(x, ...) {
  as.matrix(as.data.frame(x, ...))
}
