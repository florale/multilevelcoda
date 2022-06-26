#' @title Estimating Average Marginal Effects for Within-person Isotemporal Substitution Model.
#'
#' @description
#' Using a fitted model object, estimates the the average marginal difference 
#' in outcomes when compositional variables are substituted for a specific period
#' at `within-person` level. The resulting \code{wsubmargins} encapsulates 
#' substitution estimation across all compositional variables present
#' in the \code{\link{brmcoda}} object.
#'
#' @param object A \code{\link{brmcoda}} object. Required.
#' @param substitute A \code{data.frame} or \code{data.table} of the possible substitution of variables.
#' This dataset can be computed using function \code{possub}. Required.
#' @param minute A integer or numeric value indicating the maximum minute for which substitution model is desired.
#' Default to \code{60L} (i.e., the model loops through 1:60L minutes).
#' @param ... Additional arguments to be passed to \code{\link{describe_posterior}}.
#' 
#' @return A list containing the result of isotemporal multilevel substitution model.
#' Each elements of the list is the substitution estimation for a compositional variable.
#'
#' @importFrom data.table data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo
#' @importFrom extraoperators %snin% %sin%
#' @importFrom bayestestR describe_posterior
#' @importFrom foreach %dopar%
#' @importFrom stats fitted
#' @export
#' @examples
#' 
#' data(sbp)
#' data (mcompd)
#' 
#' ps <- possub(parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' 
#' library(doFuture)
#' registerDoFuture()
#' plan(multisession, workers = 5)
#' system.time(testwsm <- wsubmargins(object = adjm, substitute = ps, minute = 2))
wsubmargins <- function (object, substitute, minute = 60, ...) {
  
  if (isTRUE(missing(object))) {
    stop(paste(
      "'object' is a required argument and cannot be missing;",
      " it should be an object of class 'brmcoda'.", 
      " See ?wsubmargins for details.",
      sep = "\n"))
  }
  
  if (isFALSE(inherits(object, "brmcoda"))) {
    stop(sprintf(
    "Can't handle an object of class (%s) 
  It should be a fitted 'brmcoda' object
  See ?bsub for details.",
                 class(object)))
  }
  
  if (isTRUE(missing(substitute))) {
    stop(paste(
      "'substitute' is a required argument and cannot be missing;",
      " it should be a dataset of possible substitution", 
      " and can be computed using multilevelcoda::possub.", 
      " See ?wsubmargins for details.",
      sep = "\n"))
  }
  
  if(isFALSE(missing(minute))) {
    if (isFALSE(is.integer(minute))) {
      if (isFALSE(minute > 0)) {
        stop("'minute' must be an positive integer value.")
      }
    }
  } else {
    minute <- 60L
  }
  
  if (isFALSE(identical(ncol(substitute), length(object$CompIlr$parts)))) {
    stop(sprintf(
    "The number of columns in 'substitute' (%d) must be the same
  as the compositional variables in 'parts' (%d).",
                 ncol(substitute),
                 length(object$CompIlr$parts)))
  }
  
  if (isFALSE(identical(colnames(substitute), object$CompIlr$parts))) {
    stop(sprintf(
    "The names of compositional variables must be the same
  in 'substitute' (%s) and 'parts' (%s).",
                 colnames(substitute),
                 object$CompIlr$parts))
  }
  
  # between-person composition
  b <- object$CompIlr$BetweenComp
  b <- as.data.table(clo(b, total = object$CompIlr$total))
  
  min <- as.integer(minute)
  
  # model for no change
  bilr <- object$CompIlr$BetweenILR
  wilr <- as.data.table(matrix(0, nrow = nrow(bilr), ncol = ncol(bilr)))
  
  colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
  colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
  
  samed <- cbind(bilr, wilr, object$CompIlr$data)
  ysame <- fitted(object$Model, newdata = samed, re.form = NA, summary = FALSE)
  ysame <- rowMeans(ysame) # average across participants when there is no change
  
  # substitution model
  out <- .get.wsubmargins(object = object, b = b,
                          substitute = substitute,
                          ysame = ysame, min = min)
}