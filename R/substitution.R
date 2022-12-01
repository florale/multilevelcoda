#' @title Between-person Substitution Model.
#' 
#' @description 
#' Estimate the difference in outcomes
#' when compositional variables are substituted for a specific unit of change. 
#' The resulting \code{substitution} encapsulates
#' substitution estimation across all compositional variables
#' present in the \code{\link{brmcoda}} object.
#' 
#' @param object A fitted \code{\link{brmcoda}} object.
#' @param substitute A \code{data.frame} or \code{data.table} of possible substitution of variables.
#' This dataset can be computed using function \code{possub}. Required.
#' @param delta A integer, numeric value or vector indicating the amount of change in compositional parts
#' for substitution.
#' @param regrid If non-\code{NULL}, a \code{data.table} of reference grid consisting 
#' of combinations of covariates over which predictions are made.
#' Otherwise, the reference grid is constructed via \code{\link{ref_grid}}.
#' @param summary A logical value. 
#' Should the estimate at each level of the reference grid (\code{FALSE}) 
#' or their average (\code{TRUE}) be returned? Default to \code{TRUE}.
#' @param ... Additional arguments to be passed to \code{\link{describe_posterior}}.
#' 
#' @return A list containing the result of isotemporal multilevel substitution model.
#' Each elements of the list is the substitution estimation for a compositional variables, 
#' which include at least six elements.
#' \itemize{
#'   \item{\code{Mean}}{ Posterior means.}
#'   \item{\code{CI_low}} and \item{\code{CI_high}}{ 95% credible intervals.}
#'   \item{\code{Delta}}{ Amount substituted across compositional parts.}
#'   \item{\code{From}}{ Compositional variable that is substituted from.}
#'   \item{\code{To}}{ Compositional variables that is substituted to.}
#'   \item{\code{Level}}{}
#'   \item{\code{Type}}{}
#' }
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo mean.acomp
#' @importFrom extraoperators %snin% %sin%
#' @importFrom insight find_predictors
#' @importFrom emmeans ref_grid
#' @importFrom stats fitted
#' @export
#' @examples
#' \dontrun{
#' data(mcompd)
#' data(sbp)
#' data(psub)
#' cilr <- compilr(data = mcompd, sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#' 
#' # model with compositional predictor at between and between-person levels
#' m <- brmcoda(compilr = cilr, 
#'              formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + 
#'                                 wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID), 
#'              chain = 1, iter = 500)
#'              
#' subm <- substitution(object = m, substitute = psub, delta = 1,
#'                      type = "marginal", level = "between")
#' }
substitution <- function(object, substitute, delta, 
                         regrid = NULL, summary = TRUE, 
                         level = c("between", "within"),
                         type = c("conditional", "marginal"),
                         ...) {
  
  
  if (isTRUE(missing(object))) {
    stop(paste(
      "'object' is a required argument and cannot be missing;",
      "  it should be an object of class 'brmcoda'.", 
      "  See ?bsub for details.",
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
      "  it should be a dataset of possible substitution", 
      "  and can be computed using multilevelcoda::possub.", 
      "  See ?bsub for details.",
      sep = "\n"))
  }
  
  if(isFALSE(missing(delta))) {
    if (isFALSE(is.integer(delta))) {
      if (isFALSE(delta > 0)) {
        stop(" 'delta' must be an positive integer value.")
      }
    }
  } else if (isTRUE(missing(delta))){
    stop(paste(
      "'delta' is a required argument and cannot be missing;",
      "  it should be interger, numeric positive value or vector", 
      "  to specify the change in units across compositional parts", 
      sep = "\n"))
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
      "The names of compositional variables must be the 
  same in 'substitute' (%s) and 'parts' (%s).",
      colnames(substitute),
      object$CompIlr$parts))
  }
  
  if (isFALSE(is.null(regrid))) {
    if(any(c(colnames(object$CompIlr$BetweenILR), colnames(object$CompIlr$WithinILR))
           %in% c(colnames(regrid)))) {
      stop(paste(
        "'regrid' should not have any column names starting with 'bilr', 'wilr', or 'ilr'.",
        "  These variables will be calculated by substitution model.",
        "  Reference grid should contain information about the covariates used in 'brmcoda'.",
        "  Please provide a different reference grid.",
        sep = "\n"))
    }
  }

  if ("between" %in% level) {
    if("conditional" %in% type) {
      bout <- bsub(object = object, substitute = substitute, delta = delta, 
                regrid = regrid, summary = summary, 
                level = level, type = type)
    }
    if("marginal" %in% type) {
      bmout <- bsubmargins(object = object, substitute = substitute, delta = delta,
                           level = level, type = type)
    }
  }
  
  if ("within" %in% level) {
    if("conditional" %in% type) {
      wout <- wsub(object = object, substitute = substitute, delta = delta, 
                  regrid = regrid, summary = summary,
                  level = level, type = type)
    }
    if("marginal" %in% type) {
      wmout <- wsubmargins(object = object, substitute = substitute, delta = delta,
                           level = level, type = type)
    }
  }
  
  out <- list(BetweenpersonSub = if(exists("bout")) (bout) else (NULL),
              WithinpersonSub = if(exists("wout")) (wout) else (NULL),
              BetweenpersonSubMargins = if(exists("bmout")) (bmout) else (NULL),
              WithinpersonSubMargins = if(exists("wmout")) (wmout) else (NULL))
  out
}