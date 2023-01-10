#' @title Substitution Model.
#' 
#' @description 
#' Estimate the difference in an outcome
#' when compositional parts are substituted for specific unit(s). 
#' The \code{substitution} output encapsulates
#' the substitution results for all compositional parts
#' present in the \code{\link{brmcoda}} object.
#' 
#' @param object A fitted \code{\link{brmcoda}} object. Required.
#' @param delta A integer, numeric value or vector indicating the amount of substituted change between compositional parts.
#' @param basesub A \code{data.frame} or \code{data.table} of the base possible substitution of compositional parts.
#' This data set can be computed using function \code{\link{basesub}}. 
#' If \code{NULL}, all possible pairwise substitution of compositional parts are used.
#' @param regrid If non-\code{NULL}, a \code{data.table} of reference grid consisting 
#' of combinations of covariates over which predictions are made.
#' If \code{NULL}, the reference grid is constructed via \code{\link{ref_grid}}.
#' @param summary A logical value. 
#' Should the estimate at each level of the reference grid (\code{FALSE}) 
#' or their average (\code{TRUE}) be returned? Default to \code{TRUE}.
#' @param level A character string or vector. 
#' Should the estimate be at the \code{between}-person and/or \code{within}-person level? Required.
#' @param type A character string or vector. 
#' Should the estimate be \code{conditional} mean or average \code{marginal} mean? Required.
#' @param ... Additional arguments to be passed to \code{\link{describe_posterior}}.
#' 
#' @return A list containing the result of multilevel compositional substitution model.
#' Each element of the list is the estimation for a compositional part 
#' and include at least six elements.
#' \itemize{
#'   \item{\code{Mean}}{ Posterior means.}
#'   \item{\code{CI_low}} and \item{\code{CI_high}}{ 95% credible intervals.}
#'   \item{\code{Delta}}{ Amount substituted across compositional parts.}
#'   \item{\code{From}}{ Compositional part that is substituted from.}
#'   \item{\code{To}}{ Compositional parts that is substituted to.}
#'   \item{\code{Level}}{Level where changes in composition takes place. Either }
#'   \item{\code{EffectType}}{Either estimated `conditional` or average `marginal` changes.}
#' }
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo mean.acomp
#' @importFrom extraoperators %snin% %sin%
#' @importFrom insight find_predictors
#' @importFrom emmeans ref_grid
#' @importFrom stats fitted
#' @export
#' @examples
#' \donttest{
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
#' subm <- substitution(object = m, delta = c(1, 10),
#'                      type = "conditional", level = c("between", "within"))
#' }
substitution <- function(object, delta, basesub = NULL,
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
  
  if(isTRUE(missing(basesub))) {
    count <- length(object$CompIlr$parts)
    n <- count - 2
    
    subvars1 <- c(1, -1)
    subvars2 <- rep(0, n)
    subvars <- c(subvars1, subvars2)
    
    nc <- length(subvars)
    nr <- (nc - 1) * count
    
    basesub <- matrix(0, nrow = nr, ncol = nc, dimnames = list(NULL, object$CompILR$parts))
    k <- 0
    
    for(i in 1:nc)
      for(j in 1:nc)
        if(i != j) {
          k <- k + 1
          basesub[k, c(i, j)] <- c(1, -1)
        }
    
    basesub <- as.data.table(basesub)
    names(basesub) <- object$CompIlr$parts
    } else if(isFALSE(missing(basesub))) {
      
      if (isFALSE(identical(ncol(basesub), length(object$CompIlr$parts)))) {
        stop(sprintf(
        "The number of columns in 'basesub' (%d) must be the same
        as the compositional parts in 'parts' (%d).",
        ncol(basesub),
        length(object$CompIlr$parts)))
        }
      
      if (isFALSE(identical(colnames(basesub), object$CompIlr$parts))) {
        stop(sprintf(
        "The names of compositional parts must be the
        same in 'basesub' (%s) and 'parts' (%s).",
        colnames(basesub),
        object$CompIlr$parts))
      }
    }
  
  if ("between" %in% level) {
    if("conditional" %in% type) {
      bout <- bsub(object = object, basesub = basesub, delta = delta,
                   regrid = regrid, summary = summary,
                   level = "between", type = "conditional")
    }
    if("marginal" %in% type) {
      bmout <- bsubmargins(object = object, basesub = basesub, delta = delta,
                           level = "between", type = "marginal")
    }
  }
  
  if ("within" %in% level) {
    if("conditional" %in% type) {
      wout <- wsub(object = object, basesub = basesub, delta = delta,
                   regrid = regrid, summary = summary,
                   level = "within", type = "conditional")
    }
    if("marginal" %in% type) {
      wmout <- wsubmargins(object = object, basesub = basesub, delta = delta,
                           level = "within", type = "marginal")
    }
  }
  
  out <- list(BetweenpersonSub = if(exists("bout")) (bout) else (NULL),
              WithinpersonSub = if(exists("wout")) (wout) else (NULL),
              BetweenpersonSubMargins = if(exists("bmout")) (bmout) else (NULL),
              WithinpersonSubMargins = if(exists("wmout")) (wmout) else (NULL))
  out
}