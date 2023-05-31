#' @title Between-person Basic Substitution Model.
#' 
#' @description 
#' Estimate the difference in outcomes
#' when compositional parts are substituted for specific unit(s) at `between-person` level. 
#' The \code{bsub} output encapsulates
#' the substitution results for all compositional parts
#' present in the \code{\link{brmcoda}} object.
#' 
#' Notes: The reference composition for substitution model 
#' is the compositional mean of the dataset provided.
#' For average marginal effect, use \code{\link{bsubmargins}}.
#'
#' @param object A fitted \code{\link{brmcoda}} object. Required.
#' @param delta A integer, numeric value or vector indicating the amount of substituted change between compositional parts.
#' @param basesub A \code{data.frame} or \code{data.table} of the base possible substitution of compositional parts.
#' This data set can be computed using function \code{\link{basesub}}. 
#' If \code{NULL}, all possible pairwise substitution of compositional parts are used.
#' @param summary A logical value. 
#' Should the estimate at each level of the reference grid (\code{FALSE}) 
#' or their average (\code{TRUE}) be returned? Default to \code{TRUE}.
#' @param level A character string or vector. 
#' Should the estimate be at the \code{between}-person and/or \code{within}-person level? Required.
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
#'   \item{\code{Level}} { Level where changes in composition takes place.}
#' }
#' 
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
#' 
#' cilr <- compilr(data = mcompd, sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#' 
#' # model with compositional predictor at between and between-person levels
#' m <- brmcoda(compilr = cilr, 
#'              formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + 
#'                                 wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID), 
#'              chain = 1, iter = 500,
#'              backend = "cmdstanr")
#' subm <- bsub(object = m, basesub = psub, delta = 5)
#' }
bsub <- function(object,
                 basesub,
                 delta,
                 summary = TRUE,
                 ref = "unitmean",
                 level = "between",
                 weight = c("equal", "proportional"),
                 refdata = NULL,
                 ...) {
  
  # input for substitution model
  ID <- 1 # to make fitted() happy
  delta <- as.integer(delta)
  
  # error if delta out of range
  if(isTRUE(any(delta > min(refcomp)))) {
    stop(sprintf(
      "delta value should be less than or equal to %s, which is
  the amount of composition part available for pairwise substitution.",
  round(min(refcomp), 2)
    ))
  }
  
  # d0 -------------------------------
  if (isTRUE(is.null(refdata))) {
    d0 <- refdata(object = object,
                  ref = ref,
                  weight = weight,
                  build.rg = FALSE)
  } else {
    d0 <- refdata
  }
  
  # y0 --------------------------------
  y0 <- fitted(
    object$Model,
    newdata = d0,
    re_formula = NA,
    summary = FALSE)
  
  # yb ---------------------------------
  out <- get.bsub(
    object = object,
    basesub = basesub,
    refcomp = refcomp,
    delta = delta,
    y0 = y0,
    d0 = d0,
    summary = summary,
    covnames = covnames,
    refgrid = refgrid,
    level = level,
    type = type)
  
}