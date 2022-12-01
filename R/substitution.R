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
#' for which substitution model is desired.
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
#'   \item{\code{MinSubstituted}}{ Minute substituted within the composition.}
#'   \item{\code{Substitute}}{Compositional variables to be substituted from/to.}
#'   \item{\code{Predictor}}{Central compositional variable to be substituted from/to.}
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
#' subm <- substitution(object = m, substitute = psub, 
#'                      type = "conditional", level = "between")
#' }
substitution <- function(object, substitute, delta, 
                         regrid = NULL, summary = TRUE, 
                         level = c("between", "within"),
                         type = c("conditional", "marginal"),
                         ...) {
  
  if ("between" %in% level) {
    if("conditional" %in% type) {
    bout <- bsub(object = object, substitute = substitute, delta = delta, 
                regrid = regrid, summary = summary)
    }
    if("marginal" %in% type) {
      bmout <- bsubmargins(object = object, substitute = substitute, delta = delta)
    }
  }
  
  if ("within" %in% level) {
    if("conditional" %in% type) {
      wout <- wsub(object = object, substitute = substitute, delta = delta, 
                  regrid = regrid, summary = summary)
    }
    if("marginal" %in% type) {
      wmout <- wsubmargins(object = object, substitute = substitute, delta = delta)
    }
  }
  out <- list(BetweenpersonSub = if(exists("bout")) (bout) else (NULL),
              WithinpersonSub = if(exists("wout")) (wout) else (NULL),
              BetweenpersonSubMargins = if(exists("bmout")) (bmout) else (NULL),
              WithinpersonSubMargins = if(exists("wmout")) (wmout) else (NULL))
  out
}