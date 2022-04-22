#' Fit Bayesian generalised (non-)linear multilevel compositional model 
#' via full Bayesian inference using Stan,
#' when composition is the predictor.
#' 
#' This function fits a \code{brm} model to
#' between-person and within-person ILR coordinates as predictor.
#' 
#' @param formula A object of class \code{formula}, \code{brmsformula}:
#' A symbolic description of the model to be fitted. 
#' Details of the model specification can be found in \code{\link{brmsformula}}.
#' @param compilr A \code{compilr} object containing data of composition, ILR coordinates,
#' and other variables used in the model.
#' @param ... Further arguments passed to \code{\link{brm}}.
#' 
#' @return A list with two elements
#' \itemize{
#'   \item{\code{CompIlr}}{ A object of class \code{compilr} used in the \code{brm} model. }
#'   \item{\code{BrmModel}}{ An object of class \code{brmsfit}, which contains the posterior draws 
#'   along with many other useful information about the model.}
#'   
#' @importFrom brms brm
#' @export
#' @examples
#' data(mcompd)
#' data(sbp)
#'
#' ## compute compositions and ILR coordinates
#' cilr <- compilr(data = mcompd, sbp = sbp, parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#' 
#' ## inspect names of ILR coordinates bfore passing to 'brm' model
#' names(cilr$BetweenILR)
#' names(cilr$WithinILR)
#' 
#' ## run brmcoda model
#' # unadjusted model
#' brmcodatest <- brmcoda(compilr = cilr, 
#'                formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID), 
#'                core = 8, chain = 4)
#' 
#' # adjusted model
#' adjbrmcodatest <- brmcoda(compilr = cilr, 
#'                formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID), 
#'                core = 8, chain = 4)
#' print(mcm$BrmModel)
#' 
#' ## clean-up
#' rm(mcompd, sbp, cilr, mcm, mcmc)
brmcoda <- function (formula, compilr, ...) {
  if (isFALSE(inherits(compilr, "compilr"))) {
    stop("compilr must be an object of class compilr.")
  }
  tmp <- cbind(compilr$data, compilr$BetweenILR, compilr$WithinILR)
  
  m <- brm(formula,
           data = tmp,
           ...)
  
  out <- structure(
    list(
      CompIlr = compilr,
      BrmModel = m),
    class = "brmcoda")
  
  out

}