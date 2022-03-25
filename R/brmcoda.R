#' Computing brms model for between-person and within-person multilevel compositional data analysis,
#' with composition as the predictor.
#' 
#' @param formula A object of class \code{formula}, \code{brmsformula}:
#' A symbolic description of the model to be fitted.
#' @param compilr A \code{compilr} object containing composition, ILR coordinates,
#' used in the model.
#' @param ... Further arguments passed to \code{\link{brm}}.
#' 
#' @return A list with two elements
#' \itemize{
#'   \item{\code{CompIlr}}{ A object of class \code{compilr} used in the \code{brm} model. }
#'   \item{\code{Results}}{ An object of class \code{brmsfit}, which contains the posterior draws 
#'   along with many other useful information about the model.}
#'   
#' @importFrom brms brm
#' @export
#' @examples
#' data(mcompd)
#' data(sbp)
#'
#' ## compute compositions and ILR coordinates
#' cilr <- compilr(data = mcompd, sbp = sbp, composition = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#' 
#' ## inspect names of ILR coordinates bfore passing to 'brm' model
#' names(cilr$BetweenILR)
#' names(cilr$WithinILR)
#' 
#' ## run brmcoda model
#' mcm <- brmcoda(compilr = cilr, 
#'                formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + wilr1 + wilr2 + wilr3 + wilr4 + Female + Age + (1 | ID), 
#'                core = 8, chain = 4)
#' 
#' print(mcm$BrmModel)
#' 
#' ## clean-up
#' rm(mcompd, sbp, cilr, mcm)
brmcoda <- function (formula, compilr, ...) {

  tmp <- cbind(compilr$data, compilr$BetweenILR, compilr$WithinILR)
  
  m <- brm(formula,
           data = tmp,
           ...)
  
  out <- list(
    CompIlr = compilr,
    BrmModel = m)

}