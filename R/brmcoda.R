#' Fit Bayesian generalised (non-)linear multilevel compositional model 
#' via full Bayesian inference
#' 
#' Fit a \code{brm} model with multilevel ILR coordinates
#' 
#' @param formula A object of class \code{formula}, \code{brmsformula}:
#' A symbolic description of the model to be fitted. 
#' Details of the model specification can be found in \code{\link{brmsformula}}.
#' @param compilr A \code{\link{compilr}} object containing data of composition, 
#' ILR coordinates, and other variables used in the model.
#' @param ... Further arguments passed to \code{\link{brm}}.
#' 
#' @return A \code{\link{brmcoda}} with two elements
#' \itemize{
#'   \item{\code{CompILR}}{ An object of class \code{compilr} used in the \code{brm} model. }
#'   \item{\code{Model}}{ An object of class \code{brmsfit}, which contains the posterior draws 
#'   along with many other useful information about the model.}
#'   }
#' @importFrom brms brm
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   cilr <- compilr(data = mcompd, sbp = sbp,
#'                   parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#'   
#'   # inspects ILRs before passing to brmcoda
#'   names(cilr$BetweenILR)
#'   names(cilr$WithinILR)
#'   names(cilr$TotalILR)
#'   
#'   # model with compositional predictor at between and within-person levels
#'   m1 <- brmcoda(compilr = cilr,
#'                 formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                   wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   # model with compositional outcome
#'   m2 <- brmcoda(compilr = cilr,
#'                 formula = mvbind(ilr1, ilr2, ilr3, ilr4) ~ Stress + Female + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   }}
#' @export
brmcoda <- function (compilr, formula, ...) {

  if (isTRUE(missing(compilr))) {
    stop(paste(
      "'compilr' is a required argument and cannot be missing;",
      "  it should be an object of class compilr.", 
      "  See ?multilevelcoda::compilr for details.",
      sep = "\n"))
  }
  if (isFALSE(inherits(compilr, "compilr"))) {
    stop(paste(
      "compilr must be an object of class compilr.",
      "  See ?multilevelcoda::compilr for details.",
      sep = "\n"))
  }
  
  tmp <- cbind(compilr$data, compilr$BetweenILR, 
               compilr$WithinILR, compilr$TotalILR)
  
  m <- brm(formula,
           data = tmp,
           ...)
  
  structure(
    list(CompILR = compilr,
         Model = m),
    class = "brmcoda")
}
