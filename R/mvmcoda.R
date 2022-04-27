#' Fit Bayesian generalised (non-)linear multivariate multilevel compositional model 
#' via full Bayesian inference using brms,
#' when composition is the outcome.
#'
#' @param formula A object of class \code{\link{brmsformula}}, \code{\link{mvbrmsformula}}:
#' A symbolic description of the model to be fitted. 
#' Details of the model specification can be found in \code{\link{mvbrmsformula}}.
#' @param compilr A \code{\link{compilr}} object containing data of composition, 
#' ILR coordinates, and other variables used in the model.
#' @param ... Further arguments passed to \code{\link{brm}}.
#' 
#' @return 
#' @importFrom brms brm
#' @export
#' @examples
#' \donttest{
#' data(mcompd)
#' data(sbp)
#' cilr <- compilr(data = mcompd, sbp = sbp, idvar = "ID",
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' 
#' ## inspect names of ILR coordinates bfore passing to 'brm' model
#' names(cilr$TotalILR)
#' 
#' ## run mvmcoda
#' mvm <- mvmcoda(compilr = cilr, formula = mvbind(ilr1, ilr2, ilr3, ilr4) ~ STRESS + (1 | ID),
#'                cores = 8, chain = 4)
#' }
mvmcoda <- function(compilr, formula, ...) {
  
  if (isTRUE(missing(compilr))) {
    stop(paste(
      "'compilr' is a required argument and cannot be missing;",
      " it should be an object of class compilr.", 
      " See ?multilevelcoda::compilr for details.",
      sep = "\n"))
  }
  
  if (isFALSE(inherits(compilr, "compilr"))) {
    stop(paste(
      "compilr must be an object of class compilr.",
      " See ?multilevelcoda::compilr for details.",
      sep = "\n"))
  }
  message("Please be patient, 'mvmcoda' is working hard...")

  tmpd <- cbind(compilr$data, compilr$TotalILR)
  m <- brm(formula,
           data = tmpd,
           ...)
  
  out <- structure(
    list(CompIlr = compilr,
         BrmModel = m),
    class = "mvmcoda")
  out
}
