#' Fit Bayesian generalised (non-)linear multivariate multilevel compositional model 
#' via full Bayesian inference using Stan,
#' when composition is the outcome.
#'
#'
#' This function utilises 
#' @param formula A object of class \code{\link{brmsformula}}, \code{\link{mvbrmsformula}}:
#' A symbolic description of the model to be fitted. 
#' Details of the model specification can be found in \code{\link{mvbrmsformula}}.
#' @param compilr A \code{compilr} object containing data of composition, ILR coordinates,
#' and other variables used in the model.
#' @param ... Further arguments passed to \code{\link{brm}}.
#' 
#' @return 
#' @importFrom reshape2 melt
#' @importFrom brms brm
#' @export
#' @examples
#' data(mcompd)
#' cilr <- compilr(data = mcompd, sbp = sbp, idvar = "ID")
#' 
#' ## inspect names of ILR coordinates bfore passing to 'brm' model
#' names(cilr$TotalILR)
#' 
#' ## run mvmcoda
#' mvm1 <- mvmcoda(compilr = cilr, formula = mvbind(ilr1, ilr2, ilr3, ilr4) ~ STRESS + (1 | ID))
#' 
mvmcoda <- function(compilr, formula, ...) {
  
  message("Please be patient! 'mvmcoda' is working hard...")
  start_time <- Sys.time()
  
  tilr <- compilr$TotalILR

  tmpd <- cbind(data, tilr)
  
  m <- brm(eval(formula),
           data = tmpd,
           ...)
  
  out <- list(
    CompIlr = compilr,
    BrmModel = m)
  
  end_time <- Sys.time()
  
  message(paste(
  sprintf("'mvmcoda' took %d to do hard work.", end_time - start_time),
          "If you would like to run our model faster,",
          "consider using parallelisation with 'brm'",
          "See ?brms::brm for details.",
  sep = "\n"))
  
}
