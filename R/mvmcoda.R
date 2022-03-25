#' A series of functions for multilevel compositional data analysis with
#' composition as the outcome.
#'
#' This provides a few sentence description about the example function.
#'
#' @param ... Further arguments passed to \code{\link{brm}}.
#' @return 
#' @importFrom compositions acomp
#' @importFrom reshape2 melt
#' @importFrom brms brm
#' @export
#' @examples
#' data(mcompd)
#' cilr <- compilr(data = mcompd, sbp = sbp, idvar = "ID")
#' 
#' mvm1 <- mvmcoda(compilr = cilr, formula = value ~ variable + variable:STRESS - 1 + (variable -1| ID))
mvmcoda <- function(compilr, formula, ...) {
  
  message("Please be patient! 'mvmcoda' is working hard...")
  start_time <- Sys.time()
  
  tilr <- compilr$TotalILR

  tmpd <- cbind(data, tilr)
  
  vn <- colnames(tmpd) %snin% names(tilr)
  ilrn <- colnames(tmpd) %sin% names(tilr)
  
  sd <- melt(tmpd, id.vars = vn,
             measure.vars = ilrn)
  
  m <- brm(eval(formula),
           data = sd,
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
