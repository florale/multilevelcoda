#' A series of functions for single level multivariate compositional models
#' 
#' Fit Bayesian generalised (non-)linear multivariate compositional model 
#' via full Bayesian inference using Stan,
#' when composition is the outcome.
#'
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
#' names(cilr$BetweenILR)
#' 
#' ## run mvcoda
#' mv1 <- mvcoda(compilr = cilr, formula = mvbind(ilr1, ilr2, ilr3, ilr4) ~ STRESS + Age + Female)
#' 
mvcoda <- function(data, composition, sbp, formula, ...) {

  if (isFALSE(inherits(data, c("data.table", "data.frame", "matrix")))) {
    stop("'data' must be a data table, data frame or matrix.")
  }
  if (isFALSE(inherits(sbp, "matrix"))) {
    stop(sprintf("'sbp' is a '%s' but must be a matrix.",
                 paste(class(sbp), collapse = ";")))
  }
  if (isFALSE(identical(length(composition), ncol(sbp)))) {
    stop(sprintf("The number of compositional variables in 'composition' (%d) 
                 must be the same as in 'sbp' (%d).",
                 length(composition),
                 ncol(sbp)))
  }

  # compilr
  tmp <- as.data.table(data)
  psi <- gsi.buildilrBase(t(sbp))
  
  comp <- acomp(tmp[, composition, with = FALSE])
  ilr <- ilr(comp, V = psi)
  
  colnames(ilr) <- paste0("ilr", seq_len(ncol(ilr)))
  
  tmp <- cbind(tmp, ilr)
  
  # brm model
  m <- brm(eval(formula),
           data = tmp,
           ...)
  
  out <- list(BrmsResults = m,
              ILR = ilr,
              Comp = comp,
              data = tmp,
              psi = psi,
              sbp = sbp,
              composition = composition)
  
  return(out)
}
