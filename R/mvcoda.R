#' A series of functions for single level multivariate compositional models
#' 
#' Fit Bayesian generalised (non-)linear multivariate compositional model 
#' via full Bayesian inference using Stan,
#' when composition is the outcome.
#' 
#' @param data A \code{data.frame} or \code{data.table}
#' containing data of all variables used in the analysis. 
#' It must include a composition and a ID variable. Required.
#' @param formula A object of class \code{\link{brmsformula}}, \code{\link{mvbrmsformula}}:
#' A symbolic description of the model to be fitted. 
#' @param sbp A signary matrix indicating sequential binary partition. Required.
#' Details of the model specification can be found in \code{\link{mvbrmsformula}}.
#' @param parts A character vector specifying the names of compositional variables. Required.
#' @param ... Further arguments passed to \code{\link{brm}}.
#' 
#' @return A list of results. TODO
#' @importFrom compositions ilr acomp gsi.buildilrBase
#' @importFrom data.table copy as.data.table :=
#' @importFrom reshape2 melt
#' @importFrom brms brm
#' @export
#' @examples
#' data(mcompd)
#' data(sbp)
#' ## run mvcoda
#' mv <- mvcoda(mcompd, mvbind(ilr1, ilr2, ilr3, ilr4) ~ STRESS + Age + Female,
#'               sbp, c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' 
mvcoda <- function(data, formula, sbp, parts, ...) {

  if (isFALSE(inherits(data, c("data.table", "data.frame", "matrix")))) {
    stop("'data' must be a data table, data frame or matrix.")
  }
  if (isFALSE(inherits(sbp, "matrix"))) {
    stop(sprintf("'sbp' is a '%s' but must be a matrix.",
                 paste(class(sbp), collapse = ";")))
  }
  if (isFALSE(identical(length(parts), ncol(sbp)))) {
    stop(sprintf("The number of compositional variables in 'parts' (%d) 
                 must be the same as in 'sbp' (%d).",
                 length(parts),
                 ncol(sbp)))
  }

  # compilr
  tmp <- as.data.table(data)
  psi <- gsi.buildilrBase(t(sbp))
  
  comp <- acomp(tmp[, parts, with = FALSE])
  ilr <- ilr(comp, V = psi)
  
  colnames(ilr) <- paste0("ilr", seq_len(ncol(ilr)))
  
  if(any(c(colnames(ilr)) %in% colnames(tmp))) {
    stop(paste("data should not have any column names starting with 'ilr';",
               "because these variables will be used in subsequent models.",
               "Please rename them before running 'mvcoda'.",
               sep = "\n"))
  }
  
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
              parts = parts)
  
  return(out)
}
