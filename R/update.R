#' Update \code{\link{compilr}} 
#' 
#' This method allows for updating an existing \code{\link{compilr}} object.
#' 
#' @param object A \code{\link{compilr}} class object to be updated.
#' @param newdata A \code{data.frame} or \code{data.table}
#' containing data of all variables used in the analysis. 
#' It must include a composition and the same ID variable as the existing \code{\link{compilr}} object.
#' @param ... generic argument, not in use.
#' 
#' @inherit compilr return
#' 
#' @seealso \code{\link{compilr}}
#' 
#' @importFrom extraoperators %ain% %snin% %nin%
#' @method update compilr
#' 
#' @examples 
#' cilr <- compilr(data = mcompd, sbp = sbp, 
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#' 
#' # update with new data
#' newdat <- mcompd[ID != 1] # excluding ID 1
#' cilr1 <- update(object = cilr, newdata = newdat)
#' @export
update.compilr <- function(object, newdata, ...) { 
  
  if (isTRUE(missing(newdata))) {
    stop("'newdata' is required when updating a 'compilr' object.")
    } else {
      if (isFALSE(inherits(newdata, c("data.table", "data.frame", "matrix")))) {
        stop("newdata must be a data table, data frame or matrix.")
      }
    }
  
  if (isFALSE(object$parts %ain% colnames(newdata))) { # check compositional variables
    stop(sprintf(
      "The folllowing compositional parts could not be found in 'newdata'  (%s).",
      paste(object$parts %snin% colnames(newdata), collapse = ", ")))
  }
  
  if (object$idvar %nin% colnames(newdata)) { # check ID
    stop(sprintf(
      "The names of the ID variable must be the \n same in 'object' (%s) and 'newdata'.",
      object$idvar))
  }
  
  # update compilr
  sbp <- object$sbp
  parts <- object$parts
  total <- object$total
  idvar <- object$idvar
  
  compilr(newdata, sbp, parts, total, idvar)
}

#' Update \code{\link{brmcoda}} models
#' 
#' This method allows for updating an existing \code{\link{brmcoda}} object.
#' 
#' @param object A fitted \code{\link{brmcoda}} object to be updated.
#' @param formula. Changes to the formula; for details see \code{\link{update.formula}} and \code{\link{brmsformula}}.
#' @param newcilr A \code{\link{compilr}} object containing data of composition, 
#' ILR coordinates, and other variables used in the updated model.
#' @param newdata A \code{data.frame} or \code{data.table}
#' containing data of all variables used in the analysis. 
#' It must include a composition and the same ID variable as the existing \code{\link{compilr}} object.
#' @param ... Further arguments passed to \code{\link{brm}}.
#' 
#' @inherit brmcoda return
#'  
#' @seealso \code{\link{brmcoda}}
#' 
#' @method update brmcoda
#' 
#' @examples 
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#' 
#' # model with compositional predictor at between and within-person levels
#' fit <- brmcoda(compilr = compilr(data = mcompd, sbp = sbp, 
#'                                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), 
#'                                  idvar = "ID"
#'                                  ), 
#'               formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                                  wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID), 
#'               chain = 1, iter = 500,
#'               backend = "cmdstanr")
#'
#' # removing the effect of wilr1
#' fit1 <- update(fit, formula. = ~ . - wilr1)
#' 
#' # using only a subset
#' fit2 <- update(fit, newdata = mcompd[ID != 1])
#' }}
#' @export
update.brmcoda <- function(object,
                           formula. = NULL,
                           newdata = NULL,
                           newcilr = NULL,
                           ...) {
  
  if(!is.null(newcilr) && !is.null(newdata)) {
    warning(paste(
      "Either 'newcilr' or 'newdata' is required to update brmcoda,",
      "  but both were supplied.",
      "  'newcilr' will be used to update brmcoda and 'newdata' will be ignored.",
      sep = "\n"))
  }
  
  # check args
  if (is.null(formula.) && is.null(newcilr) && is.null(newdata)) {
    stop("either 'formula.', 'newdata', or 'newcilr' is required when updating a 'brmcoda' object.")
  }
  
  # updating only formula
  if (!is.null(formula.) && is.null(newcilr) && is.null(newdata)) {
    fit_new <- update(object$Model, formula. = formula., ...)
    
    structure(
      list(CompILR = object$CompILR,
           Model = fit_new),
      class = "brmcoda")
    
  } else { # updating newdata/newcilr (with and without updating formula)
    
    if(is.null(newcilr)) {
      if(!is.null(newdata)) {
        newcilr <- update(object$CompILR, newdata = newdata)
      }} else {
        if(isFALSE(inherits(newcilr, "compilr"))) {
          stop("'newcilr' must be an object of class 'compilr'.")
        }
      }
    
    newdata <- cbind(newcilr$data,
                     newcilr$BetweenILR,
                     newcilr$WithinILR,
                     newcilr$TotalILR)
    
    fit_new <- update(object$Model, formula. = formula., newdata = newdata, ...)
    
    structure(
      list(CompILR = newcilr,
           Model = fit_new),
      class = "brmcoda")
    
  }
}