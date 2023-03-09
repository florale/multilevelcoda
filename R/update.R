#' Update \code{\link{compilr}} 
#' 
#' This method allows for updating an existing \code{\link{compilr}} object.
#' 
#' @param compilr A \code{\link{compilr}} object containing data of composition, 
#' ILR coordinates, and other variables to be updated.
#' @param newdata A \code{data.frame} or \code{data.table}
#' containing data of all variables used in the analysis. 
#' It must include a composition and the same ID variable as the existing \code{\link{compilr}} object.
#'
#' @return A updated \code{\link{compilr}} object with twelve elements.
#' \itemize{
#'   \item{\code{BetweenComp}}{ A vector of class \code{acomp} representing one closed between-person composition
#'   or a matrix of class \code{acomp} representing multiple closed between-person compositions each in one row.}
#'   \item{\code{WithinComp}}{ A vector of class \code{acomp} representing one closed within-person composition
#'   or a matrix of class \code{acomp} representing multiple closed within-person compositions each in one row.}
#'   \item{\code{TotalComp}}{ A vector of class \code{acomp} representing one closed total composition
#'   or a matrix of class \code{acomp} representing multiple closed total compositions each in one row.}
#'   \item{\code{BetweenILR}}{ Isometric log ratio transform of between-person composition.}
#'   \item{\code{WithinILR}}{ Isometric log ratio transform of within-person composition.}
#'   \item{\code{TotalILR}}{ Isometric log ratio transform of total composition.}
#'   \item{\code{data}}{ The user's dataset or imputed dataset if the input data contains zeros.}
#'   \item{\code{psi}}{ A ILR matrix associated with user-defined partition structure.}
#'   \item{\code{sbp}}{ The user-defined sequential binary partition matrix.}
#'   \item{\code{parts}}{ Names of compositional variables.}
#'   \item{\code{idvar}}{ Name of the variable containing IDs.}
#'   \item{\code{total}}{ Total amount to which the compositions is closed.}
#' }
#' @exportS3Method update compilr
#' @examples 
#' data(mcompd)
#' data(sbp)
#' cilr <- compilr(data = mcompd, sbp = sbp, 
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#' 
#' newdat <- mcompd[ID != 1] # excluding ID 1
#' cilr_new <- update(compilr = cilr, newdata = newdat)
#' 
#' ## cleanup
#' rm(cilr, mcompd, sbp, cilr_new, newdat)
#' 
update.compilr <- function(compilr, newdata) { 
  # add checks
  
  sbp <- compilr$sbp
  parts <- compilr$parts
  total <- compilr$total
  idvar <- compilr$idvar
  
  compilr(newdata, sbp, parts, total, idvar)
  
}

#' Update \code{\link{brmcoda}} models
#' 
#' This method allows for updating an existing \code{\link{brmcoda}} object.
#' 
#' Details of the model specification can be found in \code{\link{brmsformula}}.
#' @param compilr A \code{\link{compilr}} object containing data of composition, 
#' ILR coordinates, and other variables used in the model.
#' @param newdata A \code{data.frame} or \code{data.table}
#' containing data of all variables used in the analysis. 
#' It must include a composition and the same ID variable as the existing \code{\link{compilr}} object.
#' @param ... Further arguments passed to \code{\link{update}}.
#' 
#' @exportS3Method update brmcoda
#' @examples 
#' \dontrun{
#' data(mcompd)
#' data(sbp)
#' cilr <- compilr(data = mcompd, sbp = sbp, 
#'         parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#' 
#' # model with compositional predictor at between and within-person levels
#' fit <- brmcoda(compilr = cilr, 
#'               formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                                  wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID), 
#'               chain = 1, iter = 500,
#'               backend = "cmdstanr")
#'
#' newdat <- mcompd[ID != 1] # excluding ID 1
#' fit_new <- update(fit, newdata = newdat)
#' }
update.brmcoda <- function(model,
                           newcilr = NULL, newdata = NULL, ...) {
  
  if (is.null(newdata)) {
    stop2("'newdata' is required when updating a 'brmcoda' object.")
  }
  
  if(!is.null(newdata)) {
    newcilr <- update(model$CompIlr, newdata = newdata)
  }
  
  newdata <- cbind(newcilr$data, newcilr$BetweenILR, 
                   newcilr$WithinILR, newcilr$TotalILR)
  
  fit_new <- update(model$Model, newdata = newdata, ...)
  
  structure(
    list(CompIlr = newcilr,
         Model = fit_new),
    class = "brmcoda")
}
