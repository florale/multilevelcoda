#' @title Estimated marginal means of composition.
#'
#' @description
#' This function obtains estimated marginal means (EMMs) of models
#' containing composition as the outcome
#'
#' @param data A fitted \code{\link{mvmcoda}} model data. Required.
#' @param specs A character vector specifying the names of the predictors over which EMMs are desired.
#' @param at A numeric vector passed to \code{\link{ref_grid}} specifying the predictor's values for EMMs.
#' @param ... Further arguments passed to \code{\link{emmeans}}.
#'
#' @return A list with four elements.
#' @importFrom compositions ilrInv
#' @importFrom emmeans emmeans
#' @importFrom data.table as.data.table
#' @importFrom reshape2 melt dcast
#' @export
#' @examples
#' 
#' emtest <- emmcoda(object = mvm1, specs = "STRESS", at = c(1:10))
#'
emmcoda <- function (object, specs, at, ...) {

  tmp <- copy(object)
  
  ls <- list(at)
  names(ls) <- specs    

  em <- emmeans(tmp$BrmModel, pairwise ~ ilr, specs, mult.name = "ilr", at = ls, ...)
  emd <- as.data.table(em$emmeans)
  
  ilr <- emd[, 1:3]
  ilr <- spread(ilr, key = "ilr", value = "emmean")
  
  ilrinv <- ilrInv(ilr[, 2:ncol(ilr)], V = tmp$CompIlr$psi)
  ilrinv <- clo(ilrinv, total = 1440)
  ilrinv <- as.data.table(ilrinv)
  ilrinv <- cbind(ilr[, 1], ilrinv)
  colnames(ilrinv) <- c(specs, tmp$CompIlr$composition)
  
  ilrinvl <- melt(ilrinv, id = specs)
  
  out <- list(
    emmILR = emd,
    emmWide = ilrinv,
    emmLong = ilrinvl,
    Predictor = specs)
  
  return(out)
}
