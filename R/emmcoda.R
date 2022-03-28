#' #' @title Estimated marginal means of composition.
#' #' 
#' #' @description 
#' #' This function obtains estimated marginal means (EMMs) of models 
#' #' containing composition as the outcome
#' #' 
#' #' @param data
#' #' 
#' #' @return
#' #' @importFrom compositions ilrInv
#' #' @importFrom emmeans ref_grid
#' #' @importFrom data.table as.data.table
#' #' @export
#' #' @examples
#' #' 
#' emmcoda <- function (data, ...) {
#' 
#' tmp <- copy(data)
#' 
#' }