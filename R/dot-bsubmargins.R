#' #' @title Margianl Effects of Between-person Isotemporal Substitution.
#' #'
#' #' @description
#' #' This function estimates the difference in outcomes
#' #' when compositional variables are substituted for a specific period
#' #' at between-person level. The model loops through all compositional variables present
#' #' in the \code{\link{brmcoda}} object.
#' #'
#' #' @param object A fitted \code{\link{brmcoda}} model data. Required.
#' #' @param substitute A data frame or data table indicating the possible substitution of variables.
#' #' This dataset can be computed using \code{possub}. Required.
#' #' @param minute A integer or numeric value indicating the maximum minute for which substitution model is desired.
#' #' Default to 60L. In this case, the model loops through 1:60L minutes.
#' #'
#' #' @return A list containing the result of isotemporal multilevel substitution model.
#' #' Each elements of the list corresspond to a compositional variable.
#' #'
#' #' @importFrom data.table as.data.table copy := setnames
#' #' @importFrom compositions acomp ilr clo
#' #' @importFrom extraoperators %snin% %sin%
#' #' @importFrom insight find_predictors
#' #' @importFrom emmeans ref_grid
#' #' @export
#' #' @examples
#' #' data(psub)
#' #'
#' #' library(doFuture)
#' #' registerDoFuture()
#' #' plan(multisession, workers = 5)
#' #' system.time(bsmtest2 <- bsubmargins(object = adjbrmcodatest, substitute = psub, minute = 5))
#' .bsubmargins <- function (object, substitute, minute = 60L) {
#' 
#'   if(isFALSE(missing(minute))) {
#'     if (isFALSE(is.integer(minute))) {
#'       if (isFALSE(is.numeric(minute))) {
#'         stop("'minute' must be an integer or a numeric value > 0.")
#'       }
#'     }
#'   } else {
#'     minute <- 60L
#'   }
#'   if (isFALSE(identical(ncol(substitute), length(object$CompIlr$parts)))) {
#'     stop(sprintf("The number of columns in 'substitute' (%d) must be the same
#'   as the compositional variables in 'parts' (%d).",
#'                  ncol(substitute),
#'                  length(object$CompIlr$parts)))
#'   }
#'   if (isFALSE(identical(colnames(substitute), object$CompIlr$parts))) {
#'     stop(sprintf("The names of compositional variables must be the same
#'   in 'substitute' (%s) and 'parts' (%s).",
#'                  colnames(substitute),
#'                  object$CompIlr$parts))
#'   }
#' 
#'   # Get between-person composition
#'   b <- object$CompIlr$BetweenComp
#'   b <- as.data.table(clo(b, total = 1440))
#' 
#'   psi <- object$CompIlr$psi
#'   min <- as.integer(minute)
#' 
#'   # model for no change
#'   ## ILRs
#'   bilr <- object$CompIlr$BetweenILR
#'   wilr <- matrix(0, nrow = nrow(bilr), ncol = ncol(bilr))
#'   wilr <- as.data.table(wilr)
#'   colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
#'   colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
#' 
#'   # prediction
#'   dsame <- cbind(bilr, wilr, object$CompIlr$data)
#'   ysame <- fitted(object$Model, newdata = dsame, re.form = NA, summary = FALSE)
#'   ysame <- rowMeans(ysame) # average across participants when there is no change
#' 
#'   iout <- .get.bsubmargins.(object = object, b = b,
#'                             substitute = substitute,
#'                             min = min, psi = psi,
#'                             ysame = ysame)
#' 
#' }
