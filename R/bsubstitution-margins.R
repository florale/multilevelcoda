#' #' @title Margianl Effects of Between-person Isotemporal Substitution.
#' #'
#' #' @description
#' #' This function estimates the difference in outcomes
#' #' when compositional variables are substituted for a specific period
#' #' at between-person level. The model loops through all compositional variables present
#' #' in the \code{\link{brmcoda}} object.
#' #'
#' #' @param data A fitted \code{\link{brmcoda}} model data. Required.
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
#' #' ps <- possub(data = mcompd, composition = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' #'
#' #' library(doFuture)
#' #' registerDoFuture()
#' #' plan(multisession, workers = 5)
#' #' system.time(bsmtest2 <- bsubmargins(data = adjbrmcodatest, substitute = ps, minute = 5))
#' bsubmargins <- function (data, substitute, minute = 60L) {
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
#'   if (isFALSE(identical(ncol(substitute), length(data$CompIlr$composition)))) {
#'     stop(sprintf("The number of columns in 'substitute' (%d) must be the same
#'   as the compositional variables in 'composition' (%d).",
#'                  ncol(substitute),
#'                  length(data$CompIlr$composition)))
#'   }
#'   if (isFALSE(identical(colnames(substitute), data$CompIlr$composition))) {
#'     stop(sprintf("The names of compositional variables must be the same
#'   in 'substitute' (%s) and 'composition' (%s).",
#'                  colnames(substitute),
#'                  data$CompIlr$composition))
#'   }
#' 
#'   tmp <- copy(data)
#' 
#'   # Get between-person composition
#'   b <- tmp$CompIlr$BetweenComp
#'   b <- as.data.table(clo(b, total = 1440))
#' 
#'   psi <- tmp$CompIlr$psi
#'   min <- as.integer(minute)
#' 
#'   # model for no change
#'   ## ILRs
#'   bilr <- tmp$CompIlr$BetweenILR
#'   wilr <- matrix(0, nrow = nrow(bilr), ncol = ncol(bilr))
#'   wilr <- as.data.table(wilr)
#'   colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
#'   colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
#' 
#'   # prediction
#'   dsame <- cbind(bilr, wilr, tmp$CompIlr$data)
#'   ysame <- fitted(tmp$BrmModel, newdata = dsame, re.form = NA, summary = FALSE)
#'   ysame <- rowMeans(ysame) # emmeans across participants when there is no change
#' 
#'   iout <- .get.bsubm2(substitute = substitute,
#'                       b = b, tmp = tmp,
#'                       min = min, psi = psi,
#'                       ysame = ysame)
#' 
#' }
