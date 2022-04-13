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
#' #' bsmtest2 <- bsubmargins(data = adjbrmcodatest, substitute = ps, minute = 10)
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
#'   # Compute between-person composition
#'   b <- tmp$CompIlr$BetweenComp[tmp$CompIlr$data[, which(!duplicated(get(tmp$CompIlr$idvar)))], ]
#'   b <- as.data.table(clo(b, total = 1440)) # fix this to allow users to change total value
#' 
#'   ID <- 1 # to make fitted() happy
#'   psi <- tmp$CompIlr$psi
#'   min <- as.integer(minute)
#' 
#'   # set up model for no change
#'   bilr <- tmp$CompIlr$BetweenILR[tmp$CompIlr$data[, which(!duplicated(get(tmp$CompIlr$idvar)))], ]
#'   wilr <- matrix(0, nrow = nrow(bilr), ncol = ncol(bilr))
#'   wilr <- as.data.table(wilr)
#'   colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
#'   colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
#' 
#'   # check covariates
#'   ilrn <- c(names(tmp$CompIlr$BetweenILR), names(tmp$CompIlr$WithinILR)) # get ilr names from model
#'   vn <- do.call(rbind, find_predictors(tmp$BrmModel)) # get all varnames from model
#' 
#'   # if there is no covariates
#'   # number of variables in the brm model = number of ilr coordinates
#'   # run unadj subsitution model
#'   if (isTRUE(identical(length(vn), length(ilrn)))) {
#'     dsame <- cbind(bilr, wilr, ID)
#'     ysame <- fitted(tmp$BrmModel, newdata = dsame, re.form = NA, summary = FALSE)
#'     ysame <- rowMeans(ysame) # emmeans across participants when there is no change
#'     
#'     iout <- .get.bsubm2(substitute = substitute,
#'                         b = b, tmp = tmp,
#'                         min = min, psi = psi,
#'                         ysame = ysame)
#'     
#'     } else { # run adjusted model
#'       # Get reference grid containing covariates
#'       refg <- as.data.table(ref_grid(tmp$BrmModel) @grid)
#'       cv <- colnames(refg) %snin% c(ilrn, ".wgt.")
#'       refg <- refg[, cv, with = FALSE] # reference grid
#'       
#'       hout <- vector("list", length = nrow(refg))
#'       if (isFALSE(nrow(refg) == 1)) {
#'         for (h in seq_len(nrow(refg))) {
#'           refgbase <- refg[h, ]
#'           refgbase <- refgbase[rep(seq_len(nrow(refgbase)), nrow(b)), ]
#'           hout[[h]] <- refgbase
#'           }
#'         refg <- do.call(rbind, hout)
#'         h <- h
#'         } else {
#'           refg <- refg
#'           h <- 1
#'           }
#'       bilr <- bilr[rep(seq_len(nrow(bilr)), h), ]
#'       wilr <- wilr[rep(seq_len(nrow(wilr)), h), ]
#'       dsame <- cbind(bilr, wilr, ID, refg)
#'       ysame <- fitted(tmp$BrmModel, newdata = dsame, re.form = NA, summary = FALSE)
#'       ysame <- rowMeans(ysame) # emmeans across participants when there is no change
#'       # look into weighted mean
#'       iout <- .get.bsubm2(substitute = substitute,
#'                           b = b, tmp = tmp,
#'                           min = min, psi = psi,
#'                           ysame = ysame, refg = refg, h = h)
#'       }
#' }
