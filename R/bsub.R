## make Rcmd check happy
utils::globalVariables(c("Mean",  "CI_low", "CI_high", "Substitute", "MinSubstituted"))

#' Between-person Substitution Model.
#'
#' Estimates the difference in outcomes
#' when compositional variables are substituted for a specific time period.
#' at between-person level.
#'
#' @param object A \code{\link{brmcoda}} object.
#' @param substitute A \code{data.frame} or \code{data.table} of the possible substitution of variables.
#' This dataset can be computed using function \code{possub}. Required.
#' @param minute A integer or numeric value indicating the maximum minute for which substitution model is desired.
#' Default to \code{60L} (i.e., the model loops through 1:60L minutes).
#' @param regrid If non-\code{NULL}, a \code{data.table} of reference grid consists 
#' of combinations of covariates over which predictions are made.
#' If \code{NULL}, the reference grid is constructed via \code{\link{ref_grid}}.
#' @param summary A logical value. 
#' Should estimated marginal means at each level of the reference grid (\code{FALSE}) 
#' or the marginal averages thereof (\code{TRUE}) be returned? Default is FALSE
#' @param ... Additional arguments to be passed to \code{\link{describe_posterior}}.
#' 
#' @return A list of results from substitution models for all compositional variables.
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo mean.acomp
#' @importFrom extraoperators %snin% %sin%
#' @importFrom insight find_predictors
#' @importFrom emmeans ref_grid
#' @importFrom foreach %dopar%
#' @export
#' @examples
#'
#' data(mcompd)
#' data(sbp)
#' ps <- possub(parts = c("TST", "WAKE", "MVPA", "LPA", "SB")
#'
#' bsubtest <- bsub(object = adjbrmcodatest, substitute = ps, minute = 10)
#'
#' ## cleanup
#' rm(bsubtest, mcompd)
bsub <- function(object, substitute, minute = 60L, 
                 regrid = NULL, summary = FALSE, 
                 ...) {

  if(isFALSE(missing(minute))) {
    if (isFALSE(is.integer(minute))) {
      if (isFALSE(is.numeric(minute))) {
        stop("'minute' must be an integer or a numeric value > 0.")
      }
    }
  } else {
    minute <- 60L
  }
  if (isFALSE(identical(ncol(substitute), length(object$CompIlr$parts)))) {
    stop(sprintf("The number of columns in 'substitute' (%d) must be the same
  as the compositional variables in 'parts' (%d).",
                 ncol(substitute),
                 length(object$CompIlr$parts)))
  }
  if (isFALSE(identical(colnames(substitute), object$CompIlr$parts))) {
    stop(sprintf("The names of compositional variables must be the same
  in 'substitute' (%s) and 'parts' (%s).",
                 colnames(substitute),
                 object$CompIlr$parts))
  }

  # Compute compositional mean
  b <- object$CompIlr$BetweenComp

  mcomp <- mean(b, robust = TRUE)
  mcomp <- clo(mcomp, total = 1440)
  mcomp <- as.data.table(t(mcomp))
  names(mcomp) <- paste0("B", names(mcomp))

  # input for substitution model
  ID <- 1 # to make fitted() happy
  psi <- object$CompIlr$psi
  min <- as.integer(minute)

  # model for no change
  bilr <- ilr(mcomp, V = psi)
  bilr <- as.data.table(t(bilr))
  wilr <- as.data.table(matrix(0, nrow = nrow(bilr), ncol = ncol(bilr)))
  colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
  colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))

  # check covariates
  ilrn <- c(names(object$CompIlr$BetweenILR), names(object$CompIlr$WithinILR)) # get ilr names in brm model
  vn <- do.call(rbind, find_predictors(object$BrmModel)) # get all varnames in brm model

  # if there is no covariates
  # number of variables in the brm model = number of ilr coordinates
  if (isTRUE(identical(length(vn), length(ilrn)))) { # unadj subsitution model
    dsame <- cbind(bilr, wilr, ID)
    ysame <- fitted(object$BrmModel, newdata = dsame, re.form = NA, summary = FALSE)
    
    # substitution model
    iout <- get.bsub(object = object,
                     substitute = substitute,
                     mcomp = mcomp,
                     min = min, psi = psi, 
                     ysame = ysame)
    } else { # adj subsitution model
      
    # reference grid containing covariates
      if (isFALSE(missing(regrid))) {
        refg <- regrid
        } else {
          refg <- as.data.table(ref_grid(object$BrmModel) @grid)
          cv <- colnames(refg) %snin% c(ilrn, ".wgt.")
          refg <- refg[, cv, with = FALSE]
        }
      
      hout <- vector("list", length = nrow(refg))
      if (isFALSE(nrow(refg) == 1)) {
        for (h in seq_len(nrow(refg))) {
          refgbase <- refg[h, ]
          dsame <- cbind(bilr, wilr, ID, refgbase)
          hout[[h]] <- dsame
          }
        dsame <- do.call(rbind, hout)
        } else {
          refg <- refg
          dsame <- cbind(bilr, wilr, ID, refg)
        }
      
      ysame <- fitted(object$BrmModel, newdata = dsame, re.form = NA, summary = FALSE)
      
      # substitution model
      iout <- get.bsub(object = object,
                       substitute = substitute,
                       mcomp = mcomp, min = min, 
                       ysame = ysame, 
                       summary = summary,
                       cv = cv, psi = psi, refg = refg)
  }
}