#' @title Between-person Substitution Model.
#' 
#' @description 
#' Estimate the difference in outcomes
#' when compositional variables are substituted for a specific time period
#' at `within-person` level. The resulting \code{bsub} encapsulates
#' substitution estimation across all compositional variables
#' present in the \code{\link{brmcoda}} object.
#' 
#' Notes: The reference composition for substitution model 
#' is the compositional mean of the dataset provided.
#' For average marginal effect, consider using \code{\link{bsubmargins}}.
#'
#' @param object A \code{\link{brmcoda}} object.
#' @param substitute A \code{data.frame} or \code{data.table} of the possible substitution of variables.
#' This dataset can be computed using function \code{possub}. Required.
#' @param minute A integer or numeric value indicating the maximum minute 
#' for which substitution model is desired.
#' Default to \code{60L} (i.e., the model loops through 1:60L minutes).
#' @param regrid If non-\code{NULL}, a \code{data.table} of reference grid consisting 
#' of combinations of covariates over which predictions are made.
#' Otherwise, the reference grid is constructed via \code{\link{ref_grid}}.
#' @param summary A logical value. 
#' Should estimated marginal means at each level of the reference grid (\code{FALSE}) 
#' or the marginal averages thereof (\code{TRUE}) be returned? Default to \code{FALSE}.
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
#' data(mcompd)
#' data(sbp)
#' data(psub)
#' \donttest{
#' cilr <- compilr(data = mcompd, sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#' 
#' m <- brmcoda(compilr = cilr, 
#'              formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + wilr1 + 
#'              wilr2 + wilr3 + wilr4 + Female + (1 | ID))
#'                
#' testbs <- bsub(object = m, substitute = psub, minute = 5)
#' }
#' ## cleanup
#' rm(bsubtest, mcompd)
bsub <- function(object, substitute, minute = 60L, 
                 regrid = NULL, summary = FALSE, 
                 ...) {

  if (isTRUE(missing(object))) {
    stop(paste(
      "'object' is a required argument and cannot be missing;",
      " it should be an object of class brmcoda.", 
      " See ?multilevelcoda::brmcoda for details.",
      sep = "\n"))
  }
  
  if (isFALSE(inherits(object, "brmcoda"))) {
    stop(paste(
      "'object' should be a fitted brmcoda object",
      " See ?multilevelcoda::brmcoda for details.",
      sep = "\n"))
  }
  
  if (isTRUE(missing(substitute))) {
    stop(paste(
      "'substitute' is a required argument and cannot be missing;",
      " it should be a dataset of possible substitution", 
      " and can be computed using multilevelcoda::possub.", 
      " See ?multilevelcoda::possub for details.",
      sep = "\n"))
  }
  
  if(isFALSE(missing(minute))) {
    if (isFALSE(is.integer(minute))) {
      if (isFALSE(minute > 0)) {
        stop("'minute' must be an positive integer value.")
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
    stop(sprintf("The names of compositional variables must be the 
    same in 'substitute' (%s) and 'parts' (%s).",
                 colnames(substitute),
                 object$CompIlr$parts))
  }
  
  if (isFALSE(is.null(regrid))) {
    if(any(c(colnames(object$CompIlr$BetweenILR), colnames(object$CompIlr$WithinILR))
           %in% c(colnames(regrid)))) {
      stop(paste(
        "'regrid' should not have any column names starting with 'bilr', 'wilr', or 'ilr'.",
        "   It should contain information about the covariates used in 'brmcoda'.",
        "   These variables are used for subsequent substitution model.",
        "   Please provide a different reference grid.",
        sep = "\n"))
    }
  }
  
  # compositional mean
  b <- object$CompIlr$BetweenComp

  mcomp <- mean(b, robust = TRUE)
  mcomp <- clo(mcomp, total = object$CompIlr$total)
  mcomp <- as.data.table(t(mcomp))

  # input for substitution model
  ID <- 1 # to make fitted() happy
  min <- as.integer(minute)

  # model for no change
  bilr <- ilr(mcomp, V = object$CompIlr$psi)
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
    if (isFALSE(is.null(regrid))) { 
      warning(paste(
        "This is an unadjusted model, but a reference grid was provided.",
        "  Please note that the covariates provided in the reference grid",
        "  need to be present in 'brmcoda' model object.",
        "  Unadjusted substitution model was estimated.",
        sep = "\n"))
    }
    
    dsame <- cbind(bilr, wilr, ID)
    ysame <- fitted(object$BrmModel, newdata = dsame, re.form = NA, summary = FALSE)
    
    # substitution model
    out <- get.bsub(object = object, b = b,
                    substitute = substitute,
                    mcomp = mcomp,
                    min = min, 
                    ysame = ysame)
    } else { # adj subsitution model
      # reference grid containing covariates
      rg <- as.data.table(ref_grid(object$BrmModel) @grid)
      cv <- colnames(rg) %snin% c(ilrn, ".wgt.")
      
      if (isFALSE(is.null(regrid))) { # check reference grid
        if(isFALSE(identical(colnames(regrid), cv))) {
          stop(paste(
            "Please provide a reference grid that contains information about",
            "  all covariates in 'brmcoda' model to estimate the substitution model.",
            sep = "\n"))
          } else {
            refg <- regrid
            }
        } else {
          refg <- rg[, cv, with = FALSE]
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
      out <- get.bsub(object = object,
                       substitute = substitute,
                       mcomp = mcomp, 
                       min = min, 
                       ysame = ysame, 
                       summary = summary,
                       cv = cv, refg = refg)
  }
}