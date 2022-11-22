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
#' @param object A fitted \code{\link{brmcoda}} object.
#' @param substitute A \code{data.frame} or \code{data.table} of possible substitution of variables.
#' This dataset can be computed using function \code{possub}. Required.
#' @param minute A integer or numeric value indicating the maximum minute 
#' for which substitution model is desired.
#' Default to \code{60L} (i.e., the model loops through 1:60L minutes).
#' @param regrid If non-\code{NULL}, a \code{data.table} of reference grid consisting 
#' of combinations of covariates over which predictions are made.
#' Otherwise, the reference grid is constructed via \code{\link{ref_grid}}.
#' @param summary A logical value. 
#' Should the estimate at each level of the reference grid (\code{FALSE}) 
#' or their average (\code{TRUE}) be returned? Default to \code{TRUE}.
#' @param ... Additional arguments to be passed to \code{\link{describe_posterior}}.
#' 
#' @return A list containing the result of isotemporal multilevel substitution model.
#' Each elements of the list is the substitution estimation for a compositional variables, 
#' which include at least six elements.
#' \itemize{
#'   \item{\code{Mean}}{ Posterior means.}
#'   \item{\code{CI_low}} and \item{\code{CI_high}}{ 95% credible intervals.}
#'   \item{\code{MinSubstituted}}{ Minute substituted within the composition.}
#'   \item{\code{Substitute}}{Compositional variables to be substituted from/to.}
#'   \item{\code{Predictor}}{Central compositional variable to be substituted from/to.}
#' }
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo mean.acomp
#' @importFrom extraoperators %snin% %sin%
#' @importFrom insight find_predictors
#' @importFrom emmeans ref_grid
#' @importFrom stats fitted
#' @export
#' @examples
#' \dontrun{
#' data(mcompd)
#' data(sbp)
#' data(psub)
#' cilr <- compilr(data = mcompd, sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#' 
#' # model with compositional predictor at between and within-person levels
#' m <- brmcoda(compilr = cilr, 
#'              formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + 
#'                                 wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID), 
#'              chain = 1, iter = 500)
#'              
#' subm <- bsub(object = m, substitute = psub, minute = 5)
#' }
bsub <- function(object, substitute, minute = 60L, 
                 regrid = NULL, summary = TRUE, 
                 ...) {

  if (isTRUE(missing(object))) {
    stop(paste(
      "'object' is a required argument and cannot be missing;",
      "  it should be an object of class 'brmcoda'.", 
      "  See ?bsub for details.",
      sep = "\n"))
  }
  
  if (isFALSE(inherits(object, "brmcoda"))) {
    stop(sprintf(
    "Can't handle an object of class (%s) 
  It should be a fitted 'brmcoda' object
  See ?bsub for details.",
                 class(object)))
  }
  
  if (isTRUE(missing(substitute))) {
    stop(paste(
      "'substitute' is a required argument and cannot be missing;",
      "  it should be a dataset of possible substitution", 
      "  and can be computed using multilevelcoda::possub.", 
      "  See ?bsub for details.",
      sep = "\n"))
  }
  
  if(isFALSE(missing(minute))) {
    if (isFALSE(is.integer(minute))) {
      if (isFALSE(minute > 0)) {
        stop(" 'minute' must be an positive integer value.")
        }
      }
    } else {
      minute <- 60L
    }
  
  if (isFALSE(identical(ncol(substitute), length(object$CompIlr$parts)))) {
    stop(sprintf(
    "The number of columns in 'substitute' (%d) must be the same
  as the compositional variables in 'parts' (%d).",
                 ncol(substitute),
                 length(object$CompIlr$parts)))
  }
  
  if (isFALSE(identical(colnames(substitute), object$CompIlr$parts))) {
    stop(sprintf(
    "The names of compositional variables must be the 
  same in 'substitute' (%s) and 'parts' (%s).",
                 colnames(substitute),
                 object$CompIlr$parts))
  }
  
  if (isFALSE(is.null(regrid))) {
    if(any(c(colnames(object$CompIlr$BetweenILR), colnames(object$CompIlr$WithinILR))
           %in% c(colnames(regrid)))) {
      stop(paste(
        "'regrid' should not have any column names starting with 'bilr', 'wilr', or 'ilr'.",
        "  These variables will be calculated by substitution model.",
        "  Reference grid should contain information about the covariates used in 'brmcoda'.",
        "  Please provide a different reference grid.",
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
  ilrn <- c(names(object$CompIlr$BetweenILR), names(object$CompIlr$WithinILR)) # get ilr names in model
  vn <- do.call(rbind, find_predictors(object$Model)) # get all varnames in model

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
    ysame <- fitted(object$Model, newdata = dsame, re_formula = NA, summary = FALSE)
    
    # substitution model
    out <- get.bsub(object = object, substitute = substitute,
                    mcomp = mcomp, min = min, ysame = ysame, summary = summary)
    
    } else { # adj subsitution model
      # reference grid containing covariates
      rg <- as.data.table(ref_grid(object$Model) @grid)
      cv <- colnames(rg) %snin% c(ilrn, ".wgt.")
      
      if (isFALSE(is.null(regrid))) { # check user's specified reference grid
        if(isFALSE(identical(colnames(regrid), cv))) { # ensure all covs are provided
          stop(paste(
            "'regrid' should contains information about",
            "  the covariates in 'brmcoda' model to estimate the substitution model.",
            "  It should not include ILR variables nor any column names starting with 'bilr', 'wilr', or 'ilr',",
            "  as these variables will be calculated by substitution model.",
            "  Please provide a different reference grid.",
            sep = "\n"))
          } else {
            refg <- regrid
          }
        } else { # use default rg
          refg <- rg[, cv, with = FALSE]
          }
      
      dsame <- cbind(bilr, wilr, ID, refg)
      ysame <- fitted(object$Model, newdata = dsame, re_formula = NA, summary = FALSE)
      
      # substitution model
      out <- get.bsub(object = object, substitute = substitute,
                      mcomp = mcomp, min = min, ysame = ysame, 
                      summary = summary, cv = cv, refg = refg)
  }
}