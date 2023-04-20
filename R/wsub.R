#' @title Within-person Basic Substitution Model
#' 
#' @description
#' Using a fitted model object, estimate the difference in outcomes
#' when compositional parts are substituted for specific unit(s) at `within-person` level. 
#' The \code{wsub} output encapsulates 
#' the substitution results for all compositional parts
#' present in the \code{\link{brmcoda}} object.
#' 
#' Notes: The reference composition for substitution model 
#' is the compositional mean of the data set provided.
#' For average marginal effect, use \code{\link{wsubmargins}}.
#'
#' @param object A fitted \code{\link{brmcoda}} object. Required.
#' @param delta A integer, numeric value or vector indicating the amount of substituted change between compositional parts.
#' @param basesub A \code{data.frame} or \code{data.table} of the base possible substitution of compositional parts.
#' This data set can be computed using function \code{\link{basesub}}. 
#' If \code{NULL}, all possible pairwise substitution of compositional parts are used.
#' @param regrid If non-\code{NULL}, a \code{data.table} of reference grid consisting 
#' of combinations of covariates over which predictions are made.
#' If \code{NULL}, the reference grid is constructed via \code{\link{ref_grid}}.
#' @param summary A logical value. 
#' Should the estimate at each level of the reference grid (\code{FALSE}) 
#' or their average (\code{TRUE}) be returned? Default to \code{TRUE}.
#' @param  recomp A numeric or integer vector used as reference composition. If \code{NULL},
#' compositional mean is used.
#' @param level A character string or vector. 
#' Should the estimate be at the \code{between}-person and/or \code{within}-person level? Required.
#' @param type A character string or vector. 
#' Should the estimate be \code{conditional} mean or average \code{marginal} mean? Required.
#' @param ... Additional arguments to be passed to \code{\link{describe_posterior}}.
#' 
#' @return A list containing the result of multilevel compositional substitution model.
#' Each element of the list is the estimation for a compositional part 
#' and include at least six elements.
#' \itemize{
#'   \item{\code{Mean}}{ Posterior means.}
#'   \item{\code{CI_low}} and \item{\code{CI_high}}{ 95% credible intervals.}
#'   \item{\code{Delta}}{ Amount substituted across compositional parts.}
#'   \item{\code{From}}{ Compositional part that is substituted from.}
#'   \item{\code{To}}{ Compositional parts that is substituted to.}
#'   \item{\code{Level}}{ Level where changes in composition takes place.}
#'   \item{\code{EffectType}}{ Either estimated `conditional` or average `marginal` changes.}
#' }
#' 
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
#'              chain = 1, iter = 500,
#'              backend = "cmdstanr")
#'              
#' subm <- wsub(object = m, basesub = psub, delta = 5)
#' }
wsub <- function(object,
                 delta,
                 basesub,
                 regrid = NULL,
                 summary = TRUE,
                 recomp = NULL,
                 level = "within",
                 type = "conditional",
                 ...) {
  
  
  # compositional mean
  b <- object$CompILR$BetweenComp

  if (is.null(recomp)) {
    mcomp <- mean(b, robust = TRUE)
    mcomp <- clo(mcomp, total = object$CompILR$total)
    mcomp <- as.data.table(t(mcomp))
    
  } else {
    if (isFALSE(identical(length(object$CompILR$parts), length(recomp)))) {
      stop(
        sprintf(
          "The number of values in parts (%d)
  must be the same as in recomp (%d).",
  length(object$CompILR$parts),
  length(recomp)
        ))
    }
    
    if (isFALSE(class(recomp) %in% c("numeric", "interger"))) {
      stop(
        sprintf(
          "recomp (%s) should be a vector of numeric or interger values.",
          class(recomp)
        ))
    }
    
    if(isFALSE(sum(recomp) == object$CompILR$total)) {
      stop(sprintf(
        "The total amount of recomp (%s) should be the same as the composition (%s).",
        sum(recomp),
        object$CompILR$total
      ))
    }
    
    if (isTRUE((any(recomp > apply(object$CompILR$data[, object$CompILR$parts, with = FALSE], 2, max)) |
                any(recomp < apply(object$CompILR$data[, object$CompILR$parts, with = FALSE], 2, min))))) {
      stop(paste(
        sprintf(
          "recomp should be numeric or interger values that are between (%s) and (%s)",
          paste0(round(apply(object$CompILR$data[, object$CompILR$parts, with = FALSE], 2, min)), collapse = ", "),
          paste0(round(apply(object$CompILR$data[, object$CompILR$parts, with = FALSE], 2, max)), collapse = ", ")),
        "\n", 
        " for",
        paste0(object$CompILR$parts, collapse = ", "),
        "respectively"
      ))
    }
    recomp <- as.integer(recomp)
    mcomp  <- clo(recomp, total = object$CompILR$total)
    mcomp  <- as.data.table(t(mcomp))
    colnames(mcomp) <- colnames(object$CompILR$BetweenComp)
  }
  
  # error if delta out of range
  if(isTRUE(any(delta > min(mcomp)))) {
    stop(sprintf(
      "delta value should be less than or equal to %s, which is
  the amount of composition part available for pairwise substitution.",
  round(min(mcomp), 2)
    ))
  }
  
  # input for substitution model
  ID <- 1 # to make fitted() happy
  delta <- as.integer(delta)
  recomp <- is.integer(recomp)
  
  # model for no change
  bilr <- ilr(mcomp, V = object$CompILR$psi)
  bilr <- as.data.table(t(bilr))
  wilr <- as.data.table(matrix(0, nrow = nrow(bilr), ncol = ncol(bilr)))
  colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
  colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
  
  # check covariates
  ilrn <- c(names(object$CompILR$BetweenILR), names(object$CompILR$WithinILR)) # get ilr names in brm model
  vn <- do.call(rbind, find_predictors(object$Model)) # get all varnames in brm model
  
  # if there is no covariates
  # number of variables in the brm model = number of ilr coordinates
  if (isTRUE(identical(length(vn), length(ilrn)))) {
    # unadj subsitution model
    if (isFALSE(is.null(regrid))) {
      warning(
        paste(
          "This is an unadjusted model, but a reference grid was provided.",
          "  Please note that the covariates provided in the reference grid",
          "  need to be present in 'brmcoda' model object.",
          "  Unadjusted substitution model was estimated.",
          sep = "\n"
        )
      )
    }
    
    dsame <- cbind(bilr, wilr, ID)
    ysame <- fitted(
      object$Model,
      newdata = dsame,
      re_formula = NA,
      summary = FALSE)
    
    # substitution model
    out <- get.wsub(
      object = object,
      basesub = basesub,
      mcomp = mcomp,
      delta = delta,
      ysame = ysame,
      summary = summary,
      level = level,
      type = type)
    
  } else {
    # adj subsitution model
    # reference grid containing covariates
    rg <- as.data.table(ref_grid(object$Model)@grid)
    cv <- colnames(rg) %snin% c(ilrn, ".wgt.")
    
    if (isFALSE(is.null(regrid))) {
      # check user's specified reference grid
      if (isFALSE(identical(colnames(regrid), cv))) {
        # ensure all covs are provided
        stop(
          paste(
            "'regrid' should contains information about",
            "  the covariates in 'brmcoda' model to estimate the substitution model.",
            "  It should not include ILR variables nor any column names starting with 'bilr', 'wilr', or 'ilr',",
            "  as these variables will be calculated by substitution model.",
            "  Please provide a different reference grid.",
            sep = "\n"))
      } else {
        refg <- regrid
      }
    } else {
      # use default rg
      refg <- rg[, cv, with = FALSE]
    }
    
    dsame <- cbind(bilr, wilr, ID, refg)
    ysame <- fitted(
      object$Model,
      newdata = dsame,
      re_formula = NA,
      summary = FALSE)
    
    # substitution model
    out <- get.wsub(
      object = object,
      basesub = basesub,
      mcomp = mcomp,
      delta = delta,
      ysame = ysame,
      summary = summary,
      cv = cv,
      refg = refg,
      level = level,
      type = type)
  }
}