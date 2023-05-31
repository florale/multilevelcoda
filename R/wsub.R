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
#' @param summary A logical value. 
#' Should the estimate at each level of the reference grid (\code{FALSE}) 
#' or their average (\code{TRUE}) be returned? Default to \code{TRUE}.
#' @param level A character string or vector. 
#' Should the estimate be at the \code{between}-person and/or \code{within}-person level? Required.
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
                 basesub,
                 delta,
                 summary = TRUE,
                 ref = "unitmean",
                 refdata,
                 level = "within",
                 weight = c("equal", "proportional"),
                 ...) {
  
  # compositional mean
  b <- object$CompILR$BetweenComp
  mcomp <- mean(b, robust = TRUE)
  mcomp <- acomp(mcomp, total = object$CompILR$total)
  mcomp <- as.data.table(t(mcomp))
  
  # refcomp --------------------
  if (is.null(refcomp)) {
    refcomp <- mcomp
    
  } else {
    if (isFALSE(identical(length(object$CompILR$parts), length(refcomp)))) {
      stop(
        sprintf(
          "The number of values in parts (%d)
  must be the same as in refcomp (%d).",
  length(object$CompILR$parts),
  length(refcomp)
        ))
    }
    
    if (isFALSE(class(refcomp) %in% c("numeric", "interger"))) {
      stop(
        sprintf(
          "refcomp (%s) should be a vector of numeric or interger values.",
          class(refcomp)
        ))
    }
    
    if(isFALSE(sum(refcomp) == object$CompILR$total)) {
      stop(sprintf(
        "The total amount of refcomp (%s) should be the same as the composition (%s).",
        sum(refcomp),
        object$CompILR$total
      ))
    }
    
    if (isTRUE((any(refcomp > apply(object$CompILR$data[, object$CompILR$parts, with = FALSE], 2, max)) |
                any(refcomp < apply(object$CompILR$data[, object$CompILR$parts, with = FALSE], 2, min))))) {
      stop(paste(
        sprintf(
          "refcomp should be numeric or interger values that are between (%s) and (%s)",
          paste0(round(apply(object$CompILR$data[, object$CompILR$parts, with = FALSE], 2, min)), collapse = ", "),
          paste0(round(apply(object$CompILR$data[, object$CompILR$parts, with = FALSE], 2, max)), collapse = ", ")),
        "\n", 
        " for",
        paste0(object$CompILR$parts, collapse = ", "),
        "respectively"
      ))
    }
    refcomp <- as.integer(refcomp)
    refcomp  <- acomp(refcomp, total = object$CompILR$total)
    refcomp  <- as.data.table(t(refcomp))
    colnames(refcomp) <- colnames(mcomp)
  }
  
  # error if delta out of range
  if(isTRUE(any(delta > min(refcomp)))) {
    stop(sprintf(
      "delta value should be less than or equal to %s, which is
  the amount of composition part available for pairwise substitution.",
  round(min(refcomp), 2)
    ))
  }
  
  # d0 ---------------------------
  # input for substitution model
  ID <- 1 # to make fitted() happy
  delta <- as.integer(delta)
  
  # bilr is between-person ilr of the ref comp (doesn't have to be compositional mean)
  bilr0 <- ilr(refcomp, V = object$CompILR$psi)
  bilr0 <- as.data.table(t(bilr0))
  
  # wcomp and wilr are the difference between the actual compositional mean of the dataset and bilr
  # is 0 if ref comp is compositional mean
  # but is different if not
  wcomp <- refcomp - mcomp
  wilr0 <- as.data.table(t(ilr(wcomp, V = object$CompILR$psi)))
  
  colnames(bilr0) <- colnames(object$CompILR$BetweenILR)
  colnames(wilr0) <- colnames(object$CompILR$WithinILR)
  
  # check covariates
  ilrnames <- c(names(object$CompILR$BetweenILR), names(object$CompILR$WithinILR)) # get ilr names in model
  varnames <- do.call(rbind, find_predictors(object$Model)) # get all varnames in model
  
  # if there is no covariates
  # number of variables in the brm model = number of ilr coordinates
  if (isTRUE(identical(length(varnames), length(ilrnames)))) { # unadj subsitution model
    if (isFALSE(is.null(refgrid))) { 
      warning(paste(
        "This is an unadjusted model, but a reference grid was provided.",
        "  Please note that the covariates provided in the reference grid",
        "  need to be present in 'brmcoda' model object.",
        "  Unadjusted substitution model was estimated.",
        sep = "\n"))
    }
    
    d0 <- cbind(bilr0, wilr0, ID)
    
  } else { # adj subsitution model
    # reference grid containing covariates
    rg <- as.data.table(ref_grid(object$Model)@grid)
    covnames <- colnames(rg) %snin% c(ilrnames, ".wgt.")
    
    if (isFALSE(is.null(refgrid))) { # check user's specified reference grid
      if(isFALSE(identical(colnames(refgrid), covnames))) { # ensure all covs are provided
        stop(paste(
          "'refgrid' should contains information about",
          "  the covariates in 'brmcoda' model to estimate the substitution model.",
          "  It should not include ILR variables nor any column names starting with 'bilr', 'wilr', or 'ilr',",
          "  as these variables will be calculated by substitution model.",
          "  Please provide a different reference grid.",
          sep = "\n"))
      } else {
        refgrid <- as.data.table(refgrid)
      }
    } else { # use default rg
      refgrid <- rg[, covnames, with = FALSE]
    }
    
    d0 <- cbind(bilr0, wilr0, ID, refgrid)
  }
  
  # y0 --------------------------------
  y0 <- fitted(
    object$Model,
    newdata = d0,
    re_formula = NA,
    summary = FALSE, 
    ...)
  
  # yw ---------------------------------
    # substitution model
    out <- get.wsub(
      object = object,
      basesub = basesub,
      recomp = recomp,
      delta = delta,
      y0 = y0,
      d0 = d0,
      summary = summary,
      covnames = covnames,
      refgrid = refgrid,
      level = level,
      type = type)

}