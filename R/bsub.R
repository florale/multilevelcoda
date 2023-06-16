#' @title Between-person Simple Substitution.
#' 
#' @description 
#' Estimate the difference in outcomes
#' when compositional parts are substituted for specific unit(s) at `between` level. 
#' The \code{bsub} output encapsulates
#' the substitution results for all compositional parts
#' present in the \code{\link{brmcoda}} object.
#' 
#' @param object A fitted \code{\link{brmcoda}} object. Required.
#' @param delta A integer, numeric value or vector indicating the amount of substituted change between compositional parts.
#' @param basesub A \code{data.frame} or \code{data.table} of the base possible substitution of compositional parts.
#' This data set can be computed using function \code{\link{basesub}}. 
#' If \code{NULL}, all possible pairwise substitution of compositional parts are used.
#' @param ref Either a character value or vector or a dataset.
#' \code{ref} can be \code{grandmean} or
#' a \code{data.frame} or \code{data.table} of user's specified reference grid consisting
#' of combinations of covariates over which predictions are made.
#' User's specified reference grid only applicable to substitution model
#' using a single reference composition value
#' (e.g., \code{clustermean} or user's specified). Default to \code{grandmean}.
#' @param summary A logical value. 
#' Should the estimate at each level of the reference grid (\code{FALSE}) 
#' or their average (\code{TRUE}) be returned? 
#' Default to \code{TRUE}.
#' Only applicable for model with covariates in addition to
#' the isometric log-ratio coordinates (i.e., adjusted model).
#' @param level A character string or vector. 
#' Should the estimate be at the \code{between} and/or \code{within} level?
#' Default to \code{between}.
#' @param weight A character value specifying the weight to use in calculation of the reference composition.
#' \code{weight} can be \code{equal} which gives equal weight across units (e.g., individuals) or
#' \code{proportional} which weights in proportion to the frequencies of units being averaged 
#' (e.g., observations across individuals)
#' Default to \code{equal}.
#' @param ... Additional arguments to be passed to \code{\link{describe_posterior}}.
#' 
#' @return A list containing the result of multilevel compositional substitution model.
#' Each element of the list is the estimation for a compositional part 
#' and include at least eight elements.
#' \itemize{
#'   \item{\code{Mean}}{ Posterior means.}
#'   \item{\code{CI_low}} and \item{\code{CI_high}}{ 95% credible intervals.}
#'   \item{\code{Delta}}{ Amount substituted across compositional parts.}
#'   \item{\code{From}}{ Compositional part that is substituted from.}
#'   \item{\code{To}}{ Compositional parts that is substituted to.}
#'   \item{\code{Level}}{ Level where changes in composition takes place. Either }
#'   \item{\code{Reference}}{ Either \code{grandmean}, \code{clustermean}, or \code{users}}
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
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#' data(mcompd)
#' data(sbp)
#' data(psub)
#' 
#' cilr <- compilr(data = mcompd, sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#' 
#' # model with compositional predictor at between and between-person levels
#' m <- brmcoda(compilr = cilr, 
#'              formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + 
#'                                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID), 
#'              chain = 1, iter = 500,
#'              backend = "cmdstanr")
#' subm <- bsub(object = m, basesub = psub, delta = 5)
#' }}
bsub <- function(object,
                 delta,
                 basesub,
                 summary = TRUE,
                 ref = "grandmean",
                 level = "between",
                 weight = NULL,
                 ...) {
  
  # d0 -------------------------------
  if (isTRUE(ref == "grandmean")) {
    d0 <- build.rg(object = object,
                   ref = ref,
                   weight = weight,
                   fill = FALSE)
  } else {
    if (isFALSE(inherits(ref, c("data.table", "data.frame", "matrix")))) {
      stop("ref must be 'grandmean' or a data table, data frame or matrix.")
    }
    if(isFALSE(
      identical(colnames(ref),
                colnames(as.data.table(ref_grid(object$Model)@grid))))) { # ensure all covs are provided
      stop(paste(
        "'ref' should contains information about",
        "  the covariates in 'brmcoda' model to estimate the substitution model.",
        "  Please provide a different reference grid or build one using `build.rg()`.",
        sep = "\n"))
    }
    d0 <- ref
    ref <- "users"
  }
  d0 <- as.data.table(d0)
  
  # error if delta out of range
  comp0 <- d0[1, colnames(object$CompILR$BetweenComp), with = FALSE]
  
  delta <- as.integer(delta)
  if(isTRUE(any(delta > min(comp0)))) {
    stop(sprintf(
      "delta value should be less than or equal to %s, which is
  the amount of composition part available for pairwise substitution.",
  round(min(comp0), 2)
    ))
  }

  # y0 --------------------------------
  y0 <- fitted(
    object$Model,
    newdata = d0,
    re_formula = NA,
    summary = FALSE)
  
  # yb ---------------------------------
  out <- get.bsub(
    object = object,
    basesub = basesub,
    comp0 = comp0,
    delta = delta,
    y0 = y0,
    d0 = d0,
    summary = summary,
    level = level,
    ref = ref)
  
}