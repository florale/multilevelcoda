#' @title Within-person Average Substitution.
#'
#' @description
#' Using a fitted model object, estimates the the average marginal difference 
#' when compositional parts are substituted for specific unit(s) at `within` level. 
#' The \code{wsubmargins} output encapsulates 
#' the substitution results for all compositional parts
#' present in the \code{\link{brmcoda}} object.
#'
#' @param object A fitted \code{\link{brmcoda}} object. Required.
#' @param delta A integer, numeric value or vector indicating the amount of substituted change between compositional parts.
#' @param basesub A \code{data.frame} or \code{data.table} of the base possible substitution of compositional parts.
#' This data set can be computed using function \code{\link{basesub}}. 
#' If \code{NULL}, all possible pairwise substitution of compositional parts are used.
#' @param ref A character string. Default to \code{clustermean}.
#' @param level A character string. Default to \code{within}.
#' @param weight A character value specifying the weight to use in calculation of the reference composition.
#' \code{weight} can be \code{equal} which gives equal weight to units (e.g., individuals) or
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
#' @importFrom compositions acomp ilr clo
#' @importFrom stats fitted
#' @export
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
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
#' subm <- wsubmargins(object = m, basesub = psub, delta = 5)
#' }}
wsubmargins <- function (object,
                         delta,
                         basesub,
                         ref = "clustermean",
                         level = "within",
                         weight = NULL,
                         ...) {
  
  d0 <- build.rg(object = object,
                 ref = ref,
                 weight = weight,
                 fill = FALSE)
  
  # error if delta out of range
  comp0 <- d0[, colnames(object$CompILR$BetweenComp), with = FALSE]
  
  delta <- as.integer(delta)
  if(isTRUE(any(all(delta) > lapply(comp0, min)))) {
    stop(sprintf(
      "delta value should be less than or equal to %s, which is
  the amount of composition part available for pairwise substitution.",
  paste0(round(min(lapply(comp0, min))), collapse = ", ")
    ))
  }

  # y0margins --------------------------------
  y0 <- fitted(
    object$Model,
    newdata = d0,
    re_formula = NA,
    summary = FALSE
  )
  y0 <- rowMeans(y0) # average across participants when there is no change
  
  # ywmargins ---------------------------------
  # substitution model
  out <- .get.wsubmargins(
    object = object,
    delta = delta,
    basesub = basesub,
    comp0 = comp0,
    d0 = d0,
    y0 = y0,
    level = level,
    ref = ref
  )
}