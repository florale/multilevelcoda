#' @title Between-person Average Marginal Substitution Model.
#'
#' @description
#' Using a fitted model object, estimates the the average marginal difference 
#' when compositional parts are substituted for specific unit(s) at `within-person` level. 
#' The \code{wsubmargins} output encapsulates 
#' the substitution results for all compositional parts
#' present in the \code{\link{brmcoda}} object.
#'
#' @param object A fitted \code{\link{brmcoda}} object. Required.
#' @param delta A integer, numeric value or vector indicating the amount of substituted change between compositional parts.
#' @param basesub A \code{data.frame} or \code{data.table} of the base possible substitution of compositional parts.
#' This data set can be computed using function \code{\link{basesub}}. 
#' If \code{NULL}, all possible pairwise substitution of compositional parts are used.
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
#' @importFrom compositions acomp ilr clo
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
#' subm <- wsubmargins(object = m, basesub = psub, delta = 5)
#' }
wsubmargins <- function (object,
                         delta,
                         basesub,
                         level = "within",
                         type = "marginal",
                         ...) {
  
  # between-person composition
  b <- object$CompILR$BetweenComp
  b <- as.data.table(clo(b, total = object$CompILR$total))
  
  # error if delta out of range
  if(isTRUE(any(delta > apply(b, 2, min)))) {
    stop(sprintf(
      "delta value should be less than or equal to %s, which is
  the amount of composition part available for pairwise substitution.",
  paste0(round(min(apply(b, 2, min))), collapse = ", ")
    ))
  }
  delta <- as.integer(delta)
  
  # model for no change
  bilr0 <- object$CompILR$BetweenILR
  wilr0 <- as.data.table(matrix(0, nrow = nrow(bilr0), ncol = ncol(bilr0)))
  
  colnames(bilr0) <- colnames(object$CompILR$BetweenILR)
  colnames(wilr0) <- colnames(object$CompILR$WithinILR)
  
  dref <- cbind(bilr0, wilr0, object$CompILR$data)
  yref <- fitted(object$Model, newdata = dref, re_formula = NA, summary = FALSE)
  yref <- rowMeans(yref) # average across participants when there is no change
  
  # substitution model
  out <- .get.wsubmargins(
    object = object,
    b = b,
    basesub = basesub,
    yref = yref,
    delta = delta,
    level = level,
    type = type)
}