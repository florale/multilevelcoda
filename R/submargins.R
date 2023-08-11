#' Average Marginal Substitution
#'
#' Using a fitted model object, estimates the the average marginal difference 
#' when compositional parts are substituted for specific unit(s). 
#' The \code{submargins} output encapsulates 
#' the substitution results for all compositional parts
#' present in the \code{\link{brmcoda}} object.
#'
#' @param object A fitted \code{\link{brmcoda}} object.
#' @param delta A integer, numeric value or vector indicating the amount of substituted change between compositional parts.
#' @param basesub A \code{data.frame} or \code{data.table} of the base possible substitution of compositional parts.
#' This data set can be computed using function \code{\link{basesub}}. 
#' If \code{NULL}, all possible pairwise substitution of compositional parts are used.
#' @param level A character string or vector. Default to \code{total}. 
#' @param type A character string or vector. 
#' Should the estimate be \code{conditional} mean or average \code{marginal} mean?
#' @param ... Additional arguments to be passed to \code{\link{describe_posterior}}.
#' 
#' @return A list containing the result of multilevel compositional substitution model.
#' Each element of the list is the estimation for a compositional part 
#' and include at least six elements.
#' \itemize{
#'   \item{\code{Mean}}{ Posterior means.}
#'   \item{\code{CI_low} and \code{CI_high}}{ 95% credible intervals.}
#'   \item{\code{Delta}}{ Amount substituted across compositional parts.}
#'   \item{\code{From}}{ Compositional part that is substituted from.}
#'   \item{\code{To}}{ Compositional parts that is substituted to.}
#'   \item{\code{Level}}{ Level where changes in composition takes place.}
#'   \item{\code{EffectType}}{ Either estimated `conditional` or average `marginal` changes.}
#' }
#'
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo
submargins <- function(object,
                       delta,
                       basesub,
                       level = "total",
                       type = "marginal",
                       ...) {
  
  # full composition
  t <- object$CompILR$TotalComp
  t <- as.data.table(clo(t, total = object$CompILR$total))
  
  # error if delta out of range
  if(isTRUE(any(delta > apply(t, 2, min)))) {
    stop(sprintf(
      "delta value should be less than or equal to %s, which is
  the amount of composition part available for pairwise substitution.",
  paste0(round(min(apply(t, 2, min))), collapse = ", ")
    ))
  }
  delta <- as.integer(delta)
  
  # model for no change
  tilr0 <- object$CompILR$TotalILR
  
  d0 <- cbind(tilr0, object$CompILR$data)
  y0 <- fitted(object$Model, newdata = d0, re_formula = NA, summary = FALSE)
  y0 <- rowMeans(y0) # average across participants when there is no change
  
  # substitution model
  out <- .get.submargins(
    object = object,
    t = t,
    basesub = basesub,
    y0 = y0,
    delta = delta,
    level = level,
    type = type)
}