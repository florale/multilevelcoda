#' Estimating Average Marginal Effects for Between-person Substitution Model.
#'
#' Using a fitted model object, estimates the average marginal difference 
#' in outcomes when compositional variables are substituted for a specific period
#' at `between-person` level. The resulting \code{bsubmargins} encapsulates 
#' substitution estimation across all compositional variables present
#' in the \code{\link{brmcoda}} object.
#'
#' @param object A \code{\link{brmcoda}} object.
#' @param substitute A \code{data.frame} or \code{data.table} of possible substitution of variables.
#' This dataset can be computed using function \code{possub}. Required.
#' @param delta A integer, numeric value or vector indicating the amount of change in compositional parts
#' for which substitution model is desired.
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
#'
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo
#' @importFrom stats fitted
#' @export
#' @examples
#' \dontrun{
#' data(psub)
#' data(mcompd)
#' data(psub)
#' cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#' 
#' m <- brmcoda(compilr = cilr, 
#'              formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + wilr1 + 
#'              wilr2 + wilr3 + wilr4 + Female + (1 | ID), chains = 1, iter = 500)
#'                
#' subm <- bsubmargins(object = m, substitute = psub, delta = 5)
#' }
bsubmargins <- function (object, substitute, delta, ...) {
  
  # between-person composition
  b <- object$CompIlr$BetweenComp
  b <- as.data.table(clo(b, total = object$CompIlr$total))

  delta <- as.integer(delta)
  
  # model for no change
  bilr <- object$CompIlr$BetweenILR
  wilr <- as.data.table(matrix(0, nrow = nrow(bilr), ncol = ncol(bilr))) # check with JW
  
  colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
  colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
  
  samed <- cbind(bilr, wilr, object$CompIlr$data)
  ysame <- fitted(object$Model, newdata = samed, re_formula = NA, summary = FALSE)
  ysame <- rowMeans(ysame) # average across participants when there is no change
  
  # substitution model
  out <- .get.bsubmargins(object = object, b = b,
                          substitute = substitute,
                          ysame = ysame, delta = delta)
  
}