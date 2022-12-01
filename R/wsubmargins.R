#' @title Between-person Marginal Substitution Model.
#'
#' @description
#' Using a fitted model object, estimates the the average marginal difference 
#' when compositional variables are substituted for a specific amount
#' at `within-person` level. The resulting \code{wsubmargins} encapsulates 
#' substitution estimation across all compositional variables present
#' in the \code{\link{brmcoda}} object.
#'
#' @param object A \code{\link{brmcoda}} object. Required.
#' @param base A \code{data.frame} or \code{data.table} of the possible substitution of variables.
#' This dataset can be computed using function \code{possub}. Required.
#' @param delta A integer, numeric value or vector indicating the amount of change in compositional parts
#' for substitution.
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
#' data(mcompd)
#' data(sbp)
#' data(psub)
#' cilr <- compilr(data = mcompd, sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#'                 
#' library(doFuture)
#' registerDoFuture()
#' plan(multisession, workers = 5)
#' 
#' # model with compositional predictor at between and within-person levels
#' m <- brmcoda(compilr = cilr, 
#'              formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + 
#'                                 wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID), 
#'              chain = 1, iter = 500)             
#' subm <- wsubmargins(object = m, base = psub, delta = 5)
#' }
wsubmargins <- function (object, base, delta, 
                         level, type,
                         ...) {
  
  # between-person composition
  b <- object$CompIlr$BetweenComp
  b <- as.data.table(clo(b, total = object$CompIlr$total))
  
  delta <- as.integer(delta)
  
  # model for no change
  bilr <- object$CompIlr$BetweenILR
  wilr <- as.data.table(matrix(0, nrow = nrow(bilr), ncol = ncol(bilr)))
  
  colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
  colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
  
  samed <- cbind(bilr, wilr, object$CompIlr$data)
  ysame <- fitted(object$Model, newdata = samed, re_formula = NA, summary = FALSE)
  ysame <- rowMeans(ysame) # average across participants when there is no change
  
  # substitution model
  out <- .get.wsubmargins(object = object, b = b,
                          base = base,
                          ysame = ysame, delta = delta, 
                          level = level, type = type)
}