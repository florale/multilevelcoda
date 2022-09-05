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
#' @param minute A integer or numeric value indicating the maximum minute 
#' for which substitution model is desired.
#' Default to \code{60L} (i.e., the model loops through 1:60L minutes).
#' @param ... Additional arguments to be passed to \code{\link{describe_posterior}}.
#' 
#' @return A list containing the result of isotemporal multilevel substitution model.
#' Each elements of the list is the substitution estimation for a compositional variables, 
#' which include at least six core items.
#' \itemize{
#'   \item{\code{Mean}}{ Posterior means.}
#'   \item{\code{CI_low}} and \item{\code{CI_high}}{ 95% credible intervals.}
#'   \item{\code{MinSubstituted}}{ Minute substituted within the composition.}
#'   \item{\code{Substitute}}{Compositional variables to be substituted from/to.}
#'   \item{\code{Predictor}}{Central compositional variable to be substituted from/to.}
#' }
#'
#' @importFrom data.table as.data.table copy := setnames rbindlist
#' @importFrom compositions acomp ilr clo
#' @importFrom extraoperators %snin% %sin%
#' @importFrom bayestestR describe_posterior
#' @importFrom foreach %dopar%
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
#' subm <- bsubmargins(object = m, substitute = psub, minute = 5)
#' }
bsubmargins <- function (object, substitute, minute = 60L, ...) {
  
  if (isTRUE(missing(object))) {
    stop(paste(
      "'object' is a required argument and cannot be missing;",
      "  it should be an object of class 'brmcoda'.", 
      "  See ?bsubmargins for details.",
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
      "  See ?bsubmargins for details.",
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
    stop(sprintf(
    "The number of columns in 'substitute' (%d) must be the same
  as the compositional variables in 'parts' (%d).",
                 ncol(substitute),
                 length(object$CompIlr$parts)))
  }
  
  if (isFALSE(identical(colnames(substitute), object$CompIlr$parts))) {
    stop(sprintf(
    "The names of compositional variables must be the same
  in 'substitute' (%s) and 'parts' (%s).",
                 colnames(substitute),
                 object$CompIlr$parts))
  }
  
  # between-person composition
  b <- object$CompIlr$BetweenComp
  b <- as.data.table(clo(b, total = object$CompIlr$total))

  min <- as.integer(minute)
  
  # model for no change
  bilr <- object$CompIlr$BetweenILR
  wilr <- as.data.table(matrix(0, nrow = nrow(bilr), ncol = ncol(bilr))) # check with JW
  
  colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
  colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
  
  samed <- cbind(bilr, wilr, object$CompIlr$data)
  ysame <- fitted(object$Model, newdata = samed, re.form = NA, summary = FALSE)
  ysame <- rowMeans(ysame) # average across participants when there is no change
  
  # substitution model
  out <- .get.bsubmargins(object = object, b = b,
                          substitute = substitute,
                          ysame = ysame, min = min)
  
}