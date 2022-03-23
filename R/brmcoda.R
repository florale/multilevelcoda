#' Computing brms model for compostion predicting an outcome
#' 
#' @param formula A object of class \code{formula}, \code{brmsformula}:
#' @param formula A symbolic description of the model to be fitted.
#' @param compilr A object of class \code{compilr} containing composition, ILR coordinates,
#' @param compilr used in the model.
#' @param data An object of class \code{data.frame} (or one that can be coerced to that class) 
#' @param data containing data of all variables other than ILR coordinates used in the model. 
#' @return A list with two elements
#' \itemize{
#'   \item{\code{CompIlr}}{ A object of class \code{compilr} used in the \code{brm} model. }
#'   \item{\code{Results}}{ An object of class \code{brmsfit}, which contains the posterior draws along with many other useful information about the model.}
#' @importFrom brms brm
#' @export
#' @examples
#' data(mcompd)
#' data(sbp)
#'
#' ## compute compositions and ILR coordinates
#' compilrtest <- compilr(data = mcompd[, 1:6], sbp = sbp, idvar = "ID")
#'
#' ## run brmcoda model
#' brmcodatest <- brmcoda(data = mcompd, compilr = compilrtest, 
#' formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID))
#'
#' ## clean-up
#' rm(mcompd, sbp, compilrtest, brmcodatest)
brmcoda <- function (compilr, formula, data) {
  
  bilr <- compilr$BetweenILR
  wilr <- compilr$WithinILR
  
  names(bilr) <- c(paste0("bilr", 1:ncol(bilr)))
  names(wilr) <- c(paste0("wilr", 1:ncol(wilr)))
  
  copyd <- cbind(data, bilr, wilr)
  
  m <- brm(eval(formula),
           data = copyd,
           chain = 4, core = 8)
  
  out <- list(
    CompIlr = compilr,
    Results = m)

}