#' Update \code{\link{brmcoda}} models
#'
#' This method allows for updating an existing \code{\link{brmcoda}} object.
#'
#' @param object A fitted \code{\link{brmcoda}} object to be updated.
#' @param formula. Changes to the formula; for details see
#' \code{\link[stats:update.formula]{update.formula}} and \code{\link[brms:brmsformula]{brmsformula}}.
#' @param newdata A \code{data.frame} or \code{data.table}
#' containing data of all variables used in the analysis.
#' It must include a composition and the same ID variable as the existing \code{\link{complr}} object.
#' @param ... Further arguments passed to \code{\link[brms:brm]{brm}}.
#'
#' @inherit brmcoda return
#'
#' @seealso \code{\link{brmcoda}}
#'
#' @method update brmcoda
#'
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'
#' # model with compositional predictor at between and within-person levels
#' fit <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID"),
#'                 formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                    wz1_1 + wz2_1 + wz3_1 + wz4_1 + Female + (1 | ID),
#'                 chain = 1, iter = 500,
#'               backend = "cmdstanr")
#'
#' # removing the effect of bz1_1
#' fit1 <- update(fit, formula. = ~ . - bz1_1)
#'
#' # using only a subset
#' fit2 <- update(fit, newdata = mcompd[ID != 1])
#' }}
#' @export
update.brmcoda <- function(object,
                           formula. = NULL,
                           newdata = NULL,
                           ...) {
  if (!is.null(newdata)) {
    update_complr   <- complr(
      data  = object$complr$datain,
      sbp   = lapply(object$complr$output, function(x)
        x$sbp),
      parts = lapply(object$complr$output, function(x)
        x$parts),
      idvar = if (!is.null(object$complr$idvar))
        (object$complr$idvar)
      else
        NULL,
      total = lapply(object$complr$output, function(x)
        x$total),
    )
    
    update_fit <- update(object$model,
                         formula. = formula.,
                         newdata = update_complr$dataout,
                         ...)
    
    structure(list(complr = update_complr, model = update_fit), class = "brmcoda")
  } else {
    update_fit <- update(object$model, formula. = formula., ...)
    
    structure(list(complr = object$complr, model = update_fit), class = "brmcoda")
  }
}
