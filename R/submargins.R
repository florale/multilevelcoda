#' Average Marginal Substitution
#'
#'
#' This function is an alias of \code{\link{substitution}} to estimates the the difference in an outcome
#' when compositional parts are substituted for specific unit(s)
#' using cluster mean (e.g., compositional mean at individual level) as reference composition. 
# #' It is recommended that users run substitution model using the \code{\link{substitution}} function.
#' 
#' @inheritParams substitution
#' 
#' @seealso \code{\link{substitution}}
#' 
#' @inherit substitution return
#'
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo
#' 
# #' @examples
# #' \donttest{
# #' if(requireNamespace("cmdstanr")){
# #' cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp, 
# #'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)
# #' 
# #' m <- brmcoda(compilr = cilr, 
# #'              formula = Stress ~ ilr1 + ilr2 + ilr3 + ilr4 + (1 | ID), 
# #'              chains = 1, iter = 500,
# #'              backend = "cmdstanr")
# #'              
# #' subm <- submargins(object = m, basesub = psub, delta = 5)
# #' }}
# #' @export
submargins <- function (object,
                        delta,
                        basesub,
                        ref = "clustermean",
                        level = "total",
                        weight = NULL,
                        ...) {
  
  # full composition
  comp0 <- object$CompILR$TotalComp
  comp0 <- as.data.table(clo(comp0, total = object$CompILR$total))
  
  # error if delta out of range
  if(isTRUE(any(delta > apply(comp0, 2, min)))) {
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
  
  # y0margins --------------------------------
  y0 <- fitted(
    object,
    newdata = d0,
    re_formula = NULL,
    summary = FALSE
  )
  y0 <- rowMeans(as.data.frame(y0)) # average across participants when there is no change
  
  # substitution model
  out <- .get.submargins(
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