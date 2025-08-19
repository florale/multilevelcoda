#' Between-person Average Substitution
#'
#' This function is an alias of \code{\link{substitution}} to estimates the the difference in an outcome
#' when compositional parts are substituted for specific unit(s) at \emph{between} level
#' using cluster mean (e.g., compositional mean at individual level) as reference composition. 
#' It is recommended that users run substitution model using the \code{\link{substitution}} function.
#' 
#' @seealso \code{\link{substitution}}
#' 
#' @inheritParams substitution
#' 
#' @inherit substitution return
#'
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#' cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp, 
#'                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)
#' 
#' m <- brmcoda(complr = cilr, 
#'              formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 + wilr1 + 
#'                                 wilr2 + wilr3 + wilr4 + Female + (1 | ID), 
#'              chains = 1, iter = 500,
#'              backend = "cmdstanr")
#'              
#' subm <- bsubmargins(object = m, base = psub, delta = 5)
#' }}
#' @export
bsubmargins <- function (object,
                         delta,
                         ref = "clustermean",
                         level = "between",
                         summary = TRUE,
                         at = NULL,
                         parts,
                         base,
                         weight = "proportional",
                         scale = c("response", "linear"),
                         cores = NULL,
                         type = "one-to-one",
                         ...) {
  
  ref   <- "clustermean"
  level <- "between"
  
  d0 <- build.rg(object = object,
                 ref = ref,
                 parts = parts,
                 level = level,
                 weight = weight,
                 fill = FALSE)
  
  # error if delta out of range
  x0 <- d0[, colnames(object$complr$between_comp), with = FALSE]
  
  delta <- as.integer(delta)
  if(isTRUE(any(all(delta) > lapply(x0, min)))) {
    stop(sprintf(
      "delta value should be less than or equal to %s, which is the amount of composition part available for pairwise substitution.",
  paste0(round(min(lapply(x0, min))), collapse = ", ")
    ))
  }
  
  # y0margins --------------------------------
  y0 <- fitted(
    object,
    newdata = d0,
    re_formula = NULL,
    scale = scale,
    summary = FALSE
  )
  y0 <- rowMeans(as.data.frame(y0)) # average across participants when there is no change
  
  # ybmargins ---------------------------------
  # substitution model
  out <- .get.bsubmargins(
    object = object,
    base = base,
    delta = delta,
    x0 = x0,
    d0 = d0,
    y0 = y0,
    level = level,
    ref = ref,
    summary = summary,
    scale = scale,
    type = type,
    cores = cores,
    ...
  )
}