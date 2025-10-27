
#' Average Substitution
#'
#' This function is an alias of \code{\link{substitution}} to estimates the difference in an outcome
#' when compositional parts are substituted for specific unit(s)
#' using cluster mean (e.g., compositional mean at individual level) as reference composition. 
#' It is recommended that users run substitution model using the \code{\link{substitution}} function.
#' 
#' @seealso \code{\link{substitution}}
#' 
#' @inheritParams substitution
#' 
#' @inherit substitution return
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#' 
#' cilr <- complr(data = mcompd, sbp = sbp, 
#'                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)
#' 
#' # model with compositional predictor at between and within-person levels
#' m <- brmcoda(complr = cilr, 
#'              formula = Stress ~ z1_1 + z2_1 + z3_1 + z4_1 + (1 | ID), 
#'              chain = 1, iter = 500,
#'              backend = "cmdstanr")
#'                      
#' subm <- submargin(object = m, base = psub, delta = 5)
#' }}
#' @export
submargin <- function (object,
                        delta,
                        ref = "clustermean",
                        level = "aggregate",
                        summary = TRUE,
                        at = NULL,
                        parts = 1,
                        base,
                        type = "one-to-one",
                        weight = "proportional",
                        scale = c("response", "linear"),
                        cores = NULL,
                        ...) {
  
  ref <- "clustermean"
  level <- "aggregate"
  
  # if parts is numeric, get_parts
  if (is.numeric(parts)) {
    parts <- .get_parts(object[["complr"]], parts)
  }
  
  # d0 -------------------------------
  d0 <- build.rg(object = object,
                 ref = ref,
                 parts = parts,
                 level = level,
                 weight = weight,
                 fill = FALSE)
  
  # error if delta out of range
  x0 <- d0[, paste0("t", parts), with = FALSE]
  
  delta <- as.integer(delta)
  if(isTRUE(any(all(delta) > lapply(x0, min)))) {
    stop(sprintf(
      "delta value should be less than or equal to %s, which is
  the amount of composition part available for pairwise substitution.",
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

  # ymargins ---------------------------------
  # substitution model
  out <- .get.submargin(
    object = object,
    base = base,
    delta = delta,
    parts = parts,
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
