#' Simple Substitution
#' 
#' This function is an alias of \code{\link{substitution}} to estimates the the difference in an outcome
#' when compositional parts are substituted for specific unit(s) 
#' using a aggregate reference composition 
#' (e.g., compositional mean at sample level, not seperated by between- and within effects).
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
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)
#' 
#' # model with compositional predictor at between and within-person levels
#' m <- brmcoda(complr = cilr, 
#'              formula = Stress ~ z1 + z2 + z3 + z4 + (1 | ID), 
#'              chain = 1, iter = 500,
#'              backend = "cmdstanr")
#'              
#' subm <- sub(object = m, base = psub, delta = 5)
#' }}
#' @export
sub <- function (object,
                 delta,
                 ref = "grandmean",
                 level = "aggregate",
                 summary = TRUE,
                 aorg = TRUE,
                 at = NULL,
                 parts,
                 base,
                 weight = "equal",
                 scale = c("response", "linear"),
                 cores = NULL,
                 type = "one-to-one",
                 ...) {
  
  level <- "aggregate"
  
  # d0 -------------------------------
  if (isTRUE(ref == "grandmean")) {
    d0 <- build.rg(object = object,
                   ref = ref,
                   parts = parts,
                   level = level,
                   weight = weight,
                   fill = FALSE)
  } else {
    if (isFALSE(inherits(ref, c("data.table", "data.frame", "matrix")))) {
      stop("ref must be \"grandmean\" or a data table, data frame or matrix.")
    }
    if(isFALSE(  # ensure all covs are provided
      (colnames(as.data.table(ref_grid(object$model)@grid)) %snin% ".wgt.") %ain% colnames(ref))) {
      stop(paste(
        "'ref' should contains information about",
        "  the covariates in 'brmcoda' model to estimate the substitution model.",
        "  Please provide a different reference grid or build one using `build.rg()`.",
        sep = "\n"))
    }
    d0 <- ref
    ref <- "users"
  }
  d0 <- as.data.table(d0)
  # 
  # error if delta out of range
  x0 <- d0[1, paste0("t", parts), with = FALSE]
  
  delta <- as.integer(delta)
  if(isTRUE(any(delta > min(x0)))) {
    stop(sprintf(
      "delta value should be less than or equal to %s, which is the amount of composition part available for pairwise substitution.",
      round(min(x0), 2)
    ))
  }
  
  # y0 --------------------------------
  y0 <- fitted(
    object,
    newdata = d0,
    re_formula = NA,
    scale = scale,
    summary = FALSE
  )
  
  # y ---------------------------------
  out <- .get.sub(
    object = object,
    base = base,
    delta = delta,
    parts = parts,
    x0 = x0,
    d0 = d0,
    y0 = y0,
    at = at,
    level = level,
    ref = ref,
    aorg = aorg,
    summary = summary,
    scale = scale,
    type = type,
    cores = cores,
    ...
  )
}

#' Average Substitution
#'
#' This function is an alias of \code{\link{substitution}} to estimates the the difference in an outcome
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
#'              formula = Stress ~ z1 + z2 + z3 + z4 + (1 | ID), 
#'              chain = 1, iter = 500,
#'              backend = "cmdstanr")
#'                      
#' subm <- submargins(object = m, base = psub, delta = 5)
#' }}
#' @export
submargins <- function (object,
                        delta,
                        ref = "clustermean",
                        level = "aggregate",
                        summary = TRUE,
                        at = NULL,
                        parts,
                        base,
                        weight = "proportional",
                        scale = c("response", "linear"),
                        cores = NULL,
                        type = "one-to-one",
                        ...) {
  
  ref <- "clustermean"
  level <- "aggregate"
  
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
  y0 <- rowMeans(as.data.frame(y0)) # average across participants when there is no change
  
  # ymargins ---------------------------------
  # substitution model
  out <- .get.submargins(
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
