#' Between-person Simple Substitution
#' 
#' This function is an alias of \code{\link{substitution}} to estimates the difference in an outcome
#' when compositional parts are substituted for specific unit(s) at \emph{between} level
#' using a single reference composition (e.g., compositional mean at sample level).
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
#' cilr <- complr(data = mcompd, sbp = sbp, 
#'                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)
#' 
#' # model with compositional predictor at between and between-person levels
#' m <- brmcoda(complr = cilr, 
#'              formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 + 
#'                                 wz1_1 + wz2_1 + wz3_1 + wz4_1 + Female + (1 | ID), 
#'              chain = 1, iter = 500,
#'              backend = "cmdstanr")
#' subm <- bsub(object = m, base = psub, delta = 5)
#' }}
#' @export
bsub <- function(object,
                 delta,
                 ref = "grandmean",
                 level = "between",
                 summary = TRUE,
                 aorg = TRUE,
                 at = NULL,
                 parts = 1,
                 base,
                 type = "one-to-one",
                 weight = "equal",
                 scale = c("response", "linear"),
                 cores = NULL,
                 ...) {
  
  # ref <- "grandmean"
  level <- "between"
  
  # if parts is numeric, get_parts
  if (is.numeric(parts)) {
    parts <- .get_parts(object[["complr"]], parts)
  }
  
  # d0 -------------------------------
  if (isTRUE(ref == "grandmean")) {
    d0 <- build.rg(object = object,
                   ref = ref,
                   parts = parts,
                   at = at,
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

  # error if delta out of range
  x0 <- d0[1, paste0("b", parts), with = FALSE]
  
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
  
  # yb ---------------------------------
  out <- .get.bsub(
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
  
  return(out)
}