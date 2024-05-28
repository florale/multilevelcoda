#' Within-person Simple Substitution
#' 
#' This function is an alias of \code{\link{substitution}} to estimates the the difference in an outcome
#' when compositional parts are substituted for specific unit(s) at \emph{within} level
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
#' 
#' cilr <- complr(data = mcompd, sbp = sbp, 
#'                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)
#' 
#' # model with compositional predictor at between and within-person levels
#' m <- brmcoda(complr = cilr, 
#'              formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 + 
#'                                 wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID), 
#'              chain = 1, iter = 500,
#'              backend = "cmdstanr")
#'              
#' subm <- wsub(object = m, basesub = psub, delta = 60)
#' }}
#' @export
wsub <- function(object,
                 basesub,
                 delta,
                 summary = TRUE,
                 ref = "grandmean",
                 level = "within",
                 weight = "equal",
                 scale = c("response", "linear"),
                 cores = NULL,
                 ...) {
  
  # ref <- "grandmean"
  level <- "within"
  
  # d0 -------------------------------
  if (isTRUE(ref == "grandmean")) {
    d0 <- build.rg(object = object,
                   ref = ref,
                   level = level,
                   weight = weight,
                   fill = FALSE)
  } else {
    if (isFALSE(inherits(ref, c("data.table", "data.frame", "matrix")))) {
      stop("ref must be \"grandmean\" or a data table, data frame or matrix.")
    }
    if(isFALSE(  # ensure all covs are provided
      (colnames(as.data.table(insight::get_datagrid(model.frame(object),
                                                    at = paste0(object$model$formula$formula[[2]]),
                                                    length = NA)))) %ain% colnames(ref))) {
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
  comp0 <- d0[1, colnames(object$complr$between_comp), with = FALSE]
  
  delta <- as.integer(delta)
  if(isTRUE(any(delta > min(comp0)))) {
    stop(sprintf(
      "delta value should be less than or equal to %s, which is the amount of composition part available for pairwise substitution.",
  round(min(comp0), 2)
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
  
  # yw ---------------------------------
    # substitution model
    out <- .get.wsub(
      object = object,
      basesub = basesub,
      delta = delta,
      comp0 = comp0,
      d0 = d0,
      y0 = y0,
      level = level,
      ref = ref,
      summary = summary,
      scale = scale,
      cores = cores
    )

}