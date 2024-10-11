#' Fit Bayesian generalised (non-)linear multilevel compositional model 
#' via full Bayesian inference
#' 
#' Fit a \code{brm} model with multilevel ILR coordinates
#' 
#' @param complr A \code{\link{complr}} object containing data of composition, 
#' ILR coordinates, and other variables used in the model.
#' @param formula A object of class \code{formula}, \code{brmsformula}:
#' A symbolic description of the model to be fitted. 
#' Details of the model specification can be found in \code{\link[brms:brmsformula]{brmsformula}}.
#' @param ... Further arguments passed to \code{\link[brms:brm]{brm}}.
#' 
#' @return A \code{\link{brmcoda}} with two elements
#'   \item{\code{complr}}{ An object of class \code{complr} used in the \code{brm} model. }
#'   \item{\code{model}}{ An object of class \code{brmsfit}, which contains the posterior draws 
#'   along with many other useful information about the model.}
#' @importFrom brms brm
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   cilr <- complr(data = mcompd, sbp = sbp,
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#'   
#'   # inspects ILRs before passing to brmcoda
#'   names(cilr$between_logratio)
#'   names(cilr$within_logratio)
#'   names(cilr$logratio)
#'   
#'   # model with compositional predictor at between and within-person levels
#'   m1 <- brmcoda(complr = cilr,
#'                 formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                                    wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   # model with compositional outcome
#'   m2 <- brmcoda(complr = cilr,
#'                 formula = mvbind(ilr1, ilr2, ilr3, ilr4) ~ Stress + Female + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   }}
#' @export
brmcoda <- function (complr, formula, ...) {
  
  if (isTRUE(missing(complr))) {
    stop(paste(
      "'complr' is a required argument and cannot be missing;",
      "  it should be an object of class complr.", 
      "  See ?multilevelcoda::complr for details.",
      sep = "\n"))
  }
  if (isFALSE(inherits(complr, "complr"))) {
    stop(paste(
      "complr must be an object of class complr.",
      "  See ?multilevelcoda::complr for details.",
      sep = "\n"))
  }
  
  tmp <- cbind(complr$data,
               complr$between_logratio,
               complr$within_logratio,
               complr$logratio)
  
  m <- brm(formula,
           data = tmp,
           ...)
  
  structure(
    list(complr = complr,
         model = m),
    class = "brmcoda")
}

#' Estimate pivot balance coordinates
#' 
#' @param object An object of class \code{brmcoda}.
#' @param ... Further arguments passed to \code{\link[brms:brm]{brm}}.
#' 
#' @return A list of \code{\link{brmcoda}} for each pivot balance coordinate.
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   cilr <- complr(data = mcompd, sbp = sbp,
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID",
#'                  total = 1440)
#'   
#'   # inspects ILRs before passing to brmcoda
#'   names(cilr$between_logratio)
#'   names(cilr$within_logratio)
#'   names(cilr$logratio)
#'   
#'   # model with compositional predictor at between and within-person levels
#'   m1 <- brmcoda(complr = cilr,
#'                 formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                                    wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   m_pb_1 <- brmcoda_pivot(m1)
#'   }}
#' @export
brmcoda_pivot <- function (object, ...) {
  
  out_d <- vector("list")
  
  # loop through parts
  for (d in object$complr$parts) {
    parts_d <- append(d, grep(d, object$complr$parts, value = T, invert = T))
    sbp_d   <- build.sbp(parts_d)
    sbp_d   <- sbp_d[, object$complr$parts]
    
    clr_d <- complr(data  = object$complr$data, 
                    sbp   = sbp_d,
                    parts = object$complr$parts,
                    idvar = object$complr$idvar,
                    total = object$complr$total)
    
    dat_d <-  cbind(clr_d$data,
                    clr_d$between_logratio,
                    clr_d$within_logratio,
                    clr_d$logratio)
    
    fit_d <- update(object$model,
                    newdata = dat_d,
                    ...)
    
    brmcoda_d <- structure(list(complr = clr_d,
                                model  = fit_d),
                           class = "brmcoda")
    
    out_d[[d]] <-   brmcoda_d
  }
  structure(out_d, class = "brmcoda_pivot")
}

