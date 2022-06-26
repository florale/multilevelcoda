#' @title Estimated marginal means of composition.
#'
#' @description
#' This function obtains estimated marginal means (EMMs) of models
#' containing composition as the outcome
#' A sample from the posterior distribution of the regression coefficients.
#'
#' @param object A fitted \code{\link{mvmcoda}} model data. Required.
#' @param x A character vector specifying the names of the predictors over which EMMs are desired.
#' If \code{NULL}, the average of all predictionu is returned. Default to \code{NULL}.
#' @param at A numeric value, vector or list passed to \code{\link{ref_grid}}
#' specifying the predictor's values for EMMs.
#' @param ... Additional arguments to be passed to \code{\link{describe_posterior}}.
#'
#' @return A list with four elements.
#' @importFrom compositions ilrInv
#' @importFrom emmeans emmeans
#' @importFrom data.table as.data.table setDT rbindlist
#' @importFrom reshape2 dcast
#' @importFrom stats predict
#' @export
emmcoda <- function (object, x = NULL, at = NULL, ...) {
  
  if(isTRUE(missing(object))) {
    stop(paste(
      "argument 'object' is missing, with no default",
      "  it should be an object of class mvmcoda.",
      "  See ?emmcoda for details.",
      sep = "\n"))
  }

  if(isFALSE(inherits(object, "mvmcoda"))) {
    stop(sprintf(
    "Can't handle an object of class (%d)
  object should be a fitted 'mvmcoda' object
  See ?emmcoda for details.",
                 class(object)))
  }
  
  # extract variable information from model
  vn <- rbindlist(find_predictors(object$Model)) # get all varnames in brm model
  vn <- unique(vn)
  rg <- as.data.table(ref_grid(object$Model) @grid)
  
  if(isFALSE(is.null(at)) && is.null(x)) {
      stop("'at' is level of x, but 'x' is missing")
    
      } else if (is.null(at) && is.null(x)) {
        # no user's specified value - use default ref grid
        cv <- colnames(rg) %snin% c("rep.meas", ".wgt.") # get cov names, including x
        rg <- rg[, cv, with = FALSE]
        rg <- unique(rg)
        newd <- rg[, ID := 1]
        
        } else if (is.null(at) && isFALSE(is.null(x))) {
          if(isFALSE(x %in% colnames(rg))) { # check name
            stop("'x' is not a predictor in 'mvmcoda' object")
          } else {
            cv <- colnames(rg) %snin% c("rep.meas", ".wgt.")
            rg <- rg[, cv, with = FALSE]
            rg <- unique(rg)
            newd <- rg[, ID := 1]
          }
          
          } else { # both x and at are provided
            if(isFALSE(x %in% colnames(rg))) { 
              stop("'x' is not a predictor in 'mvmcoda' object")
            } else if (isFALSE(is(class(at), class(x)))) {
              stop(sprintf(
                "at (%s) must be values of x (%s)",
                class(at),
                class(rg[, get(x)])))
              } else {
                cv <- colnames(rg) %snin% c("rep.meas", ".wgt.", x)
                rg <- rg[, cv, with = FALSE]
                rg <- unique(rg)
                ls <- list(at)
                names(ls) <- x
                newd <- setDT(ls)[, ID := 1]
                
                iout <- vector("list", length(at))
                if (isFALSE(nrow(vn) == 1)) {
                  for (i in seq_len(nrow(rg))) {
                    nd <- cbind(newd, rg[i, ])
                    iout[[i]] <- nd
                    }
                  newd <- do.call(rbind, iout)
                }
              }
          }
  
  yhat <- as.data.table(fitted(object$Model, newdata = newd, re.form = NA, summary = FALSE))
  yhat <- dcast(yhat, V1 + V2 ~ V3, value.var = "value") 
  # V1 - posterior draws of each row in newd, V2 - row in newd, V3 - ilr
  yhat <- yhat[order(yhat$V2),]
  # check yhat with JW
  # inverse ilr to get composition
  ilr <- yhat[, -c(1:2)]
  
  if(is.null(at)) { # average all prediction
    comp <- ilrInv(ilr, V = object$CompIlr$psi)
    comp <- as.data.table(clo(comp, total = object$CompIlr$total))
    names(comp) <- object$CompIlr$parts
    comp <- describe_posterior(comp, centrality = "mean")

    ilr <- describe_posterior(ilr, centrality = "mean")
    
  } else {
    ilr <- split(ilr, rep(1:nrow(newd), each = nrow(yhat) / nrow(newd)))
    
    comp <- lapply(ilr, function(x) {
      x <- ilrInv(x, V = object$CompIlr$psi)
      x <- as.data.table(clo(x, total = object$CompIlr$total))
      names(x) <- object$CompIlr$parts
      x <- describe_posterior(x, centrality = "mean")}) # if want to use ROPE, need to manually specify
    ilr <- lapply(ilr, function(x) {describe_posterior(x, centrality = "mean")})
    
    d <- split(newd, rep(1:nrow(newd), each = 1))
    
    dcomp <- lapply(d, function(x){x[rep(seq_len(nrow(x)), length(object$CompIlr$parts)), ]})
    dilr <- lapply(d, function(x){x[rep(seq_len(nrow(x)), ncol(object$CompIlr$TotalILR)), ]})
    
    comp <- Map(cbind, comp, dcomp)
    ilr <- Map(cbind, ilr, dilr)
    
    names(comp) <- paste0(x, at)
    names(ilr) <- paste0(x, at)
  }
  out <- list(comp, ilr)
  names(out) <- c("Composition", "ILR")
  out
}
