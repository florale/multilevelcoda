#' Draws from the Posterior Predictive Distribution
#' 
#' Compute posterior draws of the posterior predictive distribution. Can be
#' performed for the data used to fit the model (posterior predictive checks) or
#' for new data. By definition, these draws have higher variance than draws
#' of the expected value of the posterior predictive distribution computed by
#' \code{\link{posterior_epred.brmsfit}}. This is because the residual error
#' is incorporated in \code{posterior_predict}. However, the estimated means of
#' both methods averaged across draws should be very similar.
#' 
#' @param object An object of class \code{brmcoda}.
#' @param newdata An optional data.frame for which to evaluate predictions. If
#' \code{NULL} (default), the original data of the model is used.
#' \code{NA} values within factors are interpreted as if all dummy
#' variables of this factor are zero. This allows, for instance, to make
#' predictions of the grand mean when using sum coding.
#' @param re_formula formula containing group-level effects to be considered in
#' the prediction. If \code{NULL} (default), include all group-level effects;
#' if \code{NA}, include no group-level effects.
#' @param scale Either \code{"response"} or \code{"linear"} or \code{"comp"}.
#' If either \code{"response"} or \code{"linear"},
#' results are returned using the corresponding method in \code{predict.brmsfit}.
#' \code{"comp"} is only relevant to multivariate models of class \code{mvbrmsformula}
#' (i.e., when response is composition). If \code{"comp"} is specified,
#' results are returned on the compositional scale of the response variable.
#' @param summary Should summary statistics be returned instead of the raw values?
#' Default is \code{TRUE}.
#' @param ... Further arguments passed to \code{\link{predict.brmsfit}}
#' that control additional aspects of prediction.
#' 
#' @importFrom compositions ilrInv
#' @importFrom brms posterior_summary do_call
#' @importFrom abind abind
#' @method predict brmcoda
#' @export
predict.brmcoda <- function(object,
                            newdata = NULL,
                            re_formula = NULL,
                            scale = c("response", "linear", "comp"),
                            summary = TRUE,
                            ...) {
  
  if (missing(scale)) {
    out <- predict(
      object$Model,
      newdata = newdata,
      re_formula = re_formula,
      summary = summary,
      ...
    )
  } else if (scale %in% c("response", "linear")) {
    out <- predict(
      object$Model,
      newdata = newdata,
      re_formula = re_formula,
      scale = scale,
      summary = summary,
      ...
    )
  } else if (scale == "comp") {
    if (isFALSE(inherits(object$Model$formula, "mvbrmsformula"))) {
      stop(sprintf(
        "This is a %s model, but a model of class mvbrmsformula is required to
  return results on compsitional scale.",
  class(object$Model$formula)[1]))
      
    } else {
      out <- predict(
        object$Model,
        newdata = newdata,
        re_formula = re_formula,
        summary = FALSE,
        ...
      )
      out <- lapply(asplit(out, 1), function(x) {
        x <- compositions::ilrInv(x, V = gsi.buildilrBase(t(object$CompILR$sbp)))
        as.data.table(clo(x, total = object$CompILR$total))
      })
      
      out <- brms::do_call(abind::abind, c(out, along = 3))
      out <- aperm(out, c(3, 1, 2)) #draw-row-col
      
      if(isTRUE(summary)) {
        out <- brms::posterior_summary(out)
        dimnames(out)[[3]] <- object$CompILR$parts
      }
    }
  }
  out
}

#' Expected Values of the Posterior Predictive Distribution
#' 
#' Compute posterior draws of the expected value of the posterior predictive
#' distribution. Can be performed for the data used to fit the model (posterior
#' predictive checks) or for new data. By definition, these predictions have
#' smaller variance than the posterior predictions performed by the
#' \code{\link{posterior_predict.brmsfit}} method. This is because only the
#' uncertainty in the expected value of the posterior predictive distribution is
#' incorporated in the draws computed by \code{posterior_epred} while the
#' residual error is ignored there. However, the estimated means of both methods
#' averaged across draws should be very similar.
#' 
#' @inheritParams predict.brmcoda
#' @param ... Further arguments passed to \code{\link{fitted.brmsfit}}
#' that control additional aspects of prediction.
#' 
#' @importFrom compositions ilrInv
#' @importFrom brms posterior_summary do_call
#' @importFrom abind abind
#' @method fitted brmcoda
#' @export
fitted.brmcoda <- function(object,
                           newdata = NULL,
                           re_formula = NULL,
                           scale = c("response", "linear", "comp"),
                           summary = TRUE,
                           ...) {
  
  if (missing(scale)) {
    out <- fitted(
      object$Model,
      newdata = newdata,
      re_formula = re_formula,
      summary = summary,
      ...
    )
  } else if (scale %in% c("response", "linear")) {
    out <- fitted(
      object$Model,
      newdata = newdata,
      re_formula = re_formula,
      scale = scale,
      summary = summary,
      ...
    )
  } else if (scale == "comp") {
    if (isFALSE(inherits(object$Model$formula, "mvbrmsformula"))) {
      stop(sprintf(
        "This is a %s model, but a model of class mvbrmsformula is required to
  return results on compsitional scale.",
  class(object$Model$formula)[1]))
      
    } else {
      out <- fitted(
        object$Model,
        newdata = newdata,
        re_formula = re_formula,
        summary = FALSE,
        ...
      )
      out <- lapply(asplit(out, 1), function(x) {
        x <- compositions::ilrInv(x, V = gsi.buildilrBase(t(object$CompILR$sbp)))
        as.data.table(clo(x, total = object$CompILR$total))
      })
      
      out <- brms::do_call(abind::abind, c(out, along = 3))
      out <- aperm(out, c(3, 1, 2)) #draw-row-col
      
      if(isTRUE(summary)) {
        out <- brms::posterior_summary(out)
        dimnames(out)[[3]] <- object$CompILR$parts
      }
    }
  }
  out
}

#' Extract Population-Level Estimates
#' 
#' Extract the population-level ('fixed') effects
#' from a \code{brmsfit} in the \code{brmcoda} object.
#' 
#' @aliases fixef
#' 
#' @inheritParams predict.brmcoda
#' 
#' @return If \code{summary} is \code{TRUE}, a matrix returned
#' by \code{\link{posterior_summary}} for the population-level effects.
#' If \code{summary} is \code{FALSE}, a matrix with one row per
#' posterior draw and one column per population-level effect.
#' 
#' @importFrom brms fixef
#' @method fixef brmcoda
#' @export fixef
#' @export
fixef.brmcoda <- function(object, ...) {
  fixef(object$Model, ...)
}

#' Covariance and Correlation Matrix of Population-Level Effects
#'
#' Get a point estimate of the covariance or
#' correlation matrix of population-level parameters
#'
#' @inheritParams predict.brmcoda
#' 
#' @param correlation Logical; if \code{FALSE} (the default), compute
#'   the covariance matrix, if \code{TRUE}, compute the correlation matrix.
#' 
#' @return covariance or correlation matrix of population-level parameters
#' 
#' @method vcov brmcoda
#' @export
vcov.brmcoda <- function(object, ...) {
  vcov(object$Model, ...)
}

#' Extract Group-Level Estimates
#'
#' Extract the group-level ('random') effects of each level
#' 
#' @aliases ranef
#' 
#' @inheritParams predict.brmcoda
#' 
#' @return A list of 3D arrays (one per grouping factor).
#' If \code{summary} is \code{TRUE},
#' the 1st dimension contains the factor levels,
#' the 2nd dimension contains the summary statistics
#' (see \code{\link{posterior_summary}}), and
#' the 3rd dimension contains the group-level effects.
#' If \code{summary} is \code{FALSE}, the 1st dimension contains
#' the posterior draws, the 2nd dimension contains the factor levels,
#' and the 3rd dimension contains the group-level effects.
#'
#' @importFrom brms ranef
#' @method ranef brmcoda
#' @export ranef
#' @export
ranef.brmcoda <- function(object, ...) {
  ranef(object$Model, ...)
}

#' Extract Model Coefficients
#'
#' Extract model coefficients, which are the sum of population-level
#' effects and corresponding group-level effects
#' 
#' @aliases coef
#' 
#' @inheritParams predict.brmcoda
#' 
#' @return A list of 3D arrays (one per grouping factor).
#' If \code{summary} is \code{TRUE},
#' the 1st dimension contains the factor levels,
#' the 2nd dimension contains the summary statistics
#' (see \code{\link{posterior_summary}}), and
#' the 3rd dimension contains the group-level effects.
#' If \code{summary} is \code{FALSE}, the 1st dimension contains
#' the posterior draws, the 2nd dimension contains the factor levels,
#' and the 3rd dimension contains the group-level effects.
#'
#' @method coef brmcoda
#' @export
coef.brmcoda <- function(object, ...) {
  coef(object$Model, ...)
}

#' Extract Variance and Correlation Components
#'
#' This function calculates the estimated standard deviations,
#' correlations and covariances of the group-level terms
#' in a multilevel model of class \code{brmsfit}.
#' For linear models, the residual standard deviations,
#' correlations and covariances are also returned.
#'
#' @aliases VarCorr
#' 
#' @return A list of lists (one per grouping factor), each with
#' three elements: a matrix containing the standard deviations,
#' an array containing the correlation matrix, and an array
#' containing the covariance matrix with variances on the diagonal.
#'
#' @importFrom brms VarCorr
#' @method VarCorr brmcoda
#' @export VarCorr
#' @export
VarCorr.brmcoda <- function(object, ...) {
  VarCorr(object$Model, ...)
}