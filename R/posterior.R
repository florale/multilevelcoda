#' Draws from the Posterior Predictive Distribution
#' 
#' Compute posterior draws of the posterior predictive distribution
#' of a \code{brmsfit} model in the \code{brmcoda} object.
#' Can be performed for the data used to fit the model (posterior predictive checks) or
#' for new data. By definition, these draws have higher variance than draws
#' of the expected value of the posterior predictive distribution computed by
#' \code{\link{fitted.brmcoda}}. This is because the residual error
#' is incorporated in \code{posterior_predict}. However, the estimated means of
#' both methods averaged across draws should be very similar.
#' 
#' @aliases predict
#' 
#' @param object An object of class \code{brmcoda}.
#' @param acomp Should the
#' results be returned on the compositional scale of the response variable?
#' Only applicable for models with compositional response
#' @param ... Further arguments passed to \code{\link{predict.brmsfit}}
#' that control additional aspects of prediction.
#' @inheritParams brms::predict.brmsfit
#' 
#' @inherit brms::predict.brmsfit return
#' 
#' @seealso \code{\link[brms:predict.brmsfit]{predict.brmsfit}}
#' 
#' @importFrom compositions ilrInv
#' @importFrom brms posterior_summary do_call
#' @importFrom abind abind
#' @importFrom stats predict
#' @method predict brmcoda
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   ## fit a model
#'   cilr <- compilr(data = mcompd, sbp = sbp,
#'                   parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                   idvar = "ID", total = 1440)
#'   
#'   m1 <- brmcoda(compilr = cilr,
#'                 formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                   wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   ## predicted responses
#'   pred <- predict(m1)
#'   head(pred)
#'   
#'   ## fit a model with compositional outcome
#'   m2 <- brmcoda(compilr = cilr,
#'                 formula = mvbind(ilr1, ilr2, ilr3, ilr4) ~ Stress + Female + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   ## predicted responses on compositional scale
#'   predcomp <- predict(m2, acomp = TRUE)
#'   head(predcomp)
#' }}
#' @export
predict.brmcoda <- function(object,
                            acomp = FALSE,
                            summary = TRUE,
                            ...) {
  
  if (isFALSE(acomp)) {
    out <- predict(
      object$Model,
      summary = summary,
      ...
    )
  } else {
    if (isFALSE(inherits(object$Model$formula, "mvbrmsformula"))) {
      stop(sprintf(
        "This is a %s model, but a model of class mvbrmsformula is required to
  return results on compsitional scale.",
  class(object$Model$formula)[1]))
      
    } else {
      out <- predict(
        object$Model,
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
#' distribution of a \code{brmsfit} model in the \code{brmcoda} object.
#' Can be performed for the data used to fit the model (posterior
#' predictive checks) or for new data. By definition, these predictions have
#' smaller variance than the posterior predictions performed by the
#' \code{\link{predict.brmcoda}} method. This is because only the
#' uncertainty in the expected value of the posterior predictive distribution is
#' incorporated in the draws computed by \code{fitted} while the
#' residual error is ignored there. However, the estimated means of both methods
#' averaged across draws should be very similar.
#' 
#' @aliases fitted
#' 
#' @inheritParams predict.brmcoda
#' @param ... Further arguments passed to \code{\link{fitted.brmsfit}}
#' that control additional aspects of prediction.
#' 
#' @inherit brms::fitted.brmsfit return
#'   
#' @seealso \code{\link[brms:fitted.brmsfit]{fitted.brmsfit}}
#' 
#' @importFrom compositions ilrInv
#' @importFrom brms posterior_summary do_call
#' @importFrom abind abind
#' @importFrom stats fitted
#' @method fitted brmcoda
#' 
#' @examples
#' \donttest{
#' ## fit a model
#' if(requireNamespace("cmdstanr")){
#'   ## compute composition and ilr coordinates
#'   cilr <- compilr(data = mcompd, sbp = sbp,
#'                   parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                   idvar = "ID", total = 1440)
#'   
#'   ## fit a model
#'   m1 <- brmcoda(compilr = cilr,
#'                 formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                   wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   ## compute expected predictions
#'   epred <- fitted(m1)
#'   head(epred)
#'   
#'   ## fit a model with compositional outcome
#'   m2 <- brmcoda(compilr = cilr,
#'                 formula = mvbind(ilr1, ilr2, ilr3, ilr4) ~ Stress + Female + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   ## expected predictions on compositional scale
#'   epredcomp <- fitted(m2, acomp = TRUE)
#'   head(epredcomp)
#' }}
#' @export
fitted.brmcoda <- function(object,
                           acomp = FALSE,
                           summary = TRUE,
                           ...) {
  
  if (isFALSE(acomp)) {
    out <- fitted(
      object$Model,
      summary = summary,
      ...
    )
  } else {
    if (isFALSE(inherits(object$Model$formula, "mvbrmsformula"))) {
      stop(sprintf(
        "This is a %s model, but a model of class mvbrmsformula is required to
  return results on compsitional scale.",
  class(object$Model$formula)[1]))
      
    } else {
      out <- fitted(
        object$Model,
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

#' Population-Level Estimates
#' 
#' Extract the population-level ('fixed') effects
#' from the \code{brmsfit} object in a \code{brmcoda} object.
#' 
#' @aliases fixef
#' 
#' @param object An object of class \code{brmcoda}.
#' @param ... Further arguments passed to \code{\link{fixef.brmsfit}}.
#' 
#' @inherit brms::fixef.brmsfit return
#' 
#' @seealso \code{\link[brms:fixef.brmsfit]{fixef.brmsfit}}
#' 
#' @importFrom brms fixef
#' @method fixef brmcoda
#' 
#' @examples
#' \donttest{
#' ## fit a model
#' if(requireNamespace("cmdstanr")){
#'   ## fit a model
#'   m <- brmcoda(compilr = compilr(data = mcompd, sbp = sbp,
#'                                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                  idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'   
#'   ## extract population-Level coefficients
#'   fixef(m)
#' }}
#' @export fixef
#' @export
fixef.brmcoda <- function(object, ...) {
  fixef(object$Model, ...)
}

#' Covariance and Correlation Matrix of Population-Level Effects
#'
#' Get a point estimate of the covariance or
#' correlation matrix of population-level parameters
#' of the \code{brmsfit} object in a \code{brmcoda} object.
#'
#' @inheritParams fixef.brmcoda
#' @param ... Further arguments passed to \code{\link{vcov.brmsfit}}.
#' 
#' @inherit brms::vcov.brmsfit return
#' 
#' @seealso \code{\link[brms:vcov.brmsfit]{vcov.brmsfit}}
#' 
#' @importFrom stats vcov
#' @method vcov brmcoda
#' 
#' @examples
#' \donttest{
#' ## fit a model
#' if(requireNamespace("cmdstanr")){
#'   m <- brmcoda(compilr = compilr(data = mcompd, sbp = sbp,
#'                                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                  idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'   
#'   vcov(m)
#' }}
#' @export
vcov.brmcoda <- function(object, ...) {
  vcov(object$Model, ...)
}

#' Group-Level Estimates
#'
#' Extract the group-level ('random') effects of each level
#' of the \code{brmsfit} object in a \code{brmcoda} object.
#' 
#' @aliases ranef
#' 
#' @inheritParams fixef.brmcoda
#' @param ... Further arguments passed to \code{\link{ranef.brmsfit}}.
#' 
#' @inherit brms::ranef.brmsfit return
#'
#' @seealso \code{\link[brms:ranef.brmsfit]{ranef.brmsfit}}
#' 
#' @importFrom brms ranef
#' @method ranef brmcoda
#' 
#' @examples
#' \donttest{
#' ## fit a model
#' if(requireNamespace("cmdstanr")){
#'   m <- brmcoda(compilr = compilr(data = mcompd, sbp = sbp,
#'                                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                  idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'   
#'   ## extract group-level coefficients
#'   ranef(m)
#' }}
#' @export ranef
#' @export
ranef.brmcoda <- function(object, ...) {
  ranef(object$Model, ...)
}

#' Model Coefficients
#'
#' Extract model coefficients, which are the sum of population-level
#' effects and corresponding group-level effects
#' of the \code{brmsfit} object in a \code{brmcoda} object.
#' 
#' @aliases coef
#' 
#' @inheritParams fixef.brmcoda
#' @param ... Further arguments passed to \code{\link{coef.brmsfit}}.
#' 
#' @inherit brms::coef.brmsfit return
#'
#' @seealso \code{\link[brms:coef.brmsfit]{coef.brmsfit}}
#' 
#' @importFrom stats coef
#' @method coef brmcoda
#' 
#' @examples
#' \donttest{
#' ## fit a model
#' if(requireNamespace("cmdstanr")){
#'   m <- brmcoda(compilr = compilr(data = mcompd, sbp = sbp,
#'                                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                  idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'   
#'   ## extract population and group-level coefficients separately
#'   fixef(m)
#'   ranef(m)
#'   
#'   ## extract combined coefficients
#'   coef(m)
#' }}
#' @export
coef.brmcoda <- function(object, ...) {
  coef(object$Model, ...)
}

#' Extract Variance and Correlation Components
#'
#' Calculates the estimated standard deviations,
#' correlations and covariances of the group-level terms
#' of the \code{brmsfit} object in a \code{brmcoda} object.
#' 
#' @aliases VarCorr
#' 
#' @param x An object of class \code{brmcoda}.
#' @param ... Further arguments passed to \code{\link{VarCorr.brmsfit}}.
#' 
#' @inherit brms::VarCorr.brmsfit return
#' 
#' @seealso \code{\link[brms:VarCorr.brmsfit]{VarCorr.brmsfit}}
#' 
#' @importFrom brms VarCorr
#' @method VarCorr brmcoda
#' 
#' @examples
#' \donttest{
#' ## fit a model
#' if(requireNamespace("cmdstanr")){
#'   m <- brmcoda(compilr = compilr(data = mcompd, sbp = sbp,
#'                                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                  idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'   
#'   VarCorr(m)
#' }}
#' @export VarCorr
#' @export
VarCorr.brmcoda <- function(x, ...) {
  VarCorr(x$Model, ...)
}

#' Posterior Draws of Residuals/Predictive Errors
#' 
#' Compute posterior draws of residuals/predictive errors
#' 
#' @inheritParams fixef.brmcoda
#' @param ... Further arguments passed to \code{\link{residuals.brmsfit}}.
#' @inherit brms::residuals.brmsfit return
#'   
#' @seealso \code{\link[brms:residuals.brmsfit]{residuals.brmsfit}}
#' 
#' @importFrom stats residuals
#' @method residuals brmcoda
#' 
#' @examples
#' \donttest{
#' ## fit a model
#' if(requireNamespace("cmdstanr")){
#'   m <- brmcoda(compilr = compilr(data = mcompd, sbp = sbp,
#'                                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                  idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'   
#'   ## extract residuals
#'   res <- residuals(m)
#'   head(res)
#' }}
#' @export
residuals.brmcoda <- function(object, ...) {
  residuals(object$Model, ...)
}
