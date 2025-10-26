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
#' @param scale Specifically for models with compositional responses,
#' either \code{"response"} or \code{"linear"}.
#' If \code{"linear"},
#' results are returned on the log-ratio scale.
#' If \code{"response"}, results are returned on the compositional scale
#' of the response variable.
#' @param parts Only for models with compositional response
#' A optional character string specifying names of compositional parts that make up the response in \code{brmcoda} model.
#' This should correspond to a single set of names of compositional parts specified
#' in the \code{complr} object. Default to the first composition in the \code{complr} object.
#' @param ... Further arguments passed to \code{\link[brms:predict.brmsfit]{predict.brmsfit}}
#' that control additional aspects of prediction.
#'
#' @inherit brms::predict.brmsfit return
#'
#' @seealso \code{\link[brms:predict.brmsfit]{predict.brmsfit}}
#'
#' @importFrom compositions ilrInv clo
#' @importFrom brms posterior_summary do_call
#' @importFrom abind abind
#' @importFrom stats predict
#' 
#' @method predict brmcoda
#'
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   ## fit a model
#'   x <- complr(data = mcompd, sbp = sbp,
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                  idvar = "ID", total = 1440)
#'
#'   m1 <- brmcoda(complr = x,
#'                 formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                    wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'
#'   ## predicted responses
#'   pred <- predict(m1)
#'   head(pred)
#'
#'   ## fit a model with compositional outcome
#'   m2 <- brmcoda(complr = x,
#'                 formula = mvbind(z1_1, z2_1, z3_1, z4_1) ~ 
#'                           bz1_1 + bz2_1 + bz3_1 + bz4_1 + Female + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'
#'   ## predicted responses on ilr scale
#'   predilr <- predict(m2, scale = "linear")
#'   head(predilr)
#'
#'   ## predicted responses on compositional scale
#'   predcomp <- predict(m2, scale = "response")
#'   head(predcomp)
#' }}
#' @export
predict.brmcoda <- function(object,
                            scale = c("linear", "response"),
                            parts = 1,
                            ...) {
  
  if (inherits(object[["model"]][["formula"]], "mvbrmsformula")) {
    # get variables
    brmcoda_vars <- get_variables(object)
    complr_vars  <- get_variables(object$complr)
    
    # check which composition is identical used in model estimation
    idy <- which(
      vapply(complr_vars, function(y) {
        any(sapply(c("Z", "bZ", "wZ"), function(z) {
          identical(sort(y[[z]]), sort(brmcoda_vars$y))
        }))
      }, logical(1))
    )
    
    # check which type of response
    if (identical(brmcoda_vars$y, complr_vars[[idy]]$bZ)) {
      model_resp_type <- "between"
      total <- object$complr$output[[idy]]$total
    }
    if (identical(brmcoda_vars$y, complr_vars[[idy]]$wZ)) {
      model_resp_type <- "within"
      total <- 1
    }
    if (identical(brmcoda_vars$y, complr_vars[[idy]]$Z)) {
      model_resp_type <- "aggregate"
      total <- object$complr$output[[idy]]$total
    }
    if (!exists("model_resp_type")) {
      stop(
        "The response variables in the brmcoda model do not correspond to any of the compositional parts in the complr object."
      )
    }
    
    # predict
    if (identical(scale, "response")) {
      out <- predict(object$model, scale = "response", ...)
      
      # back transform
      out <- lapply(asplit(out, 1), function(x) {
        x <- ilrInv(x, V = object$complr$output[[idy]]$psi)
        as.data.table(clo(x, total = total))
      })
      out <- brms::do_call(abind::abind, c(out, along = 3))
      out <- aperm(out, c(3, 1, 2)) #draw-row-col
    } else {
      out <- predict(object$model, ...)
    }
    
  } else {
    out <- predict(object$model, ...)
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
#' @param summary Should summary statistics be returned instead of the raw values? Default is \code{TRUE}.
#' @param ... Further arguments passed to \code{\link[brms:fitted.brmsfit]{fitted.brmsfit}}
#' that control additional aspects of prediction.
#'
#' @inherit brms::fitted.brmsfit return
#'
#' @seealso \code{\link[brms:fitted.brmsfit]{fitted.brmsfit}}
#'
#' @importFrom compositions ilrInv clo
#' @importFrom brms posterior_summary do_call
#' @importFrom abind abind
#' @importFrom stats fitted
#' 
#' @method fitted brmcoda
#'
#' @examples
#' \donttest{
#' ## fit a model
#' if(requireNamespace("cmdstanr")){
#'   ## compute composition and ilr coordinates
#'   x <- complr(data = mcompd, sbp = sbp,
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                  idvar = "ID", total = 1440)
#'
#'   ## fit a model
#'   m1 <- brmcoda(complr = x,
#'                 formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                    wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'
#'   ## compute expected predictions
#'   epred <- fitted(m1)
#'   head(epred)
#'
#'   ## fit a model with compositional outcome
#'   m2 <- brmcoda(complr = x,
#'                 formula = mvbind(z1_1, z2_1, z3_1, z4_1) ~ Stress + Female + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'
#'   ## expected predictions on compositional scale
#'   epredcomp <- fitted(m2, scale = "response")
#'   head(epredcomp)
#' }}
#' @export
fitted.brmcoda <- function(object,
                           scale = c("linear", "response"),
                           parts = 1,
                           summary = TRUE,
                           ...) {
  
  if (inherits(object[["model"]][["formula"]], "mvbrmsformula")) {
    # get brmcoda variables
    brmcoda_vars <- get_variables(object)
    complr_vars  <- get_variables(object$complr)
    
    # check which composition is identical used in model estimation
    idy <- which(
      vapply(complr_vars, function(y) {
        any(sapply(c("Z", "bZ", "wZ"), function(z) {
          identical(sort(y[[z]]), sort(brmcoda_vars$y))
        }))
      }, logical(1))
    )

    # check which type of response
    if (identical(brmcoda_vars$y, complr_vars[[idy]]$bZ)) {
      model_resp_type <- "between"
      total <- object$complr$output[[idy]]$total
      Yn <- complr_vars[[idy]]$bX
    }
    if (identical(brmcoda_vars$y, complr_vars[[idy]]$wZ)) {
      model_resp_type <- "within"
      total <- 1
      Yn <- complr_vars[[idy]]$wX
    }
    if (identical(brmcoda_vars$y, complr_vars[[idy]]$Z)) {
      model_resp_type <- "aggregate"
      total <- object$complr$output[[idy]]$total
      Yn <- complr_vars[[idy]]$X
    }
    if (!exists("model_resp_type")) {
      stop(
        "The response variables in the brmcoda model do not correspond to any of the compositional parts in the complr object."
      )
    }
    
    # predict
    if (identical(scale, "response")) {
      out <- fitted(object$model, scale = "response", summary = FALSE, ...)
      
      # back transform
      out <- lapply(asplit(out, 1), function(x) {
        x <- ilrInv(x, V = object$complr$output[[idy]]$psi)
        clo(x, total = total)
      })
      out <- brms::do_call(abind::abind, c(out, along = 3)) # 1 = draw, 2 = ob, 3 = part
      out <- aperm(out, c(3, 1, 2)) #draw-row-col
      dimnames(out)[[3]] <- Yn

      # summary posteriors
      if (summary) {
        out <- posterior_summary(out)
        dimnames(out)[[3]] <- Yn
      }
    } else {
      out <- fitted(object$model, summary = summary, ...)
    }
  } else {
    out <- fitted(object$model, summary = summary, ...)
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
#' @param ... Further arguments passed to \code{\link[brms:fixef.brmsfit]{fixef.brmsfit}}.
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
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID", total = 1440),
#'                formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                                   chain = 1, iter = 500,
#'                                   backend = "cmdstanr")
#'
#'   ## extract population-Level coefficients
#'   fixef(m)
#' }}
#' @export fixef
#' @export
fixef.brmcoda <- function(object, ...) {
  fixef(object$model, ...)
}

#' Covariance and Correlation Matrix of Population-Level Effects
#'
#' Get a point estimate of the covariance or
#' correlation matrix of population-level parameters
#' of the \code{brmsfit} object in a \code{brmcoda} object.
#'
#' @inheritParams fixef.brmcoda
#' @param ... Further arguments passed to \code{\link[brms:vcov.brmsfit]{vcov.brmsfit}}.
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
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID", total = 1440),
#'                formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                                   chain = 1, iter = 500,
#'                                   backend = "cmdstanr")
#'
#'   vcov(m)
#' }}
#' @export
vcov.brmcoda <- function(object, ...) {
  vcov(object$model, ...)
}

#' Group-Level Estimates
#'
#' Extract the group-level ('random') effects of each level
#' of the \code{brmsfit} object in a \code{brmcoda} object.
#'
#' @aliases ranef
#'
#' @inheritParams fixef.brmcoda
#' @param ... Further arguments passed to \code{\link[brms:ranef.brmsfit]{ranef.brmsfit}}.
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
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID", total = 1440),
#'                formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                                   chain = 1, iter = 500,
#'                                   backend = "cmdstanr")
#'
#'   ## extract group-level coefficients
#'   ranef(m)
#' }}
#' @export ranef
#' @export
ranef.brmcoda <- function(object, ...) {
  ranef(object$model, ...)
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
#' @param ... Further arguments passed to \code{\link[brms:coef.brmsfit]{coef.brmsfit}}.
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
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID", total = 1440),
#'                formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                                   chain = 1, iter = 500,
#'                                   backend = "cmdstanr")
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
  coef(object$model, ...)
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
#' @param ... Further arguments passed to \code{\link[brms:VarCorr.brmsfit]{VarCorr.brmsfit}}.
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
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID", total = 1440),
#'                formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                                  chain = 1, iter = 500,
#'                                  backend = "cmdstanr")
#'
#'   VarCorr(m)
#' }}
#' @export VarCorr
#' @export
VarCorr.brmcoda <- function(x, ...) {
  VarCorr(x$model, ...)
}

#' Posterior Draws of Residuals/Predictive Errors
#'
#' Compute posterior draws of residuals/predictive errors
#'
#' @inheritParams fixef.brmcoda
#' @param ... Further arguments passed to \code{\link[brms:residuals.brmsfit]{residuals.brmsfit}}.
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
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID", total = 1440),
#'                formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                                   chain = 1, iter = 500,
#'                                   backend = "cmdstanr")
#'
#'   ## extract residuals
#'   res <- residuals(m)
#'   head(res)
#' }}
#' @export
residuals.brmcoda <- function(object, ...) {
  residuals(object$model, ...)
}

#' Efficient approximate leave-one-out cross-validation (LOO)
#'
#' Perform approximate leave-one-out cross-validation based
#' on the posterior likelihood using the \pkg{loo} package.
#' For more details see \code{\link[loo:loo]{loo}}.
#'
#' @aliases loo
#'
#' @param x A \code{brmcoda} object.
#' @inheritParams brms::loo.brmsfit
#' @inherit brms::loo.brmsfit return
#'
#' @seealso \code{\link[brms:loo.brmsfit]{loo.brmsfit}}
#'
#' @importFrom loo loo is.loo
#' @method loo brmcoda
#' @export
loo.brmcoda <- function(x, ...) {
  loo(x$model, ...)
}
