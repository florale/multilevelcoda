#' Checks if argument is a \code{brmcoda} object
#'
#' @param x An object of class \code{brmcoda}.
#'
#' @export
is.brmcoda <- function(x) {
  inherits(x, "brmcoda")
}


#' Extract Number of Observations from \code{brmcoda} object
#'
#' @param object A \code{brmcoda} object.
#' @param ... Further arguments to be passed to methods.
#'
#' @importFrom stats nobs
#' @method nobs brmcoda
#' @export
nobs.brmcoda <- function(object, ...) {
  nobs(object$model, ...)
}

#' Extracting the Model Frame from a Formula or Fit from \code{brmcoda} object
#'
#' @param formula A \code{brmcoda} object.
#' @param ... Further arguments to be passed to methods.
#'
#' @importFrom stats model.frame
#' @method model.frame brmcoda
#' @export
model.frame.brmcoda <- function(formula, ...) {
  model.frame(formula$model, ...)
}

#' Index \code{brmcoda} objects
#'
#' @aliases variables nvariables niterations nchains ndraws
#'
#' @param x An object of class \code{brmcoda}.
#' @param ... Arguments passed to individual methods.
#'
#' @name draws-index-brmcoda
#'
NULL

#' @rdname draws-index-brmcoda
#' @importFrom brms variables
#' @method variables brmcoda
#'
#' @seealso \code{\link[brms:variables]{variables.brmsfit}}
#'
#' @export
variables.brmcoda <- function(x, ...) {
  variables(x$model, ...)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms nvariables
#' @method nvariables brmcoda
#'
#' @seealso \code{\link[brms:nvariables.brmsfit]{nvariables.brmsfit}}
#'
#' @export
nvariables.brmcoda <- function(x, ...) {
  nvariables(x$model, ...)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms niterations
#' @method niterations brmcoda
#'
#' @seealso \code{\link[brms:niterations.brmsfit]{niterations.brmsfit}}
#'
#' @export
niterations.brmcoda <- function(x) {
  niterations(x$model)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms nchains
#' @method nchains brmcoda
#'
#' @seealso \code{\link[brms:nchains.brmsfit]{nchains.brmsfit}}
#'
#' @export
nchains.brmcoda <- function(x) {
  nchains(x$model)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms ndraws
#' @method ndraws brmcoda
#'
#' @seealso \code{\link[brms:ndraws.brmsfit]{ndraws.brmsfit}}
#'
#' @export
ndraws.brmcoda <- function(x) {
  ndraws(x$model)
}

# nwarmup.brmcoda <- function(x) {
#   if (inherits(x$model, "brmcoda")) {
#     nwarmup(x$model)
#   }
# }
#
# nthin.brmcoda <- function(x) {
#   if (inherits(x$model, "brmcoda")) {
#     nthin(x$model)
#   }
# }

#' Extract Diagnostic Quantities from \code{brmsfit} Models in \code{brmcoda}
#'
#' @name diagnostic-quantities-brmcoda
#' @aliases log_posterior nuts_params rhat neff_ratio
#'
#' @param x,object A \code{brmcoda} object or another \R object for which
#' the methods are defined.
#' @param ... Arguments passed to individual methods (if applicable).
#'
#' @return The exact form of the output depends on the method.
#'
NULL

#' @rdname diagnostic-quantities-brmcoda
#' @importFrom brms log_posterior
#' @method log_posterior brmcoda
#'
#' @seealso \code{\link[brms:log_posterior.brmsfit]{log_posterior.brmsfit}}
#'
#' @export
log_posterior.brmcoda <- function(object, ...) {
  log_posterior(object$model, ...)
}

#' @rdname diagnostic-quantities-brmcoda
#' @importFrom brms nuts_params
#' @method nuts_params brmcoda
#'
#' @seealso \code{\link[brms:nuts_params.brmsfit]{nuts_params.brmsfit}}
#'
#' @export
nuts_params.brmcoda <- function(object, ...) {
  nuts_params(object$model, ...)
}

#' @rdname diagnostic-quantities-brmcoda
#' @importFrom brms rhat
#' @method rhat brmcoda
#'
#' @seealso \code{\link[brms:rhat.brmsfit]{rhat.brmsfit}}
#'
#' @export
rhat.brmcoda <- function(x, ...) {
  rhat(x$model, ...)
}

#' @rdname diagnostic-quantities-brmcoda
#' @importFrom brms neff_ratio
#' @method neff_ratio brmcoda
#'
#' @seealso \code{\link[brms:neff_ratio.brmsfit]{neff_ratio.brmsfit}}
#'
#' @export
neff_ratio.brmcoda <- function(object, ...) {
  neff_ratio(object$model, ...)
}

#' Bayes Factors from Marginal Likelihoods
#'
#' Compute Bayes factors from marginal likelihoods
#'
#' @param x1 A \code{brmcoda} object.
#' @param x2 Another \code{brmcoda} object based on the same responses.
#' @param ... Other arguments passed to \code{\link[brms:bayes_factor.brmsfit]{bayes_factor.brmsfit}}.
#'
#' @importFrom brms bayes_factor
#' @method bayes_factor brmcoda
#'
#' @seealso \code{\link[brms:bayes_factor.brmsfit]{bayes_factor.brmsfit}}
#'
#' @export
bayes_factor.brmcoda <- function(x1, x2, ...) {
  out <- invisible(bayes_factor(x1 = x1$model, x2 = x2$model, ...))
  
  m1 <- deparse(substitute(x1))
  m2 <- deparse(substitute(x2))
  
  attr(out, "model_names") <- c(m1, m2)
  out
}

#' Extract Priors of a \code{brmsfit} from a \code{brmcoda} object
#'
#' Compute Bayes factors from marginal likelihoods
#'
#' @param object An object of class \code{brmcoda}.
#' @inheritParams brms::prior_summary.brmsfit
#'
#' @importFrom brms prior_summary
#' @method prior_summary brmcoda
#'
#' @seealso \code{\link[brms:prior_summary.brmsfit]{prior_summary.brmsfit}}
#'
#' @export
prior_summary.brmcoda <- function(object, ...) {
  prior_summary(object$model, ...)
}

#' Posterior Predictive Checks for \code{brmcoda} Objects
#'
#' Perform posterior predictive checks with the help of the \pkg{bayesplot} package.
#'
#' @aliases pp_check
#'
#' @param object An object of class \code{brmcoda}.
#' @inheritParams brms::pp_check.brmsfit
#'
#' @importFrom bayesplot pp_check
#' @method pp_check brmcoda
#'
#' @seealso \code{\link[brms:pp_check.brmsfit]{pp_check.brmsfit}}
#'
#' @export
pp_check.brmcoda <- function(object, ...) {
  pp_check(object$model, ...)
}

#' #' Posteriors Sampling Diagnostic
#' #'
#' #' Extract diagnostic metrics (Effective Sample Size (`ESS`), `Rhat` and Monte
#' #' Carlo Standard Error `MCSE`).
#' #'
#' #' @param posteriors An object of class \code{brmcoda}.
#' #' @inheritParams bayestestR::diagnostic_posterior
#' #' @param ... Other arguments passed to \code{\link{diagnostic_posterior}}.
#' #'
#' #' @importFrom bayestestR diagnostic_posterior
#' #' @method diagnostic_posterior brmcoda
#' #'
#' #' @seealso \code{\link[bayestestR:diagnostic_posterior]{diagnostic_posterior}}
#' #'
#' #' @export
#' diagnostic_posterior.brmcoda <- function(posteriors, diagnostic = c("ESS", "Rhat"), ...) {
#'   diagnostic_posterior(posterior$model, diagnostic = diagnostic, ...)
#' }
