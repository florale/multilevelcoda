#' Several methods for \code{brmcoda} objects.
#'
#' Checks if argument is a \code{brmcoda} object
#'
#' @param x An object of class \code{brmcoda}.
#'
#' @export
is.brmcoda <- function(x) {
  inherits(x, "brmcoda")
}

#' Extract Number of Observations from \pkg{brmcoda} object
#'
#' @param object A \code{brmcoda} object.
#' @param ... Further arguments to be passed to methods.
#'
#' @importFrom stats nobs
#' @method nobs brmcoda
#' @export
nobs.brmcoda <- function(object, ...) {
  nobs(object$Model, ...)
}

#' Extract Posterior Draws
#'
#' Extract posterior draws in conventional formats
#' as data.frames, matrices, or arrays.
#'
#' @inheritParams brms::as.data.frame.brmsfit
#'
#' @return A data.frame, matrix, or array containing the posterior draws.
#' 
#' @export
as.data.frame.brmcoda <- function(x, row.names = NULL, optional = TRUE, ...) {
  as.data.frame(x$Model, ...)
}

#' @rdname as.data.frame.brmcoda
#' @export
as.matrix.brmcoda <- function(x, ...) {
  as.matrix(x$Model, ...)
}

#' @rdname as.data.frame.brmcoda
#' @export
as.array.brmcoda <- function(x, ...) {
  as.array(x$Model, ...)
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
  variables(x$Model, ...)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms nvariables
#' @method nvariables brmcoda
#' 
#' @seealso \code{\link[brms:nvariables.brmsfit]{nvariables.brmsfit}}
#' 
#' @export
nvariables.brmcoda <- function(x, ...) {
  nvariables(x$Model, ...)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms niterations
#' @method niterations brmcoda
#' 
#' @seealso \code{\link[brms:niterations.brmsfit]{niterations.brmsfit}}
#' 
#' @export
niterations.brmcoda <- function(x) {
  niterations(x$Model)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms nchains
#' @method nchains brmcoda
#' 
#' @seealso \code{\link[brms:nchains.brmsfit]{nchains.brmsfit}}
#' 
#' @export 
nchains.brmcoda <- function(x) {
  nchains(x$Model)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms ndraws
#' @method ndraws brmcoda
#' 
#' @seealso \code{\link[brms:ndraws.brmsfit]{ndraws.brmsfit}}
#' 
#' @export
ndraws.brmcoda <- function(x) {
  ndraws(x$Model)
}

# nwarmup.brmcoda <- function(x) {
#   if (inherits(x$Model, "brmcoda")) {
#     nwarmup(x$Model)
#   }
# }
#
# nthin.brmcoda <- function(x) {
#   if (inherits(x$Model, "brmcoda")) {
#     nthin(x$Model)
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
  log_posterior(object$Model)
}

#' @rdname diagnostic-quantities-brmcoda
#' @importFrom brms nuts_params
#' @method nuts_params brmcoda
#' 
#' @seealso \code{\link[brms:nuts_params.brmsfit]{nuts_params.brmsfit}}
#' 
#' @export
nuts_params.brmcoda <- function(object, ...) {
  nuts_params(object$Model)
}

#' @rdname diagnostic-quantities-brmcoda
#' @importFrom brms rhat
#' @method rhat brmcoda
#' 
#' @seealso \code{\link[brms:rhat.brmsfit]{rhat.brmsfit}}
#' 
#' @export
rhat.brmcoda <- function(x, ...) {
  rhat(x$Model)
}

#' @rdname diagnostic-quantities-brmcoda
#' @importFrom brms neff_ratio
#' @method neff_ratio brmcoda
#' 
#' @seealso \code{\link[brms:neff_ratio.brmsfit]{neff_ratio.brmsfit}}
#' 
#' @export
neff_ratio.brmcoda <- function(object, ...) {
  neff_ratio(object$Model)
}
