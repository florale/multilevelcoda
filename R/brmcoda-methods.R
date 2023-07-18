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
#' Extract Number of Observations from \pkg{brmcoda} Models
#'
#' @param object A \code{brmcoda} object.
#'
#' @method nobs brmcoda
#' @export
nobs.brmcoda <- function(x) {
  nobs(x$Model)
}

#' Extract Posterior Draws
#'
#' Extract posterior draws in conventional formats
#' as data.frames, matrices, or arrays.
#'
#' @param x An object of class \code{brmcoda}.
#'
#' @export
as.data.frame.brmcoda <- function(x) {
  as.data.frame(x$Model)
}

#' @rdname as.data.frame.brmcoda
#' @export
as.matrix.brmcoda <- function(x) {
  as.matrix(x$Model)
}

#' @rdname as.data.frame.brmcoda
#' @export
as.array.brmcoda <- function(x) {
  as.array(x$Model)
}

#' Index \code{brmcoda} objects
#'
#' @param x An object of class \code{brmcoda}.
#'
#' @name draws-index-brmcoda
NULL

#' @rdname draws-index-brmcoda
#' @importFrom brms variables
#' @method variables brmcoda
#' @export variables
#' @export
variables.brmcoda <- function(x, ...) {
  variables(x$Model)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms nvariables
#' @method nvariables brmcoda
#' @export nvariables
#' @export
nvariables.brmcoda <- function(x, ...) {
  nvariables(x$Model)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms niterations
#' @method niterations brmcoda
#' @export niterations
#' @export
niterations.brmcoda <- function(x) {
  niterations(x$Model)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms nchains
#' @method nchains brmcoda
#' @export nchains
#' @export 
nchains.brmcoda <- function(x) {
  nchains(x$Model)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms ndraws
#' @method ndraws brmcoda
#' @export ndraws
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

#' Extract Diagnostic Quantities from \pkg{brmcoda} Models
#'
#' @name diagnostic-quantities
#' @aliases log_posterior nuts_params rhat neff_ratio
#'
#' @param x A \code{brmcoda} object.
#' @param ... Arguments passed to individual methods.
#'
#' @return The exact form of the output depends on the method.
#'
NULL

#' @rdname diagnostic-quantities
#' @importFrom brms log_posterior
#' @export log_posterior
#' @export
log_posterior.brmcoda <- function(x) {
  log_posterior(x$Model)
}

#' @rdname diagnostic-quantities
#' @importFrom brms nuts_params
#' @export
nuts_params.brmcoda <- function(x) {
  nuts_params(x$Model)
}

#' @rdname diagnostic-quantities
#' @importFrom brms rhat
#' @export rhat
#' @export
rhat.brmcoda <- function(x) {
  rhat(x$Model)
}

#' @rdname diagnostic-quantities
#' @importFrom brms neff_ratio
#' @export neff_ratio
#' @export
neff_ratio.brmcoda <- function(x) {
  neff_ratio(x$Model)
}

#' Interface to \pkg{shinystan}
#' 
#' Provide an interface to \pkg{shinystan} for models fitted with \pkg{brms}
#' 
#' @aliases launch_shinystan
#' 
#' @param object A fitted model object typically of class \code{brmcoda}.
#' @param ... Optional arguments to pass to 
#' \code{\link{launch_shinystan.brmsfit}} or \code{\link[shiny:runApp]{runApp}}
#'
#' @seealso \code{\link[shinystan:launch_shinystan]{launch_shinystan}}
#' 
#' @return An S4 shinystan object
#' 
#' @method launch_shinystan brmcoda
#' @importFrom shinystan launch_shinystan
#' @export launch_shinystan
#' @export
launch_shinystan.brmcoda <- function(object, ...) {
  launch_shinystan(object$Model, ...)
}
