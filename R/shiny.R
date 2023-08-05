#' Interface to \pkg{shinystan}
#' 
#' Provide an interface to \pkg{shinystan} for models fitted with \pkg{brms}
#' 
#' @aliases launch_shinystan
#' 
#' @param object A fitted model object of class \code{brmcoda}.
#' @param ... Optional arguments to pass to 
#' \code{\link{launch_shinystan.brmsfit}} or \code{\link[shiny:runApp]{runApp}}.
#'
#' @seealso \code{\link[shinystan:launch_shinystan]{launch_shinystan}}
#' 
#' @return An S4 shinystan object
#' 
#' @method launch_shinystan brmcoda
#' @importFrom brms launch_shinystan
#' @export
launch_shinystan.brmcoda <- function(object, ...) {
  launch_shinystan(object$Model, ...)
}
