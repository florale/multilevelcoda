#' Substitution Plot
#'
#' Make a plot of \code{\link{substitution}} model results.
#'
#' @param x A \code{\link{substitution}} class object.
#' @param to A character value or vector specifying the names of the compositional parts
#' that were reallocated to in the model.
#' @param ref A character value of ((\code{"grandmean"} or \code{"clustermean"} or \code{"users"}),
#' @param level A character value of (\code{"between"} or \code{"within"}).
#' @param ... Further arguments passed to \code{\link{ggplot}}.
#'
#' @return A ggplot graph object showing the estimated difference in outcome when
#' each pair of compositional variables are substituted for a specific time.
#' 
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_line geom_pointrange geom_ribbon facet_grid xlab ylab
#' @importFrom data.table copy
#' 
#' @method plot substitution
#' @export
plot.substitution <- function(x, to,
                              ref, level, ...) {
  
  if (isFALSE(any(c("grandmean", "clustermean", "users") %in% ref)) ||
      isTRUE(length(ref) > 1)) {
    stop("'ref' should be either grandmean or clustermean or users.")
  }
  ref <- as.character(ref)
  
  if (isFALSE(any(c("between", "within") %in% level)) ||
      isTRUE(length(level) > 1)) {
    stop("'level' should be either between or within.")
  }
  level <- as.character(level)
  
  # extract delta
  delta.pos <- x$delta
  delta.neg <- -1*abs(x$delta)
  delta <- c(delta.pos, delta.neg)
  
  # extract data
  tmp <- summary(object = x,
                 delta = delta,
                 to = to,
                 ref = ref,
                 level = level,
                 digits = "asis"
  )
  
  # plot
  if (isTRUE(is.sequential(delta.pos))) {
    plotsub <- ggplot(tmp, 
                      aes(x = Delta, y = Mean)) +
      geom_line(aes(colour = From), linewidth = 1) +
      geom_ribbon(
        aes(ymin = CI_low,
            ymax = CI_high, fill = From),
        alpha = 1 / 10,
        linewidth = 1 / 10) +
      geom_hline(yintercept = 0,
                 linewidth = 0.2,
                 linetype = 2) +
      geom_vline(xintercept = 0,
                 linewidth = 0.2,
                 linetype = 2) +
      facet_grid( ~ From)
    
  } else {
    plotsub <- ggplot(tmp,
                      aes(x = Delta, y = Mean)) +
      geom_line(aes(colour = From)) +
      geom_pointrange(aes(ymin = CI_low, ymax = CI_high, colour = From)) +
      geom_hline(yintercept = 0,
                 linewidth = 0.2,
                 linetype = 2) +
      geom_vline(xintercept = 0,
                 linewidth = 0.2,
                 linetype = 2) +
      facet_grid( ~ From)
    
  }
  plotsub
}

#' Trace and Density Plots for MCMC Draws plot
#'
#' Make a plot of \code{brmcoda} model results.
#'
#' @param x A \code{\link{brmcoda}} class object.
#' @param ... Further arguments passed to \code{\link{plot.brmsfit}}.
#'
#' @inherit brms::plot.brmsfit return
#'
#' @seealso \code{\link[brms:plot.brmsfit]{plot.brmsfit}}
#' 
#' @method plot brmcoda
#' @export
#' @examples
#' \dontrun{
#' cilr <- compilr(data = mcompd, sbp = sbp,
#'         parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#'
#' # model with compositional predictor at between and within-person levels
#' fit <- brmcoda(compilr = cilr,
#'                formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                                  wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                chain = 1, iter = 500)
#' plot(fit)
#' }
plot.brmcoda <- function(x, ...) {
  plot(x$Model, ...)
}

#' Create a matrix of output plots from a \code{\link{brmcoda}}'s \code{\link{brmsfit}} object
#'
#' A \code{\link[graphics:pairs]{pairs}}
#' method that is customized for MCMC output.
#'
#' @param x A \code{brmcoda} class object.
#' @param ... Further arguments passed to \code{\link{pairs.brmsfit}}.
#'
#' @inherit brms::pairs.brmsfit return
#' 
#' @seealso \code{\link[brms:pairs.brmsfit]{pairs.brmsfit}}
#' 
#' @importFrom graphics pairs
#' @method pairs brmcoda
#' @export
#' @examples
#' \dontrun{
#' cilr <- compilr(data = mcompd, sbp = sbp,
#'         parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#'
#' # model with compositional predictor at between and within-person levels
#' fit <- brmcoda(compilr = cilr,
#'                formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                                  wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                chain = 1, iter = 500)
#' pairs(fit)
#' }
pairs.brmcoda <- function(x, ...) {
  pairs(x$Model, ...)
}

#' MCMC Plots Implemented in \pkg{bayesplot}
#'
#' Call MCMC plotting functions
#' implemented in the \pkg{bayesplot} package.
#'
#' @param object A \code{brmcoda} class object.
#' @param ... Further arguments passed to \code{\link{mcmc_plot.brmsfit}}.
#' 
#' @inherit brms::mcmc_plot.brmsfit return
#' 
#' @seealso \code{\link[brms:mcmc_plot.brmsfit]{mcmc_plot.brmsfit}}
#' 
#' @importFrom brms mcmc_plot
#' @method mcmc_plot brmcoda
#' 
#' @export
#' @examples
#' \dontrun{
#' cilr <- compilr(data = mcompd, sbp = sbp,
#'         parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#'
#' # model with compositional predictor at between and within-person levels
#' fit <- brmcoda(compilr = cilr,
#'                formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                                  wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                chain = 1, iter = 500)
#' mcmc_plot(fit)
#' }
mcmc_plot.brmcoda <- function(object, ...) {
  mcmc_plot(object$Model, ...)
}