#' Substitution Plot
#'
#' Make a plot of \code{\link{substitution}} model results.
#'
#' @param x A \code{\link{substitution}} class object.
#' @param to An optional character value or vector specifying the names of the compositional parts
#' that were reallocated to in the model.
#' @param ref A character value of ((\code{"grandmean"} or \code{"clustermean"} or \code{"users"}),
#' @param level An optional character value of (\code{"between"}, \code{"within"}), or \code{"aggregate"}).
#' @param ... Further components to the plot, followed by a plus sign (+).
#'
#' @return A ggplot graph object showing the estimated difference in outcome when
#' each pair of compositional variables are substituted for a specific time.
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_line geom_pointrange geom_ribbon
#' @importFrom ggplot2 facet_wrap vars label_both label_bquote position_dodge2
#' @importFrom ggplot2 theme theme_bw scale_fill_manual scale_colour_manual
#' @importFrom data.table copy
#' @importFrom bayesplot color_scheme_get
#'
#' @method plot substitution
#' @export
plot.substitution <- function(x, to, ref, level, ...) {
  
  if (missing(ref)) {
    ref <- x$ref
  }
  if (missing(level)) {
    level <- x$level
  }
  if (missing(to)) {
    to <- x$parts
  }
  
  if (length(x$delta) > 1) {
    # extract data
    tmp <- summary(
      object = x,
      delta = sort(c(-abs(x$delta), abs(x$delta))),
      to = to,
      ref = ref,
      level = level,
      digits = "asis"
    )
    col_pal <- rev(unlist(color_scheme_get("brewer-PuBuGn")))
    names(col_pal) <- sort(x$parts)
    
    ggplot(tmp, aes(x = Delta,  y = Estimate, colour = From, fill = From)) +
      geom_hline(yintercept = 0,
                 linewidth = 0.25,
                 linetype = 2) +
      geom_vline(xintercept = 0,
                 linewidth = 0.25,
                 linetype = 2) +
      geom_line(linewidth = 0.75) +
      geom_ribbon(aes(ymin = CI_low, ymax = CI_high),
                  alpha = 3 / 10,
                  linewidth = 0) +
      scale_colour_manual(values = if (length(x$parts) < 7) col_pal else NULL) +
      scale_fill_manual(values = if (length(x$parts) < 7) col_pal else NULL) +
      facet_wrap(
        vars(From, To),
        labeller = ggplot2::label_bquote(cols = .(as.character(From)) %<-% minutes %->% .(as.character(To))),
        nrow = length(x$parts)
      ) +
      theme_bw() +
      theme(
        legend.position  = "none",
        strip.background = element_rect(fill = "transparent", color = "black", linewidth = 0.5),
        axis.ticks       = element_blank()
      )
    
  } else {
    # extract data
    tmp <- summary(
      object = x,
      delta = x$delta,
      to = to,
      ref = ref,
      level = level,
      digits = "asis"
    )
    col_pal <- unlist(color_scheme_get("blue"))
    names(col_pal) <- sort(x$parts)
    
    ggplot(tmp, aes(x = Delta, y = Estimate, colour = From)) +
      geom_hline(yintercept = 0,
                 linewidth = 0.25,
                 linetype = 2) +
      # geom_vline(xintercept = 0,
      #            linewidth = 0.2,
      #            linetype = 2) +
      geom_pointrange(aes(ymin = CI_low, ymax = CI_high), 
                      position = position_dodge2(width = 0.25),
                      size = 0.5, linewidth = 0.75
                      ) +
      scale_colour_manual(values = if (length(x$parts) < 7) col_pal else NULL) +
      facet_wrap(~ To,
                 labeller = label_both,
                 nrow = length(x$parts)) +
      # facet_wrap(
      #   vars(From, To),
      #   labeller = ggplot2::label_bquote(cols = .(as.character(From)) %->% .(as.character(To))),
      #   nrow = length(x$parts)
      # ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        strip.background = element_rect(fill = "transparent", color = "black", linewidth = 0.5),
        axis.text.x      = element_blank(),
        axis.title.x     = element_blank(),
        axis.ticks       = element_blank()
      )
  }
}

#' Trace and Density Plots for MCMC Draws plot
#'
#' Make a plot of \code{brmcoda} model results.
#'
#' @param x A \code{\link{brmcoda}} class object.
#' @param ... Further arguments passed to \code{\link[brms:plot.brmsfit]{plot.brmsfit}}.
#'
#' @inherit brms::plot.brmsfit return
#'
#' @seealso \code{\link[brms:plot.brmsfit]{plot.brmsfit}}
#'
#' @method plot brmcoda
#' @export
#' @examples
#' \dontrun{
#' cilr <- complr(data = mcompd, sbp = sbp,
#'         parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#'
#' # model with compositional predictor at between and within-person levels
#' fit <- brmcoda(complr = cilr,
#'                 formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                    wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                chain = 1, iter = 500)
#' plot(fit)
#' }
plot.brmcoda <- function(x, ...) {
  plot(x$model, ...)
}

#' Create a matrix of output plots from a \code{\link{brmcoda}}'s \code{\link[brms:brmsfit]{brmsfit}} object
#'
#' A \code{\link[graphics:pairs]{pairs}}
#' method that is customized for MCMC output.
#'
#' @param x A \code{brmcoda} class object.
#' @param ... Further arguments passed to \code{\link[brms:pairs.brmsfit]{pairs.brmsfit}}.
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
#' cilr <- complr(data = mcompd, sbp = sbp,
#'         parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#'
#' # model with compositional predictor at between and within-person levels
#' fit <- brmcoda(complr = cilr,
#'                 formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                    wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                chain = 1, iter = 500)
#' pairs(fit)
#' }
pairs.brmcoda <- function(x, ...) {
  pairs(x$model, ...)
}

#' MCMC Plots Implemented in \pkg{bayesplot}
#'
#' Call MCMC plotting functions
#' implemented in the \pkg{bayesplot} package.
#'
#' @param object A \code{brmcoda} class object.
#' @param ... Further arguments passed to \code{\link[brms:mcmc_plot.brmsfit]{mcmc_plot.brmsfit}}.
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
#' cilr <- complr(data = mcompd, sbp = sbp,
#'         parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#'
#' # model with compositional predictor at between and within-person levels
#' fit <- brmcoda(complr = cilr,
#'                 formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                    wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                chain = 1, iter = 500)
#' mcmc_plot(fit)
#' }
mcmc_plot.brmcoda <- function(object, ...) {
  mcmc_plot(object$model, ...)
}
