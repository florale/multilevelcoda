#' @title Substitution plot
#'
#' @description
#' This function is useful for visualising the
#' estimated differences in outcomes when compositional variables
#' are substituted for a specific period of time.
#'
#' @param data A dataset to use for plot.
#' It must be a component of a list resulted from one of the following functions:
#' \code{\link{wsub}}, \code{\link{bsub}}, \code{\link{wsubmargins}}, \code{\link{bsubmargins}}.
#' @param x A character string specifying name of the compostional predictor variable.
#' @param y A character string specifying the name of the outcome variable.
#' @param ... Further arguments passed to \code{\link{ggplot}}.
#'
#' @return A ggplot graph object showing the estimated difference in outcome when
#' each pair of compositional variables are substituted for a specific time.
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_line geom_ribbon facet_grid xlab ylab
#' @importFrom ggsci scale_color_simpsons
#' @importFrom ggsci scale_fill_simpsons
#' @importFrom data.table copy
#' @export
plotsub <- function(data, x, y, ...) {

  if (isFALSE(inherits(data, c("data.table", "data.frame")))) {
    stop("data must be a data table or data frame,",
         "and is an element of a wsub, bsub, wsubmargins, bsubmargins object.")
  }
  
  tmp <- copy(data)

  plot <- ggplot(tmp, aes(x = Delta, y = Mean)) +
    geom_hline(yintercept = 0, size = 0.2, linetype = 2) +
    geom_vline(xintercept = 0, size = 0.2, linetype = 2) +
    geom_line(aes(colour = To), size = 1) +
    geom_ribbon(aes(ymin = CI_low,
                    ymax = CI_high, fill = To),
                alpha = 1 / 10, size = 1 / 10) +
    facet_grid(~ To) +
    xlab(paste("Change in", eval(x), sep = " ")) +
    ylab(paste("Change in", eval(y), sep = " "))
  plot
}
