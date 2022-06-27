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
#' @param font A character string specifying the name of user's preferred font.
#' Default is \code{serif} which will choose a serif font available on the system.
#' @param ... Further arguments passed to \code{\link{ggplot}}.
#'
#' @return A ggplot graph object showing the estimated difference in outcome when
#' each pair of compositional variables are substituted for a specific time.
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_line geom_ribbon facet_grid xlab ylab
#' @importFrom cowplot theme_cowplot
#' @importFrom ggsci scale_color_simpsons
#' @importFrom data.table copy
#' @export
plotsub <- function(data, x, y, font = "serif", ...) {

  if (isFALSE(inherits(data, c("data.table", "data.frame")))) {
    stop("data must be a data table or data frame,",
         "and is an element of a wsub, bsub, wsubmargins, bsubmargins object.")
  }
  
  tmp <- copy(data)

  plot <- ggplot(tmp, aes(x = MinSubstituted, y = Mean)) +
    geom_hline(yintercept = 0, size = 0.2, linetype = 2) +
    geom_vline(xintercept = 0, size = 0.2, linetype = 2) +
    geom_line(aes(colour = Substitute), size = 1) +
    geom_ribbon(aes(ymin = CI_low,
                    ymax = CI_high, fill = Substitute),
                alpha = 1 / 10, size = 1 / 10) +
    scale_color_simpsons() + # overwrite this to specify colour
    scale_fill_simpsons() + # overwrite this to specify fill
    facet_grid(~ Substitute) +
    xlab(paste("Change in", eval(x), sep = " ")) +
    ylab(paste("Change in", eval(y), sep = " ")) +
    theme_cowplot(font_family = eval(font), font_size = 12, line_size = 0)
  plot
}

#' Marginal effects of composition plot
#'
#' This function is useful for visualising the
#' estimated changes in compositional outcome variable
#' associated with a specific change in some independent variable.
#'
#' @param data A \code{\link{emmcoda}} object to use for plot.
#' @param x A character string specifying name of the predictor variable.
#' @param y A character string specifying the name of the compositional outcome variable.
#' @param font A character string specifying the name of user's preferred font.
#' Default is \code{serif} which will choose a serif font available on the system.
#' @param ... Further arguments passed to \code{\link{ggplot}}.
#'
#' @return A ggplot graph object showing the estimated change in compositional outcome
#' associated with changes in predictor.
#' @importFrom ggplot2 ggplot aes geom_area xlab ylab guides guide_legend
#' @importFrom cowplot theme_cowplot
#' @importFrom ggsci scale_fill_simpsons
#' @importFrom data.table copy
#' @export
#' @examples
#'
#' data(mcompd)
plotemmc <- function(data, x, y, font = "serif", ...) {

  tmp <- copy(data)
  specs <- tmp$Predictor

  plot <- ggplot(tmp$emmLong, aes(get(specs), value, fill = variable)) +
    geom_area(alpha = .9, size = 1.5, colour = "white") +
    scale_fill_simpsons()+
    xlab(paste("Change in", eval(x), sep = " ")) +
    ylab(paste("Change in", eval(y), sep = " ")) +
    guides(colour = guide_legend(eval(y))) + 
    theme_cowplot(font_family = eval(font), font_size = 12, line_size = .25)

  plot
}

