#' Generate multilevel substitution plot
#'
#' This function is designed to generate a plot showing the
#' estimated differences in outcomes when compositional variables
#' are substituted for a specific period of time
#'
#' @param object A dataset resulted from the substitution model.
#' @param predictor A name of the compostional predictor variable
#' that is substituted to/from. # nolint
#' @param outcome A name of the outcome variable.
#' @return A ggplot graph object
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_line geom_ribbon facet_grid xlab ylab
#' @importFrom cowplot theme_cowplot
#' @export
#' @examples
#' ## TODO
subplot <- function(object, predictor, outcome) {
  ggplot(object, aes(x = at, y = Mean, colour = Substitute)) +
    geom_hline(yintercept = 0, size = 0.2, linetype = 2) +
    geom_vline(xintercept = 0, size = 0.2, linetype = 2) +
    geom_line(aes(colour = Substitute), size = 1) +
    geom_ribbon(aes(ymin = CI_low,
                    ymax = CI_high, colour = Substitute),
                alpha = 1 / 10, size = 1 / 10) +
    scale_color_simpsons() + # overwrite this in case people want their own colours # nolint
    facet_grid(~Substitute) +
    xlab("Change in compositional predictor") + # overwrite this
    ylab("Change in outcome") + # overwrite this
    theme_cowplot(font_family = "Times New Roman", font_size = 12, line_size = 0)
}
