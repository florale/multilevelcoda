#' Generate multilevel substitution plot
#'
#' This function is designed to generate a plot showing the
#' estimated differences in outcomes when compositional variables
#' are substituted for a specific period of time
#'
#' @param data A dataset to use for plot. 
#' @param data It must be a component of a list computed using either of the following:
#' @param data \code{wsub}, \code{bsub}, \code{wsubmargin}, \code{bsubmargin}.
#' @param iv A character string indicating name of the compostional predictor variable
#' @param dv A character string indicating the name of the outcome variable.
#' @param font A character string indicating the name of preferred font. Default is Times New Roman.
#' 
#' @return A ggplot graph object showing the estimated difference in outcome when 
#' @return each pair of compositional variables are substituted for a specific time.
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_line geom_ribbon facet_grid xlab ylab
#' @importFrom cowplot theme_cowplot
#' @importFrom ggsci scale_color_simpsons
#' @export
#' @examples
#' ## TODO
#' 
#' data(mcompd)
#' plotsub(data = bsubctest$TST, iv = "sleep", dv = "stress")
#' 
plotsub <- function(data, iv, dv, font = "Times New Roman") {
  
  niceplot <- ggplot(data, aes(x = MinSubstituted, y = Mean, colour = Substitute)) +
    geom_hline(yintercept = 0, size = 0.2, linetype = 2) +
    geom_vline(xintercept = 0, size = 0.2, linetype = 2) +
    geom_line(aes(colour = Substitute), size = 1) +
    geom_ribbon(aes(ymin = CI_low,
                    ymax = CI_high, colour = Substitute),
                alpha = 1 / 10, size = 1 / 10) +
    scale_color_simpsons() + # overwrite this to specify colour
    facet_grid(~ Substitute) +
    xlab(paste("Change in", eval(iv), sep = " ")) +
    ylab(paste("Change in", eval(dv), sep = " ")) +
    theme_cowplot(font_family = eval(font), font_size = 12, line_size = 0)
  
  return(niceplot)
}
