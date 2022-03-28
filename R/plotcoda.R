#' Generate substitution plot
#'
#' This function is useful for visualising the
#' estimated differences in outcomes when compositional variables
#' are substituted for a specific period of time.
#'
#' @param data A dataset to use for plot. 
#' It must be a component of a list resulted from one of the following functions:
#' \code{wsub}, \code{bsub}, \code{wsubmargins}, \code{bsubmargins}.
#' @param iv A character string indicating name of the compostional predictor variable.
#' @param dv A character string indicating the name of the outcome variable.
#' @param font A character string indicating the name of user's preferred font. Default is \code{Times New Roman}.
#' @param ... Further arguments passed to \code{ggplot}.
#' 
#' @return A ggplot graph object showing the estimated difference in outcome when 
#' each pair of compositional variables are substituted for a specific time.
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_line geom_ribbon facet_grid xlab ylab
#' @importFrom cowplot theme_cowplot
#' @importFrom ggsci scale_color_simpsons
#' @importFrom data.table copy
#' @export
#' @examples
#' 
#' data(mcompd)
#' plotsub(data = bsubmarginstest$TST, iv = "sleep", dv = "stress")
#' 
plotsub <- function(data, iv, dv, font = "Times New Roman", ...) {
  
  if (isFALSE(inherits(data, c("data.table", "data.frame")))) {
    stop("data must be a data table or data frame")
  }
  
  tmp <- copy(data)
  
  niceplot <- ggplot(tmp, aes(x = MinSubstituted, y = Mean, colour = Substitute)) +
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

#' #' Generate marginal effects of composition plot
#' #'
#' #' This function is useful for visualising the
#' #' estimated changes in compositional outcome variable
#' #' associated with a specific change in some independent variable.
#' #'
#' #' @param data A dataset to use for plot. 
#' #' It must be a \code{\link{mvmcoda}} or \code{\link{mvcoda}} object 
#' #' @param compilr A \code{compilr} object containing data of composition, ILR coordinates,
#' #' and other variables used for plot.
#' #' @param iv A character string indicating name of the predictor variable.
#' #' @param dv A character string indicating the name of the compositional outcome variable.
#' #' @param font A character string indicating the name of user's preferred font. Default is \code{Times New Roman}.
#' #' @param ... Further arguments passed to \code{ggplot}.
#' #' 
#' #' @return A ggplot graph object showing the estimated difference in outcome when 
#' #' each pair of compositional variables are substituted for a specific time.
#' #' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_line geom_ribbon facet_grid xlab ylab
#' #' @importFrom cowplot theme_cowplot
#' #' @importFrom ggsci scale_color_simpsons
#' #' @importFrom data.table copy
#' #' @export
#' #' @examples
#' #' 
#' #' data(mcompd)
#' #' plotsub(data = bsubctest2$TST, iv = "sleep")
#' plotcoda <- function(data, iv) {
#'   psi <- mvcoda$psi
#' 
#' }