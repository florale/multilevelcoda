#' @title Substitution plot
#'
#' @description
#' Make a plot of \code{\link{substitution}} model results.
#'
#' @param object A \code{\link{substitution}} object containing the output of substitution models.
#' @param to A character value or vector specifying the names of the compositional parts
#' that were reallocated to in the model.
#' @param ref A character value of ((\code{grandmean} or \code{clustermean} or \code{users}),
#' @param level A character value of (\code{between} or \code{within}).
#' @param ... Further arguments passed to \code{\link{ggplot}}.
#'
#' @return A ggplot graph object showing the estimated difference in outcome when
#' each pair of compositional variables are substituted for a specific time.
#' 
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_line geom_pointrange geom_ribbon facet_grid xlab ylab
#' @importFrom data.table copy
#' 
#' @exportS3Method plot substitution
plot.substitution <- function(object, to,
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
  delta.pos <- object$delta
  delta.neg <- -1*abs(object$delta)
  delta <- c(delta.pos, delta.neg)
  
  # extract data
  if (isTRUE(is.sequential(delta.pos))) {
    tmp <- summary(object = object,
                   delta = delta,
                   to = to,
                   ref = ref,
                   level = level,
                   digits = "asis"
    )
    
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
    tmp <- summary(object = object,
                   delta = delta,
                   to = to,
                   ref = ref,
                   level = level,
                   digits = "asis"
    )
    
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
