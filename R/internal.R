## make Rcmd check happy
utils::globalVariables(c("i",  "..cols", ".", "To", ".SD", "t", "head",  "fitted", 
                         "x", "object", "predict", "residuals", "tail", "vcov", "coef",
                         "Mean",  "CI_low", "CI_high", "From", "Delta", "pairs",
                         "spread", "value", "variable", "ID", "EffectType", "Level", "Reference",
                         "update", "posterior",
                         "sim", "condition",
                         "est", "lower", "upper", "JI", "N", "K", "D", "Stat", "Estimates", "MCSE",
                         "sigma", "OnTarget"
                         ))


# expand grid data frame
expand.grid.df <- function(...) Reduce(function(...) merge.data.frame(..., by = NULL, all = TRUE), list(...))

# check sequence of number
is.sequential <- function(x) {
  all(length(x) > 2 & all(abs(diff(x)) == 1))
}

#' @keywords internal
#' 
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_linerange
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_void
#' @importFrom ggplot2 unit
#' @importFrom hrbrthemes theme_ipsum 
NULL

# plot for shiny sim
.par_plot <- function(data, shiny = FALSE, d = 4, font = "Arial Narrow") {
  
  # colour palette --------------
  col_brmcoda_d3 <- 
    c("#9A5C7D", "#B98AA3", "#DCD5CE", "#8DA290", "#708885", "#5A6367", 
      "#1C1718")
  col_brmcoda_d4 <- 
    c("#9A5C7D", "#B98AA3", "#DCD5CE", "#8DA290", "#708885", "#5A6367", 
      "#456691", "#2A3E59", 
      "#1C1718")
  col_brmcoda_d5 <- 
    c("#9A5C7D", "#B98AA3", "#DCD5CE", "#8DA290", "#708885", "#5A6367", 
      "#456691", "#2A3E59", 
      "#9c8aa4", "#5E4F65", 
      "#1C1718")
  
  col_sub_d3 <- 
    c("#bf5b4b", "#A69188",
      "#EAD3BF", "#FAD899",
      "#8DA290", "#133A1B"
    )
  col_sub_d4 <-
    c("#2A3E59", "#456691",
      "#944C4C", "#C99696",
      "#bf5b4b", "#A69188",
      "#EAD3BF", "#FAD899",
      "#8DA290", "#133A1B",
      "#6d765b", "#3d251e"
    )
  col_sub_d5 <- 
    c("#1C1718", "#2A3E59",
      "#456691", "#647F9A",
      "#8CAACB", "#DCD5CE",
      "#DAA5AE", "#b6485d",
      "#944C4C", "#C99696",
      "#bf5b4b", "#bb847a",
      "#A69188", "#EAD3BF",
      "#FAD899", "#8DA290",
      "#133A1B", "#6d765b",
      "#3b4031", "#3d251e"
    )
  
  if (all(data$Stat == "bias")) {
    ylab <- "Bias"
    yintercept <- 0
    if ("Estimand" %in% colnames(data)) {
      y_lims <- c(-0.16, 0.16)
      y_breaks <- c(-0.1, 0, 0.1)
    } else {
      y_lims <- c(-0.075, 0.075)
      y_breaks <- c(-0.05, 0, 0.05)
    }
  } else if (all(data$Stat == "cover")) {
    ylab <- "Coverage"
    yintercept <- 0.95
    if ("Estimand" %in% colnames(data)) {
      y_lims <- c(0.9, 1)
      y_breaks <- c(0.9, 0.95, 1)
    } else {
      y_lims <- c(0.9, 1)
      y_breaks <- c(0.9, 0.95, 1)
    }
  } else if (all(data$Stat == "becover")) {
    ylab <- "Bias-Eliminated Coverage"
    yintercept <- 0.95
    if ("Estimand" %in% colnames(data)) {
      y_lims <- c(0.9, 1)
      y_breaks <- c(0.9, 0.95, 1)
    } else {
      y_lims <- c(0.9, 1)
      y_breaks <- c(0.9, 0.95, 1)
    }
  } else if (all(data$Stat == "mse")) {
    ylab <- "Empirical Standard Error"
    yintercept <- 0
    if ("Estimand" %in% colnames(data)) {
      y_lims <- c(0, 3.5)
      y_breaks <- c(0, 1.5, 3)
    } else {
      y_lims <- c(0, 1)
      y_breaks <- c(0, 0.5, 1)
    }
  } else if (all(data$Stat == "empse")) {
    ylab <- "Mean-squared Error"
    yintercept <- 0
    if ("Estimand" %in% colnames(data)) {
      y_lims <- c(0, 3)
      y_breaks <- c(0, 1.5, 3)
    } else {
      y_lims <- c(0, 1)
      y_breaks <- c(0, 0.5, 1)
    }
  }
  
  if (d == 4) {
    if ("Substitution" %in% colnames(data)) {
      xvar <- data$Substitution
      xtext <- 13
    } else {
      xvar <- data$Estimand
      xtext <- 10
    }
  } else if (d == 3) {
    if ("Substitution" %in% colnames(data)) {
      xvar <- data$Substitution
      xtext <- 7
    } else {
      xvar <- data$Estimand
      xtext <- 8
    }
  } else if (d == 5) {
    if ("Substitution" %in% colnames(data)) {
      xvar <- data$Substitution
      xtext <- 21
    } else {
      xvar <- data$Estimand
      xtext <- 12
    }
  }
  
  if (nlevels(xvar) == 7) {
    col <- col_brmcoda_d3
  } else if (nlevels(xvar) == 9) {
    col <- col_brmcoda_d4
  } else if (nlevels(xvar) == 11) {
    col <- col_brmcoda_d5
  } else if (nlevels(xvar) == 6) {
    col <- col_sub_d3
  } else if (nlevels(xvar) == 12) {
    col <- col_sub_d4
  } else if (nlevels(xvar) == 20) {
    col <- col_sub_d5
  }
  
  point_size <- ifelse(shiny == TRUE, 2, 2.25)
  line_size <- ifelse(shiny == TRUE, 0.75, 0.75)
  btext_size <- ifelse(shiny == TRUE, 14, 12)
  text_size <- ifelse(shiny == TRUE, 12, 13)
  yseg <- y_breaks[[1]]
  yendseg <- y_breaks[[3]]
  
  if (shiny == TRUE) {
    gg <- 
      ggplot(data, 
             aes(x = xvar, y = est, 
                 ymin = lower, ymax = upper,
                 colour = xvar)) +
      geom_hline(yintercept = yintercept, color = "#666666", linetype = "dashed", linewidth = 0.5) +
      geom_point(size = point_size) +
      geom_linerange(linewidth = line_size) +
      labs(x = "", y = ylab, colour = "Parameter") +
      scale_colour_manual(values = col) +
      scale_y_continuous(limits = y_lims,
                         breaks = y_breaks) +
      scale_x_discrete(drop = FALSE) +
      # facet_wrap(ggplot2::vars(N, K), labeller = ggplot2::label_both) +
      facet_wrap(ggplot2::vars(JI), labeller = ggplot2::label_context, strip.position = "top") +
      theme_ipsum() +
      coord_flip() +
      theme(
        axis.ticks        = element_blank(),
        panel.background  = element_rect(fill = "transparent", colour = "black", linewidth = line_size),
        panel.border      = element_rect(fill = "transparent", colour = "black", linewidth = line_size),
        # panel.grid.major  = element_blank(),
        # panel.grid.minor  = element_blank(),
        plot.background   = element_rect(fill = "transparent", colour = NA),
        axis.title.y      = element_text(size = btext_size, face = "bold"),
        axis.title.x      = element_text(size = btext_size, face = "bold"),
        axis.text.x       = element_text(size = text_size),
        axis.text.y       = element_blank(),
        title             = element_text(size = btext_size, face = "bold"),
        legend.text       = element_text(size = text_size),
        strip.text.x      = element_text(size = text_size),
        legend.position   = "none",
        panel.spacing.y   = unit(0, "lines"),
        panel.spacing.x   = unit(0.75, "lines")
        # strip.text.x      = element_blank()
      )
    plotly::ggplotly(gg, height = 1300)
    
  } else {
    gg <- 
      ggplot(data, 
             aes(x = xvar, y = est, 
                 ymin = lower, ymax = upper,
                 colour = xvar)) +
      geom_segment(aes(x = 0.5, xend = xvar, y = yintercept, yend = yintercept), color = "#666666", linetype = "dashed", linewidth = 0.5) +
      geom_segment(aes(y = yseg, yend = yendseg, x = 0.5, xend = 0.5), color = "black", linewidth = 0.5) +
      geom_text(aes(label = JI, y = yintercept, x = xtext), color = "black", family = font, vjust = "inward", hjust = "inward") +
      # geom_hline(yintercept = yintercept, color = "#666666", linetype = "dashed", linewidth = 0.5) +
      geom_point(size = point_size) +
      geom_linerange(linewidth = line_size) +
      # geom_segment(aes(x = "sigma", xend = xvar, y = yintercept, yend = yintercept), color = "#666666", linetype = "dashed", linewidth = 0.25) +
      labs(x = "", y = ylab, colour = "Parameter") +
      scale_colour_manual(values = col) +
      scale_y_continuous(limits = y_lims,
                         breaks = y_breaks) +
      scale_x_discrete(drop = FALSE, expand = c(0,1.05)) +
      # facet_wrap(ggplot2::vars(N, K), labeller = ggplot2::label_both) +
      # facet_wrap(ggplot2::vars(NK), labeller = ggplot2::label_context, strip.position = "top") +
      hrbrthemes::theme_ipsum() + theme_void() +
      coord_flip() +
      theme(
        axis.ticks        = element_blank(),
        panel.background  = element_rect(fill = "transparent", colour = NA, linewidth = line_size),
        panel.border      = element_rect(fill = "transparent", colour = NA, linewidth = line_size),
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        plot.background   = element_rect(fill = "transparent", colour = NA),
        axis.title.y      = element_text(size = btext_size, face = "bold"),
        # axis.title.x      = element_text(size = btext_size, face = "bold"),
        # axis.text.x       = element_text(size = text_size),
        # axis.text.y       = element_blank(),
        # title             = element_blank(),
        legend.text       = element_blank(),
        # strip.text.x      = element_text(size = text_size),
        legend.position   = "none",
        # panel.grid.major  = element_blank(),
        # panel.grid.minor  = element_blank(),
        axis.title.x      = element_blank(),
        axis.line.y       = element_blank(),
        axis.text.x       = element_text(size = text_size, family = font),
        axis.text.y       = element_blank()
        # strip.text.x      = element_text(size = text_size, family = font),
        # strip.background  = element_blank(),
        # strip.placement   = "outside"
        # strip.text.x      = element_blank()
      )
    gg
    
  }
}
