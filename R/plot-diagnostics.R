## clear R CMD CHECK notes
if (getRversion() >= "2.15.1")  utils::globalVariables(c(
  "obs_id", "is_extremevalue", "xvar", "yvar", "y", "y0", "xend", "yend",
  "rug_y0", "rug_y1"))

#' Diagnostics Plot for Compositional Diagnostics
#'
#' Make a scatterplot matrix of \code{\link{diagnostics}} results from a
#' \code{\link{complr}} object, with density and rug plots on the diagonal.
#'
#' @param x A \code{diagnostics} class object created by
#' \code{\link{diagnostics.complr}}.
#' @param transform Character string indicating whether to plot the raw
#' composition values in \code{x$x} or centered logratio transformed values.
#' Defaults to \code{"clr"}.
#' @param ... Currently unused.
#'
#' @return A ggplot graph object showing pairwise scatterplots for each
#' composition part, with point colour indicating extreme value status. Diagonal
#' panels show marginal densities and rug marks for observations flagged as
#' extreme values.
#'
#' @importFrom compositions clr
#' @importFrom ggplot2 ggplot aes facet_grid geom_point geom_ribbon geom_segment labs
#' @importFrom ggplot2 scale_colour_manual theme theme_bw element_blank element_rect element_text
#' @importFrom data.table .I := as.data.table data.table melt rbindlist
#' @importFrom stats density
#'
#' @method plot diagnostics
#' @export
#'
#' @examples
#' data(mcompd)
#' data(sbp)
#'
#' ids <- unique(mcompd$ID)[1:20]
#' cilr <- complr(
#'   data = mcompd[ID %in% ids, .SD[1:3], by = ID],
#'   sbp = sbp,
#'   parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'   idvar = "ID",
#'   total = 1440
#' )
#'
#' between_diag <- diagnostics(cilr, level = "between", ev.perc = .05)
#'
#' # Plot centered logratio transformed composition values
#' plot(between_diag)
#'
#' # Plot on the raw composition scale
#' plot(between_diag, transform = "raw")
plot.diagnostics <- function(x, transform = c("clr", "raw"), ...) {
  stopifnot(is.diagnostics(x))

  if (nrow(x$x) != nrow(x$distance)) {
    stop("x$x and x$distance must have the same number of rows.")
  }

  transform <- match.arg(transform)
  X_dt <- as.data.table(switch(transform,
    raw = x$x,
    clr = clr(x$x)
  ))
  D_dt <- as.data.table(x$distance)
  vars <- colnames(x$x)

  ## hack for now, just error if these variable names used
  ## todo later come up with more elegant approach to avoid name conflicts with user data
  if (any(c("is_extremevalue", "obs_id") %in% names(X_dt))) {
    stop("data must not contain a column named 'is_extremevalue' or 'obs_id'.")
  }

  X_dt[, obs_id := .I]
  X_dt[, is_extremevalue := D_dt$is_extremevalue]

  # scatter plot data
  x_long <- melt(X_dt,
    id.vars       = c("obs_id", "is_extremevalue"),
    measure.vars  = vars,
    variable.name = "xvar",
    value.name    = "x")

  y_long <- melt(X_dt,
    id.vars       = c("obs_id", "is_extremevalue"),
    measure.vars  = vars,
    variable.name = "yvar",
    value.name    = "y")

  scatter_dt <- merge(
    x_long, y_long,
    by = c("obs_id", "is_extremevalue"),
    allow.cartesian = TRUE)[xvar != yvar]

  # diagonal density / rug data
  range_dt <- melt(X_dt,
    id.vars       = c("obs_id", "is_extremevalue"),
    measure.vars  = vars,
    variable.name = "var",
    value.name    = "value")
  range_dt <- range_dt[, {
      r <- range(value, na.rm = TRUE)
      span <- diff(r)
      list(
        vmin    = r[1],
        vmax    = r[2],
        span    = span,
        dens_y0 = r[1] + 0.10 * span, # density baseline
        dens_y1 = r[1] + 0.80 * span, # density top
        rug_y0  = r[1] + 0.02 * span, # rug start
        rug_y1  = r[1] + 0.09 * span # rug end
      )
    },
    by = var]

  # calculate density for each variable for diagonal panels
  # done manually so we can scale so they fit nicely in the
  # panel, otherwise ranges can be messed up
  diag_density_dt <- rbindlist(lapply(vars, function(v) {
    rr <- range_dt[var == v]
    z <- X_dt[[v]]
    z <- z[is.finite(z)]

    ## if fewer than 5 unique values, just make a flat density
    ## otherwise use density() function
    if (length(unique(z)) < 5L) {
      d_x <- seq(rr$vmin, rr$vmax, length.out = 100)
      d_y <- rep(1, length(d_x))
    } else {
      d <- density(z, na.rm = TRUE,
        from  = rr$vmin, to = rr$vmax,
        cut = 0, n = 512)
      d_x <- d$x
      d_y <- d$y
    }

    data.table(xvar = v, yvar = v,
      x = d_x, y0 = rr$dens_y0,
      y = rr$dens_y0 + (d_y / max(d_y)) * (rr$dens_y1 - rr$dens_y0))
  }))

  ## rug for diagonals, only for extreme values
  diag_rug_dt <- melt(
    X_dt[is_extremevalue == TRUE],
    id.vars = c("obs_id", "is_extremevalue"),
    measure.vars = vars, variable.name = "var", value.name = "x")
  diag_rug_dt <- diag_rug_dt[range_dt, on = .(var)]
  diag_rug_dt <- diag_rug_dt[, .(
      xvar = var, yvar = var,
      x = x, xend = x,
      y = rug_y0, yend = rug_y1)]

  # Keep facet order fixed
  set_facets <- function(x) {
    x[, xvar := factor(xvar, levels = vars)]
    x[, yvar := factor(yvar, levels = vars)]
    x
  }

  scatter_dt <- set_facets(scatter_dt)
  diag_density_dt <- set_facets(diag_density_dt)
  diag_rug_dt <- set_facets(diag_rug_dt)

  ## make the plot
  out <- ggplot() +
    geom_point(data = scatter_dt,
      aes(x = x, y = y, colour = is_extremevalue),
      size = 1.2, shape = 16) +
    geom_ribbon(data = diag_density_dt,
      aes(x = x, ymin = y0, ymax = y),
      inherit.aes = FALSE,
      fill = "grey80", colour = NA) +
    geom_segment(data = diag_rug_dt,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE, linewidth = 0.35,
      colour = "black") +
    scale_colour_manual(
      values = c(`FALSE` = "grey70", `TRUE` = "black"),
      guide = "none") +
    facet_grid(yvar ~ xvar, scales = "free") +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.spacing = grid::unit(0.15, "lines"),
      strip.background = element_rect(fill = "grey95"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  out
}
