skip_on_cran()

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  backend <- "rstan"
  ## if using rstan backend, models can crash on Windows
  ## so skip if on windows and cannot use cmdstanr
  skip_on_os("windows")
} else {
  if (isFALSE(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))) {
    backend <- "cmdstanr"
  }
}

# Packages
library(testthat)
library(multilevelcoda)
library(brms)
library(ggplot2)

# Model
#---------------------------------------------------------------------------------------------------
data(mcompd)
data(sbp)
data(psub)

cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")

suppressWarnings(
  m <- brmcoda(compilr = cilr,
               formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
               chain = 1, iter = 500, seed = 123,
               backend = backend))
foreach::registerDoSEQ()

x <- bsub(object = m, base = psub, delta = 2)

# Testing
#---------------------------------------------------------------------------------------------------

test_that("plotsub errors for invalid input", {
  
  ## incorrect dataset
  expect_error(p <- plotsub(x, "Sleep", "Stress"))

  ## missing arguments
  expect_error(p <- plotsub(x = "Sleep", y = "Stress"))
  expect_error(p <- plotsub(data = x, x = "Sleep"))
  expect_error(p <- plotsub(data = x, y = "Stress"))
  expect_error(p <- plotsub(x$TST))

})
test_that("plot have known output", {
  
  p0 <-  plotsub(x$TST, "Sleep", "Stress")
  p1 <- ggplot(x$TST, aes(x = Delta, y = Mean)) +
    geom_hline(yintercept = 0, linewidth = 0.2, linetype = 2) +
    geom_vline(xintercept = 0, linewidth = 0.2, linetype = 2) +
    geom_line(aes(colour = From), linewidth = 1) +
    geom_ribbon(aes(ymin = CI_low,
                    ymax = CI_high, fill = From),
                alpha = 1 / 10, linewidth = 1 / 10) +
    facet_grid(~ From) +
    xlab("Change to Sleep") +
    ylab("Change in Stress")

  expect_s3_class(p0, "ggplot")
})

