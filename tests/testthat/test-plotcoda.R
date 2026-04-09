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

cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)

suppressWarnings(
  m <- brmcoda(complr = cilr,
               formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                 wz1_1 + wz2_1 + wz3_1 + wz4_1 + Female + (1 | ID),
               chain = 1, iter = 500, seed = 123,
               backend = backend,
               silent = 2, refresh = 0))
foreach::registerDoSEQ()

x_multi <- substitution(object = m, delta = 1:10,
                        ref = "grandmean",
                        level = c("between", "within"))
x_single <- substitution(object = m, delta = 5,
                         ref = "grandmean",
                         level = c("between", "within"))

# Testing
#---------------------------------------------------------------------------------------------------
test_that("plot.substitution returns a ggplot for multiple deltas", {
  tmp <- plot(x_multi)
  expect_s3_class(tmp, "ggplot")
  expect_length(tmp$layers, 4L)
})

test_that("plot.substitution returns a ggplot for a single delta", {
  tmp <- plot(x_single)
  expect_s3_class(tmp, "ggplot")
  expect_length(tmp$layers, 2L)
})

test_that("plot.brmcoda returns a non-empty list of plots", {
  tmp <- plot(m, combo = c("hist", "trace"))
  expect_type(tmp, "list")
  expect_gt(length(tmp), 0L)
})

test_that("pairs.brmcoda returns a bayesplot grid", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  suppressWarnings(tmp <- pairs(m, off_diag_args = list(size = 0.5)))
  expect_true(inherits(tmp, "bayesplot_grid"))
  expect_s3_class(tmp, "gtable")
})

test_that("mcmc_plot.brmcoda returns a ggplot", {
  tmp <- mcmc_plot(m, type = "trace")
  expect_s3_class(tmp, "ggplot")
})
