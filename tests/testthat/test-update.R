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
library(data.table)
library(multilevelcoda)
library(extraoperators)
library(brms)
library(lme4)

# model
#---------------------------------------------------------------------------------------------------
data(mcompd)
data(sbp)
data(psub)

cilr <- complr(data = mcompd[ID %in% 1:200, .SD[1:5], by = ID], sbp = sbp,
               parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")

suppressWarnings(
  fit <- brmcoda(complr = cilr,
                 formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + Female + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = backend))

parts <- colnames(psub)

# Tests update.brmcoda -----------------------------------------------------------------------------

# test_that("update.brmcoda gives errors and warnings where appropriate", {
#   
#   # warning when newdata and newcomplr provided
#   expect_warning(fit_new <- update(object = fit, 
#                                    newdata = mcompd[ID != 1:10], 
#                                    newcomplr = update(cilr, mcompd[ID != 1:10])))
#   
#   ## only newcomplr provided, incorrect
#   expect_error(fit_new <- update(object = fit, newcomplr = mcompd))
#   
#   ## missing args
#   expect_error(fit_new <- update(object = fit))
# })
# 
# test_that("update gives expected output", {
#   
#   # updating only formula
#   fit_newformula <- update(fit, formula. = ~ . - wz1_1)
#   
#   expect_true(inherits(fit_newformula, "brmcoda"))
#   expect_true(is.null(as.data.table(fit_newformula$model$fit)$b_wz1_1))
#   
#   # updating only data
#   fit_newdat <- update(fit, newdata = mcompd[ID != 1:10])
#   
#   expect_true(inherits(fit_newdat, "brmcoda"))
#   expect_true(identical(fit_newdat$complr$data, mcompd[ID != 1:10]))
#   
#   # updating both formula and data
#   fit_new <- update(fit, 
#                     formula. = ~ . - wz2_1,
#                     newdata = mcompd[ID != 1:10])
#   
#   expect_true(inherits(fit_new, "brmcoda"))
#   expect_true(is.null(as.data.table(fit_new$model$fit)$b_wz2_1))
#   expect_true(identical(fit_newdat$complr$data, mcompd[ID != 1:10]))
#   
# })
