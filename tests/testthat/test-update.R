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
                 formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
                   wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = backend))

parts <- colnames(psub)

# Tests update.complr -----------------------------------------------------------------------------

test_that("update.complr errors where appropriate", {
  
  # missing newdata
  expect_error(cilr_new <- update(object = cilr))
  
  # newdata missing comp vars
  expect_error(cilr_new <- update(object = cilr, newdata = mcompd[, -parts, with = FALSE]))
  
  # newdata missing comp vars
  expect_error(cilr_new <- update(object = cilr, newdata = mcompd[, -c("TST")]))
  
  # newdata incorrect ID
  expect_error(cilr_new <- update(object = cilr, newdata = mcompd[, -c("ID")]))
  
  # incorrect newdata
  expect_error(cilr_new <- update(object = cilr, newdata = list("a" = 1, "b" = 2)))
})

test_that("update.complr gives expected output", {
  
  newcomplr <- update(cilr, newdata = mcompd[ID != 1:10])
  
  expect_true(inherits(newcomplr, "complr"))
  
  expect_true(identical(str(cilr), str(newcomplr)))
  expect_true(identical(newcomplr$data, mcompd[ID != 1:10]))
  expect_true(identical(cilr$parts, newcomplr$parts))
  
  expect_true(identical(ncol(newcomplr$between_comp), ncol(cilr$between_comp)))
  expect_true(identical(ncol(newcomplr$within_comp), ncol(cilr$within_comp)))
  expect_true(identical(ncol(newcomplr$comp), ncol(cilr$comp)))
  expect_true(identical(ncol(newcomplr$between_logratio), ncol(cilr$between_logratio)))
  expect_true(identical(ncol(newcomplr$within_logratio), ncol(cilr$within_logratio)))
  expect_true(identical(ncol(newcomplr$logratio), ncol(cilr$logratio)))
  
  expect_true(identical(nrow(newcomplr$between_comp), nrow(mcompd[ID != 1:10])))
  expect_true(identical(nrow(newcomplr$within_comp), nrow(mcompd[ID != 1:10])))
  expect_true(identical(nrow(newcomplr$comp), nrow(mcompd[ID != 1:10])))
  expect_true(identical(nrow(newcomplr$between_logratio), nrow(mcompd[ID != 1:10])))
  expect_true(identical(nrow(newcomplr$within_logratio), nrow(mcompd[ID != 1:10])))
  expect_true(identical(nrow(newcomplr$logratio), nrow(mcompd[ID != 1:10])))
  
})

# Tests update.brmcoda -----------------------------------------------------------------------------

test_that("update.brmcoda gives errors and warnings where appropriate", {
  
  # warning when newdata and newcomplr provided
  expect_warning(fit_new <- update(object = fit, 
                                   newdata = mcompd[ID != 1:10], 
                                   newcomplr = update(cilr, mcompd[ID != 1:10])))
  
  ## only newcomplr provided, incorrect
  expect_error(fit_new <- update(object = fit, newcomplr = mcompd))
  
  ## missing args
  expect_error(fit_new <- update(object = fit))
})

test_that("update gives expected output", {
  
  # updating only formula
  fit_newformula <- update(fit, formula. = ~ . - wilr1)
  
  expect_true(inherits(fit_newformula, "brmcoda"))
  expect_true(is.null(as.data.table(fit_newformula$model$fit)$b_wilr1))
  
  # updating only data
  fit_newdat <- update(fit, newdata = mcompd[ID != 1:10])
  
  expect_true(inherits(fit_newdat, "brmcoda"))
  expect_true(identical(fit_newdat$complr$data, mcompd[ID != 1:10]))
  
  # updating both formula and data
  fit_new <- update(fit, 
                    formula. = ~ . - wilr2,
                    newdata = mcompd[ID != 1:10])
  
  expect_true(inherits(fit_new, "brmcoda"))
  expect_true(is.null(as.data.table(fit_new$model$fit)$b_wilr2))
  expect_true(identical(fit_newdat$complr$data, mcompd[ID != 1:10]))
  
})
