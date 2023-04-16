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

# Model
#---------------------------------------------------------------------------------------------------
data(mcompd)
data(sbp)
data(psub)

cilr <- compilr(data = mcompd[ID %in% 1:200, .SD[1:5], by = ID], sbp = sbp,
                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")

suppressWarnings(
  fit <- brmcoda(compilr = cilr,
               formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
               chain = 1, iter = 500, seed = 123,
               backend = backend))

parts <- colnames(psub)

# Tests update.compilr -----------------------------------------------------------------------------

test_that("update.compilr errors where appropriate", {
  
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

test_that("update.compilr gives expected output", {
  
  newcilr <- update(cilr, mcompd[ID != 1:10])
  
  expect_true(inherits(newcilr, "compilr"))
  
  expect_true(identical(str(cilr), str(newcilr)))
  expect_true(identical(newcilr$data, mcompd[ID != 1:10]))
  expect_true(identical(cilr$parts, newcilr$parts))
  
  expect_true(identical(ncol(newcilr$BetweenComp), ncol(cilr$BetweenComp)))
  expect_true(identical(ncol(newcilr$WithinComp), ncol(cilr$WithinComp)))
  expect_true(identical(ncol(newcilr$TotalComp), ncol(cilr$TotalComp)))
  expect_true(identical(ncol(newcilr$BetweenILR), ncol(cilr$BetweenILR)))
  expect_true(identical(ncol(newcilr$WithinILR), ncol(cilr$WithinILR)))
  expect_true(identical(ncol(newcilr$TotalILR), ncol(cilr$TotalILR)))

  expect_true(identical(nrow(newcilr$BetweenComp), nrow(mcompd[ID != 1:10])))
  expect_true(identical(nrow(newcilr$WithinComp), nrow(mcompd[ID != 1:10])))
  expect_true(identical(nrow(newcilr$TotalComp), nrow(mcompd[ID != 1:10])))
  expect_true(identical(nrow(newcilr$BetweenILR), nrow(mcompd[ID != 1:10])))
  expect_true(identical(nrow(newcilr$WithinILR), nrow(mcompd[ID != 1:10])))
  expect_true(identical(nrow(newcilr$TotalILR), nrow(mcompd[ID != 1:10])))
  
})
  
# Tests update.brmcoda -----------------------------------------------------------------------------

test_that("update.brmcoda gives errors and warnings where appropriate", {
  
  # warning when newdata and newcilr provided
  expect_warning(fit_new <- update(object = fit, 
                                   newdata = mcompd[ID != 1:10], 
                                   newcilr = update(cilr, mcompd[ID != 1:10])))
  
  ## only newcilr provided, incorrect
  expect_error(fit_new <- update(object = fit, newcilr = mcompd))
  
  ## missing args
  expect_error(fit_new <- update(object = fit))
})

test_that("update gives expected output", {
  
  # updating only formula
  fit_newformula <- update(fit, formula. = ~ . - wilr1)
  
  expect_true(inherits(fit_newformula, "brmcoda"))
  expect_true(is.null(as.data.table(fit_newformula$Model$fit)$b_wilr1))

  # updating only data
  fit_newdat <- update(fit, newdata = mcompd[ID != 1:10])
  
  expect_true(inherits(fit_newdat, "brmcoda"))
  expect_true(identical(fit_newdat$CompILR$data, mcompd[ID != 1:10]))
  
  # updating both formula and data
  fit_new <- update(fit, 
                    formula. = ~ . - wilr2,
                    newdata = mcompd[ID != 1:10])
  
  expect_true(inherits(fit_new, "brmcoda"))
  expect_true(is.null(as.data.table(fit_new$Model$fit)$b_wilr2))
  expect_true(identical(fit_newdat$CompILR$data, mcompd[ID != 1:10]))
  
})
