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
  m <- brmcoda(compilr = cilr,
               formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
               chain = 1, iter = 500, seed = 123,
               backend = backend))
foreach::registerDoSEQ()

x <- substitution(object = m, basesub = psub, delta = 2)

# Testing
#---------------------------------------------------------------------------------------------------

test_that("substitution errors for invalid input", {
  
  ## missing object
  expect_error(x <- substitution(basesub = psub, delta = 2))
  
  ## not brmcoda model
  m1 <- lmer(STRESS ~ 1 + (1 | ID), data = mcompd)
  expect_error(x <- substitution(object = m1, basesub = psub, delta = 2))
  
  ## incorrect delta
  expect_error(x <- substitution(object = m, basesub = psub, delta = -10))

  ## missing delta
  expect_error(x <- substitution(object = m, basesub = psub))
  
  ## reference grid has matching names with ILRs
  rg <- data.table(bilr1 = 1)
  expect_error(x <- substitution(object = m, basesub = psub, delta = 2, regrid = rg))
  
  ## basesub does not have the same components as parts in cilr
  ps <- basesub(c("WAKE", "MVPA", "LPA", "SB"))
  expect_error(x <- substitution(object = m, basesub = ps, minute = 2))
  
  ## basesub does have the same names as parts in cilr
  ps <- basesub(parts = c("Sleep", "WAKE", "MVPA", "LPA", "SB"))
  expect_error(x <- substitution(object = m, basesub = ps, minute = 2))
  
})
