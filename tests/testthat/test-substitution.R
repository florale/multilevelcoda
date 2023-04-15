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

# Tests
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
  expect_error(x <- substitution(object = m, basesub = ps, delta = 2))
  
  ## basesub does not have the same names as parts in cilr
  ps <- basesub(parts = c("Sleep", "WAKE", "MVPA", "LIPA", "SB"))
  expect_error(x <- substitution(object = m, basesub = ps, delta = 2))
  
})

test_that("substitution computes basesub if missing", {
  x1 <- substitution(object = m, delta = 2)
  expect_true(identical(x, x1))
})

test_that("substitution gives expected outputs", {
  expect_true(identical(names(x), 
                        c("BetweenpersonSub", "WithinpersonSub",
                          "BetweenpersonSubMargins", 'WithinpersonSubMargins')))
  
  ## bsub
  expect_type(x$BetweenpersonSub, "list")
  expect_equal(length(x$BetweenpersonSub), length(m$CompILR$parts))
  expect_s3_class(x$BetweenpersonSub$TST, "data.table")
  expect_s3_class(x$BetweenpersonSub$WAKE, "data.table")
  expect_s3_class(x$BetweenpersonSub$MVPA, "data.table")
  expect_s3_class(x$BetweenpersonSub$LPA, "data.table")
  expect_s3_class(x$BetweenpersonSub$SB, "data.table")
  
  expect_type(x$BetweenpersonSub$TST$Mean, "double")
  expect_type(x$BetweenpersonSub$TST$CI_low, "double")
  expect_type(x$BetweenpersonSub$TST$CI_high, "double")
  expect_type(x$BetweenpersonSub$TST$Delta, "double")
  expect_type(x$BetweenpersonSub$TST$From, "character")
  expect_type(x$BetweenpersonSub$TST$To, "character")
  
  expect_type(x$BetweenpersonSub$WAKE$Mean, "double")
  expect_type(x$BetweenpersonSub$WAKE$CI_low, "double")
  expect_type(x$BetweenpersonSub$WAKE$CI_high, "double")
  expect_type(x$BetweenpersonSub$WAKE$Delta, "double")
  expect_type(x$BetweenpersonSub$WAKE$From, "character")
  expect_type(x$BetweenpersonSub$WAKE$To, "character")
  
  expect_type(x$BetweenpersonSub$MVPA$Mean, "double")
  expect_type(x$BetweenpersonSub$MVPA$CI_low, "double")
  expect_type(x$BetweenpersonSub$MVPA$CI_high, "double")
  expect_type(x$BetweenpersonSub$MVPA$Delta, "double")
  expect_type(x$BetweenpersonSub$MVPA$From, "character")
  expect_type(x$BetweenpersonSub$MVPA$To, "character")
  
  expect_type(x$BetweenpersonSub$LPA$Mean, "double")
  expect_type(x$BetweenpersonSub$LPA$CI_low, "double")
  expect_type(x$BetweenpersonSub$LPA$CI_high, "double")
  expect_type(x$BetweenpersonSub$LPA$Delta, "double")
  expect_type(x$BetweenpersonSub$LPA$From, "character")
  expect_type(x$BetweenpersonSub$LPA$To, "character")
  
  expect_type(x$BetweenpersonSub$SB$Mean, "double")
  expect_type(x$BetweenpersonSub$SB$CI_low, "double")
  expect_type(x$BetweenpersonSub$SB$CI_high, "double")
  expect_type(x$BetweenpersonSub$SB$Delta, "double")
  expect_type(x$BetweenpersonSub$SB$From, "character")
  expect_type(x$BetweenpersonSub$SB$To, "character")
  
  expect_true(ncol(x$BetweenpersonSub$TST) >= 8)
  expect_true(ncol(x$BetweenpersonSub$WAKE) >= 8)
  expect_true(ncol(x$BetweenpersonSub$MVPA) >= 8)
  expect_true(ncol(x$BetweenpersonSub$LPA) >= 8)
  expect_true(ncol(x$BetweenpersonSub$SB) >= 8)
  
  expect_true(all(x$BetweenpersonSub$TST$To == "TST"))
  expect_true(all(x$BetweenpersonSub$WAKE$To == "WAKE"))
  expect_true(all(x$BetweenpersonSub$MVPA$To == "MVPA"))
  expect_true(all(x$BetweenpersonSub$LPA$To == "LPA"))
  expect_true(all(x$BetweenpersonSub$SB$To == "SB"))
  
  expect_true(all(x$BetweenpersonSub$TST$Level == "between"))
  expect_true(all(x$BetweenpersonSub$WAKE$Level == "between"))
  expect_true(all(x$BetweenpersonSub$MVPA$Level == "between"))
  expect_true(all(x$BetweenpersonSub$LPA$Level == "between"))
  expect_true(all(x$BetweenpersonSub$SB$Level == "between"))
  
  expect_true(all(x$BetweenpersonSub$TST$EffectType == "conditional"))
  expect_true(all(x$BetweenpersonSub$WAKE$EffectType == "conditional"))
  expect_true(all(x$BetweenpersonSub$MVPA$EffectType == "conditional"))
  expect_true(all(x$BetweenpersonSub$LPA$EffectType == "conditional"))
  expect_true(all(x$BetweenpersonSub$SB$EffectType == "conditional"))
  
  ## wsub
  expect_type(x$WithinpersonSub, "list")
  expect_equal(length(x$WithinpersonSub), length(m$CompILR$parts))
  expect_s3_class(x$WithinpersonSub$TST, "data.table")
  expect_s3_class(x$WithinpersonSub$WAKE, "data.table")
  expect_s3_class(x$WithinpersonSub$MVPA, "data.table")
  expect_s3_class(x$WithinpersonSub$LPA, "data.table")
  expect_s3_class(x$WithinpersonSub$SB, "data.table")
  
  expect_type(x$WithinpersonSub$TST$Mean, "double")
  expect_type(x$WithinpersonSub$TST$CI_low, "double")
  expect_type(x$WithinpersonSub$TST$CI_high, "double")
  expect_type(x$WithinpersonSub$TST$Delta, "double")
  expect_type(x$WithinpersonSub$TST$From, "character")
  expect_type(x$WithinpersonSub$TST$To, "character")
  
  expect_type(x$WithinpersonSub$WAKE$Mean, "double")
  expect_type(x$WithinpersonSub$WAKE$CI_low, "double")
  expect_type(x$WithinpersonSub$WAKE$CI_high, "double")
  expect_type(x$WithinpersonSub$WAKE$Delta, "double")
  expect_type(x$WithinpersonSub$WAKE$From, "character")
  expect_type(x$WithinpersonSub$WAKE$To, "character")
  
  expect_type(x$WithinpersonSub$MVPA$Mean, "double")
  expect_type(x$WithinpersonSub$MVPA$CI_low, "double")
  expect_type(x$WithinpersonSub$MVPA$CI_high, "double")
  expect_type(x$WithinpersonSub$MVPA$Delta, "double")
  expect_type(x$WithinpersonSub$MVPA$From, "character")
  expect_type(x$WithinpersonSub$MVPA$To, "character")
  
  expect_type(x$WithinpersonSub$LPA$Mean, "double")
  expect_type(x$WithinpersonSub$LPA$CI_low, "double")
  expect_type(x$WithinpersonSub$LPA$CI_high, "double")
  expect_type(x$WithinpersonSub$LPA$Delta, "double")
  expect_type(x$WithinpersonSub$LPA$From, "character")
  expect_type(x$WithinpersonSub$LPA$To, "character")
  
  expect_type(x$WithinpersonSub$SB$Mean, "double")
  expect_type(x$WithinpersonSub$SB$CI_low, "double")
  expect_type(x$WithinpersonSub$SB$CI_high, "double")
  expect_type(x$WithinpersonSub$SB$Delta, "double")
  expect_type(x$WithinpersonSub$SB$From, "character")
  expect_type(x$WithinpersonSub$SB$To, "character")
  
  expect_true(ncol(x$WithinpersonSub$TST) >= 8)
  expect_true(ncol(x$WithinpersonSub$WAKE) >= 8)
  expect_true(ncol(x$WithinpersonSub$MVPA) >= 8)
  expect_true(ncol(x$WithinpersonSub$LPA) >= 8)
  expect_true(ncol(x$WithinpersonSub$SB) >= 8)
  
  expect_true(all(x$WithinpersonSub$TST$To == "TST"))
  expect_true(all(x$WithinpersonSub$WAKE$To == "WAKE"))
  expect_true(all(x$WithinpersonSub$MVPA$To == "MVPA"))
  expect_true(all(x$WithinpersonSub$LPA$To == "LPA"))
  expect_true(all(x$WithinpersonSub$SB$To == "SB"))
  
  expect_true(all(x$WithinpersonSub$TST$Level == "within"))
  expect_true(all(x$WithinpersonSub$WAKE$Level == "within"))
  expect_true(all(x$WithinpersonSub$MVPA$Level == "within"))
  expect_true(all(x$WithinpersonSub$LPA$Level == "within"))
  expect_true(all(x$WithinpersonSub$SB$Level == "within"))
  
  expect_true(all(x$WithinpersonSub$TST$EffectType == "conditional"))
  expect_true(all(x$WithinpersonSub$WAKE$EffectType == "conditional"))
  expect_true(all(x$WithinpersonSub$MVPA$EffectType == "conditional"))
  expect_true(all(x$WithinpersonSub$LPA$EffectType == "conditional"))
  expect_true(all(x$WithinpersonSub$SB$EffectType == "conditional"))
  
  ## bsubmargins 
  expect_type(x$BetweenpersonSubMargins, "list")
  expect_equal(length(x$BetweenpersonSubMargins), length(m$CompILR$parts))
  expect_s3_class(x$BetweenpersonSubMargins$TST, "data.table")
  expect_s3_class(x$BetweenpersonSubMargins$WAKE, "data.table")
  expect_s3_class(x$BetweenpersonSubMargins$MVPA, "data.table")
  expect_s3_class(x$BetweenpersonSubMargins$LPA, "data.table")
  expect_s3_class(x$BetweenpersonSubMargins$SB, "data.table")
  
  expect_type(x$BetweenpersonSubMargins$TST$Mean, "double")
  expect_type(x$BetweenpersonSubMargins$TST$CI_low, "double")
  expect_type(x$BetweenpersonSubMargins$TST$CI_high, "double")
  expect_type(x$BetweenpersonSubMargins$TST$Delta, "double")
  expect_type(x$BetweenpersonSubMargins$TST$From, "character")
  expect_type(x$BetweenpersonSubMargins$TST$To, "character")
  
  expect_type(x$BetweenpersonSubMargins$WAKE$Mean, "double")
  expect_type(x$BetweenpersonSubMargins$WAKE$CI_low, "double")
  expect_type(x$BetweenpersonSubMargins$WAKE$CI_high, "double")
  expect_type(x$BetweenpersonSubMargins$WAKE$Delta, "double")
  expect_type(x$BetweenpersonSubMargins$WAKE$From, "character")
  expect_type(x$BetweenpersonSubMargins$WAKE$To, "character")
  
  expect_type(x$BetweenpersonSubMargins$MVPA$Mean, "double")
  expect_type(x$BetweenpersonSubMargins$MVPA$CI_low, "double")
  expect_type(x$BetweenpersonSubMargins$MVPA$CI_high, "double")
  expect_type(x$BetweenpersonSubMargins$MVPA$Delta, "double")
  expect_type(x$BetweenpersonSubMargins$MVPA$From, "character")
  expect_type(x$BetweenpersonSubMargins$MVPA$To, "character")
  
  expect_type(x$BetweenpersonSubMargins$LPA$Mean, "double")
  expect_type(x$BetweenpersonSubMargins$LPA$CI_low, "double")
  expect_type(x$BetweenpersonSubMargins$LPA$CI_high, "double")
  expect_type(x$BetweenpersonSubMargins$LPA$Delta, "double")
  expect_type(x$BetweenpersonSubMargins$LPA$From, "character")
  expect_type(x$BetweenpersonSubMargins$LPA$To, "character")
  
  expect_type(x$BetweenpersonSubMargins$SB$Mean, "double")
  expect_type(x$BetweenpersonSubMargins$SB$CI_low, "double")
  expect_type(x$BetweenpersonSubMargins$SB$CI_high, "double")
  expect_type(x$BetweenpersonSubMargins$SB$Delta, "double")
  expect_type(x$BetweenpersonSubMargins$SB$From, "character")
  expect_type(x$BetweenpersonSubMargins$SB$To, "character")
  
  expect_true(ncol(x$BetweenpersonSubMargins$TST) == 8)
  expect_true(ncol(x$BetweenpersonSubMargins$WAKE) == 8)
  expect_true(ncol(x$BetweenpersonSubMargins$MVPA) == 8)
  expect_true(ncol(x$BetweenpersonSubMargins$LPA) == 8)
  expect_true(ncol(x$BetweenpersonSubMargins$SB) == 8)
  
  expect_true(all(x$BetweenpersonSubMargins$TST$To == "TST"))
  expect_true(all(x$BetweenpersonSubMargins$WAKE$To == "WAKE"))
  expect_true(all(x$BetweenpersonSubMargins$MVPA$To == "MVPA"))
  expect_true(all(x$BetweenpersonSubMargins$LPA$To == "LPA"))
  expect_true(all(x$BetweenpersonSubMargins$SB$To == "SB"))
  
  expect_true(all(x$BetweenpersonSubMargins$TST$Level == "between"))
  expect_true(all(x$BetweenpersonSubMargins$WAKE$Level == "between"))
  expect_true(all(x$BetweenpersonSubMargins$MVPA$Level == "between"))
  expect_true(all(x$BetweenpersonSubMargins$LPA$Level == "between"))
  expect_true(all(x$BetweenpersonSubMargins$SB$Level == "between"))
  
  expect_true(all(x$BetweenpersonSubMargins$TST$EffectType == "marginal"))
  expect_true(all(x$BetweenpersonSubMargins$WAKE$EffectType == "marginal"))
  expect_true(all(x$BetweenpersonSubMargins$MVPA$EffectType == "marginal"))
  expect_true(all(x$BetweenpersonSubMargins$LPA$EffectType == "marginal"))
  expect_true(all(x$BetweenpersonSubMargins$SB$EffectType == "marginal"))
  
  ## wsubmargins
  expect_type(x$WithinpersonSubMargins, "list")
  expect_equal(length(x$WithinpersonSubMargins), length(m$CompILR$parts))
  expect_s3_class(x$WithinpersonSubMargins$TST, "data.table")
  expect_s3_class(x$WithinpersonSubMargins$WAKE, "data.table")
  expect_s3_class(x$WithinpersonSubMargins$MVPA, "data.table")
  expect_s3_class(x$WithinpersonSubMargins$LPA, "data.table")
  expect_s3_class(x$WithinpersonSubMargins$SB, "data.table")
  
  expect_type(x$WithinpersonSubMargins$TST$Mean, "double")
  expect_type(x$WithinpersonSubMargins$TST$CI_low, "double")
  expect_type(x$WithinpersonSubMargins$TST$CI_high, "double")
  expect_type(x$WithinpersonSubMargins$TST$Delta, "double")
  expect_type(x$WithinpersonSubMargins$TST$From, "character")
  expect_type(x$WithinpersonSubMargins$TST$To, "character")
  
  expect_type(x$WithinpersonSubMargins$WAKE$Mean, "double")
  expect_type(x$WithinpersonSubMargins$WAKE$CI_low, "double")
  expect_type(x$WithinpersonSubMargins$WAKE$CI_high, "double")
  expect_type(x$WithinpersonSubMargins$WAKE$Delta, "double")
  expect_type(x$WithinpersonSubMargins$WAKE$From, "character")
  expect_type(x$WithinpersonSubMargins$WAKE$To, "character")
  
  expect_type(x$WithinpersonSubMargins$MVPA$Mean, "double")
  expect_type(x$WithinpersonSubMargins$MVPA$CI_low, "double")
  expect_type(x$WithinpersonSubMargins$MVPA$CI_high, "double")
  expect_type(x$WithinpersonSubMargins$MVPA$Delta, "double")
  expect_type(x$WithinpersonSubMargins$MVPA$From, "character")
  expect_type(x$WithinpersonSubMargins$MVPA$To, "character")
  
  expect_type(x$WithinpersonSubMargins$LPA$Mean, "double")
  expect_type(x$WithinpersonSubMargins$LPA$CI_low, "double")
  expect_type(x$WithinpersonSubMargins$LPA$CI_high, "double")
  expect_type(x$WithinpersonSubMargins$LPA$Delta, "double")
  expect_type(x$WithinpersonSubMargins$LPA$From, "character")
  expect_type(x$WithinpersonSubMargins$LPA$To, "character")
  
  expect_type(x$WithinpersonSubMargins$SB$Mean, "double")
  expect_type(x$WithinpersonSubMargins$SB$CI_low, "double")
  expect_type(x$WithinpersonSubMargins$SB$CI_high, "double")
  expect_type(x$WithinpersonSubMargins$SB$Delta, "double")
  expect_type(x$WithinpersonSubMargins$SB$From, "character")
  expect_type(x$WithinpersonSubMargins$SB$To, "character")
  
  expect_true(ncol(x$WithinpersonSubMargins$TST) == 8)
  expect_true(ncol(x$WithinpersonSubMargins$WAKE) == 8)
  expect_true(ncol(x$WithinpersonSubMargins$MVPA) == 8)
  expect_true(ncol(x$WithinpersonSubMargins$LPA) == 8)
  expect_true(ncol(x$WithinpersonSubMargins$SB) == 8)
  
  expect_true(all(x$WithinpersonSubMargins$TST$To == "TST"))
  expect_true(all(x$WithinpersonSubMargins$WAKE$To == "WAKE"))
  expect_true(all(x$WithinpersonSubMargins$MVPA$To == "MVPA"))
  expect_true(all(x$WithinpersonSubMargins$LPA$To == "LPA"))
  expect_true(all(x$WithinpersonSubMargins$SB$To == "SB"))
  
  expect_true(all(x$WithinpersonSubMargins$TST$Level == "within"))
  expect_true(all(x$WithinpersonSubMargins$WAKE$Level == "within"))
  expect_true(all(x$WithinpersonSubMargins$MVPA$Level == "within"))
  expect_true(all(x$WithinpersonSubMargins$LPA$Level == "within"))
  expect_true(all(x$WithinpersonSubMargins$SB$Level == "within"))
  
  expect_true(all(x$WithinpersonSubMargins$TST$EffectType == "marginal"))
  expect_true(all(x$WithinpersonSubMargins$WAKE$EffectType == "marginal"))
  expect_true(all(x$WithinpersonSubMargins$MVPA$EffectType == "marginal"))
  expect_true(all(x$WithinpersonSubMargins$LPA$EffectType == "marginal"))
  expect_true(all(x$WithinpersonSubMargins$SB$EffectType == "marginal"))
  
})

