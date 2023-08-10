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
                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)

suppressWarnings(
  m <- brmcoda(compilr = cilr,
               formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
               chain = 1, iter = 500, seed = 123,
               backend = backend))
foreach::registerDoSEQ()

x <- substitution(object = m, delta = 2)

# Tests
#---------------------------------------------------------------------------------------------------

test_that("substitution errors for invalid input", {
  
  ## missing object
  expect_error(x <- substitution(basesub = psub, delta = 2))
  
  ## not brmcoda model
  m1 <- lmer(Stress ~ 1 + (1 | ID), data = mcompd)
  expect_error(x <- substitution(object = m1, basesub = psub, delta = 2))
  
  ## incorrect delta
  expect_error(x <- substitution(object = m, basesub = psub, delta = -10))

  ## missing delta
  expect_error(x <- substitution(object = m, basesub = psub))
  
  # ## reference grid has matching names with ILRs
  # rg <- data.table(bilr1 = 1)
  # expect_error(x <- substitution(object = m, basesub = psub, delta = 2, regrid = rg))
  
  ## basesub does not have the same components as parts in cilr
  ps <- basesub(c("WAKE", "MVPA", "LPA", "SB"))
  expect_error(x <- substitution(object = m, basesub = ps, delta = 2))
  
  ## basesub does not have the same names as parts in cilr
  ps <- basesub(parts = c("Sleep", "WAKE", "MVPA", "LIPA", "SB"))
  expect_error(x <- substitution(object = m, basesub = ps, delta = 2))
  
})

test_that("substitution gives expected outputs", {
  
  ## bsub
  expect_type(x$BetweenSub, "list")
  expect_equal(length(x$BetweenSub), length(m$CompILR$parts))
  expect_s3_class(x$BetweenSub$TST, "data.table")
  expect_s3_class(x$BetweenSub$WAKE, "data.table")
  expect_s3_class(x$BetweenSub$MVPA, "data.table")
  expect_s3_class(x$BetweenSub$LPA, "data.table")
  expect_s3_class(x$BetweenSub$SB, "data.table")
  
  expect_type(x$BetweenSub$TST$Mean, "double")
  expect_type(x$BetweenSub$TST$CI_low, "double")
  expect_type(x$BetweenSub$TST$CI_high, "double")
  expect_type(x$BetweenSub$TST$Delta, "double")
  expect_type(x$BetweenSub$TST$From, "character")
  expect_type(x$BetweenSub$TST$To, "character")
  
  expect_type(x$BetweenSub$WAKE$Mean, "double")
  expect_type(x$BetweenSub$WAKE$CI_low, "double")
  expect_type(x$BetweenSub$WAKE$CI_high, "double")
  expect_type(x$BetweenSub$WAKE$Delta, "double")
  expect_type(x$BetweenSub$WAKE$From, "character")
  expect_type(x$BetweenSub$WAKE$To, "character")
  
  expect_type(x$BetweenSub$MVPA$Mean, "double")
  expect_type(x$BetweenSub$MVPA$CI_low, "double")
  expect_type(x$BetweenSub$MVPA$CI_high, "double")
  expect_type(x$BetweenSub$MVPA$Delta, "double")
  expect_type(x$BetweenSub$MVPA$From, "character")
  expect_type(x$BetweenSub$MVPA$To, "character")
  
  expect_type(x$BetweenSub$LPA$Mean, "double")
  expect_type(x$BetweenSub$LPA$CI_low, "double")
  expect_type(x$BetweenSub$LPA$CI_high, "double")
  expect_type(x$BetweenSub$LPA$Delta, "double")
  expect_type(x$BetweenSub$LPA$From, "character")
  expect_type(x$BetweenSub$LPA$To, "character")
  
  expect_type(x$BetweenSub$SB$Mean, "double")
  expect_type(x$BetweenSub$SB$CI_low, "double")
  expect_type(x$BetweenSub$SB$CI_high, "double")
  expect_type(x$BetweenSub$SB$Delta, "double")
  expect_type(x$BetweenSub$SB$From, "character")
  expect_type(x$BetweenSub$SB$To, "character")
  
  expect_true(ncol(x$BetweenSub$TST) >= 8)
  expect_true(ncol(x$BetweenSub$WAKE) >= 8)
  expect_true(ncol(x$BetweenSub$MVPA) >= 8)
  expect_true(ncol(x$BetweenSub$LPA) >= 8)
  expect_true(ncol(x$BetweenSub$SB) >= 8)
  
  expect_true(all(x$BetweenSub$TST$To == "TST"))
  expect_true(all(x$BetweenSub$WAKE$To == "WAKE"))
  expect_true(all(x$BetweenSub$MVPA$To == "MVPA"))
  expect_true(all(x$BetweenSub$LPA$To == "LPA"))
  expect_true(all(x$BetweenSub$SB$To == "SB"))
  
  expect_true(all(x$BetweenSub$TST$Level == "between"))
  expect_true(all(x$BetweenSub$WAKE$Level == "between"))
  expect_true(all(x$BetweenSub$MVPA$Level == "between"))
  expect_true(all(x$BetweenSub$LPA$Level == "between"))
  expect_true(all(x$BetweenSub$SB$Level == "between"))
  
  ## wsub
  expect_type(x$WithinSub, "list")
  expect_equal(length(x$WithinSub), length(m$CompILR$parts))
  expect_s3_class(x$WithinSub$TST, "data.table")
  expect_s3_class(x$WithinSub$WAKE, "data.table")
  expect_s3_class(x$WithinSub$MVPA, "data.table")
  expect_s3_class(x$WithinSub$LPA, "data.table")
  expect_s3_class(x$WithinSub$SB, "data.table")
  
  expect_type(x$WithinSub$TST$Mean, "double")
  expect_type(x$WithinSub$TST$CI_low, "double")
  expect_type(x$WithinSub$TST$CI_high, "double")
  expect_type(x$WithinSub$TST$Delta, "double")
  expect_type(x$WithinSub$TST$From, "character")
  expect_type(x$WithinSub$TST$To, "character")
  
  expect_type(x$WithinSub$WAKE$Mean, "double")
  expect_type(x$WithinSub$WAKE$CI_low, "double")
  expect_type(x$WithinSub$WAKE$CI_high, "double")
  expect_type(x$WithinSub$WAKE$Delta, "double")
  expect_type(x$WithinSub$WAKE$From, "character")
  expect_type(x$WithinSub$WAKE$To, "character")
  
  expect_type(x$WithinSub$MVPA$Mean, "double")
  expect_type(x$WithinSub$MVPA$CI_low, "double")
  expect_type(x$WithinSub$MVPA$CI_high, "double")
  expect_type(x$WithinSub$MVPA$Delta, "double")
  expect_type(x$WithinSub$MVPA$From, "character")
  expect_type(x$WithinSub$MVPA$To, "character")
  
  expect_type(x$WithinSub$LPA$Mean, "double")
  expect_type(x$WithinSub$LPA$CI_low, "double")
  expect_type(x$WithinSub$LPA$CI_high, "double")
  expect_type(x$WithinSub$LPA$Delta, "double")
  expect_type(x$WithinSub$LPA$From, "character")
  expect_type(x$WithinSub$LPA$To, "character")
  
  expect_type(x$WithinSub$SB$Mean, "double")
  expect_type(x$WithinSub$SB$CI_low, "double")
  expect_type(x$WithinSub$SB$CI_high, "double")
  expect_type(x$WithinSub$SB$Delta, "double")
  expect_type(x$WithinSub$SB$From, "character")
  expect_type(x$WithinSub$SB$To, "character")
  
  expect_true(ncol(x$WithinSub$TST) >= 8)
  expect_true(ncol(x$WithinSub$WAKE) >= 8)
  expect_true(ncol(x$WithinSub$MVPA) >= 8)
  expect_true(ncol(x$WithinSub$LPA) >= 8)
  expect_true(ncol(x$WithinSub$SB) >= 8)
  
  expect_true(all(x$WithinSub$TST$To == "TST"))
  expect_true(all(x$WithinSub$WAKE$To == "WAKE"))
  expect_true(all(x$WithinSub$MVPA$To == "MVPA"))
  expect_true(all(x$WithinSub$LPA$To == "LPA"))
  expect_true(all(x$WithinSub$SB$To == "SB"))
  
  expect_true(all(x$WithinSub$TST$Level == "within"))
  expect_true(all(x$WithinSub$WAKE$Level == "within"))
  expect_true(all(x$WithinSub$MVPA$Level == "within"))
  expect_true(all(x$WithinSub$LPA$Level == "within"))
  expect_true(all(x$WithinSub$SB$Level == "within"))
  
  ## bsubmargins 
  expect_type(x$BetweenSubMargins, "list")
  expect_equal(length(x$BetweenSubMargins), length(m$CompILR$parts))
  expect_s3_class(x$BetweenSubMargins$TST, "data.table")
  expect_s3_class(x$BetweenSubMargins$WAKE, "data.table")
  expect_s3_class(x$BetweenSubMargins$MVPA, "data.table")
  expect_s3_class(x$BetweenSubMargins$LPA, "data.table")
  expect_s3_class(x$BetweenSubMargins$SB, "data.table")
  
  expect_type(x$BetweenSubMargins$TST$Mean, "double")
  expect_type(x$BetweenSubMargins$TST$CI_low, "double")
  expect_type(x$BetweenSubMargins$TST$CI_high, "double")
  expect_type(x$BetweenSubMargins$TST$Delta, "double")
  expect_type(x$BetweenSubMargins$TST$From, "character")
  expect_type(x$BetweenSubMargins$TST$To, "character")
  
  expect_type(x$BetweenSubMargins$WAKE$Mean, "double")
  expect_type(x$BetweenSubMargins$WAKE$CI_low, "double")
  expect_type(x$BetweenSubMargins$WAKE$CI_high, "double")
  expect_type(x$BetweenSubMargins$WAKE$Delta, "double")
  expect_type(x$BetweenSubMargins$WAKE$From, "character")
  expect_type(x$BetweenSubMargins$WAKE$To, "character")
  
  expect_type(x$BetweenSubMargins$MVPA$Mean, "double")
  expect_type(x$BetweenSubMargins$MVPA$CI_low, "double")
  expect_type(x$BetweenSubMargins$MVPA$CI_high, "double")
  expect_type(x$BetweenSubMargins$MVPA$Delta, "double")
  expect_type(x$BetweenSubMargins$MVPA$From, "character")
  expect_type(x$BetweenSubMargins$MVPA$To, "character")
  
  expect_type(x$BetweenSubMargins$LPA$Mean, "double")
  expect_type(x$BetweenSubMargins$LPA$CI_low, "double")
  expect_type(x$BetweenSubMargins$LPA$CI_high, "double")
  expect_type(x$BetweenSubMargins$LPA$Delta, "double")
  expect_type(x$BetweenSubMargins$LPA$From, "character")
  expect_type(x$BetweenSubMargins$LPA$To, "character")
  
  expect_type(x$BetweenSubMargins$SB$Mean, "double")
  expect_type(x$BetweenSubMargins$SB$CI_low, "double")
  expect_type(x$BetweenSubMargins$SB$CI_high, "double")
  expect_type(x$BetweenSubMargins$SB$Delta, "double")
  expect_type(x$BetweenSubMargins$SB$From, "character")
  expect_type(x$BetweenSubMargins$SB$To, "character")
  
  expect_true(ncol(x$BetweenSubMargins$TST) == 8)
  expect_true(ncol(x$BetweenSubMargins$WAKE) == 8)
  expect_true(ncol(x$BetweenSubMargins$MVPA) == 8)
  expect_true(ncol(x$BetweenSubMargins$LPA) == 8)
  expect_true(ncol(x$BetweenSubMargins$SB) == 8)
  
  expect_true(all(x$BetweenSubMargins$TST$To == "TST"))
  expect_true(all(x$BetweenSubMargins$WAKE$To == "WAKE"))
  expect_true(all(x$BetweenSubMargins$MVPA$To == "MVPA"))
  expect_true(all(x$BetweenSubMargins$LPA$To == "LPA"))
  expect_true(all(x$BetweenSubMargins$SB$To == "SB"))
  
  expect_true(all(x$BetweenSubMargins$TST$Level == "between"))
  expect_true(all(x$BetweenSubMargins$WAKE$Level == "between"))
  expect_true(all(x$BetweenSubMargins$MVPA$Level == "between"))
  expect_true(all(x$BetweenSubMargins$LPA$Level == "between"))
  expect_true(all(x$BetweenSubMargins$SB$Level == "between"))

  ## wsubmargins
  expect_type(x$WithinSubMargins, "list")
  expect_equal(length(x$WithinSubMargins), length(m$CompILR$parts))
  expect_s3_class(x$WithinSubMargins$TST, "data.table")
  expect_s3_class(x$WithinSubMargins$WAKE, "data.table")
  expect_s3_class(x$WithinSubMargins$MVPA, "data.table")
  expect_s3_class(x$WithinSubMargins$LPA, "data.table")
  expect_s3_class(x$WithinSubMargins$SB, "data.table")
  
  expect_type(x$WithinSubMargins$TST$Mean, "double")
  expect_type(x$WithinSubMargins$TST$CI_low, "double")
  expect_type(x$WithinSubMargins$TST$CI_high, "double")
  expect_type(x$WithinSubMargins$TST$Delta, "double")
  expect_type(x$WithinSubMargins$TST$From, "character")
  expect_type(x$WithinSubMargins$TST$To, "character")
  
  expect_type(x$WithinSubMargins$WAKE$Mean, "double")
  expect_type(x$WithinSubMargins$WAKE$CI_low, "double")
  expect_type(x$WithinSubMargins$WAKE$CI_high, "double")
  expect_type(x$WithinSubMargins$WAKE$Delta, "double")
  expect_type(x$WithinSubMargins$WAKE$From, "character")
  expect_type(x$WithinSubMargins$WAKE$To, "character")
  
  expect_type(x$WithinSubMargins$MVPA$Mean, "double")
  expect_type(x$WithinSubMargins$MVPA$CI_low, "double")
  expect_type(x$WithinSubMargins$MVPA$CI_high, "double")
  expect_type(x$WithinSubMargins$MVPA$Delta, "double")
  expect_type(x$WithinSubMargins$MVPA$From, "character")
  expect_type(x$WithinSubMargins$MVPA$To, "character")
  
  expect_type(x$WithinSubMargins$LPA$Mean, "double")
  expect_type(x$WithinSubMargins$LPA$CI_low, "double")
  expect_type(x$WithinSubMargins$LPA$CI_high, "double")
  expect_type(x$WithinSubMargins$LPA$Delta, "double")
  expect_type(x$WithinSubMargins$LPA$From, "character")
  expect_type(x$WithinSubMargins$LPA$To, "character")
  
  expect_type(x$WithinSubMargins$SB$Mean, "double")
  expect_type(x$WithinSubMargins$SB$CI_low, "double")
  expect_type(x$WithinSubMargins$SB$CI_high, "double")
  expect_type(x$WithinSubMargins$SB$Delta, "double")
  expect_type(x$WithinSubMargins$SB$From, "character")
  expect_type(x$WithinSubMargins$SB$To, "character")
  
  expect_true(ncol(x$WithinSubMargins$TST) == 8)
  expect_true(ncol(x$WithinSubMargins$WAKE) == 8)
  expect_true(ncol(x$WithinSubMargins$MVPA) == 8)
  expect_true(ncol(x$WithinSubMargins$LPA) == 8)
  expect_true(ncol(x$WithinSubMargins$SB) == 8)
  
  expect_true(all(x$WithinSubMargins$TST$To == "TST"))
  expect_true(all(x$WithinSubMargins$WAKE$To == "WAKE"))
  expect_true(all(x$WithinSubMargins$MVPA$To == "MVPA"))
  expect_true(all(x$WithinSubMargins$LPA$To == "LPA"))
  expect_true(all(x$WithinSubMargins$SB$To == "SB"))
  
  expect_true(all(x$WithinSubMargins$TST$Level == "within"))
  expect_true(all(x$WithinSubMargins$WAKE$Level == "within"))
  expect_true(all(x$WithinSubMargins$MVPA$Level == "within"))
  expect_true(all(x$WithinSubMargins$LPA$Level == "within"))
  expect_true(all(x$WithinSubMargins$SB$Level == "within"))
  
})

