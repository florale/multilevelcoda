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
                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)

suppressWarnings(
  m <- brmcoda(complr = cilr,
               formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                 wz1_1 + wz2_1 + wz3_1 + wz4_1 + Female + (1 | ID),
               chain = 1, iter = 500, seed = 123,
               backend = backend))
foreach::registerDoSEQ()

x <- substitution(object = m, delta = 2)

# Tests
#---------------------------------------------------------------------------------------------------

test_that("substitution errors for invalid input", {
  
  ## missing object
  expect_error(x <- substitution(delta = 2))
  
  ## not brmcoda model
  m1 <- lmer(Stress ~ 1 + (1 | ID), data = mcompd)
  expect_error(x <- substitution(object = m1, delta = 2))
  
  ## incorrect delta
  expect_error(x <- substitution(object = m, delta = -10))

  ## missing delta
  expect_error(x <- substitution(object = m))
  
  # ## reference grid has matching names with ILRs
  # rg <- data.table(bz1_1 = 1)
  # expect_error(x <- substitution(object = m, delta = 2, regrid = rg))
  
  ## base does not have the same components as parts in cilr
  ps <- build.base(c("WAKE", "MVPA", "LPA", "SB"))
  expect_error(x <- substitution(object = m, base = ps, delta = 2))
  
  ## base does not have the same names as parts in cilr
  ps <- build.base(parts = c("Sleep", "WAKE", "MVPA", "LIPA", "SB"))
  expect_error(x <- substitution(object = m, base = ps, delta = 2))
  
})

# test_that("substitution gives expected outputs", {
#   
#   ## bsub
#   expect_type(x$between_simple_sub, "list")
#   expect_equal(length(x$between_simple_sub), length(m$complr$parts))
#   expect_s3_class(x$between_simple_sub$TST, "data.table")
#   expect_s3_class(x$between_simple_sub$WAKE, "data.table")
#   expect_s3_class(x$between_simple_sub$MVPA, "data.table")
#   expect_s3_class(x$between_simple_sub$LPA, "data.table")
#   expect_s3_class(x$between_simple_sub$SB, "data.table")
#   
#   expect_type(x$between_simple_sub$TST$Mean, "double")
#   expect_type(x$between_simple_sub$TST$CI_low, "double")
#   expect_type(x$between_simple_sub$TST$CI_high, "double")
#   expect_type(x$between_simple_sub$TST$Delta, "double")
#   expect_type(x$between_simple_sub$TST$From, "character")
#   expect_type(x$between_simple_sub$TST$To, "character")
#   
#   expect_type(x$between_simple_sub$WAKE$Mean, "double")
#   expect_type(x$between_simple_sub$WAKE$CI_low, "double")
#   expect_type(x$between_simple_sub$WAKE$CI_high, "double")
#   expect_type(x$between_simple_sub$WAKE$Delta, "double")
#   expect_type(x$between_simple_sub$WAKE$From, "character")
#   expect_type(x$between_simple_sub$WAKE$To, "character")
#   
#   expect_type(x$between_simple_sub$MVPA$Mean, "double")
#   expect_type(x$between_simple_sub$MVPA$CI_low, "double")
#   expect_type(x$between_simple_sub$MVPA$CI_high, "double")
#   expect_type(x$between_simple_sub$MVPA$Delta, "double")
#   expect_type(x$between_simple_sub$MVPA$From, "character")
#   expect_type(x$between_simple_sub$MVPA$To, "character")
#   
#   expect_type(x$between_simple_sub$LPA$Mean, "double")
#   expect_type(x$between_simple_sub$LPA$CI_low, "double")
#   expect_type(x$between_simple_sub$LPA$CI_high, "double")
#   expect_type(x$between_simple_sub$LPA$Delta, "double")
#   expect_type(x$between_simple_sub$LPA$From, "character")
#   expect_type(x$between_simple_sub$LPA$To, "character")
#   
#   expect_type(x$between_simple_sub$SB$Mean, "double")
#   expect_type(x$between_simple_sub$SB$CI_low, "double")
#   expect_type(x$between_simple_sub$SB$CI_high, "double")
#   expect_type(x$between_simple_sub$SB$Delta, "double")
#   expect_type(x$between_simple_sub$SB$From, "character")
#   expect_type(x$between_simple_sub$SB$To, "character")
#   
#   expect_true(ncol(x$between_simple_sub$TST) >= 8)
#   expect_true(ncol(x$between_simple_sub$WAKE) >= 8)
#   expect_true(ncol(x$between_simple_sub$MVPA) >= 8)
#   expect_true(ncol(x$between_simple_sub$LPA) >= 8)
#   expect_true(ncol(x$between_simple_sub$SB) >= 8)
#   
#   expect_true(all(x$between_simple_sub$TST$To == "TST"))
#   expect_true(all(x$between_simple_sub$WAKE$To == "WAKE"))
#   expect_true(all(x$between_simple_sub$MVPA$To == "MVPA"))
#   expect_true(all(x$between_simple_sub$LPA$To == "LPA"))
#   expect_true(all(x$between_simple_sub$SB$To == "SB"))
#   
#   expect_true(all(x$between_simple_sub$TST$Level == "between"))
#   expect_true(all(x$between_simple_sub$WAKE$Level == "between"))
#   expect_true(all(x$between_simple_sub$MVPA$Level == "between"))
#   expect_true(all(x$between_simple_sub$LPA$Level == "between"))
#   expect_true(all(x$between_simple_sub$SB$Level == "between"))
#   
#   ## wsub
#   expect_type(x$within_simple_sub, "list")
#   expect_equal(length(x$within_simple_sub), length(m$complr$parts))
#   expect_s3_class(x$within_simple_sub$TST, "data.table")
#   expect_s3_class(x$within_simple_sub$WAKE, "data.table")
#   expect_s3_class(x$within_simple_sub$MVPA, "data.table")
#   expect_s3_class(x$within_simple_sub$LPA, "data.table")
#   expect_s3_class(x$within_simple_sub$SB, "data.table")
#   
#   expect_type(x$within_simple_sub$TST$Mean, "double")
#   expect_type(x$within_simple_sub$TST$CI_low, "double")
#   expect_type(x$within_simple_sub$TST$CI_high, "double")
#   expect_type(x$within_simple_sub$TST$Delta, "double")
#   expect_type(x$within_simple_sub$TST$From, "character")
#   expect_type(x$within_simple_sub$TST$To, "character")
#   
#   expect_type(x$within_simple_sub$WAKE$Mean, "double")
#   expect_type(x$within_simple_sub$WAKE$CI_low, "double")
#   expect_type(x$within_simple_sub$WAKE$CI_high, "double")
#   expect_type(x$within_simple_sub$WAKE$Delta, "double")
#   expect_type(x$within_simple_sub$WAKE$From, "character")
#   expect_type(x$within_simple_sub$WAKE$To, "character")
#   
#   expect_type(x$within_simple_sub$MVPA$Mean, "double")
#   expect_type(x$within_simple_sub$MVPA$CI_low, "double")
#   expect_type(x$within_simple_sub$MVPA$CI_high, "double")
#   expect_type(x$within_simple_sub$MVPA$Delta, "double")
#   expect_type(x$within_simple_sub$MVPA$From, "character")
#   expect_type(x$within_simple_sub$MVPA$To, "character")
#   
#   expect_type(x$within_simple_sub$LPA$Mean, "double")
#   expect_type(x$within_simple_sub$LPA$CI_low, "double")
#   expect_type(x$within_simple_sub$LPA$CI_high, "double")
#   expect_type(x$within_simple_sub$LPA$Delta, "double")
#   expect_type(x$within_simple_sub$LPA$From, "character")
#   expect_type(x$within_simple_sub$LPA$To, "character")
#   
#   expect_type(x$within_simple_sub$SB$Mean, "double")
#   expect_type(x$within_simple_sub$SB$CI_low, "double")
#   expect_type(x$within_simple_sub$SB$CI_high, "double")
#   expect_type(x$within_simple_sub$SB$Delta, "double")
#   expect_type(x$within_simple_sub$SB$From, "character")
#   expect_type(x$within_simple_sub$SB$To, "character")
#   
#   expect_true(ncol(x$within_simple_sub$TST) >= 8)
#   expect_true(ncol(x$within_simple_sub$WAKE) >= 8)
#   expect_true(ncol(x$within_simple_sub$MVPA) >= 8)
#   expect_true(ncol(x$within_simple_sub$LPA) >= 8)
#   expect_true(ncol(x$within_simple_sub$SB) >= 8)
#   
#   expect_true(all(x$within_simple_sub$TST$To == "TST"))
#   expect_true(all(x$within_simple_sub$WAKE$To == "WAKE"))
#   expect_true(all(x$within_simple_sub$MVPA$To == "MVPA"))
#   expect_true(all(x$within_simple_sub$LPA$To == "LPA"))
#   expect_true(all(x$within_simple_sub$SB$To == "SB"))
#   
#   expect_true(all(x$within_simple_sub$TST$Level == "within"))
#   expect_true(all(x$within_simple_sub$WAKE$Level == "within"))
#   expect_true(all(x$within_simple_sub$MVPA$Level == "within"))
#   expect_true(all(x$within_simple_sub$LPA$Level == "within"))
#   expect_true(all(x$within_simple_sub$SB$Level == "within"))
#   
#   ## bsubmargins 
#   expect_type(x$between_avg_sub, "list")
#   expect_equal(length(x$between_avg_sub), length(m$complr$parts))
#   expect_s3_class(x$between_avg_sub$TST, "data.table")
#   expect_s3_class(x$between_avg_sub$WAKE, "data.table")
#   expect_s3_class(x$between_avg_sub$MVPA, "data.table")
#   expect_s3_class(x$between_avg_sub$LPA, "data.table")
#   expect_s3_class(x$between_avg_sub$SB, "data.table")
#   
#   expect_type(x$between_avg_sub$TST$Mean, "double")
#   expect_type(x$between_avg_sub$TST$CI_low, "double")
#   expect_type(x$between_avg_sub$TST$CI_high, "double")
#   expect_type(x$between_avg_sub$TST$Delta, "double")
#   expect_type(x$between_avg_sub$TST$From, "character")
#   expect_type(x$between_avg_sub$TST$To, "character")
#   
#   expect_type(x$between_avg_sub$WAKE$Mean, "double")
#   expect_type(x$between_avg_sub$WAKE$CI_low, "double")
#   expect_type(x$between_avg_sub$WAKE$CI_high, "double")
#   expect_type(x$between_avg_sub$WAKE$Delta, "double")
#   expect_type(x$between_avg_sub$WAKE$From, "character")
#   expect_type(x$between_avg_sub$WAKE$To, "character")
#   
#   expect_type(x$between_avg_sub$MVPA$Mean, "double")
#   expect_type(x$between_avg_sub$MVPA$CI_low, "double")
#   expect_type(x$between_avg_sub$MVPA$CI_high, "double")
#   expect_type(x$between_avg_sub$MVPA$Delta, "double")
#   expect_type(x$between_avg_sub$MVPA$From, "character")
#   expect_type(x$between_avg_sub$MVPA$To, "character")
#   
#   expect_type(x$between_avg_sub$LPA$Mean, "double")
#   expect_type(x$between_avg_sub$LPA$CI_low, "double")
#   expect_type(x$between_avg_sub$LPA$CI_high, "double")
#   expect_type(x$between_avg_sub$LPA$Delta, "double")
#   expect_type(x$between_avg_sub$LPA$From, "character")
#   expect_type(x$between_avg_sub$LPA$To, "character")
#   
#   expect_type(x$between_avg_sub$SB$Mean, "double")
#   expect_type(x$between_avg_sub$SB$CI_low, "double")
#   expect_type(x$between_avg_sub$SB$CI_high, "double")
#   expect_type(x$between_avg_sub$SB$Delta, "double")
#   expect_type(x$between_avg_sub$SB$From, "character")
#   expect_type(x$between_avg_sub$SB$To, "character")
#   
#   expect_true(ncol(x$between_avg_sub$TST) == 8)
#   expect_true(ncol(x$between_avg_sub$WAKE) == 8)
#   expect_true(ncol(x$between_avg_sub$MVPA) == 8)
#   expect_true(ncol(x$between_avg_sub$LPA) == 8)
#   expect_true(ncol(x$between_avg_sub$SB) == 8)
#   
#   expect_true(all(x$between_avg_sub$TST$To == "TST"))
#   expect_true(all(x$between_avg_sub$WAKE$To == "WAKE"))
#   expect_true(all(x$between_avg_sub$MVPA$To == "MVPA"))
#   expect_true(all(x$between_avg_sub$LPA$To == "LPA"))
#   expect_true(all(x$between_avg_sub$SB$To == "SB"))
#   
#   expect_true(all(x$between_avg_sub$TST$Level == "between"))
#   expect_true(all(x$between_avg_sub$WAKE$Level == "between"))
#   expect_true(all(x$between_avg_sub$MVPA$Level == "between"))
#   expect_true(all(x$between_avg_sub$LPA$Level == "between"))
#   expect_true(all(x$between_avg_sub$SB$Level == "between"))
# 
#   ## wsubmargins
#   expect_type(x$within_avg_sub, "list")
#   expect_equal(length(x$within_avg_sub), length(m$complr$parts))
#   expect_s3_class(x$within_avg_sub$TST, "data.table")
#   expect_s3_class(x$within_avg_sub$WAKE, "data.table")
#   expect_s3_class(x$within_avg_sub$MVPA, "data.table")
#   expect_s3_class(x$within_avg_sub$LPA, "data.table")
#   expect_s3_class(x$within_avg_sub$SB, "data.table")
#   
#   expect_type(x$within_avg_sub$TST$Mean, "double")
#   expect_type(x$within_avg_sub$TST$CI_low, "double")
#   expect_type(x$within_avg_sub$TST$CI_high, "double")
#   expect_type(x$within_avg_sub$TST$Delta, "double")
#   expect_type(x$within_avg_sub$TST$From, "character")
#   expect_type(x$within_avg_sub$TST$To, "character")
#   
#   expect_type(x$within_avg_sub$WAKE$Mean, "double")
#   expect_type(x$within_avg_sub$WAKE$CI_low, "double")
#   expect_type(x$within_avg_sub$WAKE$CI_high, "double")
#   expect_type(x$within_avg_sub$WAKE$Delta, "double")
#   expect_type(x$within_avg_sub$WAKE$From, "character")
#   expect_type(x$within_avg_sub$WAKE$To, "character")
#   
#   expect_type(x$within_avg_sub$MVPA$Mean, "double")
#   expect_type(x$within_avg_sub$MVPA$CI_low, "double")
#   expect_type(x$within_avg_sub$MVPA$CI_high, "double")
#   expect_type(x$within_avg_sub$MVPA$Delta, "double")
#   expect_type(x$within_avg_sub$MVPA$From, "character")
#   expect_type(x$within_avg_sub$MVPA$To, "character")
#   
#   expect_type(x$within_avg_sub$LPA$Mean, "double")
#   expect_type(x$within_avg_sub$LPA$CI_low, "double")
#   expect_type(x$within_avg_sub$LPA$CI_high, "double")
#   expect_type(x$within_avg_sub$LPA$Delta, "double")
#   expect_type(x$within_avg_sub$LPA$From, "character")
#   expect_type(x$within_avg_sub$LPA$To, "character")
#   
#   expect_type(x$within_avg_sub$SB$Mean, "double")
#   expect_type(x$within_avg_sub$SB$CI_low, "double")
#   expect_type(x$within_avg_sub$SB$CI_high, "double")
#   expect_type(x$within_avg_sub$SB$Delta, "double")
#   expect_type(x$within_avg_sub$SB$From, "character")
#   expect_type(x$within_avg_sub$SB$To, "character")
#   
#   expect_true(ncol(x$within_avg_sub$TST) == 8)
#   expect_true(ncol(x$within_avg_sub$WAKE) == 8)
#   expect_true(ncol(x$within_avg_sub$MVPA) == 8)
#   expect_true(ncol(x$within_avg_sub$LPA) == 8)
#   expect_true(ncol(x$within_avg_sub$SB) == 8)
#   
#   expect_true(all(x$within_avg_sub$TST$To == "TST"))
#   expect_true(all(x$within_avg_sub$WAKE$To == "WAKE"))
#   expect_true(all(x$within_avg_sub$MVPA$To == "MVPA"))
#   expect_true(all(x$within_avg_sub$LPA$To == "LPA"))
#   expect_true(all(x$within_avg_sub$SB$To == "SB"))
#   
#   expect_true(all(x$within_avg_sub$TST$Level == "within"))
#   expect_true(all(x$within_avg_sub$WAKE$Level == "within"))
#   expect_true(all(x$within_avg_sub$MVPA$Level == "within"))
#   expect_true(all(x$within_avg_sub$LPA$Level == "within"))
#   expect_true(all(x$within_avg_sub$SB$Level == "within"))
#   
# })

