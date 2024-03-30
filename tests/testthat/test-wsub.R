skip_on_cran()

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  backend <- "rstan"
  ## if using rstan backend, models can crash on Windows and MAC OS
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

cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
               parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)

suppressWarnings(
  m <- brmcoda(complr = cilr,
               formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
               chain = 1, iter = 500, seed = 123,
               backend = backend))
foreach::registerDoSEQ()

x <- wsub(object = m, basesub = psub, delta = 2)

# Testing
#---------------------------------------------------------------------------------------------------

# test_that("wsub errors for invalid input", {
# 
#   ## missing object
#   expect_error(x <- wsub(basesub = psub, delta = 2))
# 
#   ## missing basesub
#   expect_error(x <- wsub(object = m, delta = 2))
# 
#   ## not brmcoda model
#   m1 <- lmer(Stress ~ 1 + (1 | ID), data = mcompd)
#   expect_error(x <- wsub(object = m1, basesub = psub, delta = 2))
# 
#   ## invalid delta
#   expect_error(x <- wsub(object = m, basesub = psub, delta = -10))
#   expect_error(x <- wsub(object = m, basesub = psub, delta = 1:10))
# 
#   ## missing delta
#   expect_error(x <- substitution(object = m1, basesub = psub))
# 
#   ## basesub does not have the same components as parts in cilr
#   ps <- basesub(c("WAKE", "MVPA", "LPA", "SB"))
#   expect_error(x <- wsub(object = m, basesub = ps, delta = 2))
# 
#   ## basesub does have the same names as parts in cilr
#   ps <- basesub(parts = c("Sleep", "WAKE", "MVPA", "LPA", "SB"))
#   expect_error(x <- wsub(object = m, basesub = ps, delta = 2))
# 
# })

# test_that("wsub works as expected for adjusted/unadjusted model", {
#   
#   ## reference grid is provided for unadjusted model
#   suppressWarnings(
#     m2 <- brmcoda(complr = cilr,
#                   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#                   chain = 1, iter = 500, seed = 123,
#                   backend = backend))
#   rg <- data.table(Age = 1)
#   expect_warning(x <- wsub(object = m2, basesub = psub, delta = 2, regrid = rg))
#   
#   ## incorect reference grid 1
#   rg <- data.table(Age = 1)
#   expect_error(x <- wsub(object = m, basesub = psub, delta = 2, regrid = rg))
#   
#   ## reference grid has matching names with ILRs
#   rg <- data.table(bilr1 = 1)
#   expect_error(x <- wsub(object = m, basesub = psub, delta = 2, regrid = rg))
#   
#   ## incorect reference grid 2
#   rg <- data.table(bilr1 = 1, Age = 1)
#   expect_error(x <- wsub(object = m, basesub = psub, delta = 2, regrid = rg))
# 
#   # delta out of range
#   expect_error(x <- wsub(object = m, basesub = psub, delta = 1000))
#   
#   ## function knows to use correct user's specified reference grid
#   rg <- data.table(Female = 1)
#   x3 <- wsub(object = m, basesub = psub, delta = 2, regrid = rg)
#   expect_true(all(x3$TST$Female == 1))
#   expect_true(all(x3$WAKE$Female == 1))
#   expect_true(all(x3$MVPA$Female == 1))
#   expect_true(all(x3$LPA$Female == 1))
#   expect_true(all(x$SB$Female == 1))
#   
#   expect_true(all(x3$TST$Female != 0))
#   expect_true(all(x3$WAKE$Female != 0))
#   expect_true(all(x3$MVPA$Female != 0))
#   expect_true(all(x3$LPA$Female != 0))
#   expect_true(all(x3$SB$Female != 0))
#   
#   ## model with unspecified reference grid works as expected
#   expect_equal(x$TST$Female, NULL)
#   expect_equal(x$WAKE$Female, NULL)
#   expect_equal(x$MVPA$Female, NULL)
#   expect_equal(x$LPA$Female, NULL)
#   expect_equal(x$SB$Female, NULL)
#   
#   ## model with unspecified reference grid works as expected
#   expect_true("Female" %nin% colnames(x$TST))
#   expect_true("Female" %nin% colnames(x$WAKE))
#   expect_true("Female" %nin% colnames(x$MVPA))
#   expect_true("Female" %nin% colnames(x$LPA))
#   expect_true("Female" %nin% colnames(x$SB))
#   
#   ## average across reference grid as default
#   x4 <- wsub(object = m, basesub = psub, delta = 2, summary = TRUE)
#   x5 <- wsub(object = m, basesub = psub, delta = 2)
#   expect_equal(x4, x5)
#   
#   ## keep prediction at each level of refrence grid 
#   cilr <- complr(data = mcompd[ID %in% c(1:5, 185:190), .SD[1:3], by = ID], sbp = sbp,
#                   parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)
#   
#   suppressWarnings(
#     m <- brmcoda(complr = cilr,
#                  formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                    wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                  chain = 1, iter = 500, seed = 123,
#                  backend = backend))
#   
#   x6 <- wsub(object = m, basesub = psub, delta = 2, summary = FALSE)
#   
#   expect_equal(nrow(x6$TST), nrow(x5$TST) * 2)
#   expect_equal(nrow(x6$WAKE), nrow(x5$WAKE) * 2)
#   expect_equal(nrow(x6$MVPA), nrow(x5$MVPA) * 2)
#   expect_equal(nrow(x6$LPA), nrow(x5$LPA) * 2)
#   expect_equal(nrow(x6$SB), nrow(x5$SB) * 2)
#   
#   expect_true("Female" %in% colnames(x6$TST))
#   expect_true("Female" %in% colnames(x6$WAKE))
#   expect_true("Female" %in% colnames(x6$MVPA))
#   expect_true("Female" %in% colnames(x6$LPA))
#   expect_true("Female" %in% colnames(x6$SB))
#   
#   expect_true(all(x6$TST$Female %in% c(0, 1)))
#   expect_true(all(x6$WAKE$Female %in% c(0, 1)))
#   expect_true(all(x6$MVPA$Female %in% c(0, 1)))
#   expect_true(all(x6$LPA$Female %in% c(0, 1)))
#   expect_true(all(x6$SB$Female %in% c(0, 1)))
#   
# })
# 
# test_that("wsub checks for user-specified reference composition", {
#   
#   # incorrect length
#   ref1 <- c(400, 60, 500, 60)
#   expect_error(wsub(object = m, basesub = psub, recomp = ref1, delta = 2))
#   
#   # incorrect class
#   ref2 <- c("400", "100", "500", "200", "200")
#   expect_error(wsub(object = m, basesub = psub, recomp = ref2, delta = 2))
#   
#   # incorrect class
#   ref3 <- c(400, 100, 500, 200, 200)
#   expect_error(x <- wsub(object = m, basesub = psub, recomp = ref3, delta = 2))
#   
#   # values outside of possible range
#   ref4 <- c(100, 100, 900, 100, 240)
#   expect_error(x <- wsub(object = m, basesub = psub, recomp = ref4, delta = 2))
#   
#   # include 0
#   ref5 <- c(100, 200, 900, 0, 240)
#   expect_error(x <- wsub(object = m, basesub = psub, recomp = ref5, delta = 2))
#   
# })

test_that("wsub outputs what expected", {
  
  ## types
  expect_type(x, "list")
  expect_equal(length(x), length(m$complr$parts))
  expect_s3_class(x$TST, "data.table")
  expect_s3_class(x$WAKE, "data.table")
  expect_s3_class(x$MVPA, "data.table")
  expect_s3_class(x$LPA, "data.table")
  expect_s3_class(x$SB, "data.table")
  
  expect_type(x$TST$Mean, "double")
  expect_type(x$TST$CI_low, "double")
  expect_type(x$TST$CI_high, "double")
  expect_type(x$TST$Delta, "double")
  expect_type(x$TST$From, "character")
  expect_type(x$TST$To, "character")
  
  expect_type(x$WAKE$Mean, "double")
  expect_type(x$WAKE$CI_low, "double")
  expect_type(x$WAKE$CI_high, "double")
  expect_type(x$WAKE$Delta, "double")
  expect_type(x$WAKE$From, "character")
  expect_type(x$WAKE$To, "character")
  
  expect_type(x$MVPA$Mean, "double")
  expect_type(x$MVPA$CI_low, "double")
  expect_type(x$MVPA$CI_high, "double")
  expect_type(x$MVPA$Delta, "double")
  expect_type(x$MVPA$From, "character")
  expect_type(x$MVPA$To, "character")
  
  expect_type(x$LPA$Mean, "double")
  expect_type(x$LPA$CI_low, "double")
  expect_type(x$LPA$CI_high, "double")
  expect_type(x$LPA$Delta, "double")
  expect_type(x$LPA$From, "character")
  expect_type(x$LPA$To, "character")
  
  expect_type(x$SB$Mean, "double")
  expect_type(x$SB$CI_low, "double")
  expect_type(x$SB$CI_high, "double")
  expect_type(x$SB$Delta, "double")
  expect_type(x$SB$From, "character")
  expect_type(x$SB$To, "character")
  
  expect_true(ncol(x$TST) >= 8)
  expect_true(ncol(x$WAKE) >= 8)
  expect_true(ncol(x$MVPA) >= 8)
  expect_true(ncol(x$LPA) >= 8)
  expect_true(ncol(x$SB) >= 8)
  
  expect_true(all(x$TST$To == "TST"))
  expect_true(all(x$WAKE$To == "WAKE"))
  expect_true(all(x$MVPA$To == "MVPA"))
  expect_true(all(x$LPA$To == "LPA"))
  expect_true(all(x$SB$To == "SB"))
  
  expect_true(all(x$TST$Level == "within"))
  expect_true(all(x$WAKE$Level == "within"))
  expect_true(all(x$MVPA$Level == "within"))
  expect_true(all(x$LPA$Level == "within"))
  expect_true(all(x$SB$Level == "within"))
  
})

test_that("wsub gives results in sensible range", {
  
  ## difference in outcome
  expect_true(x$TST$Mean %ae% "[-0.5, 0) | (0, 0.5]")
  expect_true(x$WAKE$Mean %ae% "[-0.5, 0) | (0, 0.5]")
  expect_true(x$MVPA$Mean %ae% "[-0.5, 0) | (0, 0.5]")
  expect_true(x$LPA$Mean %ae% "[-0.5, 0) | (0, 0.5]")
  expect_true(x$SB$Mean %ae% "[-0.5, 0) | (0, 0.5]")
  
  expect_true(x$TST$CI_low %ae% "[-1, 0) | (0, 1]")
  expect_true(x$WAKE$CI_low %ae% "[-1, 0) | (0, 1]")
  expect_true(x$MVPA$CI_low %ae% "[-1, 0) | (0, 1]")
  expect_true(x$LPA$CI_low %ae% "[-1, 0) | (0, 1]")
  expect_true(x$SB$CI_low %ae% "[-1, 0) | (0, 1]")
  
  expect_true(x$TST$CI_high %ae% "[-1, 0) | (0, 1]")
  expect_true(x$WAKE$CI_high %ae% "[-1, 0) | (0, 1]")
  expect_true(x$MVPA$CI_high %ae% "[-1, 0) | (0, 1]")
  expect_true(x$LPA$CI_high %ae% "[-1, 0) | (0, 1]")
  expect_true(x$SB$CI_high %ae% "[-1, 0) | (0, 1]")
  
})

test_that("wsub gives results in expected direction and magnitude", {
  
  ## values are opposite sign for opposite substitution
  for (i in seq_along(x)) {
    expect_true(all(x[[i]][, sign(Mean[sign(Delta) == 1]) 
                           %a!=% sign(Mean[sign(Delta) == -1]), by = From]$V1))
  }
  
  ## results for 1 min have smaller magnitude than 2 mins
  for (i in seq_along(x)) {
    expect_true(all(x[[i]][, abs(Mean[abs(Delta) == 1]) 
                           < abs(Mean[abs(Delta) == 2])]))
  }
})

# #---------------------------------------------------------------------------------------------------
# # Check the results from substitution model align with brm model 
# # using results from pairwise substitution
# ## Estimates should be in the direction between pairwise coordinates and  pairwise substitution 
# ## CIs should indicate consistent significance between pairwise coordinates and substitution 
# 
# sbp <- matrix(c(
#   1, 1, -1,-1, -1,
#   1, -1, 0, 0, 0,
#   0, 0, -1, -1, 1,
#   0, 0, 1, -1, 0), ncol = 5, byrow = TRUE)
# 
# cilr <- complr(data = mcompd, sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)
# 
# suppressWarnings(m <- brmcoda(complr = cilr,
#                               formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                               seed = 123))
# 
# a <- wsub(object = m, basesub = psub, delta = 2, summary = TRUE)
# 
# test_that("wsub's estimates matches with brm model's (TST vs WAKE and MVPA vs LPA)", {
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[7, 1] > 0))) { # wilr2 = more TST less WAKE
#     expect_true(all(a$TST[From == "WAKE" & Delta > 1]$Mean > 0)) # more TST less WAKE
#     expect_true(all(a$WAKE[From == "TST" & Delta > 1]$Mean < 0)) # more WAKE less TST
#     
#   } else {
#     expect_true(all(a$TST[From == "WAKE" & Delta > 1]$Mean < 0))
#     expect_true(all(a$WAKE[From == "TST" & Delta > 1]$Mean > 0))
#     
#   }
#   
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[9, 1] > 0))) { # wilr4 = more MVPA less LPA
#     expect_true(all(a$MVPA[From == "LPA" & Delta > 1]$Mean > 0)) # more MVPA less LPA
#     expect_true(all(a$LPA[From == "MVPA" & Delta > 1]$Mean < 0)) # more LPA less MVPA
#     
#   } else {
#     expect_true(all(a$MVPA[From == "LPA" & Delta > 1]$Mean < 0))
#     expect_true(all(a$LPA[From == "MVPA" & Delta > 1]$Mean > 0))
#   }
#   
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[7, 3], summary(m$model)$fixed[7, 4]))
#     == (0 %agele% c(a$TST[From == "WAKE" & Delta == 1]$CI_low,
#                     a$TST[From == "WAKE" & Delta == 1]$CI_high))))
#   
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[9, 3], summary(m$model)$fixed[9, 4]))
#     == (0 %agele% c(a$MVPA[From == "LPA" & Delta == 1]$CI_low,
#                     a$MVPA[From == "LPA" & Delta == 1]$CI_high))))
#   
# })
# 
# sbp <- matrix(c(
#   1, -1, 1,-1, -1,
#   1, 0, -1, 0, 0,
#   0, 1, 0, -1, -1,
#   0, 0, 0, 1, -1), ncol = 5, byrow = TRUE)
# 
# cilr <- complr(data = mcompd, sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)
# 
# suppressWarnings(m <- brmcoda(complr = cilr,
#                               formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                               seed = 123))
# 
# s <- wsub(object = m, basesub = psub, delta = 2, summary = TRUE)
# 
# test_that("wsub's results matches with brm model (TST vs MVPA and LPA vs SB)", { 
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[7, 1] > 0))) { # wilr2 = more TST less MVPA
#     expect_true(all(s$TST[From == "MVPA" & Delta > 1]$Mean > 0)) # more TST less MVPA
#     expect_true(all(s$MVPA[From == "TST" & Delta > 1]$Mean < 0)) # more MVPA less TST
#   } else {
#     expect_true(all(s$TST[From == "MVPA" & Delta > 1]$Mean < 0))
#     expect_true(all(s$MVPA[From == "TST" & Delta > 1]$Mean > 0))
#   }
#   
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[9, 1] > 0))) { # wilr4 = more LPA less SB
#     expect_true(all(s$LPA[From == "SB" & Delta > 1]$Mean > 0)) # more LPA less SB
#     expect_true(all(s$SB[From == "LPA" & Delta > 1]$Mean < 0)) # more SB less LPA
#   } else {
#     expect_true(all(s$LPA[From == "SB" & Delta > 1]$Mean < 0))
#     expect_true(all(s$SB[From == "LPA" & Delta > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[7, 3], summary(m$model)$fixed[7, 4]))
#     == (0 %agele% c(s$TST[From == "MVPA" & Delta == 1]$CI_low,
#                     s$TST[From == "MVPA" & Delta == 1]$CI_high))))
#   
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[9, 3], summary(m$model)$fixed[9, 4]))
#     == (0 %agele% c(s$LPA[From == "SB" & Delta == 1]$CI_low,
#                     s$LPA[From == "SB" & Delta == 1]$CI_high))))
#   
# })
# 
# sbp <- matrix(c(
#   1, -1, -1, 1, -1,
#   1, 0, 0, -1, 0,
#   0, -1, 1, 0, -1,
#   0, 1, 0, 0, -1), ncol = 5, byrow = TRUE)
# 
# cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)
# 
# suppressWarnings(m <- brmcoda(complr = cilr,
#                               formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                               chain = 1, iter = 500, seed = 123))
# 
# d <- wsub(object = m, basesub = psub, delta = 2)
# 
# test_that("wsub's results matches with brm model (TST vs LPA and WAKE vs SB)", {
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[7, 1] > 0))) { # bilr2 = more TST less LPA
#     expect_true(all(d$TST[From == "LPA" & Delta > 1]$Mean > 0)) # more TST less LPA
#     expect_true(all(d$LPA[From == "TST" & Delta > 1]$Mean < 0)) # more LPA less TST
#   } else {
#     expect_true(all(d$TST[From == "LPA" & Delta > 1]$Mean < 0))
#     expect_true(all(d$LPA[From == "TST" & Delta > 1]$Mean > 0))
#   }
#   
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[9, 1] > 0))) { # bilr4 = more WAKE less SB
#     expect_true(all(d$WAKE[From == "SB" & Delta > 1]$Mean > 0)) # more WAKE less SB
#     expect_true(all(d$SB[From == "WAKE" & Delta > 1]$Mean < 0)) # more SB less WAKE
#   } else {
#     expect_true(all(d$WAKE[From == "SB" & Delta > 1]$Mean < 0))
#     expect_true(all(d$SB[From == "WAKE" & Delta > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[7, 3], summary(m$model)$fixed[7, 4]))
#     == (0 %agele% c(d$TST[From == "LPA" & Delta == 1]$CI_low,
#                     d$TST[From == "LPA" & Delta == 1]$CI_high))))
#   
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[9, 3], summary(m$model)$fixed[9, 4]))
#     == (0 %agele% c(d$WAKE[From == "SB" & Delta == 1]$CI_low,
#                     d$WAKE[From == "SB" & Delta == 1]$CI_high))))
#   
# })
# 
# sbp <- matrix(c(
#   1, -1, -1, -1, 1,
#   1, 0, 0, 0, -1,
#   0, -1, -1, 1, 0,
#   0, 1, -1, 0, 0), ncol = 5, byrow = TRUE)
# 
# cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)
# 
# suppressWarnings(m <- brmcoda(complr = cilr,
#                               formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                               chain = 1, iter = 500, seed = 123))
# 
# f <- wsub(object = m, basesub = psub, delta = 2)
# 
# test_that("wsub's results matches with brm model (TST vs SB and WAKE vs MVPA)", {
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[7, 1] > 0))) { # bilr2 = more TST less SB
#     expect_true(all(f$TST[From == "SB" & Delta > 1]$Mean > 0)) # more TST less SB
#     expect_true(all(f$SB[From == "TST" & Delta > 1]$Mean < 0)) # more SB less TST
#   } else {
#     expect_true(all(f$TST[From == "SB" & Delta > 1]$Mean < 0))
#     expect_true(all(f$SB[From == "TST" & Delta > 1]$Mean > 0))
#   }
#   
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[9, 1] > 0))) { # bilr4 = more WAKE less MVPA
#     expect_true(all(f$WAKE[From == "MVPA" & Delta > 1]$Mean > 0)) # more WAKE less MVPA
#     expect_true(all(f$MVPA[From == "WAKE" & Delta > 1]$Mean < 0)) # more MVPA less WAKE
#   } else {
#     expect_true(all(f$WAKE[From == "MVPA" & Delta > 1]$Mean < 0))
#     expect_true(all(f$MVPA[From == "WAKE" & Delta > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[7, 3], summary(m$model)$fixed[7, 4]))
#     == (0 %agele% c(f$TST[From == "SB" & Delta == 1]$CI_low,
#                     f$TST[From == "SB" & Delta == 1]$CI_high))))
#   
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[9, 3], summary(m$model)$fixed[9, 4]))
#     == (0 %agele% c(f$WAKE[From == "MVPA" & Delta == 1]$CI_low,
#                     f$WAKE[From == "MVPA" & Delta == 1]$CI_high))))
#   
# })
# 
# sbp <- matrix(c(
#   -1, -1, 1, -1, 1,
#   0, 0, 1, 0, -1,
#   1, -1, 0, -1, 0,
#   0, 1, 0, -1, 0), ncol = 5, byrow = TRUE)
# 
# cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)
# 
# suppressWarnings(m <- brmcoda(complr = cilr,
#                               formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                               seed = 123))
# 
# g <- wsub(object = m, basesub = psub, delta = 2)
# 
# test_that("wsub's results matches with brm model (MVPA vs SB) and (WAKE vs LPA)", {
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[7, 1] > 0))) { # wilr2 = more MVPA less SB
#     expect_true(all(g$MVPA[From == "SB" & Delta > 1]$Mean > 0)) # more MVPA less SB
#     expect_true(all(g$SB[From == "MVPA" & Delta > 1]$Mean < 0)) # more SB less MVPA
#   } else {
#     expect_true(all(g$MVPA[From == "SB" & Delta > 1]$Mean < 0))
#     expect_true(all(g$SB[From == "MVPA" & Delta > 1]$Mean > 0))
#   }
#   
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[9, 1] > 0))) { # wilr4 = more WAKE less LPA
#     expect_true(all(g$WAKE[From == "LPA" & Delta > 1]$Mean > 0)) # more WAKE less LPA
#     expect_true(all(g$LPA[From == "WAKE" & Delta > 1]$Mean < 0)) # more LPA less WAKE
#   } else {
#     expect_true(all(g$WAKE[From == "LPA" & Delta > 1]$Mean < 0))
#     expect_true(all(g$LPA[From == "WAKE" & Delta > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[7, 3], summary(m$model)$fixed[7, 4]))
#     == (0 %agele% c(g$MVPA[From == "SB" & Delta == 1]$CI_low,
#                     g$MVPA[From == "SB" & Delta == 1]$CI_high))))
#   
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[9, 3], summary(m$model)$fixed[9, 4]))
#     == (0 %agele% c(g$WAKE[From == "LPA" & Delta == 1]$CI_low,
#                     g$WAKE[From == "LPA" & Delta == 1]$CI_high))))
#   
# })
# 
#---------------------------------------------------------------------------------------------------
# Test 2-component composition for consistency between brm model and substitution model
## TST vs WAKE
test_that("wsub's results matches with brm model for 2-component composition (TST vs WAKE)", {
  
  sbp <- as.matrix(data.table(1, -1))
  cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                 parts = c("TST", "WAKE"), idvar = "ID", total = 1440)
  psub <- basesub(c("TST", "WAKE"))
  suppressWarnings(
    m <- brmcoda(complr = cilr,
                 formula = Stress ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = backend))
  a <- wsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$model)$fixed[3, 1] > 0))) { 
    expect_true(all(a$TST[From == "WAKE" & Delta > 1]$Mean > 0)) 
    expect_true(all(a$WAKE[From == "TST" & Delta > 1]$Mean < 0)) 
  } else {
    expect_true(all(a$TST[From == "WAKE" & Delta > 1]$Mean < 0))
    expect_true(all(a$WAKE[From == "TST" & Delta > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$model)$fixed[3, 3], summary(m$model)$fixed[3, 4]))
    == (0 %agele% c(a$TST[From == "WAKE" & Delta == 1]$CI_low,
                    a$TST[From == "WAKE" & Delta == 1]$CI_high))))
  
})

## TST vs MVPA
test_that("wsub's results matches with brm model for 2-component composition (TST vs MVPA)", {
  
  sbp <- as.matrix(data.table(1, -1))
  cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                 parts = c("TST", "MVPA"), idvar = "ID", total = 1440)
  psub <- basesub(c("TST", "MVPA"))
  suppressWarnings(
    m <- brmcoda(complr = cilr,
                 formula = Stress ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = backend))
  b <- wsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$model)$fixed[3, 1] > 0))) { 
    expect_true(all(b$TST[From == "MVPA" & Delta > 1]$Mean > 0)) 
    expect_true(all(b$MVPA[From == "TST" & Delta > 1]$Mean < 0)) 
  } else {
    expect_true(all(b$TST[From == "MVPA" & Delta > 1]$Mean < 0))
    expect_true(all(b$MVPA[From == "TST" & Delta > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$model)$fixed[3, 3], summary(m$model)$fixed[3, 4]))
    == (0 %agele% c(b$TST[From == "MVPA" & Delta == 1]$CI_low,
                    b$TST[From == "MVPA" & Delta == 1]$CI_high))))
  
})

## TST vs LPA
test_that("wsub's results matches with brm model for 2-component composition (TST vs LPA)", {
  sbp <- as.matrix(data.table(1, -1))
  cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                 parts = c("TST", "LPA"), idvar = "ID", total = 1440)
  psub <- basesub(c("TST", "LPA"))
  suppressWarnings(
    m <- brmcoda(complr = cilr,
                 formula = Stress ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = backend))
  c <- wsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$model)$fixed[3, 1] > 0))) { 
    expect_true(all(c$TST[From == "LPA" & Delta > 1]$Mean > 0)) 
    expect_true(all(c$LPA[From == "TST" & Delta > 1]$Mean < 0)) 
  } else {
    expect_true(all(c$TST[From == "LPA" & Delta > 1]$Mean < 0))
    expect_true(all(c$LPA[From == "TST" & Delta > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$model)$fixed[3, 3], summary(m$model)$fixed[3, 4]))
    == (0 %agele% c(c$TST[From == "LPA" & Delta == 1]$CI_low,
                    c$TST[From == "LPA" & Delta == 1]$CI_high))))
  
})

## TST vs SB
test_that("wsub's results matches with brm model for 2-component composition (TST vs SB)", {
  sbp <- as.matrix(data.table(1, -1))
  cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                 parts = c("TST", "SB"), idvar = "ID", total = 1440)
  psub <- basesub(c("TST", "SB"))
  suppressWarnings(
    m <- brmcoda(complr = cilr,
                 formula = Stress ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = backend))
  d <- wsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$model)$fixed[3, 1] > 0))) { 
    expect_true(all(d$TST[From == "SB" & Delta > 1]$Mean > 0)) 
    expect_true(all(d$SB[From == "TST" & Delta > 1]$Mean < 0)) 
  } else {
    expect_true(all(d$TST[From == "SB" & Delta > 1]$Mean < 0))
    expect_true(all(d$SB[From == "TST" & Delta > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$model)$fixed[3, 3], summary(m$model)$fixed[3, 4]))
    == (0 %agele% c(d$TST[From == "SB" & Delta == 1]$CI_low,
                    d$TST[From == "SB" & Delta == 1]$CI_high))))
  
})

## WAKE vs MVPA
test_that("wsub's results matches with brm model for 2-component composition (WAKE vs MVPA)", {
  sbp <- as.matrix(data.table(1, -1))
  cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                 parts = c("WAKE", "MVPA"), idvar = "ID", total = 1440)
  psub <- basesub(c("WAKE", "MVPA"))
  
  suppressWarnings(
    m <- brmcoda(complr = cilr,
                 formula = Stress ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = backend))
  e <- wsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$model)$fixed[3, 1] > 0))) { 
    expect_true(all(e$WAKE[From == "MVPA" & Delta > 1]$Mean > 0)) 
    expect_true(all(e$MVPA[From == "WAKE" & Delta > 1]$Mean < 0)) 
  } else {
    expect_true(all(e$WAKE[From == "MVPA" & Delta > 1]$Mean < 0))
    expect_true(all(e$MVPA[From == "WAKE" & Delta > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$model)$fixed[3, 3], summary(m$model)$fixed[3, 4]))
    == (0 %agele% c(e$WAKE[From == "MVPA" & Delta == 1]$CI_low,
                    e$WAKE[From == "MVPA" & Delta == 1]$CI_high))))
  
})

## WAKE vs LPA
test_that("wsub's results matches with brm model for 2-component composition (WAKE vs LPA)", {
  sbp <- as.matrix(data.table(1, -1))
  cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                 parts = c("WAKE", "LPA"), idvar = "ID", total = 1440)
  psub <- basesub(c("WAKE", "LPA"))
  suppressWarnings(
    m <- brmcoda(complr = cilr,
                 formula = Stress ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = backend))
  f <- wsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$model)$fixed[3, 1] > 0))) { 
    expect_true(all(f$WAKE[From == "LPA" & Delta > 1]$Mean > 0)) 
    expect_true(all(f$LPA[From == "WAKE" & Delta > 1]$Mean < 0)) 
  } else {
    expect_true(all(f$WAKE[From == "LPA" & Delta > 1]$Mean < 0))
    expect_true(all(f$LPA[From == "WAKE" & Delta > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$model)$fixed[3, 3], summary(m$model)$fixed[3, 4]))
    == (0 %agele% c(f$WAKE[From == "LPA" & Delta == 1]$CI_low,
                    f$WAKE[From == "LPA" & Delta == 1]$CI_high))))
  
})

## WAKE vs SB
test_that("wsub's results matches with brm model for 2-component composition (WAKE vs SB)", {
  sbp <- as.matrix(data.table(1, -1))
  cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                 parts = c("WAKE", "SB"), idvar = "ID", total = 1440)
  psub <- basesub(c("WAKE", "SB"))
  suppressWarnings(
    m <- brmcoda(complr = cilr,
                 formula = Stress ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = backend))
  g <- wsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$model)$fixed[3, 1] > 0))) { 
    expect_true(all(g$WAKE[From == "SB" & Delta > 1]$Mean > 0)) 
    expect_true(all(g$SB[From == "WAKE" & Delta > 1]$Mean < 0)) 
  } else {
    expect_true(all(g$WAKE[From == "SB" & Delta > 1]$Mean < 0))
    expect_true(all(g$SB[From == "WAKE" & Delta > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$model)$fixed[3, 3], summary(m$model)$fixed[3, 4]))
    == (0 %agele% c(g$WAKE[From == "SB" & Delta == 1]$CI_low,
                    g$WAKE[From == "SB" & Delta == 1]$CI_high))))
  
})

## MVPA vs LPA
test_that("wsub's results matches with brm model for 2-component composition (MVPA vs LPA)", {
  
  sbp <- as.matrix(data.table(1, -1))
  cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                 parts = c("MVPA", "LPA"), idvar = "ID", total = 1440)
  psub <- basesub(c("MVPA", "LPA"))
  suppressWarnings(
    m <- brmcoda(complr = cilr,
                 formula = Stress ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = backend))
  h <- wsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$model)$fixed[3, 1] > 0))) { 
    expect_true(all(h$MVPA[From == "LPA" & Delta > 1]$Mean > 0)) 
    expect_true(all(h$LPA[From == "MVPA" & Delta > 1]$Mean < 0)) 
  } else {
    expect_true(all(h$MVPA[From == "LPA" & Delta > 1]$Mean < 0))
    expect_true(all(h$LPA[From == "MVPA" & Delta > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$model)$fixed[3, 3], summary(m$model)$fixed[3, 4]))
    == (0 %agele% c(h$MVPA[From == "LPA" & Delta == 1]$CI_low,
                    h$MVPA[From == "LPA" & Delta == 1]$CI_high))))
  
})

## MVPA vs SB
test_that("wsub's results matches with brm model for 2-component composition (MVPA vs SB)", {
  
  sbp <- as.matrix(data.table(1, -1))
  cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                 parts = c("MVPA", "SB"), idvar = "ID", total = 1440)
  psub <- basesub(c("MVPA", "SB"))
  suppressWarnings(
    m <- brmcoda(complr = cilr,
                 formula = Stress ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = backend))
  i <- wsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$model)$fixed[3, 1] > 0))) { 
    expect_true(all(i$MVPA[From == "SB" & Delta > 1]$Mean > 0)) 
    expect_true(all(i$SB[From == "MVPA" & Delta > 1]$Mean < 0)) 
  } else {
    expect_true(all(i$MVPA[From == "SB" & Delta > 1]$Mean < 0))
    expect_true(all(i$SB[From == "MVPA" & Delta > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$model)$fixed[3, 3], summary(m$model)$fixed[3, 4]))
    == (0 %agele% c(i$MVPA[From == "SB" & Delta == 1]$CI_low,
                    i$MVPA[From == "SB" & Delta == 1]$CI_high))))
  
})

## LPA vs SB
test_that("wsub's results matches with brm model for 2-component composition (LPA vs SB)", {
  
  sbp <- as.matrix(data.table(1, -1))
  cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                 parts = c("LPA", "SB"), idvar = "ID", total = 1440)
  psub <- basesub(c("LPA", "SB"))
  suppressWarnings(
    m <- brmcoda(complr = cilr,
                 formula = Stress ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = backend))
  j <- wsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$model)$fixed[3, 1] > 0))) { 
    expect_true(all(j$LPA[From == "SB" & Delta > 1]$Mean > 0)) 
    expect_true(all(j$SB[From == "LPA" & Delta > 1]$Mean < 0)) 
  } else {
    expect_true(all(j$LPA[From == "SB" & Delta > 1]$Mean < 0))
    expect_true(all(j$SB[From == "LPA" & Delta > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$model)$fixed[3, 3], summary(m$model)$fixed[3, 4]))
    == (0 %agele% c(j$LPA[From == "SB" & Delta == 1]$CI_low,
                    j$LPA[From == "SB" & Delta == 1]$CI_high))))
  
})
