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

# if (tolower(Sys.info()[["sysname"]]) == "darwin" && R.version[["arch"]] == "aarch64") {
#   skip_on_ci()
# }

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

cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")

suppressWarnings(
  m <- brmcoda(compilr = cilr,
               formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
               chain = 1, iter = 500, seed = 123,
               backend = backend))
foreach::registerDoSEQ()

x <- bsub(object = m, basesub = psub, delta = 2)

# Testing
#---------------------------------------------------------------------------------------------------

# test_that("bsub errors for invalid input", {
# 
#   ## missing object
#   expect_error(x <- bsub(basesub = psub, delta = 2))
# 
#   ## missing basesub
#   expect_error(x <- bsub(object = m, delta = 2))
# 
#   ## not brmcoda model
#   m1 <- lmer(STRESS ~ 1 + (1 | ID), data = mcompd)
#   expect_error(x <- bsub(object = m1, basesub = psub, delta = 2))
# 
#   ## invalid delta
#   expect_error(x <- bsub(object = m, basesub = psub, delta = -10))
#   expect_error(x <- bsub(object = m, basesub = psub, delta = 1:10))
# 
#   ## missing delta
#   expect_error(x <- substitution(object = m1, basesub = psub))
# 
#   ## basesub does not have the same components as parts in cilr
#   ps <- basesub(c("WAKE", "MVPA", "LPA", "SB"))
#   expect_error(x <- bsub(object = m, basesub = ps, delta = 2))
# 
#   ## basesub does have the same names as parts in cilr
#   ps <- basesub(parts = c("Sleep", "WAKE", "MVPA", "LPA", "SB"))
#   expect_error(x <- bsub(object = m, basesub = ps, delta = 2))
# 

# })

# test_that("bsub works as expected for adjusted/unadjusted model", {
#   
#   ## reference grid is provided for unadjusted model
#   suppressWarnings(
#     m2 <- brmcoda(compilr = cilr,
#                   formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#                   chain = 1, iter = 500, seed = 123,
#                   backend = "cmdstanr"))
#   rg <- data.table(Age = 1)
#   expect_warning(x <- bsub(object = m2, basesub = psub, delta = 2, regrid = rg))
#   
#   ## incorect reference grid 1
#   rg <- data.table(Age = 1)
#   expect_error(x <- bsub(object = m, basesub = psub, delta = 2, regrid = rg))
#   
#   ## reference grid has matching names with ILRs
#   rg <- data.table(bilr1 = 1)
#   expect_error(x <- bsub(object = m, basesub = psub, delta = 2, regrid = rg))
#   
#   ## incorect reference grid 2
#   rg <- data.table(bilr1 = 1, Age = 1)
#   expect_error(x <- bsub(object = m, basesub = psub, delta = 2, regrid = rg))
#   
#   # delta out of range
#   expect_error(x <- bsub(object = m, basesub = psub, delta = 1000))
#   
#   ## function knows to use correct user's specified reference grid
#   rg <- data.table(Female = 1)
#   x3 <- bsub(object = m, basesub = psub, delta = 2, regrid = rg)
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
#   x4 <- bsub(object = m, basesub = psub, delta = 2, summary = TRUE)
#   x5 <- bsub(object = m, basesub = psub, delta = 2)
#   expect_equal(x4, x5)
#   
#   ## keep prediction at each level of refrence grid 
#   cilr <- compilr(data = mcompd[ID %in% c(1:5, 185:190), .SD[1:3], by = ID], sbp = sbp,
#                   parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#   
#   suppressWarnings(
#     m <- brmcoda(compilr = cilr,
#                  formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                    wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                  chain = 1, iter = 500, seed = 123,
#                  backend = "cmdstanr"))
#   
#   x6 <- bsub(object = m, basesub = psub, delta = 2, summary = FALSE)
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
# test_that("bsub checks for user-specified reference composition", {
#   
#   # incorrect length
#   ref1 <- c(400, 60, 500, 60)
#   expect_error(bsub(object = m, basesub = psub, recomp = ref1, delta = 2))
#   
#   # incorrect class
#   ref2 <- c("400", "100", "500", "200", "200")
#   expect_error(bsub(object = m, basesub = psub, recomp = ref2, delta = 2))
#   
#   # incorrect class
#   ref3 <- c(400, 100, 500, 200, 200)
#   expect_error(x <- bsub(object = m, basesub = psub, recomp = ref3, delta = 2))
#   
#   # values outside of possible range
#   ref4 <- c(100, 100, 900, 100, 240)
#   expect_error(x <- bsub(object = m, basesub = psub, recomp = ref4, delta = 2))
#   
#   # include 0
#   ref5 <- c(100, 200, 900, 0, 240)
#   expect_error(x <- bsub(object = m, basesub = psub, recomp = ref5, delta = 2))
#   
# })

test_that("bsub outputs what expected", {
  
  ## types
  expect_type(x, "list")
  expect_equal(length(x), length(m$CompILR$parts))
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
  
  expect_true(all(x$TST$Level == "between"))
  expect_true(all(x$WAKE$Level == "between"))
  expect_true(all(x$MVPA$Level == "between"))
  expect_true(all(x$LPA$Level == "between"))
  expect_true(all(x$SB$Level == "between"))
  
  })

if (tolower(Sys.info()[["sysname"]]) != "darwin") {
  test_that("bsub gives results in sensible range", {
    
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
}

test_that("bsub gives results in expected direction and magnitude", {
    
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

#---------------------------------------------------------------------------------------------------
# Test 2-component composition for consistency between brm model and substitution model
## TST vs WAKE
test_that("bsub's results matches with brm model for 2-component composition (TST vs WAKE)", {
  
  sbp <- as.matrix(data.table(1, -1))
  cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                  parts = c("TST", "WAKE"), idvar = "ID")
  psub <- basesub(c("TST", "WAKE"))
  suppressWarnings(
    m <- brmcoda(compilr = cilr,
                 formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = "cmdstanr"))
  a <- bsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) { 
    expect_true(all(a$TST[From == "WAKE" & Delta > 1]$Mean > 0)) 
    expect_true(all(a$WAKE[From == "TST" & Delta > 1]$Mean < 0)) 
  } else {
    expect_true(all(a$TST[From == "WAKE" & Delta > 1]$Mean < 0))
    expect_true(all(a$WAKE[From == "TST" & Delta > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
    == (0 %agele% c(a$TST[From == "WAKE" & Delta == 1]$CI_low,
                    a$TST[From == "WAKE" & Delta == 1]$CI_high))))
  
})

## TST vs MVPA
test_that("bsub's results matches with brm model for 2-component composition (TST vs MVPA)", {
  
  sbp <- as.matrix(data.table(1, -1))
  cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                  parts = c("TST", "MVPA"), idvar = "ID")
  psub <- basesub(c("TST", "MVPA"))
  suppressWarnings(
    m <- brmcoda(compilr = cilr,
                 formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = "cmdstanr"))
  b <- bsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) {
    expect_true(all(b$TST[From == "MVPA" & Delta > 1]$Mean > 0))
    expect_true(all(b$MVPA[From == "TST" & Delta > 1]$Mean < 0))
  } else {
    expect_true(all(b$TST[From == "MVPA" & Delta > 1]$Mean < 0))
    expect_true(all(b$MVPA[From == "TST" & Delta > 1]$Mean > 0))
  }

  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
    == (0 %agele% c(b$TST[From == "MVPA" & Delta == 1]$CI_low,
                    b$TST[From == "MVPA" & Delta == 1]$CI_high))))

})

## TST vs LPA
test_that("bsub's results matches with brm model for 2-component composition (TST vs LPA)", {

  sbp <- as.matrix(data.table(1, -1))
  cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                  parts = c("TST", "LPA"), idvar = "ID")
  psub <- basesub(c("TST", "LPA"))
  suppressWarnings(
    m <- brmcoda(compilr = cilr,
                 formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = "cmdstanr"))
  c <- bsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) {
    expect_true(all(c$TST[From == "LPA" & Delta > 1]$Mean > 0))
    expect_true(all(c$LPA[From == "TST" & Delta > 1]$Mean < 0))
  } else {
    expect_true(all(c$TST[From == "LPA" & Delta > 1]$Mean < 0))
    expect_true(all(c$LPA[From == "TST" & Delta > 1]$Mean > 0))
  }

  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
    == (0 %agele% c(c$TST[From == "LPA" & Delta == 1]$CI_low,
                    c$TST[From == "LPA" & Delta == 1]$CI_high))))

})

## TST vs SB
test_that("bsub's results matches with brm model for 2-component composition (TST vs SB)", {

  sbp <- as.matrix(data.table(1, -1))
  cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                  parts = c("TST", "SB"), idvar = "ID")
  psub <- basesub(c("TST", "SB"))
  suppressWarnings(
    m <- brmcoda(compilr = cilr,
                 formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = "cmdstanr"))
  d <- bsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) {
    expect_true(all(d$TST[From == "SB" & Delta > 1]$Mean > 0))
    expect_true(all(d$SB[From == "TST" & Delta > 1]$Mean < 0))
  } else {
    expect_true(all(d$TST[From == "SB" & Delta > 1]$Mean < 0))
    expect_true(all(d$SB[From == "TST" & Delta > 1]$Mean > 0))
  }

  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
    == (0 %agele% c(d$TST[From == "SB" & Delta == 1]$CI_low,
                    d$TST[From == "SB" & Delta == 1]$CI_high))))

})

## WAKE vs MVPA
test_that("bsub's results matches with brm model for 2-component composition (WAKE vs MVPA)", {

  sbp <- as.matrix(data.table(1, -1))
  cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                  parts = c("WAKE", "MVPA"), idvar = "ID")
  psub <- basesub(c("WAKE", "MVPA"))
  
  suppressWarnings(
    m <- brmcoda(compilr = cilr,
                 formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = "cmdstanr"))
  e <- bsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) {
    expect_true(all(e$WAKE[From == "MVPA" & Delta > 1]$Mean > 0))
    expect_true(all(e$MVPA[From == "WAKE" & Delta > 1]$Mean < 0))
  } else {
    expect_true(all(e$WAKE[From == "MVPA" & Delta > 1]$Mean < 0))
    expect_true(all(e$MVPA[From == "WAKE" & Delta > 1]$Mean > 0))
  }

  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
    == (0 %agele% c(e$WAKE[From == "MVPA" & Delta == 1]$CI_low,
                    e$WAKE[From == "MVPA" & Delta == 1]$CI_high))))

})

## WAKE vs LPA
test_that("bsub's results matches with brm model for 2-component composition (WAKE vs LPA)", {

  sbp <- as.matrix(data.table(1, -1))
  cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                  parts = c("WAKE", "LPA"), idvar = "ID")
  psub <- basesub(c("WAKE", "LPA"))
  suppressWarnings(
    m <- brmcoda(compilr = cilr,
                 formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = "cmdstanr"))
  f <- bsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) {
    expect_true(all(f$WAKE[From == "LPA" & Delta > 1]$Mean > 0))
    expect_true(all(f$LPA[From == "WAKE" & Delta > 1]$Mean < 0))
  } else {
    expect_true(all(f$WAKE[From == "LPA" & Delta > 1]$Mean < 0))
    expect_true(all(f$LPA[From == "WAKE" & Delta > 1]$Mean > 0))
  }

  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
    == (0 %agele% c(f$WAKE[From == "LPA" & Delta == 1]$CI_low,
                    f$WAKE[From == "LPA" & Delta == 1]$CI_high))))

})

## WAKE vs SB
test_that("bsub's results matches with brm model for 2-component composition (WAKE vs SB)", {

  sbp <- as.matrix(data.table(1, -1))
  cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                  parts = c("WAKE", "SB"), idvar = "ID")
  psub <- basesub(c("WAKE", "SB"))
  suppressWarnings(
    m <- brmcoda(compilr = cilr,
                 formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = "cmdstanr"))
  g <- bsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) {
    expect_true(all(g$WAKE[From == "SB" & Delta > 1]$Mean > 0))
    expect_true(all(g$SB[From == "WAKE" & Delta > 1]$Mean < 0))
  } else {
    expect_true(all(g$WAKE[From == "SB" & Delta > 1]$Mean < 0))
    expect_true(all(g$SB[From == "WAKE" & Delta > 1]$Mean > 0))
  }

  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
    == (0 %agele% c(g$WAKE[From == "SB" & Delta == 1]$CI_low,
                    g$WAKE[From == "SB" & Delta == 1]$CI_high))))

})

## MVPA vs LPA
test_that("bsub's results matches with brm model for 2-component composition (MVPA vs LPA)", {

  sbp <- as.matrix(data.table(1, -1))
  cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                  parts = c("MVPA", "LPA"), idvar = "ID")
  psub <- basesub(c("MVPA", "LPA"))
  suppressWarnings(m <- brmcoda(compilr = cilr,
                                formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                                chain = 1, iter = 500, seed = 123,
                                backend = "cmdstanr"))
  h <- bsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) {
    expect_true(all(h$MVPA[From == "LPA" & Delta > 1]$Mean > 0))
    expect_true(all(h$LPA[From == "MVPA" & Delta > 1]$Mean < 0))
  } else {
    expect_true(all(h$MVPA[From == "LPA" & Delta > 1]$Mean < 0))
    expect_true(all(h$LPA[From == "MVPA" & Delta > 1]$Mean > 0))
  }

  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
    == (0 %agele% c(h$MVPA[From == "LPA" & Delta == 1]$CI_low,
                    h$MVPA[From == "LPA" & Delta == 1]$CI_high))))

})

## MVPA vs SB
test_that("bsub's results matches with brm model for 2-component composition (MVPA vs SB)", {

  sbp <- as.matrix(data.table(1, -1))
  cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                  parts = c("MVPA", "SB"), idvar = "ID")
  psub <- basesub(c("MVPA", "SB"))
  suppressWarnings(
    m <- brmcoda(compilr = cilr,
                 formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = "cmdstanr"))
  i <- bsub(object = m, basesub = psub, delta = 1:2)
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) {
    expect_true(all(i$MVPA[From == "SB" & Delta > 1]$Mean > 0))
    expect_true(all(i$SB[From == "MVPA" & Delta > 1]$Mean < 0))
  } else {
    expect_true(all(i$MVPA[From == "SB" & Delta > 1]$Mean < 0))
    expect_true(all(i$SB[From == "MVPA" & Delta > 1]$Mean > 0))
  }

  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
    == (0 %agele% c(i$MVPA[From == "SB" & Delta == 1]$CI_low,
                    i$MVPA[From == "SB" & Delta == 1]$CI_high))))

})

## LPA vs SB
test_that("bsub's results matches with brm model for 2-component composition (LPA vs SB)", {

  sbp <- as.matrix(data.table(1, -1))
  cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                  parts = c("LPA", "SB"), idvar = "ID")
  psub <- basesub(c("LPA", "SB"))
  suppressWarnings(
    m <- brmcoda(compilr = cilr,
                 formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                 chain = 1, iter = 500, seed = 123,
                 backend = "cmdstanr"))
  j <- bsub(object = m, basesub = psub, delta = 1:2)
  
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) {
    expect_true(all(j$LPA[From == "SB" & Delta > 1]$Mean > 0))
    expect_true(all(j$SB[From == "LPA" & Delta > 1]$Mean < 0))
  } else {
    expect_true(all(j$LPA[From == "SB" & Delta > 1]$Mean < 0))
    expect_true(all(j$SB[From == "LPA" & Delta > 1]$Mean > 0))
  }

  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
    == (0 %agele% c(j$LPA[From == "SB" & Delta == 1]$CI_low,
                    j$LPA[From == "SB" & Delta == 1]$CI_high))))

})

# #---------------------------------------------------------------------------------------------------
# ## test that results from different sbp are (nearly) identical
# # model sub dataset, chain = 1, iter = 500,
# all.equal(x$TST$Mean, f$TST$Mean, tolerance = .15) # "Mean relative difference: 0.09914661"
# all.equal(x$WAKE$Mean, f$WAKE$Mean) # "Mean relative difference: 0.0698194"
# all.equal(x$MVPA$Mean, f$MVPA$Mean) # "Mean relative difference: 0.0781148"
# all.equal(x$LPA$Mean, f$LPA$Mean) # "Mean relative difference: 0.08996186"
# all.equal(x$SB$Mean, f$SB$Mean) # "Mean relative difference: 0.07883451"
# 
# # model as usual
# data("sbp")
# cilr <- compilr(data = mcompd, sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
# 
# suppressWarnings(m <- brmcoda(compilr = cilr,
#                               formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                               chain = 10, iter = 10000, cores = 8, seed = 123))
# 
# h <- bsub(object = m, basesub = psub, delta = 1:2)
# 
# sbp <- matrix(c(
#   1, -1, -1, -1, 1,
#   1, 0, 0, 0, -1,
#   0, 1, 1, -1, 0,
#   0, 1, -1, 0, 0), ncol = 5, byrow = TRUE)
# 
# cilr <- compilr(data = mcompd, sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
# 
# suppressWarnings(m <- brmcoda(compilr = cilr,
#                               formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                               chain = 10, iter = 10000, cores = 8, seed = 123))
# 
# j <- bsub(object = m, basesub = psub, delta = 1:2)
# 
# # chain = 4, iter = 4000
# all.equal(h$TST$Mean, j$TST$Mean) # "Mean relative difference: 0.008437258"
# all.equal(h$WAKE$Mean, j$WAKE$Mean) # "Mean relative difference: 0.009946997"
# all.equal(h$MVPA$Mean, j$MVPA$Mean) # "Mean relative difference: 0.01363061"
# all.equal(h$LPA$Mean, j$LPA$Mean) # "Mean relative difference: 0.01120686"
# all.equal(h$SB$Mean, j$SB$Mean) # "Mean relative difference: 0.01221668"
# 
# # chain = 4, iter = 10000
# all.equal(h$TST$Mean, j$TST$Mean) # "Mean relative difference: 0.009482279"
# all.equal(h$WAKE$Mean, j$WAKE$Mean) # "Mean relative difference: 0.003160447"
# all.equal(h$MVPA$Mean, j$MVPA$Mean) # "Mean relative difference: 0.01092816"
# all.equal(h$LPA$Mean, j$LPA$Mean) # "Mean relative difference: 0.01207495"
# all.equal(h$SB$Mean, j$SB$Mean) # "Mean relative difference: 0.008510161"
# 
# # chain = 8, iter = 4000 - 2nd BEST PERF
# all.equal(h$TST$Mean, j$TST$Mean) # "Mean relative difference: 0.006502826"
# all.equal(h$WAKE$Mean, j$WAKE$Mean) # "Mean relative difference: 0.005147993"
# all.equal(h$MVPA$Mean, j$MVPA$Mean) # "Mean relative difference: 0.007032578"
# all.equal(h$LPA$Mean, j$LPA$Mean) # "Mean relative difference: 0.009487211"
# all.equal(h$SB$Mean, j$SB$Mean) # "Mean relative difference: 0.006800187"
# 
# # chain = 10, iter = 4000 - BEST PERF
# all.equal(h$TST$Mean, j$TST$Mean) # "Mean relative difference: 0.006383796"
# all.equal(h$WAKE$Mean, j$WAKE$Mean) # "Mean relative difference: 0.003114632"
# all.equal(h$MVPA$Mean, j$MVPA$Mean) # "Mean relative difference: 0.005826673"
# all.equal(h$LPA$Mean, j$LPA$Mean) # "Mean relative difference: 0.009096224"
# all.equal(h$SB$Mean, j$SB$Mean) # "Mean relative difference: 0.005906579"
# 
# # chain = 10, iter = 10000
# all.equal(h$TST$Mean, j$TST$Mean) # "Mean relative difference: 0.007480993"
# all.equal(h$WAKE$Mean, j$WAKE$Mean) # "Mean relative difference: 0.005886063"
# all.equal(h$MVPA$Mean, j$MVPA$Mean) # "Mean relative difference: 0.01124038"
# all.equal(h$LPA$Mean, j$LPA$Mean) # "Mean relative difference: 0.009384262"
# all.equal(h$SB$Mean, j$SB$Mean) # "Mean relative difference: 0.008712406"
# 
# ## NOTE: iteration slows down bsub
# ## number of chains probably matters most
# #---------------------------------------------------------------------------------------------------
