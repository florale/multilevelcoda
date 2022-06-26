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

#Compiling Stan program...
#error: PCH file uses an older PCH format that is no longer supported
#1 error generated.
#make: *** [/var/folders/fd/yjc3115x11v279rv07jxd0f40000gn/T/RtmpcRAlDr/model-15dc2cb19f7a] Error 1
#Error: An error occured during compilation! See the message above for more information.

# Model
#---------------------------------------------------------------------------------------------------
data(mcompd)
data(sbp)
data(psub)

cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")

suppressWarnings(m <- brmcoda(compilr = cilr,
                              formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
                                wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
                              chain = 1, iter = 500, seed = 123))
foreach::registerDoSEQ()

x <- wsub(object = m, substitute = psub, minute = 2)

# Testing
#---------------------------------------------------------------------------------------------------
test_that("wsub errors for invalid input", {
  
  ## check errors for missing object
  expect_error(x <- wsub(substitute = psub, minute = 2))
  
  ## check errors for missing substitute
  expect_error(x <- wsub(objetc = m, minute = 2))
  
  ## check errors when reference grid has matching names with ILRs
  rg <- data.table(bilr1 = 1)
  expect_error(x <- wsub(object = m, substitute = psub, minute = 2, regrid = rg))
  
})

test_that("wsub outputs what expected", {
  
  ## check types
  expect_type(x, "list")
  expect_equal(length(x), length(m$CompIlr$parts))
  expect_s3_class(x$TST, "data.table")
  expect_s3_class(x$WAKE, "data.table")
  expect_s3_class(x$MVPA, "data.table")
  expect_s3_class(x$LPA, "data.table")
  expect_s3_class(x$SB, "data.table")
  
  expect_type(x$TST$Mean, "double")
  expect_type(x$TST$CI_low, "double")
  expect_type(x$TST$CI_high, "double")
  expect_type(x$TST$MinSubstituted, "double")
  expect_type(x$TST$Substitute, "character")
  expect_type(x$TST$Predictor, "character")
  
  expect_type(x$WAKE$Mean, "double")
  expect_type(x$WAKE$CI_low, "double")
  expect_type(x$WAKE$CI_high, "double")
  expect_type(x$WAKE$MinSubstituted, "double")
  expect_type(x$WAKE$Substitute, "character")
  expect_type(x$WAKE$Predictor, "character")
  
  expect_type(x$MVPA$Mean, "double")
  expect_type(x$MVPA$CI_low, "double")
  expect_type(x$MVPA$CI_high, "double")
  expect_type(x$MVPA$MinSubstituted, "double")
  expect_type(x$MVPA$Substitute, "character")
  expect_type(x$MVPA$Predictor, "character")
  
  expect_type(x$LPA$Mean, "double")
  expect_type(x$LPA$CI_low, "double")
  expect_type(x$LPA$CI_high, "double")
  expect_type(x$LPA$MinSubstituted, "double")
  expect_type(x$LPA$Substitute, "character")
  expect_type(x$LPA$Predictor, "character")
  
  expect_type(x$SB$Mean, "double")
  expect_type(x$SB$CI_low, "double")
  expect_type(x$SB$CI_high, "double")
  expect_type(x$SB$MinSubstituted, "double")
  expect_type(x$SB$Substitute, "character")
  expect_type(x$SB$Predictor, "character")
  
  expect_true(ncol(x$TST) >= 6)
  expect_true(ncol(x$WAKE) >= 6)
  expect_true(ncol(x$MVPA) >= 6)
  expect_true(ncol(x$LPA) >= 6)
  expect_true(ncol(x$SB) >= 6)
  
  expect_true(all(x$TST$Predictor == "TST"))
  expect_true(all(x$WAKE$Predictor == "WAKE"))
  expect_true(all(x$MVPA$Predictor == "MVPA"))
  expect_true(all(x$LPA$Predictor == "LPA"))
  expect_true(all(x$SB$Predictor == "SB"))
  
})

test_that("wsub gives results in sensible range", {
  
  ## check values of difference in outcome
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
  
  ## check that values are opposite sign for opposite substitution
  for (i in seq_along(x)) {
    expect_true(all(x[[i]][, sign(Mean[sign(MinSubstituted) == 1]) 
                           %a!=% sign(Mean[sign(MinSubstituted) == -1]), by = Substitute]$V1))
  }
  
  ## check that results for 1 min have smaller magnitude than 2 mins
  for (i in seq_along(x)) {
    expect_true(all(x[[i]][, abs(Mean[abs(MinSubstituted) == 1]) 
                           < abs(Mean[abs(MinSubstituted) == 2])]))
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
# cilr <- compilr(data = mcompd, sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
# 
# suppressWarnings(m <- brmcoda(compilr = cilr,
#                               formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                               seed = 123))
# 
# a <- wsub(object = m, substitute = psub, minute = 2, summary = TRUE)
# 
# test_that("wsub's estimates matches with brm model's (TST vs WAKE and MVPA vs LPA)", {
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$Model)$fixed[7, 1] > 0))) { # wilr2 = more TST less WAKE
#     expect_true(all(a$TST[Substitute == "WAKE" & MinSubstituted > 1]$Mean > 0)) # more TST less WAKE
#     expect_true(all(a$WAKE[Substitute == "TST" & MinSubstituted > 1]$Mean < 0)) # more WAKE less TST
#     
#   } else {
#     expect_true(all(a$TST[Substitute == "WAKE" & MinSubstituted > 1]$Mean < 0))
#     expect_true(all(a$WAKE[Substitute == "TST" & MinSubstituted > 1]$Mean > 0))
#     
#   }
#   
#   if (isTRUE(suppressWarnings(summary(m$Model)$fixed[9, 1] > 0))) { # wilr4 = more MVPA less LPA
#     expect_true(all(a$MVPA[Substitute == "LPA" & MinSubstituted > 1]$Mean > 0)) # more MVPA less LPA
#     expect_true(all(a$LPA[Substitute == "MVPA" & MinSubstituted > 1]$Mean < 0)) # more LPA less MVPA
#     
#   } else {
#     expect_true(all(a$MVPA[Substitute == "LPA" & MinSubstituted > 1]$Mean < 0))
#     expect_true(all(a$LPA[Substitute == "MVPA" & MinSubstituted > 1]$Mean > 0))
#   }
#   
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$Model)$fixed[7, 3], summary(m$Model)$fixed[7, 4]))
#     == (0 %agele% c(a$TST[Substitute == "WAKE" & MinSubstituted == 1]$CI_low,
#                     a$TST[Substitute == "WAKE" & MinSubstituted == 1]$CI_high))))
#   
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$Model)$fixed[9, 3], summary(m$Model)$fixed[9, 4]))
#     == (0 %agele% c(a$MVPA[Substitute == "LPA" & MinSubstituted == 1]$CI_low,
#                     a$MVPA[Substitute == "LPA" & MinSubstituted == 1]$CI_high))))
#   
# })
# 
# sbp <- matrix(c(
#   1, -1, 1,-1, -1,
#   1, 0, -1, 0, 0,
#   0, 1, 0, -1, -1,
#   0, 0, 0, 1, -1), ncol = 5, byrow = TRUE)
# 
# cilr <- compilr(data = mcompd, sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
# 
# suppressWarnings(m <- brmcoda(compilr = cilr,
#                               formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                               seed = 123))
# 
# s <- wsub(object = m, substitute = psub, minute = 2, summary = TRUE)
# 
# test_that("wsub's results matches with brm model (TST vs MVPA and LPA vs SB)", { 
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$Model)$fixed[7, 1] > 0))) { # wilr2 = more TST less MVPA
#     expect_true(all(s$TST[Substitute == "MVPA" & MinSubstituted > 1]$Mean > 0)) # more TST less MVPA
#     expect_true(all(s$MVPA[Substitute == "TST" & MinSubstituted > 1]$Mean < 0)) # more MVPA less TST
#   } else {
#     expect_true(all(s$TST[Substitute == "MVPA" & MinSubstituted > 1]$Mean < 0))
#     expect_true(all(s$MVPA[Substitute == "TST" & MinSubstituted > 1]$Mean > 0))
#   }
#   
#   if (isTRUE(suppressWarnings(summary(m$Model)$fixed[9, 1] > 0))) { # wilr4 = more LPA less SB
#     expect_true(all(s$LPA[Substitute == "SB" & MinSubstituted > 1]$Mean > 0)) # more LPA less SB
#     expect_true(all(s$SB[Substitute == "LPA" & MinSubstituted > 1]$Mean < 0)) # more SB less LPA
#   } else {
#     expect_true(all(s$LPA[Substitute == "SB" & MinSubstituted > 1]$Mean < 0))
#     expect_true(all(s$SB[Substitute == "LPA" & MinSubstituted > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$Model)$fixed[7, 3], summary(m$Model)$fixed[7, 4]))
#     == (0 %agele% c(s$TST[Substitute == "MVPA" & MinSubstituted == 1]$CI_low,
#                     s$TST[Substitute == "MVPA" & MinSubstituted == 1]$CI_high))))
#   
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$Model)$fixed[9, 3], summary(m$Model)$fixed[9, 4]))
#     == (0 %agele% c(s$LPA[Substitute == "SB" & MinSubstituted == 1]$CI_low,
#                     s$LPA[Substitute == "SB" & MinSubstituted == 1]$CI_high))))
#   
# })
# 
# sbp <- matrix(c(
#   1, -1, -1, 1, -1,
#   1, 0, 0, -1, 0,
#   0, -1, 1, 0, -1,
#   0, 1, 0, 0, -1), ncol = 5, byrow = TRUE)
# 
# cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
# 
# suppressWarnings(m <- brmcoda(compilr = cilr,
#                               formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                               chain = 1, iter = 500, seed = 123))
# 
# d <- wsub(object = m, substitute = psub, minute = 2)
# 
# test_that("wsub's results matches with brm model (TST vs LPA and WAKE vs SB)", {
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$Model)$fixed[7, 1] > 0))) { # bilr2 = more TST less LPA
#     expect_true(all(d$TST[Substitute == "LPA" & MinSubstituted > 1]$Mean > 0)) # more TST less LPA
#     expect_true(all(d$LPA[Substitute == "TST" & MinSubstituted > 1]$Mean < 0)) # more LPA less TST
#   } else {
#     expect_true(all(d$TST[Substitute == "LPA" & MinSubstituted > 1]$Mean < 0))
#     expect_true(all(d$LPA[Substitute == "TST" & MinSubstituted > 1]$Mean > 0))
#   }
#   
#   if (isTRUE(suppressWarnings(summary(m$Model)$fixed[9, 1] > 0))) { # bilr4 = more WAKE less SB
#     expect_true(all(d$WAKE[Substitute == "SB" & MinSubstituted > 1]$Mean > 0)) # more WAKE less SB
#     expect_true(all(d$SB[Substitute == "WAKE" & MinSubstituted > 1]$Mean < 0)) # more SB less WAKE
#   } else {
#     expect_true(all(d$WAKE[Substitute == "SB" & MinSubstituted > 1]$Mean < 0))
#     expect_true(all(d$SB[Substitute == "WAKE" & MinSubstituted > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$Model)$fixed[7, 3], summary(m$Model)$fixed[7, 4]))
#     == (0 %agele% c(d$TST[Substitute == "LPA" & MinSubstituted == 1]$CI_low,
#                     d$TST[Substitute == "LPA" & MinSubstituted == 1]$CI_high))))
#   
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$Model)$fixed[9, 3], summary(m$Model)$fixed[9, 4]))
#     == (0 %agele% c(d$WAKE[Substitute == "SB" & MinSubstituted == 1]$CI_low,
#                     d$WAKE[Substitute == "SB" & MinSubstituted == 1]$CI_high))))
#   
# })
# 
# sbp <- matrix(c(
#   1, -1, -1, -1, 1,
#   1, 0, 0, 0, -1,
#   0, -1, -1, 1, 0,
#   0, 1, -1, 0, 0), ncol = 5, byrow = TRUE)
# 
# cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
# 
# suppressWarnings(m <- brmcoda(compilr = cilr,
#                               formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                               chain = 1, iter = 500, seed = 123))
# 
# f <- wsub(object = m, substitute = psub, minute = 2)
# 
# test_that("wsub's results matches with brm model (TST vs SB and WAKE vs MVPA)", {
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$Model)$fixed[7, 1] > 0))) { # bilr2 = more TST less SB
#     expect_true(all(f$TST[Substitute == "SB" & MinSubstituted > 1]$Mean > 0)) # more TST less SB
#     expect_true(all(f$SB[Substitute == "TST" & MinSubstituted > 1]$Mean < 0)) # more SB less TST
#   } else {
#     expect_true(all(f$TST[Substitute == "SB" & MinSubstituted > 1]$Mean < 0))
#     expect_true(all(f$SB[Substitute == "TST" & MinSubstituted > 1]$Mean > 0))
#   }
#   
#   if (isTRUE(suppressWarnings(summary(m$Model)$fixed[9, 1] > 0))) { # bilr4 = more WAKE less MVPA
#     expect_true(all(f$WAKE[Substitute == "MVPA" & MinSubstituted > 1]$Mean > 0)) # more WAKE less MVPA
#     expect_true(all(f$MVPA[Substitute == "WAKE" & MinSubstituted > 1]$Mean < 0)) # more MVPA less WAKE
#   } else {
#     expect_true(all(f$WAKE[Substitute == "MVPA" & MinSubstituted > 1]$Mean < 0))
#     expect_true(all(f$MVPA[Substitute == "WAKE" & MinSubstituted > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$Model)$fixed[7, 3], summary(m$Model)$fixed[7, 4]))
#     == (0 %agele% c(f$TST[Substitute == "SB" & MinSubstituted == 1]$CI_low,
#                     f$TST[Substitute == "SB" & MinSubstituted == 1]$CI_high))))
#   
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$Model)$fixed[9, 3], summary(m$Model)$fixed[9, 4]))
#     == (0 %agele% c(f$WAKE[Substitute == "MVPA" & MinSubstituted == 1]$CI_low,
#                     f$WAKE[Substitute == "MVPA" & MinSubstituted == 1]$CI_high))))
#   
# })
# 
# sbp <- matrix(c(
#   -1, -1, 1, -1, 1,
#   0, 0, 1, 0, -1,
#   1, -1, 0, -1, 0,
#   0, 1, 0, -1, 0), ncol = 5, byrow = TRUE)
# 
# cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
# 
# suppressWarnings(m <- brmcoda(compilr = cilr,
#                               formula = STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                                 wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                               seed = 123))
# 
# g <- wsub(object = m, substitute = psub, minute = 2)
# 
# test_that("wsub's results matches with brm model (MVPA vs SB) and (WAKE vs LPA)", {
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$Model)$fixed[7, 1] > 0))) { # wilr2 = more MVPA less SB
#     expect_true(all(g$MVPA[Substitute == "SB" & MinSubstituted > 1]$Mean > 0)) # more MVPA less SB
#     expect_true(all(g$SB[Substitute == "MVPA" & MinSubstituted > 1]$Mean < 0)) # more SB less MVPA
#   } else {
#     expect_true(all(g$MVPA[Substitute == "SB" & MinSubstituted > 1]$Mean < 0))
#     expect_true(all(g$SB[Substitute == "MVPA" & MinSubstituted > 1]$Mean > 0))
#   }
#   
#   if (isTRUE(suppressWarnings(summary(m$Model)$fixed[9, 1] > 0))) { # wilr4 = more WAKE less LPA
#     expect_true(all(g$WAKE[Substitute == "LPA" & MinSubstituted > 1]$Mean > 0)) # more WAKE less LPA
#     expect_true(all(g$LPA[Substitute == "WAKE" & MinSubstituted > 1]$Mean < 0)) # more LPA less WAKE
#   } else {
#     expect_true(all(g$WAKE[Substitute == "LPA" & MinSubstituted > 1]$Mean < 0))
#     expect_true(all(g$LPA[Substitute == "WAKE" & MinSubstituted > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$Model)$fixed[7, 3], summary(m$Model)$fixed[7, 4]))
#     == (0 %agele% c(g$MVPA[Substitute == "SB" & MinSubstituted == 1]$CI_low,
#                     g$MVPA[Substitute == "SB" & MinSubstituted == 1]$CI_high))))
#   
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$Model)$fixed[9, 3], summary(m$Model)$fixed[9, 4]))
#     == (0 %agele% c(g$WAKE[Substitute == "LPA" & MinSubstituted == 1]$CI_low,
#                     g$WAKE[Substitute == "LPA" & MinSubstituted == 1]$CI_high))))
#   
# })
# 
#---------------------------------------------------------------------------------------------------
# Test 2-component composition for consistency between brm model and substitution model
## TST vs WAKE
sbp <- as.matrix(data.table(1, -1))
cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("TST", "WAKE"), idvar = "ID")
psub <- possub(c("TST", "WAKE"))
suppressWarnings(m <- brmcoda(compilr = cilr,
                               formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                               chain = 1, iter = 500, seed = 123))
a <- wsub(object = m, substitute = psub, minute = 2)

test_that("wsub's results matches with brm model for 2-component composition (TST vs WAKE)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[3, 1] > 0))) { 
    expect_true(all(a$TST[Substitute == "WAKE" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(a$WAKE[Substitute == "TST" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(a$TST[Substitute == "WAKE" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(a$WAKE[Substitute == "TST" & MinSubstituted > 1]$Mean > 0))
  }

  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[3, 3], summary(m$Model)$fixed[3, 4]))
    == (0 %agele% c(a$TST[Substitute == "WAKE" & MinSubstituted == 1]$CI_low,
                    a$TST[Substitute == "WAKE" & MinSubstituted == 1]$CI_high))))
  
})

## TST vs MVPA
sbp <- as.matrix(data.table(1, -1))
cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("TST", "MVPA"), idvar = "ID")
psub <- possub(c("TST", "MVPA"))
suppressWarnings(m <- brmcoda(compilr = cilr,
                              formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                              chain = 1, iter = 500, seed = 123))
b <- wsub(object = m, substitute = psub, minute = 2)

test_that("wsub's results matches with brm model for 2-component composition (TST vs MVPA)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[3, 1] > 0))) { 
    expect_true(all(b$TST[Substitute == "MVPA" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(b$MVPA[Substitute == "TST" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(b$TST[Substitute == "MVPA" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(b$MVPA[Substitute == "TST" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[3, 3], summary(m$Model)$fixed[3, 4]))
    == (0 %agele% c(b$TST[Substitute == "MVPA" & MinSubstituted == 1]$CI_low,
                    b$TST[Substitute == "MVPA" & MinSubstituted == 1]$CI_high))))
  
})

## TST vs LPA
sbp <- as.matrix(data.table(1, -1))
cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("TST", "LPA"), idvar = "ID")
psub <- possub(c("TST", "LPA"))
suppressWarnings(m <- brmcoda(compilr = cilr,
                              formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                              chain = 1, iter = 500, seed = 123))
c <- wsub(object = m, substitute = psub, minute = 2)

test_that("wsub's results matches with brm model for 2-component composition (TST vs LPA)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[3, 1] > 0))) { 
    expect_true(all(c$TST[Substitute == "LPA" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(c$LPA[Substitute == "TST" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(c$TST[Substitute == "LPA" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(c$LPA[Substitute == "TST" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[3, 3], summary(m$Model)$fixed[3, 4]))
    == (0 %agele% c(c$TST[Substitute == "LPA" & MinSubstituted == 1]$CI_low,
                    c$TST[Substitute == "LPA" & MinSubstituted == 1]$CI_high))))
  
})

## TST vs SB
sbp <- as.matrix(data.table(1, -1))
cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("TST", "SB"), idvar = "ID")
psub <- possub(c("TST", "SB"))
suppressWarnings(m <- brmcoda(compilr = cilr,
                              formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                              chain = 1, iter = 500, seed = 123))
d <- wsub(object = m, substitute = psub, minute = 2)

test_that("wsub's results matches with brm model for 2-component composition (TST vs SB)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[3, 1] > 0))) { 
    expect_true(all(d$TST[Substitute == "SB" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(d$SB[Substitute == "TST" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(d$TST[Substitute == "SB" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(d$SB[Substitute == "TST" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[3, 3], summary(m$Model)$fixed[3, 4]))
    == (0 %agele% c(d$TST[Substitute == "SB" & MinSubstituted == 1]$CI_low,
                    d$TST[Substitute == "SB" & MinSubstituted == 1]$CI_high))))
  
})

## WAKE vs MVPA
sbp <- as.matrix(data.table(1, -1))
cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("WAKE", "MVPA"), idvar = "ID")
psub <- possub(c("WAKE", "MVPA"))

suppressWarnings(m <- brmcoda(compilr = cilr,
                              formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                              chain = 1, iter = 500, seed = 123))
e <- wsub(object = m, substitute = psub, minute = 2)

test_that("wsub's results matches with brm model for 2-component composition (WAKE vs MVPA)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[3, 1] > 0))) { 
    expect_true(all(e$WAKE[Substitute == "MVPA" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(e$MVPA[Substitute == "WAKE" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(e$WAKE[Substitute == "MVPA" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(e$MVPA[Substitute == "WAKE" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[3, 3], summary(m$Model)$fixed[3, 4]))
    == (0 %agele% c(e$WAKE[Substitute == "MVPA" & MinSubstituted == 1]$CI_low,
                    e$WAKE[Substitute == "MVPA" & MinSubstituted == 1]$CI_high))))
  
})

## WAKE vs LPA
sbp <- as.matrix(data.table(1, -1))
cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("WAKE", "LPA"), idvar = "ID")
psub <- possub(c("WAKE", "LPA"))
suppressWarnings(m <- brmcoda(compilr = cilr,
                              formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                              chain = 1, iter = 500, seed = 123))
f <- wsub(object = m, substitute = psub, minute = 2)

test_that("wsub's results matches with brm model for 2-component composition (WAKE vs LPA)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[3, 1] > 0))) { 
    expect_true(all(f$WAKE[Substitute == "LPA" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(f$LPA[Substitute == "WAKE" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(f$WAKE[Substitute == "LPA" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(f$LPA[Substitute == "WAKE" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[3, 3], summary(m$Model)$fixed[3, 4]))
    == (0 %agele% c(f$WAKE[Substitute == "LPA" & MinSubstituted == 1]$CI_low,
                    f$WAKE[Substitute == "LPA" & MinSubstituted == 1]$CI_high))))
  
})

## WAKE vs SB
sbp <- as.matrix(data.table(1, -1))
cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("WAKE", "SB"), idvar = "ID")
psub <- possub(c("WAKE", "SB"))
suppressWarnings(m <- brmcoda(compilr = cilr,
                              formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                              chain = 1, iter = 500, seed = 123))
g <- wsub(object = m, substitute = psub, minute = 2)

test_that("wsub's results matches with brm model for 2-component composition (WAKE vs SB)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[3, 1] > 0))) { 
    expect_true(all(g$WAKE[Substitute == "SB" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(g$SB[Substitute == "WAKE" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(g$WAKE[Substitute == "SB" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(g$SB[Substitute == "WAKE" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[3, 3], summary(m$Model)$fixed[3, 4]))
    == (0 %agele% c(g$WAKE[Substitute == "SB" & MinSubstituted == 1]$CI_low,
                    g$WAKE[Substitute == "SB" & MinSubstituted == 1]$CI_high))))
  
})

## MVPA vs LPA
sbp <- as.matrix(data.table(1, -1))
cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("MVPA", "LPA"), idvar = "ID")
psub <- possub(c("MVPA", "LPA"))
suppressWarnings(m <- brmcoda(compilr = cilr,
                              formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                              chain = 1, iter = 500, seed = 123))
h <- wsub(object = m, substitute = psub, minute = 2)

test_that("wsub's results matches with brm model for 2-component composition (MVPA vs LPA)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[3, 1] > 0))) { 
    expect_true(all(h$MVPA[Substitute == "LPA" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(h$LPA[Substitute == "MVPA" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(h$MVPA[Substitute == "LPA" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(h$LPA[Substitute == "MVPA" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[3, 3], summary(m$Model)$fixed[3, 4]))
    == (0 %agele% c(h$MVPA[Substitute == "LPA" & MinSubstituted == 1]$CI_low,
                    h$MVPA[Substitute == "LPA" & MinSubstituted == 1]$CI_high))))
  
})

## MVPA vs SB
sbp <- as.matrix(data.table(1, -1))
cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("MVPA", "SB"), idvar = "ID")
psub <- possub(c("MVPA", "SB"))
suppressWarnings(m <- brmcoda(compilr = cilr,
                              formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                              chain = 1, iter = 500, seed = 123))
i <- wsub(object = m, substitute = psub, minute = 2)

test_that("wsub's results matches with brm model for 2-component composition (MVPA vs SB)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[3, 1] > 0))) { 
    expect_true(all(i$MVPA[Substitute == "SB" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(i$SB[Substitute == "MVPA" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(i$MVPA[Substitute == "SB" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(i$SB[Substitute == "MVPA" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[3, 3], summary(m$Model)$fixed[3, 4]))
    == (0 %agele% c(i$MVPA[Substitute == "SB" & MinSubstituted == 1]$CI_low,
                    i$MVPA[Substitute == "SB" & MinSubstituted == 1]$CI_high))))
  
})

## LPA vs SB
sbp <- as.matrix(data.table(1, -1))
cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("LPA", "SB"), idvar = "ID")
psub <- possub(c("LPA", "SB"))
suppressWarnings(m <- brmcoda(compilr = cilr,
                              formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                              chain = 1, iter = 500, seed = 123))
j <- wsub(object = m, substitute = psub, minute = 2)

test_that("wsub's results matches with brm model for 2-component composition (LPA vs SB)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[3, 1] > 0))) { 
    expect_true(all(j$LPA[Substitute == "SB" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(j$SB[Substitute == "LPA" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(j$LPA[Substitute == "SB" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(j$SB[Substitute == "LPA" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[3, 3], summary(m$Model)$fixed[3, 4]))
    == (0 %agele% c(j$LPA[Substitute == "SB" & MinSubstituted == 1]$CI_low,
                    j$LPA[Substitute == "SB" & MinSubstituted == 1]$CI_high))))
  
})
