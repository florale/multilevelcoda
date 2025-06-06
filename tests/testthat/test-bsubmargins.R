# skip_on_cran()
# 
# if (!requireNamespace("cmdstanr", quietly = TRUE)) {
#   backend <- "rstan"
#   ## if using rstan backend, models can crash on Windows
#   ## so skip if on windows and cannot use cmdstanr
#   skip_on_os("windows")
# } else {
#   if (isFALSE(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))) {
#     backend <- "cmdstanr"
#   }
# }
# 
# # Packages
# library(testthat)
# library(data.table)
# library(multilevelcoda)
# library(extraoperators)
# library(brms)
# library(lme4)
# 
# # model
# #---------------------------------------------------------------------------------------------------
# data(mcompd)
# data(sbp)
# data(psub)
# 
# cilr <- complr(data = mcompd[ID %in% 1:200, .SD[1:5], by = ID], sbp = sbp,
#                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID", total = 1440)
# 
# suppressWarnings(
#   m <- brmcoda(complr = cilr,
#                formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#                  wilr1 + wilr2 + wilr3 + wilr4 + Female + (1 | ID),
#                chain = 1, iter = 500, seed = 123,
#                backend = backend))
# foreach::registerDoSEQ()
# 
# x <- bsubmargins(object = m, basesub = psub, delta = 2)
# 
# # Testing
# #---------------------------------------------------------------------------------------------------
# 
# # test_that("bsubmargins errors for invalid input", {
# # 
# #   ## missing object
# #   expect_error(x <- bsubmargins(basesub = psub, delta = 2))
# # 
# #   ## missing basesub
# #   expect_error(x <- bsubmargins(object = m, delta = 2))
# # 
# #   ## not brmcoda model
# #   m1 <- lmer(Stress ~ 1 + (1 | ID), data = mcompd)
# #   expect_error(x <- bsubmargins(object = m1, basesub = psub, delta = 2))
# # 
# #   ## not brmcoda model
# #   expect_error(x <- bsubmargins(object = m, basesub = psub, delta = -10))
# #   expect_error(x <- bsubmargins(object = m, basesub = psub, delta = 1:10))
# # 
# #   ## missing delta
# #   expect_error(x <- substitution(object = m1, basesub = psub))
# # 
# #   ## basesub has the same components as parts in cilr
# #   ps <- build.basesub(c("WAKE", "MVPA", "LPA", "SB"))
# #   expect_error(x <- bsubmargins(object = m, basesub = ps, delta = 2))
# # 
# #   ## basesub has the same names as parts in cilr
# #   ps <- build.basesub(parts = c("Sleep", "WAKE", "MVPA", "LPA", "SB"))
# #   expect_error(x <- bsubmargins(object = m, basesub = ps, delta = 2))
# # })
# 
# test_that("bsubmargins errors when delta out of range", {
#   expect_error(x <- bsubmargins(object = m, basesub = psub, delta = 1000))
# })
# test_that("bsubmargins outputs what expected", {
#   
#   ## types
#   expect_type(x, "list")
#   expect_equal(length(x), length(m$complr$parts))
#   expect_s3_class(x$TST, "data.table")
#   expect_s3_class(x$WAKE, "data.table")
#   expect_s3_class(x$MVPA, "data.table")
#   expect_s3_class(x$LPA, "data.table")
#   expect_s3_class(x$SB, "data.table")
#   
#   expect_type(x$TST$Mean, "double")
#   expect_type(x$TST$CI_low, "double")
#   expect_type(x$TST$CI_high, "double")
#   expect_type(x$TST$Delta, "double")
#   expect_type(x$TST$From, "character")
#   expect_type(x$TST$To, "character")
#   
#   expect_type(x$WAKE$Mean, "double")
#   expect_type(x$WAKE$CI_low, "double")
#   expect_type(x$WAKE$CI_high, "double")
#   expect_type(x$WAKE$Delta, "double")
#   expect_type(x$WAKE$From, "character")
#   expect_type(x$WAKE$To, "character")
#   
#   expect_type(x$MVPA$Mean, "double")
#   expect_type(x$MVPA$CI_low, "double")
#   expect_type(x$MVPA$CI_high, "double")
#   expect_type(x$MVPA$Delta, "double")
#   expect_type(x$MVPA$From, "character")
#   expect_type(x$MVPA$To, "character")
#   
#   expect_type(x$LPA$Mean, "double")
#   expect_type(x$LPA$CI_low, "double")
#   expect_type(x$LPA$CI_high, "double")
#   expect_type(x$LPA$Delta, "double")
#   expect_type(x$LPA$From, "character")
#   expect_type(x$LPA$To, "character")
#   
#   expect_type(x$SB$Mean, "double")
#   expect_type(x$SB$CI_low, "double")
#   expect_type(x$SB$CI_high, "double")
#   expect_type(x$SB$Delta, "double")
#   expect_type(x$SB$From, "character")
#   expect_type(x$SB$To, "character")
#   
#   expect_true(ncol(x$TST) == 8)
#   expect_true(ncol(x$WAKE) == 8)
#   expect_true(ncol(x$MVPA) == 8)
#   expect_true(ncol(x$LPA) == 8)
#   expect_true(ncol(x$SB) == 8)
#   
#   expect_true(all(x$TST$To == "TST"))
#   expect_true(all(x$WAKE$To == "WAKE"))
#   expect_true(all(x$MVPA$To == "MVPA"))
#   expect_true(all(x$LPA$To == "LPA"))
#   expect_true(all(x$SB$To == "SB"))
#   
#   expect_true(all(x$TST$Level == "between"))
#   expect_true(all(x$WAKE$Level == "between"))
#   expect_true(all(x$MVPA$Level == "between"))
#   expect_true(all(x$LPA$Level == "between"))
#   expect_true(all(x$SB$Level == "between"))
#   
# })
# 
# test_that("bsubmargins gives results in sensible range", {
#   
#   ## difference in outcome
#   expect_true(x$TST$Mean %ae% "[-0.5, 0) | (0, 0.5]")
#   expect_true(x$WAKE$Mean %ae% "[-0.5, 0) | (0, 0.5]")
#   expect_true(x$MVPA$Mean %ae% "[-0.5, 0) | (0, 0.5]")
#   expect_true(x$LPA$Mean %ae% "[-0.5, 0) | (0, 0.5]")
#   expect_true(x$SB$Mean %ae% "[-0.5, 0) | (0, 0.5]")
#   
#   expect_true(x$TST$CI_low %ae% "[-1, 0) | (0, 1]")
#   expect_true(x$WAKE$CI_low %ae% "[-1, 0) | (0, 1]")
#   expect_true(x$MVPA$CI_low %ae% "[-1, 0) | (0, 1]")
#   expect_true(x$LPA$CI_low %ae% "[-1, 0) | (0, 1]")
#   expect_true(x$SB$CI_low %ae% "[-1, 0) | (0, 1]")
#   
#   expect_true(x$TST$CI_high %ae% "[-1, 0) | (0, 1]")
#   expect_true(x$WAKE$CI_high %ae% "[-1, 0) | (0, 1]")
#   expect_true(x$MVPA$CI_high %ae% "[-1, 0) | (0, 1]")
#   expect_true(x$LPA$CI_high %ae% "[-1, 0) | (0, 1]")
#   expect_true(x$SB$CI_high %ae% "[-1, 0) | (0, 1]")
#   
# })
# 
# test_that("bsubmargins gives results in expected direction and magnitude", {
#   
#   ## values are opposite sign for opposite substitution
#   # for (i in seq_along(x)) {
#   #   expect_true(all(x[[i]][, sign(Mean[sign(Delta) == 1]) 
#   #                          %a!=% sign(Mean[sign(Delta) == -1]), by = From]$V1))
#   # }
#   
#   ## results for 1 min have smaller magnitude than 2 mins
#   for (i in seq_along(x)) {
#     expect_true(all(x[[i]][, abs(Mean[abs(Delta) == 1]) 
#                            < abs(Mean[abs(Delta) == 2])]))
#   }
# })
# 
# #---------------------------------------------------------------------------------------------------
# # Test 2-component composition for consistency between brm model and substitution model
# # using results from pairwise substitution
# ## Estimates should be in the direction between pairwise coordinates and  pairwise substitution 
# ## CIs should indicate consistent significance between pairwise coordinates and substitution 
# 
# ## TST vs WAKE
# test_that("bsubmargins's results matches with brm for 2-component composition (TST vs WAKE)", {
#   
#   sbp <- as.matrix(data.table(1, -1))
#   cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                   parts = c("TST", "WAKE"), idvar = "ID", total = 1440)
#   psub <- build.basesub(c("TST", "WAKE"))
#   suppressWarnings(
#     m <- brmcoda(complr = cilr,
#                  formula = Stress ~ bilr1 + wilr1 + (1 | ID),
#                  chain = 1, iter = 500, seed = 123,
#                  backend = backend))
#   a <- bsubmargins(object = m, basesub = psub, delta = 1:2)
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[2, 1] > 0))) { 
#     expect_true(all(a$TST[From == "WAKE" & Delta > 1]$Mean > 0)) 
#     expect_true(all(a$WAKE[From == "TST" & Delta > 1]$Mean < 0)) 
#   } else {
#     expect_true(all(a$TST[From == "WAKE" & Delta > 1]$Mean < 0))
#     expect_true(all(a$WAKE[From == "TST" & Delta > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[2, 3], summary(m$model)$fixed[2, 4]))
#     == (0 %agele% c(a$TST[From == "WAKE" & Delta == 1]$CI_low,
#                     a$TST[From == "WAKE" & Delta == 1]$CI_high))))
#   
# })
# 
# ## TST vs MVPA
# test_that("bsubmargins's results matches with brm for 2-component composition (TST vs MVPA)", {
#   
#   sbp <- as.matrix(data.table(1, -1))
#   cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                   parts = c("TST", "MVPA"), idvar = "ID", total = 1440)
#   psub <- build.basesub(c("TST", "MVPA"))
#   suppressWarnings(
#     m <- brmcoda(complr = cilr,
#                  formula = Stress ~ bilr1 + wilr1 + (1 | ID),
#                  chain = 1, iter = 500, seed = 123,
#                  backend = backend))
#   b <- bsubmargins(object = m, basesub = psub, delta = 1:2)
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[2, 1] > 0))) { 
#     expect_true(all(b$TST[From == "MVPA" & Delta > 1]$Mean > 0)) 
#     expect_true(all(b$MVPA[From == "TST" & Delta > 1]$Mean < 0)) 
#   } else {
#     expect_true(all(b$TST[From == "MVPA" & Delta > 1]$Mean < 0))
#     expect_true(all(b$MVPA[From == "TST" & Delta > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[2, 3], summary(m$model)$fixed[2, 4]))
#     == (0 %agele% c(b$TST[From == "MVPA" & Delta == 1]$CI_low,
#                     b$TST[From == "MVPA" & Delta == 1]$CI_high))))
#   
# })
# 
# ## TST vs LPA
# test_that("bsubmargins's results matches with brm model for 2-component composition (TST vs LPA)", {
# 
#   sbp <- as.matrix(data.table(1, -1))
#   cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                   parts = c("TST", "LPA"), idvar = "ID", total = 1440)
#   psub <- build.basesub(c("TST", "LPA"))
#   suppressWarnings(
#     m <- brmcoda(complr = cilr,
#                  formula = Stress ~ bilr1 + wilr1 + (1 | ID),
#                  chain = 1, iter = 500, seed = 123,
#                  backend = backend))
#   c <- bsubmargins(object = m, basesub = psub, delta = 1:2)
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[2, 1] > 0))) { 
#     expect_true(all(c$TST[From == "LPA" & Delta > 1]$Mean > 0)) 
#     expect_true(all(c$LPA[From == "TST" & Delta > 1]$Mean < 0)) 
#   } else {
#     expect_true(all(c$TST[From == "LPA" & Delta > 1]$Mean < 0))
#     expect_true(all(c$LPA[From == "TST" & Delta > 1]$Mean > 0))
#   }
# 
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[2, 3], summary(m$model)$fixed[2, 4]))
#     == (0 %agele% c(c$TST[From == "LPA" & Delta == 1]$CI_low,
#                     c$TST[From == "LPA" & Delta == 1]$CI_high))))
# 
# })
# 
# ## TST vs SB
# test_that("bsubmargins's results matches with brm model for 2-component composition (TST vs SB)", {
#   
#   sbp <- as.matrix(data.table(1, -1))
#   cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                   parts = c("TST", "SB"), idvar = "ID", total = 1440)
#   psub <- build.basesub(c("TST", "SB"))
#   suppressWarnings(
#     m <- brmcoda(complr = cilr,
#                  formula = Stress ~ bilr1 + wilr1 + (1 | ID),
#                  chain = 1, iter = 500, seed = 123,
#                  backend = backend))
#   d <- bsubmargins(object = m, basesub = psub, delta = 1:2)
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[2, 1] > 0))) { 
#     expect_true(all(d$TST[From == "SB" & Delta > 1]$Mean > 0)) 
#     expect_true(all(d$SB[From == "TST" & Delta > 1]$Mean < 0)) 
#   } else {
#     expect_true(all(d$TST[From == "SB" & Delta > 1]$Mean < 0))
#     expect_true(all(d$SB[From == "TST" & Delta > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[2, 3], summary(m$model)$fixed[2, 4]))
#     == (0 %agele% c(d$TST[From == "SB" & Delta == 1]$CI_low,
#                     d$TST[From == "SB" & Delta == 1]$CI_high))))
#   
# })
# 
# ## WAKE vs MVPA
# test_that("bsubmargins's results matches with brm for 2-component composition (WAKE vs MVPA)", {
#   
#   sbp <- as.matrix(data.table(1, -1))
#   cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                   parts = c("WAKE", "MVPA"), idvar = "ID", total = 1440)
#   psub <- build.basesub(c("WAKE", "MVPA"))
#   
#   suppressWarnings(
#     m <- brmcoda(complr = cilr,
#                  formula = Stress ~ bilr1 + wilr1 + (1 | ID),
#                  chain = 1, iter = 500, seed = 123,
#                  backend = backend))
#   e <- bsubmargins(object = m, basesub = psub, delta = 1:2)
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[2, 1] > 0))) { 
#     expect_true(all(e$WAKE[From == "MVPA" & Delta > 1]$Mean > 0)) 
#     expect_true(all(e$MVPA[From == "WAKE" & Delta > 1]$Mean < 0)) 
#   } else {
#     expect_true(all(e$WAKE[From == "MVPA" & Delta > 1]$Mean < 0))
#     expect_true(all(e$MVPA[From == "WAKE" & Delta > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[2, 3], summary(m$model)$fixed[2, 4]))
#     == (0 %agele% c(e$WAKE[From == "MVPA" & Delta == 1]$CI_low,
#                     e$WAKE[From == "MVPA" & Delta == 1]$CI_high))))
#   
# })
# 
# ## WAKE vs LPA
# test_that("bsubmargins's results matches with brm for 2-component composition (WAKE vs LPA)", {
#   
#   sbp <- as.matrix(data.table(1, -1))
#   cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                   parts = c("WAKE", "LPA"), idvar = "ID", total = 1440)
#   psub <- build.basesub(c("WAKE", "LPA"))
#   suppressWarnings(
#     m <- brmcoda(complr = cilr,
#                  formula = Stress ~ bilr1 + wilr1 + (1 | ID),
#                  chain = 1, iter = 500, seed = 123,
#                  backend = backend))
#   f <- bsubmargins(object = m, basesub = psub, delta = 1:2)
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[2, 1] > 0))) { 
#     expect_true(all(f$WAKE[From == "LPA" & Delta > 1]$Mean > 0)) 
#     expect_true(all(f$LPA[From == "WAKE" & Delta > 1]$Mean < 0)) 
#   } else {
#     expect_true(all(f$WAKE[From == "LPA" & Delta > 1]$Mean < 0))
#     expect_true(all(f$LPA[From == "WAKE" & Delta > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[2, 3], summary(m$model)$fixed[2, 4]))
#     == (0 %agele% c(f$WAKE[From == "LPA" & Delta == 1]$CI_low,
#                     f$WAKE[From == "LPA" & Delta == 1]$CI_high))))
#   
# })
# 
# ## WAKE vs SB
# test_that("bsubmargins's results matches with brm model for 2-component composition (WAKE vs SB)", {
#   
#   sbp <- as.matrix(data.table(1, -1))
#   cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                   parts = c("WAKE", "SB"), idvar = "ID", total = 1440)
#   psub <- build.basesub(c("WAKE", "SB"))
#   suppressWarnings(
#     m <- brmcoda(complr = cilr,
#                  formula = Stress ~ bilr1 + wilr1 + (1 | ID),
#                  chain = 1, iter = 500, seed = 123,
#                  backend = backend))
#   g <- bsubmargins(object = m, basesub = psub, delta = 1:2)
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[2, 1] > 0))) { 
#     expect_true(all(g$WAKE[From == "SB" & Delta > 1]$Mean > 0)) 
#     expect_true(all(g$SB[From == "WAKE" & Delta > 1]$Mean < 0)) 
#   } else {
#     expect_true(all(g$WAKE[From == "SB" & Delta > 1]$Mean < 0))
#     expect_true(all(g$SB[From == "WAKE" & Delta > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[2, 3], summary(m$model)$fixed[2, 4]))
#     == (0 %agele% c(g$WAKE[From == "SB" & Delta == 1]$CI_low,
#                     g$WAKE[From == "SB" & Delta == 1]$CI_high))))
#   
# })
# 
# ## MVPA vs LPA
# test_that("bsubmargins's results matches with brm for 2-component composition (MVPA vs LPA)", {
#   
#   sbp <- as.matrix(data.table(1, -1))
#   cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                   parts = c("MVPA", "LPA"), idvar = "ID", total = 1440)
#   psub <- build.basesub(c("MVPA", "LPA"))
#   suppressWarnings(
#     m <- brmcoda(complr = cilr,
#                  formula = Stress ~ bilr1 + wilr1 + (1 | ID),
#                  chain = 1, iter = 500, seed = 123,
#                  backend = backend))
#   h <- bsubmargins(object = m, basesub = psub, delta = 1:2)
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[2, 1] > 0))) { 
#     expect_true(all(h$MVPA[From == "LPA" & Delta > 1]$Mean > 0)) 
#     expect_true(all(h$LPA[From == "MVPA" & Delta > 1]$Mean < 0)) 
#   } else {
#     expect_true(all(h$MVPA[From == "LPA" & Delta > 1]$Mean < 0))
#     expect_true(all(h$LPA[From == "MVPA" & Delta > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[2, 3], summary(m$model)$fixed[2, 4]))
#     == (0 %agele% c(h$MVPA[From == "LPA" & Delta == 1]$CI_low,
#                     h$MVPA[From == "LPA" & Delta == 1]$CI_high))))
#   
# })
# 
# ## MVPA vs SB
# test_that("bsubmargins's results matches with brm model for 2-component composition (MVPA vs SB)", {
#   
#   sbp <- as.matrix(data.table(1, -1))
#   cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                   parts = c("MVPA", "SB"), idvar = "ID", total = 1440)
#   psub <- build.basesub(c("MVPA", "SB"))
#   suppressWarnings(
#     m <- brmcoda(complr = cilr,
#                  formula = Stress ~ bilr1 + wilr1 + (1 | ID),
#                  chain = 1, iter = 500, seed = 123,
#                  backend = backend))
#   i <- bsubmargins(object = m, basesub = psub, delta = 1:2)
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[2, 1] > 0))) { 
#     expect_true(all(i$MVPA[From == "SB" & Delta > 1]$Mean > 0)) 
#     expect_true(all(i$SB[From == "MVPA" & Delta > 1]$Mean < 0)) 
#   } else {
#     expect_true(all(i$MVPA[From == "SB" & Delta > 1]$Mean < 0))
#     expect_true(all(i$SB[From == "MVPA" & Delta > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[2, 3], summary(m$model)$fixed[2, 4]))
#     == (0 %agele% c(i$MVPA[From == "SB" & Delta == 1]$CI_low,
#                     i$MVPA[From == "SB" & Delta == 1]$CI_high))))
#   
# })
# 
# ## LPA vs SB
# test_that("bsubmargins's results matches with brm model for 2-component composition (LPA vs SB)", {
#   
#   sbp <- as.matrix(data.table(1, -1))
#   cilr <- complr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
#                   parts = c("LPA", "SB"), idvar = "ID", total = 1440)
#   psub <- build.basesub(c("LPA", "SB"))
#   suppressWarnings(
#     m <- brmcoda(complr = cilr,
#                  formula = Stress ~ bilr1 + wilr1 + (1 | ID),
#                  chain = 1, iter = 500, seed = 123,
#                  backend = backend))
#   j <- bsubmargins(object = m, basesub = psub, delta = 1:2)
#   
#   ## Estimates
#   if (isTRUE(suppressWarnings(summary(m$model)$fixed[2, 1] > 0))) { 
#     expect_true(all(j$LPA[From == "SB" & Delta > 1]$Mean > 0)) 
#     expect_true(all(j$SB[From == "LPA" & Delta > 1]$Mean < 0)) 
#   } else {
#     expect_true(all(j$LPA[From == "SB" & Delta > 1]$Mean < 0))
#     expect_true(all(j$SB[From == "LPA" & Delta > 1]$Mean > 0))
#   }
#   
#   # CIs
#   suppressWarnings(expect_true(
#     (0 %gele% c(summary(m$model)$fixed[2, 3], summary(m$model)$fixed[2, 4]))
#     == (0 %agele% c(j$LPA[From == "SB" & Delta == 1]$CI_low,
#                     j$LPA[From == "SB" & Delta == 1]$CI_high))))
#   
# })
