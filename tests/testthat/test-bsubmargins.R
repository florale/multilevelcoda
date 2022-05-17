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

####################################################################################################
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

x <- bsubmargins(object = m, substitute = psub, minute = 2)

test_that("bsubmargins errors for invalid input", {
  
  ## check errors for missing object
  expect_error(x <- bsubmargins(substitute = psub, minute = 2))
  
  ## check errors for missing substitute
  expect_error(x <- bsubmargins(objetc = m, minute = 2))
  
})

test_that("bsubmargins outputs what expected", {
  
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

test_that("bsubmargins gives results in sensible range", {
  
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

test_that("bsubmargins gives results in expected direction and magnitude", {
  
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

####################################################################################################
# Test 2-component composition for consistency between brm model and substitution model
# using results from pairwise substitution
## Estimates should be in the direction between pairwise coordinates and  pairwise substitution 
## CIs should indicate consistent significance between pairwise coordinates and  substitution 

## TST vs WAKE
sbp <- as.matrix(data.table(1, -1))
cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("TST", "WAKE"), idvar = "ID")
psub <- possub(c("TST", "WAKE"))
suppressWarnings(m <- brmcoda(compilr = cilr,
                              formula = STRESS ~ bilr1 + wilr1 + (1 | ID),
                              chain = 1, iter = 500, seed = 123))
a <- bsubmargins(object = m, substitute = psub, minute = 2)

test_that("bsubmargins's results matches with brm model for 2-component composition (TST vs WAKE)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) { 
    expect_true(all(a$TST[Substitute == "WAKE" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(a$WAKE[Substitute == "TST" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(a$TST[Substitute == "WAKE" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(a$WAKE[Substitute == "TST" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
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
b <- bsubmargins(object = m, substitute = psub, minute = 2)

test_that("bsubmargins's results matches with brm model for 2-component composition (TST vs MVPA)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) { 
    expect_true(all(b$TST[Substitute == "MVPA" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(b$MVPA[Substitute == "TST" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(b$TST[Substitute == "MVPA" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(b$MVPA[Substitute == "TST" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
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
c <- bsubmargins(object = m, substitute = psub, minute = 2)

test_that("bsubmargins's results matches with brm model for 2-component composition (TST vs LPA)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) { 
    expect_true(all(c$TST[Substitute == "LPA" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(c$LPA[Substitute == "TST" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(c$TST[Substitute == "LPA" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(c$LPA[Substitute == "TST" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
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
d <- bsubmargins(object = m, substitute = psub, minute = 2)

test_that("bsubmargins's results matches with brm model for 2-component composition (TST vs SB)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) { 
    expect_true(all(d$TST[Substitute == "SB" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(d$SB[Substitute == "TST" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(d$TST[Substitute == "SB" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(d$SB[Substitute == "TST" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
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
e <- bsubmargins(object = m, substitute = psub, minute = 2)

test_that("bsubmargins's results matches with brm model for 2-component composition (WAKE vs MVPA)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) { 
    expect_true(all(e$WAKE[Substitute == "MVPA" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(e$MVPA[Substitute == "WAKE" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(e$WAKE[Substitute == "MVPA" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(e$MVPA[Substitute == "WAKE" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
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
f <- bsubmargins(object = m, substitute = psub, minute = 2)

test_that("bsubmargins's results matches with brm model for 2-component composition (WAKE vs LPA)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) { 
    expect_true(all(f$WAKE[Substitute == "LPA" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(f$LPA[Substitute == "WAKE" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(f$WAKE[Substitute == "LPA" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(f$LPA[Substitute == "WAKE" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
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
g <- bsubmargins(object = m, substitute = psub, minute = 2)

test_that("bsubmargins's results matches with brm model for 2-component composition (WAKE vs SB)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) { 
    expect_true(all(g$WAKE[Substitute == "SB" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(g$SB[Substitute == "WAKE" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(g$WAKE[Substitute == "SB" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(g$SB[Substitute == "WAKE" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
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
h <- bsubmargins(object = m, substitute = psub, minute = 2)

test_that("bsubmargins's results matches with brm model for 2-component composition (MVPA vs LPA)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) { 
    expect_true(all(h$MVPA[Substitute == "LPA" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(h$LPA[Substitute == "MVPA" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(h$MVPA[Substitute == "LPA" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(h$LPA[Substitute == "MVPA" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
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
i <- bsubmargins(object = m, substitute = psub, minute = 2)

test_that("bsubmargins's results matches with brm model for 2-component composition (MVPA vs SB)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) { 
    expect_true(all(i$MVPA[Substitute == "SB" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(i$SB[Substitute == "MVPA" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(i$MVPA[Substitute == "SB" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(i$SB[Substitute == "MVPA" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
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
j <- bsubmargins(object = m, substitute = psub, minute = 2)

test_that("bsubmargins's results matches with brm model for 2-component composition (LPA vs SB)", {
  
  ## Estimates
  if (isTRUE(suppressWarnings(summary(m$Model)$fixed[2, 1] > 0))) { 
    expect_true(all(j$LPA[Substitute == "SB" & MinSubstituted > 1]$Mean > 0)) 
    expect_true(all(j$SB[Substitute == "LPA" & MinSubstituted > 1]$Mean < 0)) 
  } else {
    expect_true(all(j$LPA[Substitute == "SB" & MinSubstituted > 1]$Mean < 0))
    expect_true(all(j$SB[Substitute == "LPA" & MinSubstituted > 1]$Mean > 0))
  }
  
  # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[2, 3], summary(m$Model)$fixed[2, 4]))
    == (0 %agele% c(j$LPA[Substitute == "SB" & MinSubstituted == 1]$CI_low,
                    j$LPA[Substitute == "SB" & MinSubstituted == 1]$CI_high))))
  
})
