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

data(mcompd)
data(sbp)
data(psub)

cilr <- compilr(data = mcompd[ID %in% 1:10, .SD[1:3], by = ID], sbp = sbp,
                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")

suppressWarnings(
  m <- brmcoda(compilr = cilr, 
               formula = mvbind(ilr1, ilr2, ilr3, ilr4) ~ STRESS + (1 | ID),
               chain = 1, iter = 500, cores = 8,
               backend = backend))

test_that("emmcoda gives errors for invalid input", {

   ## missing object
   expect_error(x <- emmcoda(x = "STRESS"))

   ## invalid object
   expect_error(x <- emmcoda(object = cilr))

   ## only at, not x was provided
   expect_error(x <- emmcoda(object = m, at = 2))

   ## wrong x name
   expect_error(x <- emmcoda(object = m, x = "x"))

})

test_that("emmcoda uses and average across default reference grid when specified", {
   e <- emmcoda(object = m, x = "STRESS")
   expect_equal(e$ILR$STRESS$STRESS, as.data.table(ref_grid(m$Model) @grid)$STRESS)
   expect_equal(unique(e$Composition$STRESS$STRESS), 
                unique(as.data.table(ref_grid(m$Model) @grid)$STRESS))

})

test_that("emmcoda uses user's ref grid when specified", {
   e <- emmcoda(object = m, x = "STRESS", at = c(1:2))
   
   expect_true(all(c(e$ILR$STRESS1$STRESS, e$ILR$STRESS2$STRESS) %in% c(1, 2)))
   expect_true(all(c(e$Composition$STRESS1$STRESS, e$Composition$STRESS2$STRESS) %in% c(1,2)))
})

test_that("emmcoda outputs results in expected direction and signficance", {

   e <- emmcoda(object = m, x = "STRESS")
   
   # CIs
  suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[5, 3], summary(m$Model)$fixed[5, 4]))
    == (0 %agele% c(e$ILR$STRESS[1, 4], e$ILR$STRESS[1, 5]))))
  
    suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[6, 3], summary(m$Model)$fixed[6, 4]))
    == (0 %agele% c(e$ILR$STRESS[2, 4], e$ILR$STRESS[2, 5]))))
  
    suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[7, 3], summary(m$Model)$fixed[7, 4]))
    == (0 %agele% c(e$ILR$STRESS[3, 4], e$ILR$STRESS[3, 5]))))
  
    suppressWarnings(expect_true(
    (0 %gele% c(summary(m$Model)$fixed[8, 3], summary(m$Model)$fixed[8, 4]))
    == (0 %agele% c(e$ILR$STRESS[4, 4], e$ILR$STRESS[4, 5]))))
  
})
