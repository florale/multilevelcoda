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

library(testthat)
library(data.table)
library(multilevelcoda)
library(brms)

data(mcompd)
data(sbp)

cilr <- complr(
  data = mcompd[ID %in% 1:10, .SD[1:3], by = ID],
  sbp = sbp,
  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
  idvar = "ID"
)

suppressWarnings(
  fit <- brmcoda(
    complr = cilr,
    formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
      wz1_1 + wz2_1 + wz3_1 + wz4_1 + Female + (1 | ID),
    chain = 1,
    iter = 500,
    seed = 123,
    backend = backend,
    silent = 2,
    refresh = 0
  )
)
foreach::registerDoSEQ()

newdata_subset <- mcompd[ID %in% 1:5, .SD[1:2], by = ID]

test_that("update.brmcoda allows a plain refresh update", {
  fit_updated <- update(fit, silent = 2, refresh = 0)

  expect_s3_class(fit_updated, "brmcoda")
  expect_equal(nobs(fit_updated), nobs(fit))
  expect_equal(variables(fit_updated), variables(fit))
})

test_that("update.brmcoda updates the formula", {
  fit_newformula <- update(fit, formula. = ~ . - wz1_1, silent = 2, refresh = 0)

  expect_s3_class(fit_newformula, "brmcoda")
  expect_false(any(grepl("b_wz1_1$", variables(fit_newformula))))
  expect_equal(nobs(fit_newformula), nobs(fit))
})

test_that("update.brmcoda updates the data", {
  fit_newdata <- update(fit, newdata = newdata_subset, silent = 2, refresh = 0)

  expect_s3_class(fit_newdata, "brmcoda")
  expect_equal(fit_newdata$complr$datain, as.data.table(newdata_subset))
  expect_equal(nrow(fit_newdata$complr$dataout), nrow(newdata_subset))
})

test_that("update.brmcoda updates both formula and data", {
  fit_new <- update(
    fit,
    formula. = ~ . - wz2_1,
    newdata = newdata_subset,
    silent = 2,
    refresh = 0
  )

  expect_s3_class(fit_new, "brmcoda")
  expect_false(any(grepl("b_wz2_1$", variables(fit_new))))
  expect_equal(fit_new$complr$datain, as.data.table(newdata_subset))
  expect_equal(nrow(fit_new$complr$dataout), nrow(newdata_subset))
})
