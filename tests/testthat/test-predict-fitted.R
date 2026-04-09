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
library(withr)

data(mcompd)
data(sbp)

parts <- c("TST", "WAKE", "MVPA", "LPA", "SB")
total_time <- 1440

cilr <- complr(
  data = mcompd[ID %in% 1:10, .SD[1:3], by = ID],
  sbp = sbp,
  parts = parts,
  idvar = "ID",
  total = total_time
)

suppressWarnings(
  fit_scalar <- brmcoda(
    complr = cilr,
    formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
      wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
    chain = 1,
    iter = 500,
    seed = 123,
    backend = backend,
    silent = 2,
    refresh = 0
  )
)

suppressWarnings(
  fit_comp <- brmcoda(
    complr = cilr,
    formula = mvbind(z1_1, z2_1, z3_1, z4_1) ~ Stress + Female + (1 | ID),
    chain = 1,
    iter = 500,
    seed = 123,
    backend = backend,
    silent = 2,
    refresh = 0
  )
)

foreach::registerDoSEQ()

newdata_subset <- cilr$dataout[ID %in% 1:5, .SD[1:2], by = ID]

test_that("predict.brmcoda matches brms predict for scalar responses", {
  out <- with_seed(123, predict(fit_scalar))
  expected <- with_seed(123, predict(fit_scalar$model))

  expect_equal(out, expected)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(nrow(cilr$dataout), 4L))
  expect_identical(colnames(out), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
})

test_that("fitted.brmcoda matches brms fitted for scalar responses", {
  out <- fitted(fit_scalar)
  expected <- fitted(fit_scalar$model)

  expect_equal(out, expected)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(nrow(cilr$dataout), 4L))
  expect_identical(colnames(out), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
})

test_that("predict.brmcoda and fitted.brmcoda respect newdata for scalar responses", {
  pred <- predict(fit_scalar, newdata = newdata_subset)
  fitd <- fitted(fit_scalar, newdata = newdata_subset)

  expect_equal(nrow(pred), nrow(newdata_subset))
  expect_equal(nrow(fitd), nrow(newdata_subset))
  expect_identical(colnames(pred), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  expect_identical(colnames(fitd), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
})

test_that("linear-scale predict and fitted match brms for compositional responses", {
  pred <- with_seed(123, predict(fit_comp, scale = "linear"))
  pred_expected <- with_seed(123, predict(fit_comp$model))
  fitd <- fitted(fit_comp, scale = "linear")
  fitd_expected <- fitted(fit_comp$model)

  expect_equal(pred, pred_expected)
  expect_equal(fitd, fitd_expected)
  expect_equal(dim(pred), c(nrow(cilr$dataout), 4L, 4L))
  expect_equal(dim(fitd), c(nrow(cilr$dataout), 4L, 4L))
  expect_identical(dimnames(pred)[[2]], c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  expect_identical(dimnames(fitd)[[2]], c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  expect_identical(dimnames(pred)[[3]], c("z11", "z21", "z31", "z41"))
  expect_identical(dimnames(fitd)[[3]], c("z11", "z21", "z31", "z41"))
})

test_that("predict.brmcoda back-transforms compositional responses to closed compositions", {
  out <- predict(fit_comp, scale = "response")
  closure <- apply(out, c(1, 2), sum)

  expect_type(out, "double")
  expect_equal(dim(out), c(nrow(cilr$dataout), 4L, 5L))
  expect_identical(dimnames(out)[[3]], parts)
  expect_equal(
    unname(closure),
    matrix(total_time, nrow = nrow(closure), ncol = ncol(closure)),
    tolerance = 1e-6
  )
})

test_that("fitted.brmcoda back-transforms compositional responses on response scale", {
  out <- fitted(fit_comp, scale = "response")

  expect_type(out, "double")
  expect_equal(dim(out), c(nrow(cilr$dataout), 4L, 5L))
  expect_identical(dimnames(out)[[2]], c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  expect_identical(dimnames(out)[[3]], paste0("t", parts))
})

test_that("fitted.brmcoda summary = FALSE returns raw compositional draws", {
  out <- fitted(fit_comp, scale = "response", summary = FALSE)
  closure <- apply(out, c(1, 2), sum)

  expect_type(out, "double")
  expect_equal(dim(out), c(250L, nrow(cilr$dataout), 5L))
  expect_identical(dimnames(out)[[3]], paste0("t", parts))
  expect_equal(
    unname(closure),
    matrix(total_time, nrow = nrow(closure), ncol = ncol(closure)),
    tolerance = 1e-6
  )
})

test_that("predict.brmcoda and fitted.brmcoda respect newdata for compositional responses", {
  pred <- predict(fit_comp, scale = "response", newdata = newdata_subset)
  fitd <- fitted(fit_comp, scale = "response", newdata = newdata_subset)
  fitd_raw <- fitted(fit_comp, scale = "response", newdata = newdata_subset, summary = FALSE)
  pred_closure <- apply(pred, c(1, 2), sum)
  fitd_raw_closure <- apply(fitd_raw, c(1, 2), sum)

  expect_equal(dim(pred), c(nrow(newdata_subset), 4L, 5L))
  expect_equal(dim(fitd), c(nrow(newdata_subset), 4L, 5L))
  expect_equal(dim(fitd_raw), c(250L, nrow(newdata_subset), 5L))
  expect_equal(
    unname(pred_closure),
    matrix(total_time, nrow = nrow(pred_closure), ncol = ncol(pred_closure)),
    tolerance = 1e-6
  )
  expect_equal(
    unname(fitd_raw_closure),
    matrix(total_time, nrow = nrow(fitd_raw_closure), ncol = ncol(fitd_raw_closure)),
    tolerance = 1e-6
  )
})
