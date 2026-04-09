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

parts_full <- c("TST", "WAKE", "MVPA", "LPA", "SB")
parts_pair <- c("TST", "WAKE")

cilr_full <- complr(
  data = mcompd[ID %in% 1:10, .SD[1:3], by = ID],
  sbp = sbp,
  parts = parts_full,
  idvar = "ID",
  total = 1440
)

suppressWarnings(
  fit_full <- brmcoda(
    complr = cilr_full,
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

sbp_pair <- as.matrix(data.table(1, -1))
cilr_pair <- complr(
  data = mcompd[ID %in% 1:10, .SD[1:3], by = ID],
  sbp = sbp_pair,
  parts = parts_pair,
  idvar = "ID",
  total = 1440
)

suppressWarnings(
  fit_pair <- brmcoda(
    complr = cilr_pair,
    formula = Stress ~ z1_1 + (1 | ID),
    chain = 1,
    iter = 500,
    seed = 123,
    backend = backend,
    silent = 2,
    refresh = 0
  )
)

foreach::registerDoSEQ()

pc_rotate <- pivot_coord_rotate(fit_full)
pc_rotate_raw <- pivot_coord_rotate(fit_full, summary = FALSE)
pc_rotate_generic <- pivot_coord(fit_full)

pc_refit <- pivot_coord_refit(fit_pair, silent = 2, refresh = 0)
pc_refit_raw <- pivot_coord_refit(
  fit_pair,
  summary = FALSE,
  silent = 2,
  refresh = 0
)
pc_refit_generic <- pivot_coord(
  fit_pair,
  method = "refit",
  summary = FALSE,
  silent = 2,
  refresh = 0
)

test_that("pivot_coord defaults to rotate", {
  expect_s3_class(pc_rotate_generic, "pivot_coord")
  expect_identical(pc_rotate_generic$method, "rotate")
  expect_equal(pc_rotate_generic$output, pc_rotate$output)
})

test_that("pivot_coord uses rotate when both methods are supplied", {
  tmp <- pivot_coord(fit_full, method = c("rotate", "refit"))

  expect_s3_class(tmp, "pivot_coord")
  expect_identical(tmp$method, "rotate")
  expect_equal(tmp$output, pc_rotate$output)
})

test_that("pivot_coord_rotate returns summarized pivot coordinates", {
  expect_s3_class(pc_rotate, "pivot_coord")
  expect_identical(pc_rotate$method, "rotate")
  expect_true(pc_rotate$summary)
  expect_named(pc_rotate$output, parts_full)

  for (part in parts_full) {
    tmp <- pc_rotate$output[[part]]

    expect_type(tmp, "double")
    expect_equal(dim(tmp), c(4L, 2L))
    expect_identical(
      rownames(tmp),
      c("Estimate", "Est.Error", "CI_low", "CI_high")
    )
    expect_identical(colnames(tmp), c("bz1_1", "wz1_1"))
  }
})

test_that("pivot_coord_rotate summary = FALSE returns draws", {
  expect_s3_class(pc_rotate_raw, "pivot_coord")
  expect_identical(pc_rotate_raw$method, "rotate")
  expect_false(pc_rotate_raw$summary)
  expect_named(pc_rotate_raw$output, parts_full)

  for (part in parts_full) {
    tmp <- pc_rotate_raw$output[[part]]

    expect_type(tmp, "double")
    expect_equal(nrow(tmp), 250L)
    expect_identical(colnames(tmp), c("bz1_1", "wz1_1"))
  }
})

test_that("summary.pivot_coord returns summarized output for rotate results", {
  tmp <- summary(pc_rotate)
  tmp_asis <- summary(pc_rotate, digits = "asis")

  expect_type(tmp, "list")
  expect_named(tmp, parts_full)
  expect_equal(tmp_asis, pc_rotate$output)

  for (part in parts_full) {
    expect_equal(dim(tmp[[part]]), c(4L, 2L))
    expect_identical(colnames(tmp[[part]]), c("bz1_1", "wz1_1"))
    expect_identical(
      rownames(tmp[[part]]),
      c("Estimate", "Est.Error", "CI_low", "CI_high")
    )
  }
})

test_that("summary.pivot_coord errors for raw draws", {
  expect_error(
    summary(pc_rotate_raw),
    "summary' method is currently not available for draws"
  )
})

test_that("pivot_coord_refit returns summarized pivot coordinates", {
  expect_s3_class(pc_refit, "pivot_coord")
  expect_identical(pc_refit$method, "refit")
  expect_true(pc_refit$summary)
  expect_named(pc_refit$output, parts_pair)

  for (part in parts_pair) {
    tmp <- pc_refit$output[[part]]

    expect_type(tmp, "double")
    expect_equal(dim(tmp), c(4L, 1L))
    expect_identical(
      rownames(tmp),
      c("Estimate", "Est.Error", "CI_low", "CI_high")
    )
    expect_identical(colnames(tmp), "z1_1")
  }
})

test_that("pivot_coord_refit summary = FALSE returns draws", {
  expect_s3_class(pc_refit_raw, "pivot_coord")
  expect_identical(pc_refit_raw$method, "refit")
  expect_false(pc_refit_raw$summary)
  expect_named(pc_refit_raw$output, parts_pair)

  for (part in parts_pair) {
    tmp <- pc_refit_raw$output[[part]]

    expect_type(tmp, "double")
    expect_equal(nrow(tmp), 250L)
    expect_identical(colnames(tmp), "z1_1")
  }
})

test_that("pivot_coord forwards summary to refit path", {
  expect_s3_class(pc_refit_generic, "pivot_coord")
  expect_identical(pc_refit_generic$method, "refit")
  expect_false(pc_refit_generic$summary)
  expect_named(pc_refit_generic$output, parts_pair)

  for (part in parts_pair) {
    tmp <- pc_refit_generic$output[[part]]

    expect_type(tmp, "double")
    expect_equal(nrow(tmp), 250L)
    expect_identical(colnames(tmp), "z1_1")
  }
})

test_that("summary.pivot_coord works for refit results", {
  tmp <- summary(pc_refit)

  expect_type(tmp, "list")
  expect_named(tmp, parts_pair)

  for (part in parts_pair) {
    expect_equal(dim(tmp[[part]]), c(4L, 1L))
    expect_identical(colnames(tmp[[part]]), "z1_1")
    expect_identical(
      rownames(tmp[[part]]),
      c("Estimate", "Est.Error", "CI_low", "CI_high")
    )
  }
})
