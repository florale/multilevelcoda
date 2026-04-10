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

data(mcompd, package = "multilevelcoda")
data(sbp, package = "multilevelcoda")

x1 <- complr(
    data = mcompd, sbp = sbp,
    parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")

test_that("is.complr", {
  expect_true(is.complr(x1))
})

test_that("mean.complr", {
  expect_type(mean(x1), "list")
})

test_that("var.complr", {
  expect_type(compositions::var(x1), "list")
})

test_that("as.data.frame.complr", {
  expect_s3_class(as.data.frame(x1), "data.frame")
})

test_that("as.matrix.complr", {
  expect_type(as.matrix(x1), "double")
})

test_that("summary.complr", {
  output <- capture.output(x1_summary <- summary(x1))

  expect_s3_class(x1_summary, "summary.complr")
  expect_named(x1_summary, c("meta", "transform", "geometry", "output"))
  expect_named(x1_summary$meta, c("nobs", "ngrps", "idvar", "transform_type"))
  expect_named(x1_summary$transform, c("total", "psi", "sbp"))
  expect_named(
    x1_summary$geometry,
    c("composition_geometry", "logratio_class")
  )
  expect_s3_class(x1_summary$output, "data.frame")
  expect_identical(colnames(x1_summary$output), "Composition 1")
  expect_identical(
    rownames(x1_summary$output),
    c("X", "bX", "wX", "Z", "bZ", "wZ", "total")
  )
  expect_identical(x1_summary$meta$nobs, nrow(mcompd))
  expect_identical(x1_summary$meta$ngrps, length(unique(mcompd$ID)))
  expect_identical(x1_summary$meta$idvar, "ID")
  expect_identical(x1_summary$meta$transform_type, x1$transform)
  expect_equal(as.numeric(x1_summary$output["total", 1]), 1)
  expect_match(
    paste(output, collapse = "\n"),
    "Summary of Composition and Logratio Variables",
    fixed = TRUE
  )
  expect_match(
    paste(output, collapse = "\n"),
    "Number of groups \\(IDs\\): 266"
  )
})

test_that("print.complr", {
  summary_output <- capture.output(summary(x1))
  print_output <- capture.output(x1_print <- print(x1))

  expect_s3_class(x1_print, "summary.complr")
  expect_equal(print_output, summary_output)
})

## model with compositional predictor at between and within-person levels
m1 <- brmcoda(
    complr = x1,
    formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
        wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
    chain = 1, warmup = 250, iter = 500, seed = 123,
    backend = "cmdstanr", silent = 2, refresh = 0)

test_that("is.brmcoda", {
  expect_true(is.brmcoda(m1))
})

test_that("nobs.brmcoda", {
  expect_equal(nobs(m1), nrow(mcompd))
})

test_that("model.frame.brmcoda", {
  expect_s3_class(model.frame(m1), "data.frame")
})

test_that("variables.brmcoda", {
  expect_type(variables(m1), "character")
})

test_that("nvariables.brmcoda", {
  expect_equal(nvariables(m1), 280L)
})

test_that("niterations.brmcoda", {
  expect_equal(niterations(m1), 250L)
})

test_that("nchains.brmcoda", {
  expect_equal(nchains(m1), 1L)
})

test_that("ndraws.brmcoda", {
  expect_equal(ndraws(m1), 250L)
})

test_that("log_posterior.brmcoda", {
  expect_s3_class(log_posterior(m1), "data.frame")
})

test_that("nuts_params.brmcoda", {
  expect_s3_class(nuts_params(m1), "data.frame")
})

test_that("rhat.brmcoda", {
  tmpx <- rhat(m1)
  expect_true(all(tmpx > .9 & tmpx < 1.1))
})

test_that("neff_ratio.brmcoda", {
  suppressWarnings(tmpx <- neff_ratio(m1))
  expect_true(all(tmpx < 5 & tmpx > 0))
})

test_that("prior_summary.brmcoda", {
  expect_s3_class(prior_summary(m1), "data.frame")
})

test_that("summary.brmcoda", {
  m1_summary <- summary(m1)
  output <- capture.output(print(m1_summary))

  expect_s3_class(m1_summary, "brmssummary")
  expect_match(paste(output, collapse = "\n"), "Family: gaussian", fixed = TRUE)
  expect_match(
    paste(output, collapse = "\n"),
    "Number of observations: 3540",
    fixed = TRUE
  )
  expect_match(
    paste(output, collapse = "\n"),
    "total post-warmup draws = 250",
    fixed = TRUE
  )
})

test_that("print.brmcoda", {
  summary_output <- capture.output(print(summary(m1)))
  print_output <- capture.output(m1_print <- print(m1))

  expect_s3_class(m1_print, "brmssummary")
  expect_equal(print_output, summary_output)
})

test_that("pp_check.brmcoda", {
  tmpx <- pp_check(m1, type = "error_hist", ndraws = 10)
  expect_s3_class(tmpx, "ggplot")
})

test_that("fixef.brmcoda", {
  expect_type(fixef(m1), "double")
})

test_that("vcov.brmcoda", {
  expect_type(vcov(m1), "double")
})

test_that("ranef.brmcoda", {
  expect_type(ranef(m1), "list")
})

test_that("coef.brmcoda", {
  expect_type(VarCorr(m1), "list")
})

test_that("residuals.brmcoda", {
  expect_type(residuals(m1), "double")
})

## note loo and bayes factor methods not tested due to run length
