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

parts <- c("TST", "WAKE", "MVPA", "LPA", "SB")

data(mcompd)
data(sbp)
data(psub)

cilr <- complr(
  data = mcompd[ID %in% 1:10, .SD[1:3], by = ID],
  sbp = sbp,
  parts = parts,
  idvar = "ID",
  total = 1440
)

suppressWarnings(
  m <- brmcoda(
    complr = cilr,
    formula = Stress ~ z1_1 + z2_1 + z3_1 + z4_1 + (1 | ID),
    chain = 1,
    iter = 500,
    seed = 123,
    backend = backend,
    silent = 2,
    refresh = 0
  )
)
foreach::registerDoSEQ()

x <- submargin(object = m, base = psub, delta = 1:2)

stress_result <- function(x, part) {
  x[[part]][["Stress"]]
}

fit_pairwise_model <- function(pair) {
  sbp_pair <- as.matrix(data.table(1, -1))
  cilr_pair <- complr(
    data = mcompd[ID %in% 1:10, .SD[1:3], by = ID],
    sbp = sbp_pair,
    parts = pair,
    idvar = "ID",
    total = 1440
  )

  suppressWarnings(
    brmcoda(
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
}

expect_pairwise_direction_and_ci <- function(pair) {
  fit_pair <- fit_pairwise_model(pair)
  pair_result <- submargin(
    object = fit_pair,
    base = build.base(pair),
    delta = 1:2
  )

  to_first <- pair_result[[pair[1]]][["Stress"]]
  to_second <- pair_result[[pair[2]]][["Stress"]]
  coef_summary <- fixef(fit_pair$model, summary = TRUE)
  coef_estimate <- coef_summary["z1_1", "Estimate"]
  coef_has_zero <- coef_summary["z1_1", "Q2.5"] <= 0 &&
    0 <= coef_summary["z1_1", "Q97.5"]

  first_ci <- to_first[From == pair[2] & Delta == 1, .(CI_low, CI_high)]
  first_has_zero <- first_ci$CI_low <= 0 && 0 <= first_ci$CI_high

  if (isTRUE(coef_estimate > 0)) {
    expect_true(all(to_first[From == pair[2] & Delta > 1]$Estimate > 0))
    expect_true(all(to_second[From == pair[1] & Delta > 1]$Estimate < 0))
  } else {
    expect_true(all(to_first[From == pair[2] & Delta > 1]$Estimate < 0))
    expect_true(all(to_second[From == pair[1] & Delta > 1]$Estimate > 0))
  }

  expect_identical(coef_has_zero, first_has_zero)
}

test_that("submargin errors for missing required arguments", {
  expect_error(submargin(base = psub, delta = 2))
  expect_error(submargin(object = m, delta = 2))
  expect_error(submargin(object = m, base = psub))
})

test_that("submargin returns the expected nested result structure", {
  expect_type(x, "list")
  expect_named(x, parts)

  for (part in parts) {
    tmp <- stress_result(x, part)

    expect_type(x[[part]], "list")
    expect_named(x[[part]], "Stress")
    expect_s3_class(tmp, "data.table")
    expect_named(
      tmp,
      c("Estimate", "Est.Error", "CI_low", "CI_high",
        "Delta", "From", "To", "Level", "Reference")
    )
    expect_equal(nrow(tmp), 16L)
    expect_true(all(tmp$To == part))
    expect_true(all(tmp$Level == "aggregate"))
    expect_true(all(tmp$Reference == "clustermean"))
    expect_setequal(unique(tmp$From), setdiff(parts, part))
    expect_setequal(unique(tmp$Delta), c(-2, -1, 1, 2))
    expect_type(tmp$Estimate, "double")
    expect_type(tmp$Est.Error, "double")
    expect_type(tmp$CI_low, "double")
    expect_type(tmp$CI_high, "double")
    expect_type(tmp$Delta, "double")
    expect_type(tmp$From, "character")
    expect_type(tmp$To, "character")
  }
})

test_that("submargin returns estimates in a sensible range", {
  for (part in parts) {
    tmp <- stress_result(x, part)

    expect_true(all(abs(tmp$Estimate) <= 0.5))
    expect_true(all(abs(tmp$CI_low) <= 1))
    expect_true(all(abs(tmp$CI_high) <= 1))
  }
})

test_that("submargin estimates grow in magnitude from 1 to 2 minutes", {
  for (part in parts) {
    tmp <- stress_result(x, part)

    opposite_check <- tmp[
      ,
      .(opposite = sign(Estimate[Delta == 1]) != sign(Estimate[Delta == -1])),
      by = From
    ]
    magnitude_check <- tmp[
      ,
      .(smaller = abs(Estimate[abs(Delta) == 1]) < abs(Estimate[abs(Delta) == 2])),
      by = From
    ]

    expect_true(all(opposite_check$opposite))
    expect_true(all(magnitude_check$smaller))
  }
})

pairwise_parts <- combn(parts, 2, simplify = FALSE)

for (pair in pairwise_parts) {
  test_that(
    sprintf(
      "submargin matches brm for 2-component composition (%s vs %s)",
      pair[1], pair[2]
    ),
    {
      expect_pairwise_direction_and_ci(pair)
    }
  )
}
