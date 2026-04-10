skip_on_cran()

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  backend <- "rstan"
  skip_on_os("windows")
} else {
  if (isFALSE(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))) {
    backend <- "cmdstanr"
  } else {
    backend <- "rstan"
  }
}

library(testthat)
library(data.table)
library(multilevelcoda)

parts <- c("TST", "WAKE", "MVPA", "LPA", "SB")

data(mcompd)
data(sbp)
data(psub)

fit_brmcoda <- function(complr_obj, formula) {
  suppressWarnings(
    brmcoda(
      complr = complr_obj,
      formula = formula,
      chain = 1,
      iter = 300,
      warmup = 150,
      seed = 123,
      backend = backend,
      silent = 2,
      refresh = 0
    )
  )
}

expect_standard_sub_output <- function(out, level, ref) {
  expect_type(out, "list")
  expect_named(out, parts)

  for (part in parts) {
    tmp <- out[[part]][["Stress"]]
    expect_type(out[[part]], "list")
    expect_named(out[[part]], "Stress")
    expect_s3_class(tmp, "data.table")
    expect_named(
      tmp,
      c("Estimate", "Est.Error", "CI_low", "CI_high",
        "Delta", "From", "To", "Level", "Reference")
    )
    expect_equal(nrow(tmp), 16L)
    expect_setequal(unique(tmp$Delta), c(-2, -1, 1, 2))
    expect_true(all(tmp$To == part))
    expect_true(all(tmp$Level == level))
    expect_true(all(tmp$Reference == ref))
  }
}

ml_data <- mcompd[ID %in% 1:10, .SD[1:3], by = ID]
wide_data <- unique(mcompd[, .(ID, Stress, Female, Age, TST, WAKE, MVPA, LPA, SB)])

cilr_ml <- complr(
  data = ml_data,
  sbp = sbp,
  parts = parts,
  idvar = "ID",
  total = 1440
)

cilr_wide <- complr(
  data = wide_data[, .(Stress, Female, Age, TST, WAKE, MVPA, LPA, SB)],
  sbp = sbp,
  parts = parts,
  total = 1440
)

fit_ml <- fit_brmcoda(
  cilr_ml,
  Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
    wz1_1 + wz2_1 + wz3_1 + wz4_1 + Female + Age + (1 | ID)
)

fit_nocomp <- fit_brmcoda(
  cilr_ml,
  Stress ~ Female + Age + (1 | ID)
)

fit_partial_bw <- fit_brmcoda(
  cilr_ml,
  Stress ~ bz1_1 + wz1_1 + Female + Age + (1 | ID)
)

fit_agg_ml <- fit_brmcoda(
  cilr_ml,
  Stress ~ z1_1 + z2_1 + z3_1 + z4_1 + Female + Age + (1 | ID)
)

fit_partial_agg <- fit_brmcoda(
  cilr_ml,
  Stress ~ z1_1 + z2_1 + Female + Age + (1 | ID)
)

fit_single <- fit_brmcoda(
  cilr_wide,
  Stress ~ z1_1 + z2_1 + z3_1 + z4_1 + Female + Age
)

foreach::registerDoSEQ()

test_that("simple substitution wrappers cover standard summary paths", {
  b_out <- bsub(fit_ml, delta = 1:2, base = psub)
  w_out <- wsub(fit_ml, delta = 1:2, base = psub)
  s_out <- sub(fit_agg_ml, delta = 1:2, base = psub)

  expect_standard_sub_output(b_out, "between", "grandmean")
  expect_standard_sub_output(w_out, "within", "grandmean")
  expect_standard_sub_output(s_out, "aggregate", "grandmean")
})

test_that("wrapper-specific ref, range, and raw-output branches are exercised", {
  rg_between <- build.rg(
    object = fit_ml,
    ref = "grandmean",
    at = NULL,
    parts = parts,
    level = "between",
    weight = "equal"
  )
  rg_within <- build.rg(
    object = fit_ml,
    ref = "grandmean",
    at = NULL,
    parts = parts,
    level = "within",
    weight = "equal"
  )
  rg_agg <- build.rg(
    object = fit_agg_ml,
    ref = "grandmean",
    at = NULL,
    parts = parts,
    level = "aggregate",
    weight = "equal"
  )

  user_b <- bsub(fit_ml, delta = 1, base = psub, ref = rg_between[1])
  user_w <- wsub(fit_ml, delta = 1, base = psub, ref = rg_within[1])
  user_s <- sub(fit_agg_ml, delta = 1, base = psub, ref = rg_agg[1])

  expect_true(all(user_b[[1]][[1]]$Reference == "users"))
  expect_true(all(user_w[[1]][[1]]$Reference == "users"))
  expect_true(all(user_s[[1]][[1]]$Reference == "users"))

  expect_error(bsub(fit_ml, delta = 5000, base = psub), "less than or equal")
  expect_error(wsub(fit_ml, delta = 5000, base = psub), "less than or equal")
  expect_message(
    sub(fit_agg_ml, delta = 5000, base = psub),
    "unrealistic substitution"
  )

  expect_error(bsub(fit_ml, delta = 1, base = psub, ref = 1), "ref must be")
  expect_error(wsub(fit_ml, delta = 1, base = psub, ref = 1), "ref must be")
  expect_error(sub(fit_agg_ml, delta = 1, base = psub, ref = 1), "ref must be")

  expect_error(
    bsub(fit_ml, delta = 1, base = psub, ref = data.table(Female = 1, Age = mean(ml_data$Age))),
    "covariates"
  )
  expect_error(
    wsub(fit_ml, delta = 1, base = psub, ref = data.table(Female = 1, Age = mean(ml_data$Age))),
    "covariates"
  )
  expect_error(
    sub(fit_agg_ml, delta = 1, base = psub, ref = data.table(Female = 1, Age = mean(ml_data$Age))),
    "covariates"
  )

  raw_b <- bsub(
    fit_ml,
    delta = 1,
    base = psub,
    at = list(Female = c(0, 1)),
    aorg = FALSE,
    summary = FALSE
  )
  expect_named(raw_b[[1]][["Stress"]], c("posterior", "grid"))
  expect_equal(sort(raw_b[[1]][["Stress"]]$grid$Female), c(0, 1))
  expect_true(all(c("Female0", "Female1") %in% names(raw_b[[1]][["Stress"]]$posterior)))
})

test_that("substitution object and summary cover filter, print, and one-to-all branches", {
  sub_obj <- substitution(
    object = fit_ml,
    delta = 1:2,
    ref = c("grandmean", "clustermean"),
    level = c("between", "within")
  )
  one_to_all <- substitution(
    object = fit_ml,
    delta = 1:2,
    ref = "grandmean",
    level = c("between", "within"),
    type = "one-to-all"
  )
  draws_obj <- substitution(
    object = fit_ml,
    delta = 1,
    ref = "grandmean",
    level = "between",
    summary = FALSE
  )

  expect_true(is.substitution(sub_obj))
  expect_s3_class(sub_obj, "substitution")
  expect_named(
    sub_obj,
    c(
      "between_simple_sub", "within_simple_sub", "simple_sub",
      "between_avg_sub", "within_avg_sub", "avg_sub",
      "brmsformula", "delta", "ref", "level", "weight",
      "parts", "at", "summary", "type"
    )
  )

  summary_default <- summary(sub_obj)
  summary_filtered <- summary(
    sub_obj,
    delta = 1,
    to = "TST",
    from = "WAKE",
    ref = "grandmean",
    level = "between"
  )
  summary_one_to_all <- summary(one_to_all, delta = 1, to = "TST")

  expect_s3_class(summary_default, "data.table")
  expect_true(all(c("Estimate", "Est.Error", "CI_low", "CI_high", "Resp") %in% names(summary_default)))
  expect_equal(nrow(summary_filtered), 1L)
  expect_setequal(unique(summary_one_to_all$From), c("remaining", "TST"))
  expect_setequal(unique(summary_one_to_all$To), c("remaining", "TST"))

  expect_error(summary(sub_obj, ref = "bad"), "'ref' should be either one")
  expect_error(summary(sub_obj, level = "bad"), "'level' should be either one")
  expect_error(summary(sub_obj, to = "bad"), "'to' should be names")
  expect_error(summary(sub_obj, from = "bad"), "'from' should be names")
  expect_error(summary(sub_obj, delta = 999), "empty data.table")
  expect_error(summary(draws_obj), "currently not available")

  expect_output(print(sub_obj), "Estimate")
})

test_that("substitution covers explicit ref, weight, and type branches", {
  rg_between <- build.rg(
    object = fit_ml,
    ref = "grandmean",
    at = NULL,
    parts = parts,
    level = "between",
    weight = "equal"
  )
  rg_within <- build.rg(
    object = fit_ml,
    ref = "grandmean",
    at = NULL,
    parts = parts,
    level = "within",
    weight = "equal"
  )
  rg_agg <- build.rg(
    object = fit_agg_ml,
    ref = "grandmean",
    at = NULL,
    parts = parts,
    level = "aggregate",
    weight = "equal"
  )

  explicit_type <- substitution(
    object = fit_ml,
    delta = 1,
    type = "one-to-one"
  )
  between_user <- substitution(
    object = fit_ml,
    delta = 1,
    ref = rg_between[1],
    level = "between",
    at = list(Female = c(0, 1))
  )
  within_user <- substitution(
    object = fit_ml,
    delta = 1,
    ref = rg_within[1],
    level = "within",
    at = list(Female = c(0, 1))
  )
  aggregate_user <- substitution(
    object = fit_agg_ml,
    delta = 1,
    ref = rg_agg[1],
    level = "aggregate",
    at = list(Female = c(0, 1))
  )
  between_cluster <- substitution(
    object = fit_ml,
    delta = 1,
    ref = "clustermean",
    level = "between",
    weight = "proportional"
  )
  within_cluster <- substitution(
    object = fit_ml,
    delta = 1,
    ref = "clustermean",
    level = "within",
    weight = "proportional"
  )
  aggregate_cluster <- substitution(
    object = fit_agg_ml,
    delta = 1,
    ref = "clustermean",
    level = "aggregate",
    weight = "proportional"
  )

  expect_true(is.substitution(explicit_type))
  expect_true(is.substitution(between_user))
  expect_true(is.substitution(within_user))
  expect_true(is.substitution(aggregate_user))
  expect_true(is.substitution(between_cluster))
  expect_true(is.substitution(within_cluster))
  expect_true(is.substitution(aggregate_cluster))

  expect_true(is.list(between_user$between_simple_sub))
  expect_true(is.list(within_user$within_simple_sub))
  expect_true(is.list(aggregate_user$simple_sub))
  expect_true(all(between_cluster$between_avg_sub[[1]][[1]]$Reference == "clustermean"))
  expect_true(all(within_cluster$within_avg_sub[[1]][[1]]$Reference == "clustermean"))
  expect_true(all(aggregate_cluster$avg_sub[[1]][[1]]$Reference == "clustermean"))
})

test_that("substitution validates early error branches and fallback behavior", {
  fake_non_ilr <- structure(list(complr = list(transform = "alr")), class = "brmcoda")

  expect_error(substitution(delta = 1), "required argument")
  expect_error(substitution(object = list(), delta = 1), "fitted 'brmcoda' object")
  expect_error(substitution(object = fake_non_ilr, delta = 1), "should be fitted with ilr transform")
  expect_error(substitution(object = fit_ml), "'delta' is a required argument")
  expect_error(substitution(object = fit_ml, delta = -1), "positive integer value")
  expect_error(substitution(object = fit_ml, delta = 1, base = list()), "data frame or data table")
  expect_error(substitution(object = fit_nocomp, delta = 1), "No fixed effects of composition")
  expect_error(
    substitution(object = fit_partial_bw, delta = 1, level = "aggregate"),
    "aggregate' substitution analysis cannot be computed"
  )
  expect_true(is.substitution(substitution(object = fit_partial_bw, delta = 1, level = "between")))
  expect_true(is.substitution(substitution(object = fit_partial_bw, delta = 1, level = "within")))
  expect_error(
    substitution(object = fit_partial_agg, delta = 1, level = "aggregate"),
    "cannot be recycled"
  )
  expect_true(is.substitution(substitution(object = fit_single, delta = 1, level = "between")))

  expect_warning(
    single_level <- substitution(
      object = fit_single,
      delta = 1,
      ref = "clustermean",
      level = c("between", "within", "aggregate")
    ),
    "Can only use grandmean for single level model"
  )
  expect_identical(single_level$ref, "grandmean")
  expect_identical(single_level$level, "aggregate")
  expect_true(is.list(single_level$simple_sub))
})

test_that("brmcoda warning path for legacy ilr names is exercised", {
  legacy_result <- NULL

  expect_warning(
    legacy_result <- tryCatch(
      brmcoda(
        complr = cilr_ml,
        formula = Stress ~ ilr1 + Female + Age + (1 | ID),
        chain = 1,
        iter = 100,
        warmup = 50,
        seed = 123,
        backend = backend,
        silent = 2,
        refresh = 0
      ),
      error = identity
    ),
    "updated naming convention"
  )

  expect_s3_class(legacy_result, "error")
})
