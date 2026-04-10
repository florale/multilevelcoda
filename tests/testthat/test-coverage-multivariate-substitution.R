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

parts_x <- c("TST", "WAKE")
parts_y <- c("MVPA", "LPA", "SB")

data(mcompd)

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

mv_data <- mcompd[ID %in% 1:8, .SD[1:3], by = ID]
mv_complr <- complr(
  data = mv_data,
  parts = list(parts_x, parts_y),
  total = list(480, 960),
  idvar = "ID",
  transform = "ilr"
)

pair_complr <- complr(
  data = mv_data,
  sbp = as.matrix(data.table(1, -1)),
  parts = parts_x,
  idvar = "ID",
  total = 480
)

fit_mv_between <- fit_brmcoda(
  mv_complr,
  mvbind(bz1_2, bz2_2) ~ bz1_1 + Female + Age + (1 | ID)
)

fit_mv_within <- fit_brmcoda(
  mv_complr,
  mvbind(wz1_2, wz2_2) ~ wz1_1 + Female + Age + (1 | ID)
)

fit_mv_agg <- fit_brmcoda(
  mv_complr,
  mvbind(z1_2, z2_2) ~ z1_1 + Female + Age + (1 | ID)
)

fit_pair_bw <- fit_brmcoda(
  pair_complr,
  Stress ~ bz1_1 + wz1_1 + Female + Age + (1 | ID)
)

fit_pair_agg <- fit_brmcoda(
  pair_complr,
  Stress ~ z1_1 + Female + Age + (1 | ID)
)

pair_base <- build.base(parts_x)
full_ref <- unique(mv_data[, .(TST, WAKE, Female, Age)])[1]
partial_ref <- full_ref[, .(TST, WAKE, Female)]

foreach::registerDoSEQ()

expect_response_closure <- function(out, total) {
  closure <- apply(out, c(1, 2), sum)
  expect_equal(
    unname(closure),
    matrix(total, nrow = nrow(closure), ncol = ncol(closure)),
    tolerance = 1e-6
  )
}

expect_raw_sub_output <- function(out) {
  expect_named(out[[1]][[1]], c("posterior", "grid"))
  expect_s3_class(out[[1]][[1]]$grid, "data.table")
  expect_equal(sort(out[[1]][[1]]$grid$Female), c(0, 1))
  expect_true(all(vapply(out[[1]][[1]]$posterior, is.data.frame, logical(1))))
  expect_true(all(c("Delta", "From", "To", "Level", "Reference") %in%
    names(out[[1]][[1]]$posterior[[1]])))
}

test_that("multivariate brmcoda methods cover response-type branches", {
  vars_between <- get_variables(fit_mv_between)
  vars_within <- get_variables(fit_mv_within)
  vars_agg <- get_variables(fit_mv_agg)

  expect_identical(vars_between$resp_type, "multivariate-non-compositional")
  expect_identical(vars_between$fixef_type, "between")
  expect_identical(vars_within$resp_type, "multivariate-non-compositional")
  expect_identical(vars_within$fixef_type, "within")
  expect_identical(vars_agg$resp_type, "multivariate-non-compositional")
  expect_identical(vars_agg$fixef_type, "aggregate")
  expect_identical(vars_agg$ranef_type, "multilevel")

  between_pred <- predict(fit_mv_between, scale = "response")
  between_fit <- fitted(fit_mv_between, scale = "response", summary = FALSE)
  within_pred <- predict(fit_mv_within, scale = "response")
  within_fit <- fitted(fit_mv_within, scale = "response", summary = FALSE)
  agg_pred <- predict(fit_mv_agg, scale = "response")
  agg_fit <- fitted(fit_mv_agg, scale = "response", summary = FALSE)

  expect_equal(dim(between_pred), c(nrow(mv_complr$dataout), 4L, length(parts_y)))
  expect_equal(dim(between_fit), c(150L, nrow(mv_complr$dataout), length(parts_y)))
  expect_equal(dim(within_pred), c(nrow(mv_complr$dataout), 4L, length(parts_y)))
  expect_equal(dim(within_fit), c(150L, nrow(mv_complr$dataout), length(parts_y)))
  expect_equal(dim(agg_pred), c(nrow(mv_complr$dataout), 4L, length(parts_y)))
  expect_equal(dim(agg_fit), c(150L, nrow(mv_complr$dataout), length(parts_y)))

  expect_response_closure(between_pred, 960)
  expect_response_closure(between_fit, 960)
  expect_response_closure(within_pred, 1)
  expect_response_closure(within_fit, 1)
  expect_response_closure(agg_pred, 960)
  expect_response_closure(agg_fit, 960)

  expect_type(coef(fit_mv_between), "list")
})

test_that("build.rg user-reference validation branches are exercised", {
  expect_error(
    build.rg(
      object = fit_mv_agg,
      ref = full_ref,
      at = NULL,
      parts = parts_x,
      level = "aggregate",
      weight = "equal",
      fill = FALSE
    ),
    "same as the composition"
  )
  expect_error(
    build.rg(
      object = fit_mv_between,
      ref = partial_ref,
      at = list(Age = c(mean(mv_data$Age), mean(mv_data$Age) + 1)),
      parts = parts_x,
      level = "between",
      weight = "equal",
      fill = TRUE
    ),
    "same as the composition"
  )
  expect_error(
    build.rg(
      object = fit_mv_agg,
      ref = mv_complr$dataout[1, .(TST = tTST, WAKE = tWAKE, Female, Age)],
      at = NULL,
      parts = parts_x,
      level = "aggregate",
      weight = "equal",
      fill = FALSE
    ),
    "list of length"
  )
})

test_that("multivariate simple substitution wrappers cover raw response and linear paths", {
  raw_between <- bsub(
    fit_mv_between,
    delta = 1:2,
    parts = parts_x,
    base = pair_base,
    type = "one-to-all",
    summary = FALSE,
    aorg = FALSE,
    at = list(Female = c(0, 1)),
    scale = "response",
    cores = 1
  )
  raw_within <- wsub(
    fit_mv_within,
    delta = 1:2,
    parts = parts_x,
    base = pair_base,
    type = "one-to-all",
    summary = FALSE,
    aorg = FALSE,
    at = list(Female = c(0, 1)),
    scale = "response",
    cores = 1
  )
  raw_agg <- sub(
    fit_mv_agg,
    delta = 1:2,
    parts = parts_x,
    base = pair_base,
    summary = FALSE,
    aorg = FALSE,
    at = list(Female = c(0, 1)),
    scale = "response",
    cores = 1
  )
  linear_between <- bsub(
    fit_mv_between,
    delta = 1,
    parts = parts_x,
    base = pair_base,
    scale = "linear"
  )

  expect_named(raw_between, parts_x)
  expect_named(raw_within, parts_x)
  expect_named(raw_agg, parts_x)
  expect_length(raw_between[[1]], length(parts_y))
  expect_length(raw_within[[1]], length(parts_y))
  expect_length(raw_agg[[1]], length(parts_y))
  expect_true(all(names(raw_between[[1]]) != ""))
  expect_true(all(names(raw_within[[1]]) != ""))
  expect_true(all(names(raw_agg[[1]]) != ""))

  expect_raw_sub_output(raw_between)
  expect_raw_sub_output(raw_within)
  expect_raw_sub_output(raw_agg)

  expect_length(linear_between[[1]], 2L)
  expect_true(all(names(linear_between[[1]]) %in% c("z1_2", "z2_2")))
  suppressWarnings(
    expect_error(
      sub(
        fit_mv_agg,
        delta = 1:2,
        parts = parts_x,
        base = pair_base,
        type = "one-to-all",
        summary = FALSE,
        aorg = FALSE,
        at = list(Female = c(0, 1)),
        scale = "response",
        cores = 1
      ),
      "Supplied 2 items"
    )
  )
})

test_that("multivariate margin wrappers cover raw outputs and one-to-all failures", {
  raw_bmargin <- bsubmargin(
    fit_pair_bw,
    delta = 1:2,
    parts = parts_x,
    base = pair_base,
    summary = FALSE,
    scale = "response",
    cores = 1
  )
  raw_wmargin <- wsubmargin(
    fit_pair_bw,
    delta = 1:2,
    parts = parts_x,
    base = pair_base,
    summary = FALSE,
    scale = "response",
    cores = 1
  )
  raw_smargin <- submargin(
    fit_pair_agg,
    delta = 1:2,
    parts = parts_x,
    base = pair_base,
    summary = FALSE,
    scale = "response",
    cores = 1
  )

  expect_s3_class(raw_bmargin[[1]][[1]], "data.table")
  expect_s3_class(raw_wmargin[[1]][[1]], "data.table")
  expect_s3_class(raw_smargin[[1]][[1]], "data.table")
  expect_gt(nrow(raw_bmargin[[1]][[1]]), 0L)
  expect_gt(nrow(raw_wmargin[[1]][[1]]), 0L)
  expect_gt(nrow(raw_smargin[[1]][[1]]), 0L)
  expect_true(all(c("Delta", "From", "To", "Level", "Reference") %in% names(raw_bmargin[[1]][[1]])))
  expect_true(all(c("Delta", "From", "To", "Level", "Reference") %in% names(raw_wmargin[[1]][[1]])))
  expect_true(all(c("Delta", "From", "To", "Level", "Reference") %in% names(raw_smargin[[1]][[1]])))

  suppressWarnings(
    expect_error(
      bsubmargin(
        fit_pair_bw,
        delta = 1:2,
        parts = parts_x,
        base = pair_base,
        type = "one-to-all",
        summary = FALSE,
        scale = "response",
        cores = 1
      ),
      "Supplied 2 items"
    )
  )
  suppressWarnings(
    expect_error(
      wsubmargin(
        fit_pair_bw,
        delta = 1:2,
        parts = parts_x,
        base = pair_base,
        type = "one-to-all",
        summary = FALSE,
        scale = "response",
        cores = 1
      ),
      "Supplied 2 items"
    )
  )
  suppressWarnings(
    expect_error(
      submargin(
        fit_pair_agg,
        delta = 1:2,
        parts = parts_x,
        base = pair_base,
        type = "one-to-all",
        summary = FALSE,
        scale = "response",
        cores = 1
      ),
      "Supplied 2 items"
    )
  )
})
