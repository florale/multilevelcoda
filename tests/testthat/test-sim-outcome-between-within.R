library(testthat)
library(multilevelcoda)

expect_helper_decomposition <- function(data, var, expected_between) {
  expect_equal(data[[paste0(var, "_between")]], expected_between, tolerance = 1e-8)
  expect_equal(data[[paste0(var, "_within")]], data[[var]] - expected_between, tolerance = 1e-8)
}

lag_by_group <- function(x, group) {
  out <- rep(NA_real_, length(x))
  for (idx in split(seq_along(x), group)) {
    if (length(idx) > 1L) {
      out[idx[-1L]] <- x[idx[-length(idx)]]
    }
  }
  out
}

test_that("normal multilevel helpers use fixed plus retained random intercept", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    seed = 301,
    generators = list(
      x = gen_mvn(
        "x",
        level = "multilevel",
        fixed_intercept = 2,
        random_cov = 0.25,
        residual_cov = 0.75
      ),
      y = gen_outcome(y ~ between(x) + within(x), residual_cov = 0)
    )
  )

  expected_between <- 2 + sim$data$.mlsim_x_random_intercept
  expect_helper_decomposition(sim$data, "x", expected_between)
})

test_that("normal multilevel helpers ignore retained scale random intercepts", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    seed = 311,
    generators = list(
      x = gen_mvn(
        "x",
        level = "multilevel",
        fixed_intercept = 1,
        scale_fixed_intercept = log(0.5),
        random_cov = diag(c(0.25, 0.10))
      ),
      y = gen_outcome(y ~ between(x) + within(x), residual_cov = 0)
    )
  )

  expected_between <- 1 + sim$data$.mlsim_x_random_intercept
  expect_true(".mlsim_x_scale_random_intercept" %in% names(sim$data))
  expect_helper_decomposition(sim$data, "x", expected_between)
})

test_that("predictor lag1 uses lagged within component for grouped data", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    seed = 312,
    generators = list(
      x = gen_mvn(
        "x",
        level = "multilevel",
        fixed_intercept = 2,
        random_cov = 0.25,
        residual_cov = 0.75
      ),
      y = gen_outcome(y ~ lag1(x), residual_cov = 0)
    )
  )

  x_within <- sim$data$x - (2 + sim$data$.mlsim_x_random_intercept)
  expect_equal(sim$data$lag_x, lag_by_group(x_within, sim$data$group_id), tolerance = 1e-8)
})

test_that("predictor lag1 falls back to lagged observed group-mean deviations", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    seed = 313,
    generators = list(
      x = gen_mvn("x", fixed_intercept = 2, residual_cov = 1),
      y = gen_outcome(y ~ lag1(x), residual_cov = 0)
    )
  )

  x_between <- stats::ave(
    sim$data$x,
    sim$data$group_id,
    FUN = function(z) mean(z, na.rm = TRUE)
  )
  expect_equal(sim$data$lag_x, lag_by_group(sim$data$x - x_between, sim$data$group_id), tolerance = 1e-8)
})

test_that("single-level predictor lag1 keeps previous-row behavior", {
  sim <- simulate_data(
    n = 6,
    seed = 314,
    generators = list(
      x = gen_mvn("x", fixed_intercept = 2, residual_cov = 1),
      y = gen_outcome(y ~ lag1(x), residual_cov = 0)
    )
  )

  expect_equal(sim$data$lag_x, c(NA, sim$data$x[-nrow(sim$data)]), tolerance = 1e-8)
})

test_that("binary categorical predictor lag1 uses lagged within indicator", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    seed = 315,
    generators = list(
      arm = gen_categorical(
        "arm",
        level = "multilevel",
        categories = c("control", "treat"),
        fixed_intercept = c(treat = stats::qlogis(0.6)),
        random_cov = 0.15
      ),
      y = gen_outcome(y ~ lag1(arm), residual_cov = 0)
    )
  )

  arm_between <- stats::plogis(
    stats::qlogis(0.6) + sim$data$.mlsim_arm_random_intercept_treat
  )
  arm_within <- as.numeric(as.character(sim$data$arm) == "treat") - arm_between
  expect_equal(sim$data$lag_arm, lag_by_group(arm_within, sim$data$group_id), tolerance = 1e-8)
})

test_that("helpers fall back to observed group means without retained random effects", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    seed = 302,
    generators = list(
      x = gen_mvn("x", fixed_intercept = 2, residual_cov = 1),
      y = gen_outcome(y ~ between(x) + within(x), residual_cov = 0)
    )
  )

  expected_between <- stats::ave(
    sim$data$x,
    sim$data$group_id,
    FUN = function(z) mean(z, na.rm = TRUE)
  )
  expect_helper_decomposition(sim$data, "x", expected_between)
})

test_that("non-normal multilevel helpers use observed-scale expected means", {
  poisson <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    seed = 303,
    generators = list(
      p = gen_poisson("p", level = "multilevel", fixed_intercept = log(2.5), random_cov = 0.15),
      y = gen_outcome(y ~ between(p) + within(p), residual_cov = 0)
    )
  )
  expect_helper_decomposition(
    poisson$data,
    "p",
    exp(log(2.5) + poisson$data$.mlsim_p_random_intercept)
  )

  negbin <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    seed = 304,
    generators = list(
      nb = gen_negbin(
        "nb",
        level = "multilevel",
        fixed_intercept = log(3),
        random_cov = 0.15,
        scale_fixed_intercept = log(4)
      ),
      y = gen_outcome(y ~ between(nb) + within(nb), residual_cov = 0)
    )
  )
  expect_helper_decomposition(
    negbin$data,
    "nb",
    exp(log(3) + negbin$data$.mlsim_nb_random_intercept)
  )

  gamma <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    seed = 305,
    generators = list(
      g = gen_gamma(
        "g",
        level = "multilevel",
        fixed_intercept = log(1.5),
        random_cov = 0.15,
        scale_fixed_intercept = log(3)
      ),
      y = gen_outcome(y ~ between(g) + within(g), residual_cov = 0)
    )
  )
  expect_helper_decomposition(
    gamma$data,
    "g",
    exp(log(1.5) + gamma$data$.mlsim_g_random_intercept)
  )

  beta <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    seed = 306,
    generators = list(
      b = gen_beta(
        "b",
        level = "multilevel",
        fixed_intercept = stats::qlogis(0.4),
        random_cov = 0.15,
        scale_fixed_intercept = log(12)
      ),
      y = gen_outcome(y ~ between(b) + within(b), residual_cov = 0)
    )
  )
  expect_helper_decomposition(
    beta$data,
    "b",
    stats::plogis(stats::qlogis(0.4) + beta$data$.mlsim_b_random_intercept)
  )

  binomial <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    seed = 307,
    generators = list(
      count = gen_binomial(
        "count",
        level = "multilevel",
        size = 6,
        fixed_intercept = stats::qlogis(0.35),
        random_cov = 0.15
      ),
      y = gen_outcome(y ~ between(count) + within(count), residual_cov = 0)
    )
  )
  expect_helper_decomposition(
    binomial$data,
    "count",
    6 * stats::plogis(stats::qlogis(0.35) + binomial$data$.mlsim_count_random_intercept)
  )
})

test_that("binary categorical helpers use probability and realised indicator", {
  factor_sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    seed = 308,
    generators = list(
      arm = gen_categorical(
        "arm",
        level = "multilevel",
        categories = c("control", "treat"),
        fixed_intercept = c(treat = stats::qlogis(0.6)),
        random_cov = 0.15
      ),
      y = gen_outcome(y ~ between(arm) + within(arm), residual_cov = 0)
    )
  )
  factor_between <- stats::plogis(
    stats::qlogis(0.6) + factor_sim$data$.mlsim_arm_random_intercept_treat
  )
  expect_equal(factor_sim$data$arm_between, factor_between, tolerance = 1e-8)
  expect_equal(
    factor_sim$data$arm_within,
    as.numeric(as.character(factor_sim$data$arm) == "treat") - factor_between,
    tolerance = 1e-8
  )

  integer_sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    seed = 309,
    generators = list(
      arm = gen_categorical(
        "arm",
        level = "multilevel",
        categories = c("control", "treat"),
        fixed_intercept = c(treat = stats::qlogis(0.6)),
        random_cov = 0.15,
        output = "integer"
      ),
      y = gen_outcome(y ~ between(arm) + within(arm), residual_cov = 0)
    )
  )
  integer_between <- stats::plogis(
    stats::qlogis(0.6) + integer_sim$data$.mlsim_arm_random_intercept_treat
  )
  expect_equal(integer_sim$data$arm_between, integer_between, tolerance = 1e-8)
  expect_equal(
    integer_sim$data$arm_within,
    as.numeric(integer_sim$data$arm == 1L) - integer_between,
    tolerance = 1e-8
  )
})

test_that("multicategory categorical helpers fail clearly", {
  expect_error(
    simulate_data(
      n_groups = 3,
      n_per_group = 4,
      seed = 310,
      generators = list(
        arm = gen_categorical(
          "arm",
          level = "multilevel",
          categories = c("control", "low", "high"),
          fixed_intercept = c(low = 0.2, high = -0.1),
          random_cov = diag(2)
        ),
        y = gen_outcome(y ~ between(arm), residual_cov = 0)
      )
    ),
    "only support binary multilevel categorical predictors"
  )
})
