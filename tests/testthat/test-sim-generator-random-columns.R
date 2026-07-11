library(testthat)
library(multilevelcoda)

expect_group_constant <- function(data, col, group = "group_id") {
  by_group <- split(data[[col]], data[[group]])
  expect_true(all(vapply(by_group, function(x) length(unique(x)) == 1L, logical(1))))
}

expect_no_realized_random_metadata <- function(metadata) {
  expect_null(metadata$random_intercepts)
  expect_null(metadata$scale_random_intercepts)
  expect_null(metadata$scale_linear_predictor)

  group_parameter_names <- names(metadata$group_parameters)
  if (is.null(group_parameter_names)) {
    group_parameter_names <- character()
  }
  expect_false("random_intercepts" %in% group_parameter_names)
  expect_false("scale_random_intercepts" %in% group_parameter_names)
}

expect_level2_matches_single_direct <- function(generator, column, seed = 198) {
  group_sizes <- c(2L, 3L, 1L, 2L)
  single <- simulate_data(
    n = length(group_sizes),
    seed = seed,
    generators = list(x = generator("single"))
  )
  level2 <- simulate_data(
    n_groups = length(group_sizes),
    n_per_group = group_sizes,
    seed = seed,
    generators = list(x = generator("level2"))
  )

  expect_equal(level2$data[[column]], rep(single$data[[column]], group_sizes))
}

test_that("level2 direct generators share single-level resolver before expansion", {
  expect_level2_matches_single_direct(
    function(level) gen_mvn("x", level = level, fixed_intercept = 2, residual_cov = 4),
    "x"
  )
  expect_level2_matches_single_direct(
    function(level) gen_categorical("x", level = level, fixed_intercept = stats::qlogis(0.35)),
    "x"
  )
  expect_level2_matches_single_direct(
    function(level) gen_binomial("x", level = level, size = 5, fixed_intercept = stats::qlogis(0.4)),
    "x"
  )
  expect_level2_matches_single_direct(
    function(level) gen_poisson("x", level = level, fixed_intercept = log(2)),
    "x"
  )
  expect_level2_matches_single_direct(
    function(level) gen_negbin("x", level = level, fixed_intercept = log(2), scale_fixed_intercept = log(3)),
    "x"
  )
  expect_level2_matches_single_direct(
    function(level) gen_gamma("x", level = level, fixed_intercept = log(20), scale_fixed_intercept = log(4)),
    "x"
  )
  expect_level2_matches_single_direct(
    function(level) gen_beta("x", level = level, fixed_intercept = stats::qlogis(0.7), scale_fixed_intercept = log(12)),
    "x"
  )
})

test_that("categorical generators require fixed intercepts", {
  expect_error(
    gen_categorical("x"),
    "Categorical generator requires `fixed_intercept`."
  )
  expect_error(
    gen_categorical("x", level = "level2"),
    "Categorical generator requires `fixed_intercept`."
  )
  expect_error(
    gen_categorical("x", fixed_intercept = 0, random_cov = 1),
    "`random_cov` requires `level = \"multilevel\"`."
  )
})

test_that("generators require explicit link-scale parameters", {
  expect_error(gen_mvn("x"), "MVN generator requires `fixed_intercept`.")
  expect_error(gen_mvn("x", fixed_intercept = 0), "requires `residual_cov` or `scale_fixed_intercept`")
  expect_error(gen_binomial("x", fixed_intercept = stats::qlogis(0.5)), "Binomial generator requires `size`.")
  expect_error(gen_poisson("x"), "Poisson generator requires `fixed_intercept`.")
  expect_error(gen_negbin("x", fixed_intercept = log(2)), "Negative binomial generator requires `scale_fixed_intercept`.")
  expect_error(gen_gamma("x", fixed_intercept = log(2)), "Gamma generator requires `scale_fixed_intercept`.")
  expect_error(gen_beta("x", fixed_intercept = stats::qlogis(0.5)), "Beta generator requires `scale_fixed_intercept`.")
})

test_that("removed direct distribution arguments are rejected", {
  expect_error(gen_mvn("x", mean = 0), "unused argument")
  expect_error(gen_mvn("x", cov = 1), "unused argument")
  expect_error(gen_categorical("x", prob = 0.5), "unused argument")
  expect_error(gen_binomial("x", size = 4, prob = 0.5), "unused argument")
  expect_error(gen_poisson("x", lambda = 2), "unused argument")
  expect_error(gen_negbin("x", size = 3), "unused argument")
  expect_error(gen_negbin("x", mu = 2), "unused argument")
  expect_error(gen_gamma("x", shape = 2), "unused argument")
  expect_error(gen_gamma("x", mean = 1), "unused argument")
  expect_error(gen_gamma("x", rate = 1), "unused argument")
  expect_error(gen_gamma("x", scale = 1), "unused argument")
  expect_error(gen_beta("x", mean = 0.5), "unused argument")
  expect_error(gen_beta("x", precision = 10), "unused argument")
})

test_that("univariate MVN records explicit scalar covariance", {
  sim <- simulate_data(
    n = 4,
    seed = 199,
    generators = list(x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1))
  )

  expect_true("x" %in% names(sim$data))
  expect_identical(sim$generator_specs$x$type, "mvn")
  expect_identical(sim$generator_metadata$x$distribution, "mvn")
  expect_equal(sim$generator_metadata$x$fixed_intercept, c(x = 0))
  expect_equal(sim$generator_metadata$x$row_parameters$mean, matrix(0, 4, 1, dimnames = list(NULL, "x")))
  expect_equal(sim$generator_metadata$x$row_parameters$residual_cov, matrix(1, 1, 1, dimnames = list("x", "x")))
})

test_that("univariate MVN scalar covariance matches univariate normal draws", {
  sim <- simulate_data(
    n = 5,
    seed = 200,
    generators = list(x = gen_mvn("x", fixed_intercept = 2, residual_cov = 4))
  )

  expect_equal(sim$data$x, withr::with_seed(200, stats::rnorm(5, mean = 2, sd = 2)))
  expect_equal(sim$generator_metadata$x$row_parameters$residual_cov, matrix(4, 1, 1, dimnames = list("x", "x")))
})

test_that("univariate MVN level2 output is constant within group", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 3,
    seed = 200,
    generators = list(x = gen_mvn("x", level = "level2", fixed_intercept = 5, residual_cov = 1))
  )

  expect_group_constant(sim$data, "x")
  expect_identical(sim$generator_metadata$x$parameter_level, "group")
  expect_identical(sim$generator_metadata$x$distribution, "mvn")
})

test_that("MVN generators record column roles at single and level2 designs", {
  single <- simulate_data(
    n = 3,
    seed = 200,
    generators = list(x = gen_mvn("x", fixed_intercept = 1, residual_cov = 0))
  )
  level2 <- simulate_data(
    n_groups = 3,
    n_per_group = c(1, 2, 3),
    seed = 201,
    generators = list(x = gen_mvn("x", level = "level2", fixed_intercept = 2, residual_cov = 0))
  )

  expect_equal(
    as.data.frame(single$generator_metadata$x$column_roles),
    data.frame(column = "x", variable = "x", component = "observed", level = "row")
  )
  expect_equal(
    as.data.frame(level2$generator_metadata$x$column_roles),
    data.frame(column = "x", variable = "x", component = "between", level = "group")
  )
  expect_group_constant(level2$data, "x")
})

test_that("univariate MVN rejects removed normal variance aliases", {
  expect_error(gen_mvn("x", sd = 1), "unused argument")
  expect_error(
    gen_mvn(
      "x",
      level = "multilevel",
      fixed_intercept = 0,
      random_cov = 1,
      residual_var = 1
    ),
    "unused argument"
  )
})

test_that("multilevel generators require explicit random and residual parameters", {
  expect_error(
    gen_poisson("x", level = "multilevel", fixed_intercept = 0),
    "Multilevel poisson generator requires `random_cov`."
  )
  expect_error(
    gen_mvn("x", level = "multilevel", fixed_intercept = 0, random_cov = 0),
    "requires `residual_cov` or `scale_fixed_intercept`"
  )
  expect_no_error(
    gen_mvn(
      c("x", "z"),
      level = "multilevel",
      fixed_intercept = c(0, 0),
      random_cov = diag(4),
      scale_fixed_intercept = c(0, 0)
    )
  )
})

test_that("multilevel univariate MVN writes location random intercept columns", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 3,
    seed = 201,
    generators = list(
      x = gen_mvn(
        "x",
        level = "multilevel",
        fixed_intercept = 2,
        random_cov = 0.5,
        residual_cov = 0
      )
    )
  )

  random_col <- ".mlsim_x_random_intercept"

  expect_true(all(c("x", random_col) %in% names(sim$data)))
  expect_identical(sim$generator_specs$x$vars, "x")
  expect_identical(sim$generator_metadata$x$vars, "x")
  expect_equal(sim$data$x, 2 + sim$data[[random_col]], tolerance = 1e-8)
  expect_group_constant(sim$data, random_col)
  expect_gt(length(unique(sim$data[[random_col]])), 1L)
  expect_true(summary(sim)$generators$has_random_effects)
  expect_equal(sim$generator_metadata$x$residual_cov, matrix(0, 1, 1, dimnames = list("x", "x")))
  expect_no_realized_random_metadata(sim$generator_metadata$x)
})

test_that("multilevel univariate MVN writes scale random intercept columns", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 3,
    seed = 202,
    generators = list(
      x = gen_mvn(
        "x",
        level = "multilevel",
        fixed_intercept = 0,
        random_cov = diag(c(0.5, 0.2)),
        scale_fixed_intercept = 0
      )
    )
  )

  random_cols <- c(
    ".mlsim_x_random_intercept",
    ".mlsim_x_scale_random_intercept"
  )

  expect_true(all(c("x", random_cols) %in% names(sim$data)))
  lapply(random_cols, function(col) expect_group_constant(sim$data, col))
  expect_no_realized_random_metadata(sim$generator_metadata$x)
})

test_that("shared multilevel families write location and scale random columns", {
  sim <- simulate_data(
    n_groups = 3,
    n_per_group = 2,
    seed = 203,
    generators = list(
      p = gen_poisson("p", level = "multilevel", fixed_intercept = 0.3, random_cov = 0.2),
      nb = gen_negbin(
        "nb",
        level = "multilevel",
        fixed_intercept = 0.2,
        scale_fixed_intercept = log(3),
        random_cov = diag(c(0.2, 0.1))
      ),
      g = gen_gamma(
        "g",
        level = "multilevel",
        fixed_intercept = 0.1,
        scale_fixed_intercept = log(2),
        random_cov = diag(c(0.2, 0.1))
      ),
      b = gen_beta(
        "b",
        level = "multilevel",
        fixed_intercept = 0.1,
        scale_fixed_intercept = log(10),
        random_cov = diag(c(0.2, 0.1))
      )
    )
  )

  expected_cols <- c(
    ".mlsim_p_random_intercept",
    ".mlsim_nb_random_intercept",
    ".mlsim_nb_scale_random_intercept",
    ".mlsim_g_random_intercept",
    ".mlsim_g_scale_random_intercept",
    ".mlsim_b_random_intercept",
    ".mlsim_b_scale_random_intercept"
  )

  expect_true(all(expected_cols %in% names(sim$data)))
  lapply(expected_cols, function(col) expect_group_constant(sim$data, col))
  lapply(sim$generator_metadata, expect_no_realized_random_metadata)
})

test_that("multilevel MVN writes ILR-based random columns for compositional output", {
  sim <- simulate_data(
    n_groups = 3,
    n_per_group = 2,
    seed = 204,
    generators = list(
      comp = gen_mvn(
        c("z1", "z2"),
        level = "multilevel",
        fixed_intercept = c(0, 0),
        random_cov = diag(c(0.2, 0.3, 0.1, 0.15)),
        scale_fixed_intercept = c(0, 0),
        residual_cor = diag(2),
        compositional = TRUE,
        parts = c("sleep", "active", "sedentary"),
        keep_ilr = FALSE
      )
    )
  )

  random_cols <- c(
    ".mlsim_z1_random_intercept",
    ".mlsim_z2_random_intercept",
    ".mlsim_z1_scale_random_intercept",
    ".mlsim_z2_scale_random_intercept"
  )

  decomposition_cols <- c(
    "z1_between", "z2_between", "z1_within", "z2_within",
    "sleep_between", "active_between", "sedentary_between"
  )
  expect_false(any(c("z1", "z2") %in% names(sim$data)))
  expect_true(all(c("sleep", "active", "sedentary", decomposition_cols, random_cols) %in% names(sim$data)))
  expect_identical(
    sim$generator_specs$comp$vars,
    c("sleep", "active", "sedentary", decomposition_cols)
  )
  expect_identical(
    sim$generator_metadata$comp$vars,
    c("sleep", "active", "sedentary", decomposition_cols)
  )
  lapply(random_cols, function(col) expect_group_constant(sim$data, col))
  expect_no_realized_random_metadata(sim$generator_metadata$comp)
})

test_that("multilevel categorical writes one random column per non-reference category", {
  sim <- simulate_data(
    n_groups = 3,
    n_per_group = 2,
    seed = 205,
    generators = list(
      arm = gen_categorical(
        "arm",
        level = "multilevel",
        categories = c("control", "treat", "high dose"),
        fixed_intercept = setNames(c(0.2, -0.1), c("treat", "high dose")),
        random_cov = diag(2),
        reference = "control",
        output = "character"
      )
    )
  )

  random_cols <- c(
    ".mlsim_arm_random_intercept_treat",
    ".mlsim_arm_random_intercept_high.dose"
  )

  expect_true(all(c("arm", random_cols) %in% names(sim$data)))
  lapply(random_cols, function(col) expect_group_constant(sim$data, col))
  expect_no_realized_random_metadata(sim$generator_metadata$arm)
})

test_that("internal random column names cannot collide with existing columns", {
  expect_error(
    simulate_data(
      n_groups = 2,
      n_per_group = 2,
      seed = 206,
      generators = list(
        existing = gen_mvn(".mlsim_x_random_intercept", fixed_intercept = 0, residual_cov = 1),
        x = gen_mvn(
          "x",
          level = "multilevel",
          fixed_intercept = 0,
          random_cov = 0.2,
          residual_cov = 1
        )
      )
    ),
    "would overwrite existing columns"
  )
})

test_that("custom generator column roles are validated and preserved", {
  custom_with_roles <- function(context, vars, level) {
    list(
      data = data.frame(x = seq_len(context$n_rows)),
      names = vars,
      metadata = list(
        column_roles = data.frame(
          column = vars,
          variable = vars,
          component = "observed",
          level = "row"
        )
      )
    )
  }

  sim <- simulate_data(
    n = 3,
    generators = list(x = gen_custom("x", custom_with_roles))
  )

  expect_equal(
    as.data.frame(sim$generator_metadata$x$column_roles),
    data.frame(column = "x", variable = "x", component = "observed", level = "row")
  )

  expect_error(
    simulate_data(
      n = 3,
      generators = list(
        x = gen_custom("x", function(context, vars, level) {
          list(
            data = data.frame(x = seq_len(context$n_rows)),
            names = vars,
            metadata = list(
              column_roles = data.frame(
                column = "missing",
                variable = vars,
                component = "observed",
                level = "row"
              )
            )
          )
        })
      )
    ),
    "must refer to generated columns"
  )

  expect_error(
    simulate_data(
      n = 3,
      generators = list(
        x = gen_custom("x", function(context, vars, level) {
          list(
            data = data.frame(x = seq_len(context$n_rows)),
            names = vars,
            metadata = list(
              column_roles = data.frame(
                column = c(vars, vars),
                variable = c(vars, vars),
                component = c("observed", "observed"),
                level = c("row", "row")
              )
            )
          )
        })
      )
    ),
    "duplicate variable/component"
  )
})

test_that("compositional multilevel MVN emits between/within columns and roles", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    group_id = "ID",
    seed = 61,
    generators = list(
      comp = gen_mvn(
        c("ilr1", "ilr2"),
        level = "multilevel",
        fixed_intercept = c(0.2, -0.1),
        random_cov = diag(2) * 0.2,
        residual_cov = diag(2) * 0.4,
        compositional = TRUE,
        parts = c("sleep", "active", "sedentary"),
        total = 24
      )
    )
  )

  data <- sim$data
  expect_true(all(c(
    "ilr1_between", "ilr2_between", "ilr1_within", "ilr2_within",
    "sleep_between", "active_between", "sedentary_between"
  ) %in% names(data)))

  expect_equal(data$ilr1, data$ilr1_between + data$ilr1_within)
  expect_equal(data$ilr2, data$ilr2_between + data$ilr2_within)
  expect_equal(
    rowSums(as.matrix(data[, c("sleep_between", "active_between", "sedentary_between"), with = FALSE])),
    rep(24, nrow(data))
  )
  back <- multilevelcoda:::.mlsim_ilr_inverse(
    as.matrix(data[, c("ilr1_between", "ilr2_between"), with = FALSE]),
    parts = c("sleep", "active", "sedentary"),
    total = 24,
    sbp = sim$generator_metadata$comp$sbp
  )$values
  expect_equal(unname(back), unname(as.matrix(
    data[, c("sleep_between", "active_between", "sedentary_between"), with = FALSE]
  )))

  roles <- sim$generator_metadata$comp$column_roles
  expect_identical(sort(roles$column[roles$component == "between"]), c("ilr1_between", "ilr2_between"))
  expect_identical(sort(roles$column[roles$component == "within"]), c("ilr1_within", "ilr2_within"))
  expect_identical(sort(roles$column[roles$component == "observed"]), c("ilr1", "ilr2"))
  expect_false(any(c("sleep", "active", "sedentary") %in% roles$variable))
  expect_identical(
    sim$generator_metadata$comp$between_parts_vars,
    c("sleep_between", "active_between", "sedentary_between")
  )

  group_constant <- sim$data[, lapply(.SD, function(x) max(x) - min(x)), by = ID,
                             .SDcols = c("ilr1_between", "sleep_between")]
  expect_true(all(unlist(group_constant[, -1]) <= 1e-12))
})

test_that("compositional MVN roles adapt to keep_ilr and level", {
  no_ilr <- simulate_data(
    n_groups = 3,
    n_per_group = 4,
    group_id = "ID",
    seed = 62,
    generators = list(
      comp = gen_mvn(
        c("ilr1", "ilr2"),
        level = "multilevel",
        fixed_intercept = c(0, 0),
        random_cov = diag(2) * 0.2,
        residual_cov = diag(2) * 0.4,
        compositional = TRUE,
        parts = c("a", "b", "c"),
        keep_ilr = FALSE
      )
    )
  )
  expect_false(any(c("ilr1", "ilr2") %in% names(no_ilr$data)))
  expect_true(all(c("ilr1_between", "ilr1_within") %in% names(no_ilr$data)))
  roles <- no_ilr$generator_metadata$comp$column_roles
  expect_false("observed" %in% roles$component)
  expect_setequal(roles$component, c("between", "within"))

  level2 <- simulate_data(
    n_groups = 3,
    n_per_group = 4,
    group_id = "ID",
    seed = 63,
    generators = list(
      comp = gen_mvn(
        c("ilr1", "ilr2"),
        level = "level2",
        fixed_intercept = c(0, 0),
        residual_cov = diag(2),
        compositional = TRUE,
        parts = c("a", "b", "c")
      )
    )
  )
  roles2 <- level2$generator_metadata$comp$column_roles
  expect_identical(unique(roles2$component), "between")
  expect_identical(unique(roles2$level), "group")

  single <- simulate_data(
    n = 6,
    seed = 64,
    generators = list(
      comp = gen_mvn(
        c("ilr1", "ilr2"),
        fixed_intercept = c(0, 0),
        residual_cov = diag(2),
        compositional = TRUE,
        parts = c("a", "b", "c")
      )
    )
  )
  roles3 <- single$generator_metadata$comp$column_roles
  expect_identical(unique(roles3$component), "observed")
})
