library(testthat)
library(multilevelcoda)

outcome_acceptance_params <- function() {
  outcomes <- c("ilr1", "ilr2")
  beta_location <- matrix(
    0,
    nrow = 4,
    ncol = 2,
    dimnames = list(
      c("(Intercept)", "treatmenttreatment", "between(stress)", "within(stress)"),
      outcomes
    )
  )
  beta_location["between(stress)", ] <- c(0.15, -0.05)
  beta_location["within(stress)", ] <- c(0.05, 0.03)

  beta_scale <- matrix(
    log(c(0.4, 0.35, 0.45, 0.4)),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("(Intercept)", "treatmenttreatment"), outcomes)
  )
  correlation <- diag(2)
  dimnames(correlation) <- list(outcomes, outcomes)

  beta_ar <- array(
    0,
    dim = c(2, 2, 2),
    dimnames = list(c("ar1()", "treatmenttreatment:ar1()"), outcomes, outcomes)
  )
  beta_ar["ar1()", , ] <- matrix(c(0.25, 0.02, -0.01, 0.20), 2, 2, byrow = TRUE)
  beta_ar["treatmenttreatment:ar1()", , ] <- matrix(c(0.05, 0, 0, 0.03), 2, 2, byrow = TRUE)

  random_names <- c(
    sprintf("location|outcome=%s|term=(Intercept)", outcomes),
    sprintf(
      "ar|term=ar1()|to=%s|from=%s",
      rep(outcomes, each = length(outcomes)),
      rep(outcomes, times = length(outcomes))
    ),
    sprintf("scale|outcome=%s|term=(Intercept)", outcomes)
  )
  random_cov <- matrix(0, length(random_names), length(random_names))
  dimnames(random_cov) <- list(random_names, random_names)

  list(
    location = list(beta = beta_location),
    scale = list(beta = beta_scale, correlation = correlation),
    ar = list(beta = beta_ar),
    random = list(ID = list(covariance = random_cov))
  )
}

test_that("multilevel MVN predictors emit visible role columns", {
  sim <- simulate_data(
    n_groups = 3,
    n_per_group = 4,
    seed = 100,
    generators = list(
      stress = gen_mvn(
        "stress",
        level = "multilevel",
        fixed_intercept = 0,
        random_cov = 0.2,
        residual_cov = 0.8
      )
    )
  )

  expect_true(all(c("stress", "stress_between", "stress_within") %in% names(sim$data)))
  expect_equal(sim$data$stress, sim$data$stress_between + sim$data$stress_within)
  expect_true(all(vapply(split(sim$data$stress_between, sim$data$group_id), function(x) {
    length(unique(x)) == 1L
  }, logical(1))))
  expect_equal(
    sim$generator_metadata$stress$column_roles$component,
    c("observed", "between", "within")
  )
})

test_that("gen_outcome simulates non-AR Gaussian outcomes", {
  beta_location <- matrix(
    c(0, 0.5),
    nrow = 2,
    dimnames = list(c("(Intercept)", "between(stress)"), "y")
  )
  beta_scale <- matrix(log(0.1), nrow = 1, dimnames = list("(Intercept)", "y"))

  sim <- simulate_data(
    n_groups = 3,
    n_per_group = 4,
    seed = 101,
    generators = list(
      stress = gen_mvn("stress", level = "multilevel", fixed_intercept = 0, random_cov = 0.2, residual_cov = 0.8),
      y = gen_outcome(
        y ~ between(stress),
        scale = sigma ~ 1,
        params = list(location = list(beta = beta_location), scale = list(beta = beta_scale)),
        burnin = 0
      )
    )
  )

  expect_true("y" %in% names(sim$data))
  expect_equal(dim(sim$generator_metadata$y$mu), c(nrow(sim$data), 1L))
  expect_equal(sim$generator_metadata$y$selected_column_roles$column, "stress_between")
})

test_that("gen_outcome supports residual VAR and compositional output", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    group_id = "ID",
    time_id = "day",
    seed = 102,
    generators = list(
      treatment = gen_categorical(
        "treatment",
        level = "level2",
        categories = c("control", "treatment"),
        fixed_intercept = stats::qlogis(0.5),
        output = "factor"
      ),
      stress = gen_mvn("stress", level = "multilevel", fixed_intercept = 0, random_cov = 0.2, residual_cov = 0.8),
      outcome = gen_outcome(
        mvbind(ilr1, ilr2) ~ treatment + between(stress) + within(stress) +
          ar1() + treatment:ar1() + (1 + ar1() | ID),
        scale = sigma ~ treatment + (1 | ID),
        params = outcome_acceptance_params(),
        burnin = 10,
        composition = list(parts = c("sleep", "sedentary", "activity"), total = 24)
      )
    )
  )

  expect_true(all(c("ilr1", "ilr2", "sleep", "sedentary", "activity") %in% names(sim$data)))
  expect_true(all(sim$data$sleep > 0 & sim$data$sedentary > 0 & sim$data$activity > 0))
  expect_equal(rowSums(as.matrix(sim$data[, c("sleep", "sedentary", "activity"), with = FALSE])), rep(24, nrow(sim$data)))
  expect_equal(dim(sim$generator_metadata$outcome$ar$phi_by_row), c(nrow(sim$data), 2L, 2L))
  expect_true("treatmenttreatment:ar1()" %in% sim$generator_metadata$outcome$expected_parameter_names$ar)
  expect_false(is.null(sim$generator_metadata$outcome$ar$stability$stability_attempts_by_group_level))

  attempts <- sim$generator_metadata$outcome$ar$stability$stability_attempts_by_group_level
  expect_equal(
    sim$generator_metadata$outcome$ar$stability$stability_acceptance_rate,
    length(attempts) / sum(attempts)
  )

  phi_by_term <- sim$generator_metadata$outcome$ar$phi_by_group_and_term
  expect_equal(names(phi_by_term), sim$generator_metadata$outcome$expected_parameter_names$ar)
  expect_equal(dim(phi_by_term[["ar1()"]]), c(4L, 2L, 2L))
  u <- sim$generator_metadata$outcome$group_level_effects
  expect_equal(
    phi_by_term[["ar1()"]][1L, "ilr1", "ilr2"],
    outcome_acceptance_params()$ar$beta["ar1()", "ilr1", "ilr2"] +
      u[1L, "ar|term=ar1()|to=ilr1|from=ilr2"]
  )

  expect_identical(
    sim$generator_metadata$outcome$basis_package_version,
    as.character(utils::packageVersion("compositions"))
  )
})

test_that("gen_outcome validates dynamic formula restrictions", {
  beta_location <- matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))
  beta_scale <- matrix(log(0.1), nrow = 1, dimnames = list("(Intercept)", "y"))
  params <- list(location = list(beta = beta_location), scale = list(beta = beta_scale))

  expect_error(
    simulate_data(
      n = 3,
      generators = list(y = gen_outcome(y ~ 1, scale = sigma ~ ar1(), params = params, burnin = 0))
    ),
    "ar1\\(\\) is not allowed in the scale formula"
  )

  beta_ar <- array(1.1, dim = c(1, 1, 1), dimnames = list("ar1()", "y", "y"))
  expect_error(
    simulate_data(
      n = 5,
      time_id = "time",
      generators = list(
        y = gen_outcome(
          y ~ ar1(),
          scale = sigma ~ 1,
          params = list(location = list(beta = beta_location), scale = list(beta = beta_scale), ar = list(beta = beta_ar)),
          burnin = 1
        )
      )
    ),
    "Population-level AR matrix is unstable"
  )

  expect_error(
    simulate_data(
      n = 5,
      time_id = "time",
      generators = list(
        outcome_template = gen_template(
          y ~ ar1():ar1(),
          scale = sigma ~ 1,
          burnin = 0
        )
      )
    ),
    "Terms containing more than one `ar1\\(\\)` are not supported"
  )

  expect_error(
    simulate_data(
      n_groups = 2,
      n_per_group = 3,
      group_id = "ID",
      time_id = "time",
      generators = list(
        treatment = gen_categorical(
          "treatment",
          level = "level2",
          categories = c("control", "treatment"),
          fixed_intercept = stats::qlogis(0.5),
          output = "factor"
        ),
        outcome_template = gen_template(
          y ~ treatment:ar1():ar1(),
          scale = sigma ~ 1,
          burnin = 0
        )
      )
    ),
    "Terms containing more than one `ar1\\(\\)` are not supported"
  )

  valid_interaction <- simulate_data(
    n_groups = 2,
    n_per_group = 3,
    group_id = "ID",
    time_id = "time",
    generators = list(
      treatment = gen_categorical(
        "treatment",
        level = "level2",
        categories = c("control", "treatment"),
        fixed_intercept = stats::qlogis(0.5),
        output = "factor"
      ),
      outcome_template = gen_template(
        y ~ treatment * ar1(),
        scale = sigma ~ 1,
        burnin = 0
      )
    )
  )
  expect_equal(
    valid_interaction$generator_metadata$outcome_template$expected_parameter_names$ar,
    c("ar1()", "treatmenttreatment:ar1()")
  )
})

test_that("gen_outcome rejects between columns that vary within groups", {
  params <- list(
    location = list(beta = matrix(c(0, 1), nrow = 2, dimnames = list(c("(Intercept)", "between(x)"), "y"))),
    scale = list(beta = matrix(log(0.1), nrow = 1, dimnames = list("(Intercept)", "y")))
  )
  varying_between <- gen_custom("x", function(context, vars, level) {
    list(
      data = data.frame(x = seq_len(context$n_rows)),
      names = "x",
      metadata = list(column_roles = data.frame(
        column = "x",
        variable = "x",
        component = "between",
        level = "group"
      ))
    )
  })

  expect_error(
    simulate_data(
      n_groups = 3,
      n_per_group = 4,
      seed = 401,
      generators = list(
        x = varying_between,
        y = gen_outcome(y ~ between(x), scale = sigma ~ 1, params = params, burnin = 0)
      )
    ),
    "constant within each active group"
  )
})

test_that("gen_outcome rejects helper terms inside unsupported calls", {
  params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.1), nrow = 1, dimnames = list("(Intercept)", "y")))
  )

  expect_error(
    simulate_data(
      n_groups = 3,
      n_per_group = 4,
      seed = 402,
      generators = list(
        stress = gen_mvn("stress", level = "multilevel", fixed_intercept = 0, random_cov = 0.2, residual_cov = 0.8),
        y = gen_outcome(y ~ I(within(stress)^2), scale = sigma ~ 1, params = params, burnin = 0)
      )
    ),
    "`within\\(\\)` terms may only appear as main effects or inside interactions"
  )
})

test_that("group-level ar1() requires a matching population-level ar1() term", {
  params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.1), nrow = 1, dimnames = list("(Intercept)", "y")))
  )

  expect_error(
    simulate_data(
      n_groups = 3,
      n_per_group = 4,
      group_id = "ID",
      time_id = "time",
      seed = 403,
      generators = list(
        y = gen_outcome(y ~ (1 + ar1() | ID), scale = sigma ~ 1, params = params, burnin = 0)
      )
    ),
    "matching population-level terms in the same component; missing: ar1\\(\\)\\. Add the matching `ar1\\(\\)` term"
  )
})

test_that("gen_outcome records high-persistence warnings without random effects", {
  make_params <- function(phi) {
    list(
      location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
      scale = list(beta = matrix(log(0.1), nrow = 1, dimnames = list("(Intercept)", "y"))),
      ar = list(beta = array(phi, dim = c(1, 1, 1), dimnames = list("ar1()", "y", "y")))
    )
  }

  high_persistence <- NULL
  expect_warning(
    high_persistence <- simulate_data(
      n = 8,
      time_id = "time",
      seed = 144,
      generators = list(
        outcome = gen_outcome(
          y ~ ar1(),
          scale = sigma ~ 1,
          params = make_params(0.96),
          burnin = 2
        )
      )
    ),
    "spectral radius above 0.95"
  )
  high_stability <- high_persistence$generator_metadata$outcome$ar$stability
  expect_true(high_stability$high_persistence_warning)
  expect_equal(high_stability$max_spectral_radius_overall, 0.96)

  low_persistence <- NULL
  expect_warning(
    low_persistence <- simulate_data(
      n = 8,
      time_id = "time",
      seed = 145,
      generators = list(
        outcome = gen_outcome(
          y ~ ar1(),
          scale = sigma ~ 1,
          params = make_params(0.2),
          burnin = 2
        )
      )
    ),
    NA
  )
  low_stability <- low_persistence$generator_metadata$outcome$ar$stability
  expect_false(low_stability$high_persistence_warning)
  expect_equal(low_stability$max_spectral_radius_overall, 0.2)
})

test_that("gen_outcome validates AR time indices before simulation", {
  beta_location <- matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))
  beta_scale <- matrix(log(0.1), nrow = 1, dimnames = list("(Intercept)", "y"))
  beta_ar <- array(0.2, dim = c(1, 1, 1), dimnames = list("ar1()", "y", "y"))
  params <- list(
    location = list(beta = beta_location),
    scale = list(beta = beta_scale),
    ar = list(beta = beta_ar)
  )
  generator <- gen_outcome(
    y ~ ar1(),
    scale = sigma ~ 1,
    params = params,
    burnin = 1
  )

  expect_error(
    simulate_data(n = 4, generators = list(outcome = generator)),
    "requires `simulate_data\\(\\)` to supply `time_id`"
  )
  expect_error(
    simulate_data(
      n = 4,
      time_id = "time",
      time_values = letters[1:4],
      generators = list(outcome = generator)
    ),
    "requires numeric, Date, or POSIXct time values"
  )
  expect_error(
    simulate_data(
      n = 4,
      time_id = "time",
      time_values = c(1, 2, 4, 5),
      generators = list(outcome = generator)
    ),
    "requires complete equally spaced time"
  )
})

test_that("gen_outcome can shrink unstable group-level AR draws", {
  random_names <- c(
    "location|outcome=y|term=(Intercept)",
    "ar|term=ar1()|to=y|from=y"
  )
  random_cov <- diag(c(0, 1))
  dimnames(random_cov) <- list(random_names, random_names)
  params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.1), nrow = 1, dimnames = list("(Intercept)", "y"))),
    ar = list(beta = array(0, dim = c(1, 1, 1), dimnames = list("ar1()", "y", "y"))),
    random = list(ID = list(covariance = random_cov))
  )

  sim <- simulate_data(
    n_groups = 3,
    n_per_group = 4,
    group_id = "ID",
    time_id = "day",
    seed = 104,
    generators = list(
      outcome = gen_outcome(
        y ~ ar1() + (ar1() | ID),
        scale = sigma ~ 1,
        params = params,
        burnin = 1,
        ar_stability = "shrink",
        max_stability_attempts = 500,
        shrink_target_radius = 0.05
      )
    )
  )

  stability <- sim$generator_metadata$outcome$ar$stability
  expect_identical(stability$stability_rule, "shrink")
  expect_true(any(stability$shrink_factor_by_group_level < 1))
  expect_true(all(stability$spectral_radius_after_shrinkage < 0.05))
  expect_lt(stability$max_spectral_radius_overall, 0.05)
})

test_that("gen_outcome requires exact parameter order", {
  beta_location <- matrix(
    0,
    nrow = 2,
    dimnames = list(c("between(stress)", "(Intercept)"), "y")
  )
  beta_scale <- matrix(log(0.1), nrow = 1, dimnames = list("(Intercept)", "y"))

  expect_error(
    simulate_data(
      n_groups = 2,
      n_per_group = 3,
      seed = 103,
      generators = list(
        stress = gen_mvn("stress", level = "multilevel", fixed_intercept = 0, random_cov = 0.2, residual_cov = 0.8),
        y = gen_outcome(
          y ~ between(stress),
          scale = sigma ~ 1,
          params = list(location = list(beta = beta_location), scale = list(beta = beta_scale)),
          burnin = 0
        )
      )
    ),
    "names must match exactly"
  )
})

test_that("gen_outcome supports univariate outcome and design combinations", {
  make_params <- function(outcome, grouped) {
    params <- list(
      location = list(beta = matrix(
        c(0.1, 0.2),
        nrow = 2,
        dimnames = list(c("(Intercept)", "x"), outcome)
      )),
      scale = list(beta = matrix(
        log(0.15),
        nrow = 1,
        dimnames = list("(Intercept)", outcome)
      ))
    )
    if (isTRUE(grouped)) {
      random_name <- sprintf("location|outcome=%s|term=(Intercept)", outcome)
      random_cov <- matrix(0.02, nrow = 1, dimnames = list(random_name, random_name))
      params$random <- list(ID = list(covariance = random_cov))
    }
    params
  }

  cases <- list(
    list(label = "single non-compositional", grouped = FALSE, compositional = FALSE, outcome = "y"),
    list(label = "single two-part composition", grouped = FALSE, compositional = TRUE, outcome = "ilr"),
    list(label = "multilevel non-compositional", grouped = TRUE, compositional = FALSE, outcome = "y"),
    list(label = "multilevel two-part composition", grouped = TRUE, compositional = TRUE, outcome = "ilr")
  )

  for (i in seq_along(cases)) {
    case <- cases[[i]]
    params <- make_params(case$outcome, case$grouped)
    formula <- if (isTRUE(case$grouped)) {
      stats::as.formula(sprintf("%s ~ x + (1 | ID)", case$outcome))
    } else {
      stats::as.formula(sprintf("%s ~ x", case$outcome))
    }
    outcome_generator <- if (isTRUE(case$compositional)) {
      gen_outcome(
        formula,
        scale = sigma ~ 1,
        params = params,
        burnin = 0,
        composition = list(parts = c("part_a", "part_b"), total = 24)
      )
    } else {
      gen_outcome(
        formula,
        scale = sigma ~ 1,
        params = params,
        burnin = 0
      )
    }

    sim <- if (isTRUE(case$grouped)) {
      simulate_data(
        n_groups = 3,
        n_per_group = 4,
        group_id = "ID",
        seed = 200 + i,
        generators = list(
          x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
          outcome = outcome_generator
        )
      )
    } else {
      simulate_data(
        n = 8,
        seed = 200 + i,
        generators = list(
          x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
          outcome = outcome_generator
        )
      )
    }

    metadata <- sim$generator_metadata$outcome
    expect_identical(metadata$parsed$outcomes, case$outcome, info = case$label)
    expect_equal(metadata$expected_parameter_names$location, c("(Intercept)", "x"), info = case$label)
    expect_identical(metadata$level, if (isTRUE(case$grouped)) "multilevel" else "single", info = case$label)

    if (isTRUE(case$grouped)) {
      expect_equal(
        metadata$parsed$random_effect_names,
        sprintf("location|outcome=%s|term=(Intercept)", case$outcome),
        info = case$label
      )
      expect_equal(dim(metadata$group_level_effects), c(3L, 1L), info = case$label)
    } else {
      expect_null(metadata$group_level_effects, info = case$label)
    }

    if (isTRUE(case$compositional)) {
      expect_true(all(c(case$outcome, "part_a", "part_b") %in% names(sim$data)), info = case$label)
      expect_true(metadata$compositional, info = case$label)
      expect_true(all(sim$data$part_a > 0 & sim$data$part_b > 0), info = case$label)
      expect_equal(
        rowSums(as.matrix(sim$data[, c("part_a", "part_b"), with = FALSE])),
        rep(24, nrow(sim$data)),
        tolerance = 1e-8,
        info = case$label
      )
    } else {
      expect_true(case$outcome %in% names(sim$data), info = case$label)
      expect_false(any(c("part_a", "part_b") %in% names(sim$data)), info = case$label)
      expect_false(metadata$compositional, info = case$label)
    }
  }
})

test_that("gen_template emits no columns and creates univariate defaults", {
  sim <- simulate_data(
    n = 4,
    generators = list(
      outcome_template = gen_template(
        y ~ 1,
        scale = sigma ~ 1,
        burnin = 0
      )
    )
  )

  metadata <- sim$generator_metadata$outcome_template
  expect_false("y" %in% names(sim$data))
  expect_identical(metadata$distribution, "outcome_template")
  expect_identical(metadata$vars, character())
  expect_equal(
    metadata$params$location$beta,
    matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))
  )
  expect_equal(
    metadata$params$scale$beta,
    matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))
  )
  expect_equal(
    metadata$params$scale$correlation,
    matrix(1, nrow = 1, dimnames = list("y", "y"))
  )
  expect_null(metadata$params$ar)
  expect_null(metadata$params$random)
})

test_that("gen_template creates exact AR and random parameter names", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    group_id = "ID",
    time_id = "day",
    seed = 104,
    generators = list(
      treatment = gen_categorical(
        "treatment",
        level = "level2",
        categories = c("control", "treatment"),
        fixed_intercept = stats::qlogis(0.5),
        output = "factor"
      ),
      stress = gen_mvn("stress", level = "multilevel", fixed_intercept = 0, random_cov = 0.2, residual_cov = 0.8),
      outcome_template = gen_template(
        mvbind(ilr1, ilr2) ~ treatment + between(stress) + within(stress) +
          ar1() + treatment:ar1() + (1 + ar1() | ID),
        scale = sigma ~ treatment + (1 | ID),
        burnin = 10,
        composition = list(parts = c("sleep", "sedentary", "activity"), total = 24)
      )
    )
  )

  params <- sim$generator_metadata$outcome_template$params
  expected <- sim$generator_metadata$outcome_template$expected_parameter_names
  expect_false(any(c("ilr1", "ilr2", "sleep", "sedentary", "activity") %in% names(sim$data)))
  expect_identical(rownames(params$location$beta), expected$location)
  expect_identical(colnames(params$location$beta), c("ilr1", "ilr2"))
  expect_identical(rownames(params$scale$beta), expected$scale)
  expect_equal(params$scale$correlation, diag(2), ignore_attr = TRUE)
  expect_identical(dimnames(params$scale$correlation), list(c("ilr1", "ilr2"), c("ilr1", "ilr2")))
  expect_identical(dimnames(params$ar$beta), list(expected$ar, c("ilr1", "ilr2"), c("ilr1", "ilr2")))
  expect_true(all(params$ar$beta == 0))
  expect_identical(
    rownames(params$random$ID$covariance),
    rownames(outcome_acceptance_params()$random$ID$covariance)
  )
  expect_true(all(params$random$ID$covariance == 0))
  expect_equal(
    sim$generator_metadata$outcome_template$selected_column_roles$column,
    c("stress_between", "stress_within")
  )
  expect_equal(
    sim$generator_metadata$outcome_template$composition$resolved$parts,
    c("sleep", "sedentary", "activity")
  )
})

test_that("gen_template params can be passed to gen_outcome", {
  template <- simulate_data(
    n_groups = 3,
    n_per_group = 4,
    group_id = "ID",
    time_id = "day",
    seed = 105,
    generators = list(
      treatment = gen_categorical(
        "treatment",
        level = "level2",
        categories = c("control", "treatment"),
        fixed_intercept = stats::qlogis(0.5),
        output = "factor"
      ),
      stress = gen_mvn("stress", level = "multilevel", fixed_intercept = 0, random_cov = 0.2, residual_cov = 0.8),
      outcome_template = gen_template(
        mvbind(ilr1, ilr2) ~ treatment + between(stress) + within(stress) + ar1(),
        scale = sigma ~ treatment,
        burnin = 3,
        composition = list(parts = c("sleep", "sedentary", "activity"), total = 24)
      )
    )
  )
  params <- template$generator_metadata$outcome_template$params
  params$location$beta["between(stress)", "ilr1"] <- 0.2
  params$scale$beta["(Intercept)", ] <- log(c(0.2, 0.3))
  params$ar$beta["ar1()", , ] <- diag(c(0.1, 0.2))

  sim <- simulate_data(
    n_groups = 3,
    n_per_group = 4,
    group_id = "ID",
    time_id = "day",
    seed = 106,
    generators = list(
      treatment = gen_categorical(
        "treatment",
        level = "level2",
        categories = c("control", "treatment"),
        fixed_intercept = stats::qlogis(0.5),
        output = "factor"
      ),
      stress = gen_mvn("stress", level = "multilevel", fixed_intercept = 0, random_cov = 0.2, residual_cov = 0.8),
      outcome = gen_outcome(
        mvbind(ilr1, ilr2) ~ treatment + between(stress) + within(stress) + ar1(),
        scale = sigma ~ treatment,
        params = params,
        burnin = 3,
        composition = list(parts = c("sleep", "sedentary", "activity"), total = 24)
      )
    )
  )

  expect_true(all(c("ilr1", "ilr2", "sleep", "sedentary", "activity") %in% names(sim$data)))
  expect_equal(rowSums(as.matrix(sim$data[, c("sleep", "sedentary", "activity"), with = FALSE])), rep(24, nrow(sim$data)))
})

test_that("gen_template uses prior column roles for between and within terms", {
  expect_error(
    simulate_data(
      n_groups = 2,
      n_per_group = 3,
      generators = list(
        outcome_template = gen_template(
          y ~ between(stress),
          scale = sigma ~ 1,
          burnin = 0
        )
      )
    ),
    "requires exactly one prior generator column role"
  )
})

test_that("prep_sim_analysis recomputes observed between and within predictors", {
  params <- list(
    location = list(beta = matrix(
      c(0, 0.3, 0.2),
      nrow = 3,
      dimnames = list(c("(Intercept)", "between(x)", "within(x)"), "y")
    )),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y")))
  )
  sim <- simulate_data(
    n_groups = 3,
    n_per_group = 4,
    group_id = "ID",
    seed = 300,
    generators = list(
      x = gen_mvn("x", level = "multilevel", fixed_intercept = 0, random_cov = 0.3, residual_cov = 0.7),
      outcome = gen_outcome(
        y ~ between(x) + within(x),
        scale = sigma ~ 1,
        params = params,
        burnin = 0
      )
    )
  )
  sim$data$x_between <- 999
  sim$data$x_within <- -999

  analysis <- prep_sim_analysis(sim)
  expected <- data.table::copy(sim$data)
  expected[, x_between := mean(x), by = ID]
  expected[, x_within := x - x_between]

  expect_true(inherits(analysis, "mlsim_analysis"))
  expect_null(analysis$complr)
  expect_true(inherits(analysis$formula, "brmsformula"))
  expect_equal(analysis$data$x_between, expected$x_between)
  expect_equal(analysis$data$x_within, expected$x_within)
  expect_true(all(analysis$metadata$derived_roles$overwritten))
})

test_that("prep_sim_analysis creates complr objects for two-part compositions", {
  params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "ilr"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "ilr")))
  )
  sim <- simulate_data(
    n = 6,
    seed = 301,
    generators = list(
      outcome = gen_outcome(
        ilr ~ 1,
        scale = sigma ~ 1,
        params = params,
        burnin = 0,
        composition = list(parts = c("part_a", "part_b"), total = 24)
      )
    )
  )

  analysis <- prep_sim_analysis(sim)

  expect_true(is.complr(analysis$complr))
  expect_true("z1_1" %in% names(analysis$data))
  expect_identical(unname(analysis$metadata$response_map), "z1_1")
  expect_match(paste(deparse(analysis$formula), collapse = " "), "z1_1", fixed = TRUE)
})

test_that("prep_sim_analysis creates within lags for single-level AR outcomes", {
  params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y"))),
    ar = list(beta = array(0.2, dim = c(1, 1, 1), dimnames = list("ar1()", "y", "y")))
  )
  sim <- simulate_data(
    n = 6,
    time_id = "time",
    seed = 302,
    generators = list(
      outcome = gen_outcome(
        y ~ ar1(),
        scale = sigma ~ 1,
        params = params,
        burnin = 2
      )
    )
  )

  analysis <- prep_sim_analysis(sim)
  dropped <- prep_sim_analysis(sim, drop_lag_na = TRUE)
  formula_text <- paste(deparse(analysis$formula), collapse = " ")

  expect_equal(nrow(analysis$data), nrow(sim$data))
  expect_true("lag_y_within" %in% names(analysis$data))
  expect_true(is.na(analysis$data$lag_y_within[[1L]]))
  expect_equal(
    analysis$data$lag_y_within[-1L],
    head(sim$data$y, -1L) - mean(sim$data$y)
  )
  expect_equal(nrow(dropped$data), nrow(sim$data) - 1L)
  expect_match(formula_text, "lag_y_within", fixed = TRUE)
})

test_that("prep_sim_analysis translates dynamic compositional formulas", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    group_id = "ID",
    time_id = "day",
    seed = 303,
    generators = list(
      treatment = gen_categorical(
        "treatment",
        level = "level2",
        categories = c("control", "treatment"),
        fixed_intercept = stats::qlogis(0.5),
        output = "factor"
      ),
      stress = gen_mvn("stress", level = "multilevel", fixed_intercept = 0, random_cov = 0.2, residual_cov = 0.8),
      outcome = gen_outcome(
        mvbind(ilr1, ilr2) ~ treatment + between(stress) + within(stress) +
          ar1() + treatment:ar1() + (1 + ar1() | ID),
        scale = sigma ~ treatment + (1 | ID),
        params = outcome_acceptance_params(),
        burnin = 4,
        composition = list(parts = c("sleep", "sedentary", "activity"), total = 24)
      )
    )
  )

  analysis <- prep_sim_analysis(sim)
  dropped <- prep_sim_analysis(sim, drop_lag_na = TRUE)
  formula_text <- paste(deparse(analysis$formula), collapse = " ")
  first_rows <- analysis$data[, .I[1L], by = ID]$V1

  expect_true(inherits(analysis$formula, "mvbrmsformula"))
  expect_true(is.complr(analysis$complr))
  expect_equal(nrow(analysis$data), nrow(sim$data))
  expect_equal(nrow(dropped$data), nrow(sim$data) - length(unique(sim$data$ID)))
  expect_true(all(c("lag_z1_1_within", "lag_z2_1_within") %in% names(analysis$complr$dataout)))
  expect_true(all(is.na(analysis$data$lag_z1_1_within[first_rows])))
  expect_match(formula_text, "stress_between", fixed = TRUE)
  expect_match(formula_text, "stress_within", fixed = TRUE)
  expect_match(formula_text, "treatment:lag_z1_1_within", fixed = TRUE)
  expect_match(formula_text, "1 + lag_z1_1_within + lag_z2_1_within | ID", fixed = TRUE)
})

test_that("prep_sim_analysis reports missing or ambiguous outcome context", {
  expect_error(
    prep_sim_analysis(simulate_data(
      n = 3,
      generators = list(x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1))
    )),
    "Could not infer the outcome generator"
  )

  params_y1 <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y1"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y1")))
  )
  params_y2 <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y2"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y2")))
  )
  ambiguous <- simulate_data(
    n = 3,
    generators = list(
      y1 = gen_outcome(y1 ~ 1, scale = sigma ~ 1, params = params_y1, burnin = 0),
      y2 = gen_outcome(y2 ~ 1, scale = sigma ~ 1, params = params_y2, burnin = 0)
    )
  )
  expect_error(prep_sim_analysis(ambiguous), "found 2 outcome generators")
  expect_true(inherits(prep_sim_analysis(ambiguous, outcome = "y1"), "mlsim_analysis"))

  params_between <- list(
    location = list(beta = matrix(0, nrow = 2, dimnames = list(c("(Intercept)", "between(x)"), "y"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y")))
  )
  grouped <- simulate_data(
    n_groups = 2,
    n_per_group = 3,
    group_id = "ID",
    seed = 304,
    generators = list(
      x = gen_mvn("x", level = "multilevel", fixed_intercept = 0, random_cov = 0.2, residual_cov = 0.8),
      outcome = gen_outcome(y ~ between(x), scale = sigma ~ 1, params = params_between, burnin = 0)
    )
  )
  grouped$data[, ID := NULL]
  expect_error(prep_sim_analysis(grouped), "grouping column")

  params_ar <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y"))),
    ar = list(beta = array(0.1, dim = c(1, 1, 1), dimnames = list("ar1()", "y", "y")))
  )
  ar_sim <- simulate_data(
    n = 5,
    time_id = "time",
    seed = 305,
    generators = list(
      outcome = gen_outcome(y ~ ar1(), scale = sigma ~ 1, params = params_ar, burnin = 1)
    )
  )
  ar_sim$metadata$time_id <- NULL
  expect_error(prep_sim_analysis(ar_sim), "time_id")
})

test_that(".mlsim_phi_radii matches direct eigendecomposition and handles duplicates", {
  m1 <- matrix(c(0.5, 0.1, -0.2, 0.3), 2, 2)
  m2 <- matrix(c(0.9, 0, 0.4, -0.4), 2, 2)
  phi <- array(NA_real_, c(5, 2, 2))
  phi[1, , ] <- m1
  phi[2, , ] <- m2
  phi[3, , ] <- m1
  phi[4, , ] <- m2
  phi[5, , ] <- m1
  direct <- vapply(seq_len(5), function(i) {
    max(Mod(eigen(phi[i, , ], only.values = TRUE)$values))
  }, numeric(1))
  expect_equal(multilevelcoda:::.mlsim_phi_radii(phi), direct)

  set.seed(42)
  uni <- array(rnorm(7), c(7, 1, 1))
  expect_equal(multilevelcoda:::.mlsim_phi_radii(uni), abs(as.vector(uni)))
})

test_that(".mlsim_outcome_phi row subsetting matches full assembly", {
  set.seed(7)
  n <- 12L
  outcomes <- c("y1", "y2")
  Q <- cbind(rep(1, n), rnorm(n))
  colnames(Q) <- c("ar1()", "x:ar1()")
  group_index <- rep(1:3, each = 4L)
  spec <- list(
    X = matrix(1, n, 1L),
    Q = Q,
    has_ar = TRUE,
    outcomes = outcomes,
    random = list(
      ar_terms = "ar1()",
      ar_Z = matrix(1, n, 1L, dimnames = list(NULL, "ar1()"))
    ),
    series = list(group_index = group_index)
  )
  beta <- array(
    rnorm(8, sd = 0.1),
    dim = c(2, 2, 2),
    dimnames = list(colnames(Q), outcomes, outcomes)
  )
  checked <- list(ar = list(beta = beta))
  ar_names <- multilevelcoda:::.mlsim_random_effect_names("ar", "ar1()", outcomes)
  draws <- matrix(
    rnorm(3L * length(ar_names), sd = 0.05),
    nrow = 3L,
    ncol = length(ar_names),
    dimnames = list(NULL, ar_names)
  )

  full <- multilevelcoda:::.mlsim_outcome_phi(spec, checked, draws)
  rows <- which(group_index == 2L)
  sub <- multilevelcoda:::.mlsim_outcome_phi(spec, checked, draws, rows = rows)
  expect_equal(sub$phi, full$phi[rows, , , drop = FALSE])
  expect_equal(sub$fixed_phi, full$fixed_phi[rows, , , drop = FALSE])
})

test_that("simulated values satisfy the exact residual VAR recursion", {
  params <- outcome_acceptance_params()
  params$scale$correlation[1L, 2L] <- 0.4
  params$scale$correlation[2L, 1L] <- 0.4
  re <- params$random$ID$covariance
  diag(re) <- c(0.2, 0.15, 0.01, 0, 0, 0.01, 0.05, 0.05)
  params$random$ID$covariance <- re

  sim <- simulate_data(
    n_groups = 8,
    n_per_group = 60,
    group_id = "ID",
    time_id = "day",
    seed = 606,
    generators = list(
      treatment = gen_categorical(
        "treatment",
        level = "level2",
        categories = c("control", "treatment"),
        fixed_intercept = stats::qlogis(0.5),
        output = "factor"
      ),
      stress = gen_mvn("stress", level = "multilevel", fixed_intercept = 0, random_cov = 0.2, residual_cov = 0.8),
      outcome = gen_outcome(
        mvbind(ilr1, ilr2) ~ treatment + between(stress) + within(stress) +
          ar1() + treatment:ar1() + (1 + ar1() | ID),
        scale = sigma ~ treatment + (1 | ID),
        params = params,
        burnin = 5,
        composition = list(parts = c("sleep", "sedentary", "activity"), total = 24)
      )
    )
  )

  md <- sim$generator_metadata$outcome
  z <- as.matrix(sim$data[, c("ilr1", "ilr2"), with = FALSE])

  # Generated outcomes decompose exactly into location plus residual state.
  expect_equal(unname(z), unname(md$mu + md$residual))

  # Residuals follow the defining recursion e_t = Phi_t e_{t-1} + eps_t exactly.
  phi <- md$ar$phi_by_row
  groups <- as.integer(sim$data$ID)
  for (g in unique(groups)) {
    rows <- which(groups == g)
    for (i in seq_along(rows)[-1L]) {
      r <- rows[[i]]
      p <- rows[[i - 1L]]
      expect_equal(
        unname(md$residual[r, ]),
        unname(as.vector(phi[r, , ] %*% md$residual[p, ]) + md$innovation[r, ]),
        tolerance = 1e-12
      )
    }
  }

  # Innovations divided by their row-wise conditional SDs are standardized
  # draws with the requested conditional correlation.
  scaled <- md$innovation / md$sigma
  expect_lt(max(abs(apply(scaled, 2L, stats::sd) - 1)), 0.1)
  expect_lt(abs(stats::cor(scaled)[1L, 2L] - 0.4), 0.1)

  # Stored stability radii match direct per-row eigendecompositions.
  radii_direct <- vapply(seq_len(nrow(z)), function(i) {
    max(Mod(eigen(phi[i, , ], only.values = TRUE)$values))
  }, numeric(1))
  tab <- md$ar$stability$spectral_radius_by_group_level_and_row
  expect_equal(tab$spectral_radius, radii_direct[tab$row])
  expect_equal(
    md$ar$stability$max_spectral_radius_by_group_level,
    unname(as.numeric(tapply(radii_direct[tab$row], tab$group, max)))
  )
})
