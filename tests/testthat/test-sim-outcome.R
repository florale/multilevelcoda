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
  expect_equal(nrow(dropped$complr$datain), nrow(dropped$data))
  expect_equal(nrow(dropped$complr$output[[1L]]$Z), nrow(dropped$data))
  expect_equal(nrow(dropped$complr$output[[1L]]$dataout), nrow(dropped$data))
  expect_true(all(c("lag_z1_1_within", "lag_z2_1_within") %in% names(analysis$complr$dataout)))
  expect_true(all(is.na(analysis$data$lag_z1_1_within[first_rows])))
  expect_match(formula_text, "stress_between", fixed = TRUE)
  expect_match(formula_text, "stress_within", fixed = TRUE)
  expect_match(formula_text, "treatment:lag_z1_1_within", fixed = TRUE)
  # the grouping factor appears in both the mean and sigma formulas, so the
  # analysis model links the random effects across formulas to estimate their
  # correlation (the simulator draws them from one joint covariance)
  expect_match(formula_text, "1 + lag_z1_1_within + lag_z2_1_within | p1 | ID", fixed = TRUE)
  expect_match(formula_text, "sigma = sigma ~ treatment + (1 | p1 | ID)", fixed = TRUE)
})

test_that("prep_sim_analysis lag_center = none uses raw lags", {
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

  analysis <- prep_sim_analysis(sim, lag_center = "none")
  default_analysis <- prep_sim_analysis(sim)
  formula_text <- paste(deparse(analysis$formula), collapse = " ")

  expect_true("lag_y" %in% names(analysis$data))
  expect_false("lag_y_within" %in% names(analysis$data))
  expect_true(is.na(analysis$data$lag_y[[1L]]))
  expect_equal(analysis$data$lag_y[-1L], head(sim$data$y, -1L))
  expect_match(formula_text, "lag_y", fixed = TRUE)
  expect_false(grepl("lag_y_within", formula_text, fixed = TRUE))
  expect_identical(analysis$metadata$lag_center, "none")
  expect_identical(analysis$metadata$lag_columns, "lag_y")
  expect_identical(default_analysis$metadata$lag_center, "within")
  expect_error(prep_sim_analysis(sim, lag_center = "person"), "arg")
})

test_that("prep_sim_analysis errors when lag columns collide with existing columns", {
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
  sim$data[, lag_y := 0]

  expect_error(prep_sim_analysis(sim, lag_center = "none"), "lag_y")
  expect_no_error(prep_sim_analysis(sim))

  sim$data[, lag_y := NULL]
  sim$data[, lag_y_within := 0]
  expect_error(prep_sim_analysis(sim), "lag_y_within")
  expect_no_error(prep_sim_analysis(sim, lag_center = "none"))
})

test_that("prep_sim_analysis lag_center = none translates dynamic compositional formulas", {
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

  analysis <- prep_sim_analysis(sim, lag_center = "none")
  formula_text <- paste(deparse(analysis$formula), collapse = " ")
  first_rows <- analysis$data[, .I[1L], by = ID]$V1
  shifted <- analysis$data[, shift(z1_1), by = ID]$V1

  expect_true(all(c("lag_z1_1", "lag_z2_1") %in% names(analysis$complr$dataout)))
  expect_false(any(c("lag_z1_1_within", "lag_z2_1_within") %in% names(analysis$data)))
  expect_true(all(is.na(analysis$data$lag_z1_1[first_rows])))
  expect_equal(analysis$data$lag_z1_1, shifted)
  expect_match(formula_text, "treatment:lag_z1_1", fixed = TRUE)
  expect_match(formula_text, "1 + lag_z1_1 + lag_z2_1 | p1 | ID", fixed = TRUE)
})

test_that("prep_sim_analysis lags by time so gaps yield NA in grouped data", {
  params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y"))),
    ar = list(beta = array(0.2, dim = c(1, 1, 1), dimnames = list("ar1()", "y", "y")))
  )
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    group_id = "ID",
    time_id = "day",
    seed = 306,
    generators = list(
      outcome = gen_outcome(y ~ ar1(), scale = sigma ~ 1, params = params, burnin = 2)
    )
  )
  gap_id <- sim$data$ID[[1L]]
  # a skipped day, as in real intensive longitudinal data
  sim$data <- sim$data[!(ID == gap_id & day == 3)]

  analysis <- prep_sim_analysis(sim)
  dropped <- prep_sim_analysis(sim, drop_lag_na = TRUE)
  gap_rows <- analysis$data[ID == gap_id]
  other <- analysis$data[ID != gap_id]

  # the day after the gap must not silently lag across it
  expect_true(is.na(gap_rows[day == 4, lag_y_within]))
  expect_equal(
    gap_rows[day == 5, lag_y_within],
    gap_rows[day == 4, y] - mean(gap_rows$y)
  )
  expect_equal(
    other[, lag_y_within[-1L] - (head(y, -1L) - mean(y)), by = ID]$V1,
    rep(0, 3L * 4L)
  )
  expect_true(all(is.na(analysis$data[, lag_y_within[1L], by = ID]$V1)))
  # one first row per group plus the row after the gap
  expect_equal(nrow(dropped$data), nrow(sim$data) - 5L)
  expect_identical(analysis$metadata$time_step, 1)
  expect_identical(analysis$metadata$lag_gaps, 1L)
  expect_identical(analysis$metadata$lag_irregular, 0L)
})

test_that("prep_sim_analysis lags by time so gaps yield NA in single-level data", {
  params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y"))),
    ar = list(beta = array(0.2, dim = c(1, 1, 1), dimnames = list("ar1()", "y", "y")))
  )
  sim <- simulate_data(
    n = 6,
    time_id = "time",
    seed = 307,
    generators = list(
      outcome = gen_outcome(y ~ ar1(), scale = sigma ~ 1, params = params, burnin = 2)
    )
  )
  intact <- prep_sim_analysis(sim)
  explicit <- prep_sim_analysis(sim, time_step = 1)
  expect_equal(explicit$data, intact$data)
  expect_identical(intact$metadata$time_step, 1)
  expect_identical(intact$metadata$lag_gaps, 0L)

  sim$data <- sim$data[time != 3]
  analysis <- prep_sim_analysis(sim)

  expect_true(is.na(analysis$data[time == 4, lag_y_within]))
  expect_equal(
    analysis$data[time == 5, lag_y_within],
    analysis$data[time == 4, y] - mean(analysis$data$y)
  )
  expect_identical(analysis$metadata$lag_gaps, 1L)
})

test_that("prep_sim_analysis keeps complr aligned when gaps drop lag rows", {
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    group_id = "ID",
    time_id = "day",
    seed = 308,
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
  gap_id <- sim$data$ID[[1L]]
  sim$data <- sim$data[!(ID == gap_id & day == 3)]

  analysis <- prep_sim_analysis(sim)
  dropped <- prep_sim_analysis(sim, drop_lag_na = TRUE)

  expect_true(all(is.na(
    analysis$data[ID == gap_id & day == 4, c("lag_z1_1_within", "lag_z2_1_within"), with = FALSE]
  )))
  # one first row per group plus the post-gap row
  expect_equal(nrow(dropped$data), nrow(sim$data) - 5L)
  expect_equal(nrow(dropped$complr$datain), nrow(dropped$data))
  expect_equal(nrow(dropped$complr$output[[1L]]$Z), nrow(dropped$data))
  expect_equal(nrow(dropped$complr$output[[1L]]$dataout), nrow(dropped$data))
  expect_identical(dropped$metadata$lag_gaps, 1L)
  expect_identical(dropped$metadata$time_step, 1)
})

heterogeneous_step_sim <- function(seed = 309, n_groups = 3, n_per_group = 5) {
  params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y"))),
    ar = list(beta = array(0.2, dim = c(1, 1, 1), dimnames = list("ar1()", "y", "y")))
  )
  simulate_data(
    n_groups = n_groups,
    n_per_group = n_per_group,
    group_id = "ID",
    time_id = "day",
    seed = seed,
    generators = list(
      outcome = gen_outcome(y ~ ar1(), scale = sigma ~ 1, params = params, burnin = 2)
    )
  )
}

test_that("prep_sim_analysis warns when group time steps are heterogeneous", {
  sim <- heterogeneous_step_sim()
  wide_id <- sim$data$ID[[1L]]
  # one participant observed every two days among daily participants
  sim$data[ID == wide_id, day := day * 2L]

  w <- testthat::capture_warnings(analysis <- prep_sim_analysis(sim))
  expect_length(w, 2L)
  expect_match(w[[1L]], "spaced more widely")
  expect_match(w[[1L]], as.character(wide_id), fixed = TRUE)
  expect_match(w[[2L]], "All `ar1\\(\\)` lag values are `NA`")
  expect_match(w[[2L]], as.character(wide_id), fixed = TRUE)

  expect_true(analysis$metadata$time_step_heterogeneous)
  expect_identical(analysis$metadata$time_step, 1)
  expect_equal(
    unname(analysis$metadata$time_step_by_group[as.character(wide_id)]),
    2
  )
  expect_identical(analysis$metadata$lag_na_groups, as.character(wide_id))
  expect_true(all(is.na(analysis$data[ID == wide_id, lag_y_within])))
  # the daily groups keep their usual single first-row NA
  expect_true(all(
    analysis$data[ID != wide_id, sum(is.na(lag_y_within)), by = ID]$V1 == 1L
  ))

  # an explicit step is a deliberate choice: only the consequence is reported
  w_explicit <- testthat::capture_warnings(
    explicit <- prep_sim_analysis(sim, time_step = 1)
  )
  expect_length(w_explicit, 1L)
  expect_match(w_explicit, "All `ar1\\(\\)` lag values are `NA`")
  expect_true(explicit$metadata$time_step_heterogeneous)
})

test_that("prep_sim_analysis stays silent when groups share a step but start apart", {
  sim <- heterogeneous_step_sim(seed = 310)
  shift_id <- sim$data$ID[[1L]]
  sim$data[ID == shift_id, day := day + 10L]

  expect_no_warning(analysis <- prep_sim_analysis(sim))
  expect_false(analysis$metadata$time_step_heterogeneous)
  expect_identical(analysis$metadata$lag_na_groups, character())
  expect_identical(analysis$metadata$dropped_groups, character())
})

test_that("prep_sim_analysis warns when drop_lag_na removes whole groups", {
  sim <- heterogeneous_step_sim(seed = 311)
  gid <- sim$data$ID[[1L]]
  first_day <- min(sim$data[ID == gid, day])
  # a group reduced to a single observation has only NA lags
  sim$data <- sim$data[ID != gid | day == first_day]

  w <- testthat::capture_warnings(analysis <- prep_sim_analysis(sim))
  expect_length(w, 1L)
  expect_match(w, "contribute no usable lagged rows")
  expect_match(w, as.character(gid), fixed = TRUE)
  expect_true(gid %in% analysis$data$ID)
  expect_identical(analysis$metadata$lag_na_groups, as.character(gid))
  expect_identical(analysis$metadata$dropped_groups, character())

  w_drop <- testthat::capture_warnings(
    dropped <- prep_sim_analysis(sim, drop_lag_na = TRUE)
  )
  expect_length(w_drop, 1L)
  expect_match(w_drop, "removed every row")
  expect_match(w_drop, as.character(gid), fixed = TRUE)
  expect_false(gid %in% dropped$data$ID)
  expect_identical(dropped$metadata$dropped_groups, as.character(gid))
})

test_that("prep_sim_analysis group step checks skip single-level data", {
  params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y"))),
    ar = list(beta = array(0.2, dim = c(1, 1, 1), dimnames = list("ar1()", "y", "y")))
  )
  sim <- simulate_data(
    n = 6,
    time_id = "time",
    seed = 312,
    generators = list(
      outcome = gen_outcome(y ~ ar1(), scale = sigma ~ 1, params = params, burnin = 2)
    )
  )
  sim$data <- sim$data[time != 3]

  expect_no_warning(analysis <- prep_sim_analysis(sim))
  expect_no_warning(dropped <- prep_sim_analysis(sim, drop_lag_na = TRUE))
  expect_null(analysis$metadata$time_step_by_group)
  expect_false(analysis$metadata$time_step_heterogeneous)
  expect_identical(dropped$metadata$dropped_groups, character())
})

test_that("prep_sim_analysis step checks handle Date and POSIXct times", {
  sim <- heterogeneous_step_sim(seed = 313)
  wide_id <- sim$data$ID[[1L]]
  sim$data[ID == wide_id, day := day * 2L]
  sim$data[, day := as.Date("2026-01-01") + day]

  w <- testthat::capture_warnings(analysis <- prep_sim_analysis(sim))
  expect_match(w[[1L]], "spaced more widely")
  expect_true(analysis$metadata$time_step_heterogeneous)
  expect_equal(
    unname(analysis$metadata$time_step_by_group[as.character(wide_id)]),
    2
  )

  hourly <- heterogeneous_step_sim(seed = 314)
  hourly$data[, day := as.POSIXct("2026-01-01 00:00:00", tz = "UTC") + day * 3600]
  expect_no_warning(hourly_analysis <- prep_sim_analysis(hourly))
  expect_false(hourly_analysis$metadata$time_step_heterogeneous)
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

test_that("gen_outcome simulates poisson outcomes", {
  beta_location <- matrix(
    c(1, 0.3),
    nrow = 2,
    dimnames = list(c("(Intercept)", "between(stress)"), "y")
  )
  random_cov <- matrix(
    0.1,
    dimnames = list("location|outcome=y|term=(Intercept)", "location|outcome=y|term=(Intercept)")
  )

  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 5,
    group_id = "ID",
    seed = 401,
    generators = list(
      stress = gen_mvn("stress", level = "multilevel", fixed_intercept = 0, random_cov = 0.2, residual_cov = 0.8),
      y = gen_outcome(
        y ~ between(stress) + (1 | ID),
        params = list(
          location = list(beta = beta_location),
          random = list(ID = list(covariance = random_cov))
        ),
        burnin = 0,
        family = "poisson"
      )
    )
  )
  md <- sim$generator_metadata$y

  expect_true(all(sim$data$y >= 0))
  expect_equal(sim$data$y, round(sim$data$y))
  expect_identical(md$family, "poisson")
  expect_identical(md$parsed$family, "poisson")
  expect_null(md$parsed$dpar)
  expect_equal(md$mu, exp(md$eta))
  expect_identical(md$expected_parameter_names$scale, character())
  expect_null(md$sigma)
  expect_null(md$residual)
  expect_null(md$innovation)
  expect_null(md$ar)
  expect_null(md$Sigma_epsilon)
  expect_null(md$ilr_conditional_correlation)
  expect_true("location|outcome=y|term=(Intercept)" %in% md$parsed$random_effect_names)
})

test_that("gen_outcome resolves binomial trials in all shapes", {
  params <- list(
    location = list(beta = matrix(0.4, nrow = 1, dimnames = list("(Intercept)", "y")))
  )
  make_sim <- function(trials, seed) {
    simulate_data(
      n = 40,
      seed = seed,
      generators = list(
        y = gen_outcome(y ~ 1, params = params, burnin = 0, family = "binomial", trials = trials)
      )
    )
  }

  scalar_sim <- make_sim(10, 402)
  fn_sim <- make_sim(function(n) rep(5L, n), 403)
  dist_sim <- make_sim(list(distribution = "uniform", min = 3, max = 9), 404)

  for (sim in list(scalar_sim, fn_sim, dist_sim)) {
    md <- sim$generator_metadata$y
    expect_true("y_trials" %in% names(sim$data))
    expect_equal(sim$data$y_trials, as.integer(sim$data$y_trials))
    expect_true(all(sim$data$y >= 0))
    expect_true(all(sim$data$y <= sim$data$y_trials))
    expect_identical(md$family, "binomial")
    expect_identical(md$trials_column, "y_trials")
    expect_identical(md$parsed$trials_column, "y_trials")
    expect_true("y_trials" %in% md$vars)
    expect_equal(md$prob, stats::plogis(md$eta))
    expect_equal(as.vector(md$mu), as.vector(md$prob) * sim$data$y_trials)
  }
  expect_true(all(scalar_sim$data$y_trials == 10L))
  expect_true(all(fn_sim$data$y_trials == 5L))
  expect_true(all(dist_sim$data$y_trials >= 3L & dist_sim$data$y_trials <= 9L))
})

test_that("gen_outcome simulates negbin, gamma, and beta with dpar scale models", {
  intercept_beta <- function(value) {
    matrix(value, nrow = 1, dimnames = list("(Intercept)", "y"))
  }
  simple_params <- list(
    location = list(beta = intercept_beta(1)),
    scale = list(beta = intercept_beta(log(2)))
  )

  nb_sim <- simulate_data(
    n = 60,
    seed = 405,
    generators = list(
      y = gen_outcome(y ~ 1, scale = shape ~ 1, params = simple_params, burnin = 0, family = "negbin")
    )
  )
  nb_md <- nb_sim$generator_metadata$y
  expect_true(all(nb_sim$data$y >= 0))
  expect_equal(nb_sim$data$y, round(nb_sim$data$y))
  expect_identical(nb_md$dpar$name, "shape")
  expect_identical(nb_md$parsed$dpar, "shape")
  expect_equal(nb_md$dpar$value, exp(nb_md$dpar$eta))

  gamma_sim <- simulate_data(
    n = 60,
    seed = 406,
    generators = list(
      y = gen_outcome(y ~ 1, scale = shape ~ 1, params = simple_params, burnin = 0, family = "gamma")
    )
  )
  expect_true(all(gamma_sim$data$y > 0))
  expect_identical(gamma_sim$generator_metadata$y$dpar$name, "shape")

  beta_names <- c(
    "location|outcome=y|term=(Intercept)",
    "phi|outcome=y|term=(Intercept)"
  )
  beta_cov <- diag(c(0.2, 0.1))
  dimnames(beta_cov) <- list(beta_names, beta_names)
  beta_params <- list(
    location = list(beta = intercept_beta(0)),
    scale = list(beta = intercept_beta(log(15))),
    random = list(ID = list(covariance = beta_cov))
  )
  beta_sim <- simulate_data(
    n_groups = 6,
    n_per_group = 5,
    group_id = "ID",
    seed = 407,
    generators = list(
      y = gen_outcome(
        y ~ 1 + (1 | ID),
        scale = phi ~ 1 + (1 | ID),
        params = beta_params,
        burnin = 0,
        family = "beta"
      )
    )
  )
  beta_md <- beta_sim$generator_metadata$y
  expect_true(all(beta_sim$data$y > 0))
  expect_true(all(beta_sim$data$y < 1))
  expect_identical(beta_md$dpar$name, "phi")
  expect_identical(sort(beta_md$parsed$random_effect_names), sort(beta_names))
  expect_identical(
    colnames(beta_md$group_level_effects),
    c("location|outcome=y|term=(Intercept)", "phi|outcome=y|term=(Intercept)")
  )
})

test_that("gen_outcome validates non-Gaussian restrictions", {
  params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y")))
  )

  expect_error(
    gen_outcome(y ~ 1, params = params, burnin = 0, family = "student"),
    "`family` must be one of"
  )
  expect_error(
    gen_outcome(mvbind(y1, y2) ~ 1, params = params, burnin = 0, family = "poisson"),
    "only supported for `family = \"gaussian\"`"
  )
  expect_error(
    gen_outcome(
      y ~ 1,
      params = params,
      burnin = 0,
      family = "poisson",
      composition = list(parts = c("sleep", "wake"))
    ),
    "only supported for `family = \"gaussian\"`"
  )
  expect_error(
    gen_outcome(y ~ ar1(), params = params, burnin = 0, family = "poisson"),
    "`ar1\\(\\)` is currently only supported"
  )
  expect_error(
    gen_outcome(y ~ (1 + ar1() | ID), scale = shape ~ 1, params = params, burnin = 0, family = "negbin"),
    "`ar1\\(\\)` is currently only supported"
  )
  expect_error(
    gen_outcome(y ~ 1, scale = shape ~ ar1(), params = params, burnin = 0, family = "negbin"),
    "`ar1\\(\\)` is currently only supported"
  )
  expect_error(
    gen_outcome(y ~ 1, scale = sigma ~ 1, params = params, burnin = 0, family = "poisson"),
    "`scale` must not be supplied"
  )
  expect_error(
    gen_outcome(y ~ 1, params = params, burnin = 0, family = "negbin"),
    "Use `scale = shape ~ 1`"
  )
  expect_error(
    gen_outcome(y ~ 1, scale = sigma ~ 1, params = params, burnin = 0, family = "negbin"),
    "must have `shape` on the left-hand side"
  )
  expect_error(
    gen_outcome(y ~ 1, params = params, burnin = 0, family = "binomial"),
    "`trials` is required"
  )
  expect_error(
    gen_outcome(y ~ 1, scale = shape ~ 1, params = params, burnin = 0, family = "gamma", trials = 5),
    "must only be supplied"
  )
  expect_error(
    simulate_data(
      n = 5,
      seed = 408,
      generators = list(
        y = gen_outcome(y ~ 1, params = params, burnin = 0, family = "binomial", trials = 0)
      )
    ),
    "at least 1"
  )

  params_with_scale <- c(
    params,
    list(scale = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))))
  )
  expect_error(
    simulate_data(
      n = 5,
      seed = 409,
      generators = list(
        y = gen_outcome(y ~ 1, params = params_with_scale, burnin = 0, family = "poisson")
      )
    ),
    "`params\\$scale` must not be supplied"
  )

  params_with_correlation <- list(
    location = params$location,
    scale = list(
      beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y")),
      correlation = matrix(1, dimnames = list("y", "y"))
    )
  )
  expect_error(
    simulate_data(
      n = 5,
      seed = 410,
      generators = list(
        y = gen_outcome(y ~ 1, scale = shape ~ 1, params = params_with_correlation, burnin = 0, family = "negbin")
      )
    ),
    "only used for `family = \"gaussian\"`"
  )
})

test_that("gen_outcome guards against numerically extreme family parameters", {
  extreme_location <- list(
    location = list(beta = matrix(1000, nrow = 1, dimnames = list("(Intercept)", "y")))
  )
  expect_error(
    simulate_data(
      n = 5,
      seed = 411,
      generators = list(
        y = gen_outcome(y ~ 1, params = extreme_location, burnin = 0, family = "poisson")
      )
    ),
    "Non-finite poisson mean values"
  )

  extreme_scale <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(1000, nrow = 1, dimnames = list("(Intercept)", "y")))
  )
  expect_error(
    simulate_data(
      n = 5,
      seed = 412,
      generators = list(
        y = gen_outcome(y ~ 1, scale = shape ~ 1, params = extreme_scale, burnin = 0, family = "negbin")
      )
    ),
    "Non-finite negbin auxiliary scale parameter"
  )

  boundary_beta <- list(
    location = list(beta = matrix(20, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.01), nrow = 1, dimnames = list("(Intercept)", "y")))
  )
  expect_warning(
    simulate_data(
      n = 40,
      seed = 413,
      generators = list(
        y = gen_outcome(y ~ 1, scale = phi ~ 1, params = boundary_beta, burnin = 0, family = "beta")
      )
    ),
    "0/1 boundary"
  )
})

test_that("gen_template emits family-aware parameter templates", {
  poisson_sim <- simulate_data(
    n = 5,
    seed = 414,
    generators = list(
      template = gen_template(y ~ 1, burnin = 0, family = "poisson")
    )
  )
  poisson_md <- poisson_sim$generator_metadata$template
  expect_null(poisson_md$params$scale)
  expect_null(poisson_md$params$ar)
  expect_equal(
    poisson_md$params$location$beta,
    matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))
  )
  expect_identical(poisson_md$family, "poisson")
  expect_identical(poisson_md$parsed$family, "poisson")
  expect_null(poisson_md$parsed$dpar)

  binomial_sim <- simulate_data(
    n = 5,
    seed = 415,
    generators = list(
      template = gen_template(y ~ 1, burnin = 0, family = "binomial")
    )
  )
  expect_identical(binomial_sim$generator_metadata$template$parsed$trials_column, "y_trials")

  negbin_sim <- simulate_data(
    n = 5,
    seed = 416,
    generators = list(
      template = gen_template(y ~ 1, scale = shape ~ 1, burnin = 0, family = "negbin")
    )
  )
  negbin_params <- negbin_sim$generator_metadata$template$params
  expect_equal(
    negbin_params$scale$beta,
    matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))
  )
  expect_null(negbin_params$scale$correlation)

  beta_sim <- simulate_data(
    n_groups = 3,
    n_per_group = 2,
    group_id = "ID",
    seed = 417,
    generators = list(
      template = gen_template(
        y ~ 1 + (1 | ID),
        scale = phi ~ 1 + (1 | ID),
        burnin = 0,
        family = "beta"
      )
    )
  )
  beta_cov <- beta_sim$generator_metadata$template$params$random$ID$covariance
  expect_identical(
    rownames(beta_cov),
    c("location|outcome=y|term=(Intercept)", "phi|outcome=y|term=(Intercept)")
  )
})

test_that("gen_template params round-trip into gen_outcome for non-Gaussian families", {
  template_sim <- simulate_data(
    n = 30,
    seed = 418,
    generators = list(
      x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
      template = gen_template(y ~ x, scale = shape ~ 1, burnin = 0, family = "negbin")
    )
  )
  params <- template_sim$generator_metadata$template$params
  params$location$beta["(Intercept)", "y"] <- 1
  params$location$beta["x", "y"] <- 0.4
  params$scale$beta["(Intercept)", "y"] <- log(3)

  sim <- simulate_data(
    n = 30,
    seed = 419,
    generators = list(
      x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
      y = gen_outcome(y ~ x, scale = shape ~ 1, params = params, burnin = 0, family = "negbin")
    )
  )
  expect_true(all(sim$data$y >= 0))
  expect_identical(sim$generator_metadata$y$family, "negbin")
})

test_that("prep_sim_analysis emits family-correct brms formulas", {
  intercept_beta <- function(value) {
    matrix(value, nrow = 1, dimnames = list("(Intercept)", "y"))
  }

  poisson_sim <- simulate_data(
    n = 20,
    seed = 420,
    generators = list(
      y = gen_outcome(y ~ 1, params = list(location = list(beta = intercept_beta(1))), burnin = 0, family = "poisson")
    )
  )
  poisson_analysis <- prep_sim_analysis(poisson_sim)
  expect_identical(poisson_analysis$formula$family$family, "poisson")
  expect_length(poisson_analysis$formula$pforms, 0L)
  expect_identical(poisson_analysis$metadata$family, "poisson")

  binomial_sim <- simulate_data(
    n = 20,
    seed = 421,
    generators = list(
      y = gen_outcome(
        y ~ 1,
        params = list(location = list(beta = intercept_beta(0))),
        burnin = 0,
        family = "binomial",
        trials = 8
      )
    )
  )
  binomial_analysis <- prep_sim_analysis(binomial_sim)
  expect_identical(binomial_analysis$formula$family$family, "binomial")
  expect_match(
    deparse(binomial_analysis$formula$formula[[2L]]),
    "y | trials(y_trials)",
    fixed = TRUE
  )
  expect_true("y_trials" %in% names(binomial_analysis$data))

  dpar_params <- list(
    location = list(beta = intercept_beta(1)),
    scale = list(beta = intercept_beta(log(2)))
  )
  negbin_sim <- simulate_data(
    n = 20,
    seed = 422,
    generators = list(
      y = gen_outcome(y ~ 1, scale = shape ~ 1, params = dpar_params, burnin = 0, family = "negbin")
    )
  )
  negbin_analysis <- prep_sim_analysis(negbin_sim)
  expect_identical(negbin_analysis$formula$family$family, "negbinomial")
  expect_identical(deparse(negbin_analysis$formula$pforms$shape), "shape ~ 1")

  gamma_sim <- simulate_data(
    n = 20,
    seed = 423,
    generators = list(
      y = gen_outcome(y ~ 1, scale = shape ~ 1, params = dpar_params, burnin = 0, family = "gamma")
    )
  )
  gamma_analysis <- prep_sim_analysis(gamma_sim)
  expect_identical(gamma_analysis$formula$family$family, "gamma")
  expect_identical(gamma_analysis$formula$family$link, "log")

  beta_params <- list(
    location = list(beta = intercept_beta(0)),
    scale = list(beta = intercept_beta(log(10)))
  )
  beta_sim <- simulate_data(
    n = 20,
    seed = 424,
    generators = list(
      y = gen_outcome(y ~ 1, scale = phi ~ 1, params = beta_params, burnin = 0, family = "beta")
    )
  )
  beta_analysis <- prep_sim_analysis(beta_sim)
  expect_identical(beta_analysis$formula$family$family, "beta")
  expect_identical(deparse(beta_analysis$formula$pforms$phi), "phi ~ 1")
})

test_that("prep_sim_analysis links random effects shared across mean and scale formulas", {
  names_linked <- c(
    "location|outcome=y|term=(Intercept)",
    "phi|outcome=y|term=(Intercept)"
  )
  cov_linked <- diag(c(0.2, 0.1))
  dimnames(cov_linked) <- list(names_linked, names_linked)
  linked_params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(15), nrow = 1, dimnames = list("(Intercept)", "y"))),
    random = list(ID = list(covariance = cov_linked))
  )
  linked_sim <- simulate_data(
    n_groups = 5,
    n_per_group = 4,
    group_id = "ID",
    seed = 425,
    generators = list(
      y = gen_outcome(
        y ~ 1 + (1 | ID),
        scale = phi ~ 1 + (1 | ID),
        params = linked_params,
        burnin = 0,
        family = "beta"
      )
    )
  )
  linked_analysis <- prep_sim_analysis(linked_sim)
  expect_identical(deparse(linked_analysis$formula$formula), "y ~ 1 + (1 | p1 | ID)")
  expect_identical(deparse(linked_analysis$formula$pforms$phi), "phi ~ 1 + (1 | p1 | ID)")
  expect_true(linked_analysis$metadata$link_random)

  unlinked_analysis <- prep_sim_analysis(linked_sim, link_random = FALSE)
  expect_identical(deparse(unlinked_analysis$formula$formula), "y ~ 1 + (1 | ID)")
  expect_identical(deparse(unlinked_analysis$formula$pforms$phi), "phi ~ 1 + (1 | ID)")
  expect_false(unlinked_analysis$metadata$link_random)

  name_single <- "location|outcome=y|term=(Intercept)"
  cov_single <- matrix(0.2, dimnames = list(name_single, name_single))
  single_params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.5), nrow = 1, dimnames = list("(Intercept)", "y"))),
    random = list(ID = list(covariance = cov_single))
  )
  single_sim <- simulate_data(
    n_groups = 5,
    n_per_group = 4,
    group_id = "ID",
    seed = 426,
    generators = list(
      y = gen_outcome(
        y ~ 1 + (1 | ID),
        scale = sigma ~ 1,
        params = single_params,
        burnin = 0
      )
    )
  )
  single_analysis <- prep_sim_analysis(single_sim)
  expect_identical(deparse(single_analysis$formula$formula), "y ~ 1 + (1 | ID)")
  expect_identical(deparse(single_analysis$formula$pforms$sigma), "sigma ~ 1")
})

test_that("prep_sim_analysis falls back to gaussian for legacy metadata", {
  params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y")))
  )
  sim <- simulate_data(
    n = 10,
    seed = 427,
    generators = list(
      y = gen_outcome(y ~ 1, scale = sigma ~ 1, params = params, burnin = 0)
    )
  )
  sim$generator_metadata$y$parsed$family <- NULL
  sim$generator_metadata$y$parsed$dpar <- NULL

  analysis <- prep_sim_analysis(sim)
  expect_identical(analysis$formula$family$family, "gaussian")
  expect_identical(deparse(analysis$formula$pforms$sigma), "sigma ~ 1")
  expect_identical(analysis$metadata$family, "gaussian")
})

test_that("burnin is only required when ar1() appears in the location formula", {
  params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.5), nrow = 1, dimnames = list("(Intercept)", "y")))
  )
  expect_error(
    gen_outcome(y ~ ar1(), scale = sigma ~ 1, params = params),
    "`burnin` is required when `ar1\\(\\)` appears"
  )
  expect_error(
    gen_template(y ~ ar1(), scale = sigma ~ 1),
    "`burnin` is required when `ar1\\(\\)` appears"
  )

  sim <- simulate_data(
    n = 6,
    seed = 511,
    generators = list(
      y = gen_outcome(y ~ 1, scale = sigma ~ 1, params = params)
    )
  )
  expect_true("y" %in% names(sim$data))
  expect_identical(sim$generator_metadata$y$burnin$burnin, 0L)

  poisson_params <- list(
    location = list(beta = matrix(log(2), nrow = 1, dimnames = list("(Intercept)", "y")))
  )
  poisson_sim <- simulate_data(
    n = 6,
    seed = 512,
    generators = list(
      y = gen_outcome(y ~ 1, params = poisson_params, family = "poisson")
    )
  )
  expect_true("y" %in% names(poisson_sim$data))

  template_sim <- simulate_data(
    n = 6,
    seed = 513,
    generators = list(
      tmpl = gen_template(y ~ 1, scale = sigma ~ 1)
    )
  )
  expect_named(
    template_sim$generator_metadata$tmpl$params,
    c("location", "scale")
  )
})

test_that("stability metadata records the time-varying spectral norm guarantee", {
  ar_params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.5), nrow = 1, dimnames = list("(Intercept)", "y"))),
    ar = list(beta = array(0.4, dim = c(1, 1, 1), dimnames = list("ar1()", "y", "y")))
  )
  sim <- simulate_data(
    n_groups = 3,
    n_per_group = 5,
    group_id = "ID",
    time_id = "day",
    seed = 514,
    generators = list(
      y = gen_outcome(y ~ ar1(), scale = sigma ~ 1, params = ar_params, burnin = 10)
    )
  )
  stability <- sim$generator_metadata$y$ar$stability
  expect_false(stability$ar_matrices_vary_within_group)
  expect_true(stability$time_varying_stability_guaranteed)
  expect_null(stability$max_spectral_norm_overall)
})

test_that("gen_outcome resolves between/within of compositional ILR predictors", {
  params <- list(
    location = list(beta = matrix(
      c(0, 0.5, -0.3, 0.2, 0.1), 5, 1,
      dimnames = list(
        c("(Intercept)", "between(ilr1)", "between(ilr2)", "within(ilr1)", "within(ilr2)"),
        "y"
      )
    )),
    scale = list(beta = matrix(log(0.3), 1, 1, dimnames = list("(Intercept)", "y"))),
    random = list(ID = list(covariance = matrix(
      0.1, 1, 1,
      dimnames = list(
        "location|outcome=y|term=(Intercept)",
        "location|outcome=y|term=(Intercept)"
      )
    )))
  )
  sim <- simulate_data(
    n_groups = 5,
    n_per_group = 6,
    group_id = "ID",
    seed = 71,
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
      ),
      y = gen_outcome(
        y ~ between(ilr1) + between(ilr2) + within(ilr1) + within(ilr2) + (1 | ID),
        scale = sigma ~ 1,
        params = params
      )
    )
  )
  metadata <- sim$generator_metadata$y
  X <- cbind(
    1,
    sim$data$ilr1_between, sim$data$ilr2_between,
    sim$data$ilr1_within, sim$data$ilr2_within
  )
  mu_manual <- as.vector(
    X %*% params$location$beta +
      metadata$group_level_effects[as.integer(sim$data$ID), 1L]
  )
  expect_equal(as.vector(metadata$mu), mu_manual)
  selected <- metadata$selected_column_roles
  expect_setequal(selected$variable, c("ilr1", "ilr2"))
  expect_identical(unique(selected$generator), "comp")
})

test_that("prep_sim_analysis maps ILR between/within terms to complr coordinates", {
  params <- list(
    location = list(beta = matrix(
      c(0, 0.5, 0.2, 0.3), 4, 1,
      dimnames = list(
        c("(Intercept)", "between(ilr1)", "within(ilr1)", "between(x)"),
        "y"
      )
    )),
    scale = list(beta = matrix(log(0.3), 1, 1, dimnames = list("(Intercept)", "y")))
  )
  sim <- simulate_data(
    n_groups = 6,
    n_per_group = 5,
    group_id = "ID",
    seed = 72,
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
      ),
      x = gen_mvn(
        "x",
        level = "multilevel",
        fixed_intercept = 0,
        random_cov = 0.3,
        residual_cov = 1
      ),
      y = gen_outcome(
        y ~ between(ilr1) + within(ilr1) + between(x),
        scale = sigma ~ 1,
        params = params
      )
    )
  )
  analysis <- prep_sim_analysis(sim)

  expect_s3_class(analysis$complr, "complr")
  expect_length(analysis$complr$output, 1L)
  expect_null(analysis$metadata$outcome_composition_index)
  expect_identical(analysis$metadata$composition_generators, "comp")
  expect_identical(
    analysis$metadata$special_term_map,
    c("between(ilr1)" = "bz1_1", "within(ilr1)" = "wz1_1")
  )
  expect_identical(
    deparse(analysis$formula$formula),
    "y ~ bz1_1 + wz1_1 + x_between"
  )
  expect_true(all(c("bz1_1", "wz1_1", "x_between") %in% names(analysis$data)))
  # ILR variables are not manifest-centered; scalar predictors still are
  expect_false("ilr1" %in% analysis$metadata$derived_roles$variable)
  expect_true("x" %in% analysis$metadata$derived_roles$variable)

  # complr between coordinate equals the ILR of the closed arithmetic-mean composition
  parts <- as.matrix(sim$data[, c("sleep", "active", "sedentary"), with = FALSE])
  person_mean <- rowsum(parts, group = as.integer(sim$data$ID)) /
    as.vector(table(as.integer(sim$data$ID)))
  closed <- person_mean / rowSums(person_mean) * 24
  bz_manual <- multilevelcoda:::.mlsim_ilr_inverse
  basis <- compositions::gsi.buildilrBase(t(sim$generator_metadata$comp$sbp))
  bz_expected <- compositions::ilr(compositions::acomp(closed), V = basis)
  expect_equal(
    unname(analysis$data$bz1_1),
    unname(as.numeric(bz_expected[, 1L])[as.integer(sim$data$ID)]),
    tolerance = 1e-8
  )
})

test_that("prep_sim_analysis supports compositional predictors and outcomes together", {
  correlation <- diag(2)
  dimnames(correlation) <- list(c("o1", "o2"), c("o1", "o2"))
  params <- list(
    location = list(beta = matrix(
      c(0, 0.4, 0, -0.2), 2, 2,
      dimnames = list(c("(Intercept)", "between(ilr1)"), c("o1", "o2"))
    )),
    scale = list(
      beta = matrix(log(c(0.4, 0.35)), 1, 2, dimnames = list("(Intercept)", c("o1", "o2"))),
      correlation = correlation
    )
  )
  sim <- simulate_data(
    n_groups = 6,
    n_per_group = 5,
    group_id = "ID",
    seed = 73,
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
      ),
      out = gen_outcome(
        mvbind(o1, o2) ~ between(ilr1),
        scale = sigma ~ 1,
        params = params,
        composition = list(parts = c("rest", "move", "sit"), total = 24)
      )
    )
  )
  analysis <- prep_sim_analysis(sim)

  expect_length(analysis$complr$output, 2L)
  expect_identical(analysis$metadata$outcome_composition_index, 2L)
  expect_identical(analysis$metadata$composition_generators, "comp")
  expect_identical(
    analysis$metadata$response_map,
    c(o1 = "z1_2", o2 = "z2_2")
  )
  expect_identical(
    analysis$metadata$special_term_map,
    c("between(ilr1)" = "bz1_1")
  )
  expect_true(all(c("z1_2", "z2_2", "bz1_1") %in% names(analysis$data)))
  formulas <- analysis$formula$forms
  expect_identical(deparse(formulas[[1L]]$formula), "z1_2 ~ bz1_1")
  expect_identical(deparse(formulas[[2L]]$formula), "z2_2 ~ bz1_1")
})
