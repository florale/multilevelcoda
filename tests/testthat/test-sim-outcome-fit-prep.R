library(testthat)
library(multilevelcoda)

test_that("brmcoda helper columns remap scale helpers", {
  helper_map <- data.frame(
    column = c("lag_y", "lag_x"),
    brmcoda_column = c("lag_z1_1", "lag_x"),
    brmcoda_source = c("z1_1", "x"),
    response_derived = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  expect_identical(
    multilevelcoda:::.mlsim_brmcoda_helper_columns(c("lag_y", "lag_x"), helper_map),
    c("lag_z1_1", "lag_x")
  )
  expect_identical(
    multilevelcoda:::.mlsim_brmcoda_helper_columns(character(), helper_map),
    character()
  )
  expect_error(
    multilevelcoda:::.mlsim_brmcoda_helper_columns(c("lag_missing"), helper_map),
    "unknown helper columns"
  )
})

test_that("brmcoda prep returns scale helper columns in data namespace", {
  template <- simulate_data(
    n_groups = 2,
    n_per_group = 4,
    seed = 122,
    generators = list(
      x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
      y_template = outcome_template(
        cbind(y1, y2) ~ x + (1 | group_id),
        compositional = TRUE,
        parts = c("p1", "p2", "p3"),
        scale_formula = sigma ~ within(x)
      )
    )
  )$generator_metadata$y_template

  sim <- simulate_data(
    n_groups = 2,
    n_per_group = 4,
    seed = 123,
    generators = list(
      x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
      y = gen_outcome(
        cbind(y1, y2) ~ x + (1 | group_id),
        compositional = TRUE,
        parts = c("p1", "p2", "p3"),
        scale_formula = sigma ~ within(x),
        random_cov = template$random_cov
      )
    )
  )

  prep <- prepare_outcome_fit(sim, target = "brmcoda")

  expect_type(prep$scale_helper_columns, "character")
  expect_identical(prep$scale_helper_columns, "x_within")
  expect_true(all(prep$scale_helper_columns %in% names(prep$data)))
  expect_true(all(prep$scale_helper_columns %in% all.vars(prep$scale_formula)))
})

test_that("outcome random effects require explicit random_cov", {
  expect_error(
    simulate_data(
      n_groups = 2,
      n_per_group = 3,
      generators = list(
        x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
        y = gen_outcome(y ~ x + (1 | group_id))
      )
    ),
    "random_cov.*required"
  )

  expect_no_error(
    simulate_data(
      n_groups = 2,
      n_per_group = 3,
      generators = list(
        x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
        y_template = outcome_template(y ~ x + (1 | group_id))
      )
    )
  )
})

test_that("random_cov block lists are normalized to template order", {
  template <- simulate_data(
    n_groups = 3,
    n_per_group = 2,
    seed = 124,
    generators = list(
      x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
      y_template = outcome_template(cbind(y1, y2) ~ x + (1 | group_id))
    )
  )$generator_metadata$y_template

  reversed_random_cov <- rev(template$random_cov)
  sim <- simulate_data(
    n_groups = 3,
    n_per_group = 2,
    seed = 125,
    generators = list(
      x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
      y = gen_outcome(
        cbind(y1, y2) ~ x + (1 | group_id),
        coefficients = template$coefficients,
        residual_cov = template$residual_cov,
        random_cov = reversed_random_cov
      )
    )
  )

  expect_identical(
    names(sim$generator_metadata$y$random_cov),
    names(template$random_cov)
  )
})

test_that("random_cov block lists still require exact names", {
  template <- simulate_data(
    n_groups = 3,
    n_per_group = 2,
    seed = 126,
    generators = list(
      x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
      y_template = outcome_template(cbind(y1, y2) ~ x + (1 | group_id))
    )
  )$generator_metadata$y_template

  expect_error(
    simulate_data(
      n_groups = 3,
      n_per_group = 2,
      generators = list(
        x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
        y = gen_outcome(
          cbind(y1, y2) ~ x + (1 | group_id),
          random_cov = template$random_cov[-1L]
        )
      )
    ),
    "named list matching"
  )

  extra_random_cov <- template$random_cov
  extra_random_cov$extra <- template$random_cov[[1L]]
  expect_error(
    simulate_data(
      n_groups = 3,
      n_per_group = 2,
      generators = list(
        x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
        y = gen_outcome(
          cbind(y1, y2) ~ x + (1 | group_id),
          random_cov = extra_random_cov
        )
      )
    ),
    "named list matching"
  )

  duplicated_random_cov <- template$random_cov
  names(duplicated_random_cov)[[2L]] <- names(duplicated_random_cov)[[1L]]
  expect_error(
    simulate_data(
      n_groups = 3,
      n_per_group = 2,
      generators = list(
        x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
        y = gen_outcome(
          cbind(y1, y2) ~ x + (1 | group_id),
          random_cov = duplicated_random_cov
        )
      )
    ),
    "named list matching"
  )

  unnamed_random_cov <- unname(template$random_cov)
  expect_error(
    simulate_data(
      n_groups = 3,
      n_per_group = 2,
      generators = list(
        x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
        y = gen_outcome(
          cbind(y1, y2) ~ x + (1 | group_id),
          random_cov = unnamed_random_cov
        )
      )
    ),
    "named list matching"
  )
})
