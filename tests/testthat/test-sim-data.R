library(testthat)
library(multilevelcoda)

test_that("simulate_data rejects conflicting design size arguments", {
  gens <- list(x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1))

  expect_error(
    simulate_data(n = 100, n_groups = 3, n_per_group = 2, generators = gens),
    "either `n` or `n_per_group`"
  )
  expect_error(
    simulate_data(n = 4, n_per_group = 99, generators = gens),
    "`n_per_group` requires `n_groups`"
  )

  # unambiguous designs still work
  expect_equal(nrow(simulate_data(n = 4, generators = gens, seed = 1)$data), 4L)
  expect_equal(
    nrow(simulate_data(n_groups = 3, n_per_group = 2, generators = gens, seed = 1)$data),
    6L
  )
  expect_equal(
    nrow(simulate_data(n = 7, n_groups = 3, generators = gens, seed = 1)$data),
    7L
  )
})

test_that(".mlsim_validate_sbp requires a nested (orthonormal) partition", {
  nested <- rbind(c(1, -1, -1), c(0, 1, -1))
  colnames(nested) <- c("a", "b", "c")
  expect_true(multilevelcoda:::.mlsim_validate_sbp(nested, 2L))

  non_nested <- rbind(c(1, -1, 0), c(0, 1, -1))
  colnames(non_nested) <- c("a", "b", "c")
  expect_error(
    multilevelcoda:::.mlsim_validate_sbp(non_nested, 2L),
    "not orthonormal"
  )

  # the same check must trigger through the composition resolver
  expect_error(
    multilevelcoda:::.mlsim_resolve_composition(sbp = non_nested, d = 2L),
    "not orthonormal"
  )
})
