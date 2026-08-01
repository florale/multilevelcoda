library(testthat)
library(multilevelcoda)

dynamic_recovery_params <- function() {
  outcomes <- c("ilr1", "ilr2")
  location <- matrix(
    c(0.1, -0.2),
    nrow = 1,
    dimnames = list("(Intercept)", outcomes)
  )
  scale <- matrix(
    log(c(0.4, 0.35)),
    nrow = 1,
    dimnames = list("(Intercept)", outcomes)
  )
  correlation <- matrix(c(1, 0.3, 0.3, 1), 2, dimnames = list(outcomes, outcomes))
  ar <- array(
    0,
    dim = c(1, 2, 2),
    dimnames = list("ar1()", outcomes, outcomes)
  )
  ar["ar1()", , ] <- matrix(c(0.25, 0.02, -0.01, 0.20), 2, 2, byrow = TRUE)

  random_names <- c(
    sprintf("location|outcome=%s|term=(Intercept)", outcomes),
    sprintf(
      "ar|term=ar1()|to=%s|from=%s",
      rep(outcomes, each = length(outcomes)),
      rep(outcomes, times = length(outcomes))
    ),
    sprintf("scale|outcome=%s|term=(Intercept)", outcomes)
  )
  covariance <- diag(c(0.20, 0.18, 0.12, 0.06, 0.06, 0.12, 0.16, 0.14)^2)
  dimnames(covariance) <- list(random_names, random_names)
  covariance[1L, 3L] <- covariance[3L, 1L] <- 0.20 * 0.20 * 0.12

  list(
    location = list(beta = location),
    scale = list(beta = scale, correlation = correlation),
    ar = list(beta = ar),
    random = list(ID = list(covariance = covariance))
  )
}

dynamic_recovery_sim <- function() {
  simulate_data(
    n_groups = 4,
    n_per_group = 6,
    group_id = "ID",
    time_id = "day",
    seed = 401,
    generators = list(
      outcome = gen_outcome(
        mvbind(ilr1, ilr2) ~ ar1() + (1 + ar1() | ID),
        scale = sigma ~ 1 + (1 | ID),
        params = dynamic_recovery_params(),
        burnin = 4,
        composition = list(parts = c("sleep", "sedentary", "activity"), total = 24)
      )
    )
  )
}

test_that("truth table aligns dynamic multivariate parameters with brms names", {
  params <- dynamic_recovery_params()
  analysis <- prep_sim_analysis(dynamic_recovery_sim())
  truth <- as.data.frame(analysis$truth)

  fixed <- truth[truth$type == "fixed", ]
  expect_setequal(
    fixed$parameter,
    c(
      "z11_Intercept", "z21_Intercept",
      "sigma_z11_Intercept", "sigma_z21_Intercept",
      "z11_lag_z1_1_within", "z11_lag_z2_1_within",
      "z21_lag_z1_1_within", "z21_lag_z2_1_within"
    )
  )
  expect_identical(
    fixed$truth[fixed$parameter == "z11_Intercept"],
    params$location$beta["(Intercept)", "ilr1"]
  )
  expect_identical(
    fixed$truth[fixed$parameter == "sigma_z21_Intercept"],
    params$scale$beta["(Intercept)", "ilr2"]
  )
  expect_identical(
    fixed$truth[fixed$parameter == "z11_lag_z2_1_within"],
    params$ar$beta["ar1()", "ilr1", "ilr2"]
  )

  random_sd <- truth[truth$type == "random_sd", ]
  expect_identical(unique(random_sd$group), "ID")
  expect_identical(
    random_sd$parameter,
    c(
      "z11_Intercept", "z21_Intercept",
      "z11_lag_z1_1_within", "z11_lag_z2_1_within",
      "z21_lag_z1_1_within", "z21_lag_z2_1_within",
      "sigma_z11_Intercept", "sigma_z21_Intercept"
    )
  )
  expect_equal(
    random_sd$truth,
    sqrt(diag(params$random$ID$covariance)),
    ignore_attr = TRUE
  )
  expect_identical(random_sd$simulator_name, rownames(params$random$ID$covariance))

  random_cor <- truth[truth$type == "random_cor", ]
  expect_identical(nrow(random_cor), 28L)
  linked <- random_cor[random_cor$coef1 == "z11_Intercept" &
                         random_cor$coef2 == "z11_lag_z1_1_within", ]
  expect_equal(linked$truth, 0.20, tolerance = 1e-12)

  rescor <- truth[truth$type == "rescor", ]
  expect_identical(rescor$parameter, "rescor(z11,z21)")
  expect_identical(rescor$truth, 0.3)
})

test_that("lag_center = none renames AR truth", {
  sim <- dynamic_recovery_sim()
  analysis <- prep_sim_analysis(sim, lag_center = "none")
  default_analysis <- prep_sim_analysis(sim)
  truth <- as.data.frame(analysis$truth)

  fixed <- truth[truth$type == "fixed", ]
  expect_setequal(
    fixed$parameter,
    c(
      "z11_Intercept", "z21_Intercept",
      "sigma_z11_Intercept", "sigma_z21_Intercept",
      "z11_lag_z1_1", "z11_lag_z2_1",
      "z21_lag_z1_1", "z21_lag_z2_1"
    )
  )

  # simulator names are the stable cross-parametrization key
  expect_identical(analysis$truth$simulator_name, default_analysis$truth$simulator_name)
  expect_identical(analysis$truth$truth, default_analysis$truth$truth)
})

test_that("link_random = FALSE drops cross-block random correlations from truth", {
  sim <- dynamic_recovery_sim()
  linked <- prep_sim_analysis(sim)
  unlinked <- prep_sim_analysis(sim, link_random = FALSE)
  linked_truth <- as.data.frame(linked$truth)
  unlinked_truth <- as.data.frame(unlinked$truth)

  expect_identical(sum(linked_truth$type == "random_cor"), 28L)
  # 15 mu-block pairs (6 coefficients) plus 1 sigma-block pair (2 coefficients)
  expect_identical(sum(unlinked_truth$type == "random_cor"), 16L)
  cross <- unlinked_truth[unlinked_truth$type == "random_cor" &
                            grepl("^sigma_", unlinked_truth$coef2) &
                            !grepl("^sigma_", unlinked_truth$coef1), ]
  expect_identical(nrow(cross), 0L)
  expect_identical(
    sum(unlinked_truth$type == "random_sd"),
    sum(linked_truth$type == "random_sd")
  )
})

test_that("univariate outcomes use unprefixed parameter names", {
  params <- list(
    location = list(beta = matrix(
      c(0.4, 0.3, 0.2),
      nrow = 3,
      dimnames = list(c("(Intercept)", "between(x)", "within(x)"), "y")
    )),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y"))),
    random = list(ID = list(covariance = matrix(
      0.25,
      dimnames = list(
        "location|outcome=y|term=(Intercept)",
        "location|outcome=y|term=(Intercept)"
      )
    )))
  )
  sim <- simulate_data(
    n_groups = 3,
    n_per_group = 4,
    group_id = "ID",
    seed = 402,
    generators = list(
      x = gen_mvn("x", level = "multilevel", fixed_intercept = 0, random_cov = 0.3, residual_cov = 0.7),
      outcome = gen_outcome(
        y ~ between(x) + within(x) + (1 | ID),
        scale = sigma ~ 1,
        params = params,
        burnin = 0
      )
    )
  )
  truth <- as.data.frame(prep_sim_analysis(sim)$truth)

  expect_setequal(
    truth$parameter[truth$type == "fixed"],
    c("Intercept", "x_between", "x_within", "sigma_Intercept")
  )
  random_sd <- truth[truth$type == "random_sd", ]
  expect_identical(random_sd$parameter, "Intercept")
  expect_identical(random_sd$group, "ID")
  expect_identical(random_sd$truth, 0.5)
  expect_false(any(truth$type == "rescor"))
})

test_that("moderated AR terms produce interaction names for every to/from pair", {
  outcomes <- c("y1", "y2")
  ar <- array(
    0.1,
    dim = c(2, 2, 2),
    dimnames = list(c("ar1()", "treatmenttreatment:ar1()"), outcomes, outcomes)
  )
  params <- list(
    location = list(beta = matrix(
      0,
      nrow = 2,
      ncol = 2,
      dimnames = list(c("(Intercept)", "treatmenttreatment"), outcomes)
    )),
    scale = list(
      beta = matrix(log(0.3), nrow = 1, ncol = 2, dimnames = list("(Intercept)", outcomes)),
      correlation = matrix(c(1, 0, 0, 1), 2, dimnames = list(outcomes, outcomes))
    ),
    ar = list(beta = ar)
  )
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 6,
    group_id = "ID",
    time_id = "day",
    seed = 403,
    generators = list(
      treatment = gen_categorical(
        "treatment",
        level = "level2",
        categories = c("control", "treatment"),
        fixed_intercept = stats::qlogis(0.5),
        output = "factor"
      ),
      outcome = gen_outcome(
        mvbind(y1, y2) ~ treatment + ar1() + treatment:ar1(),
        scale = sigma ~ 1,
        params = params,
        burnin = 4
      )
    )
  )
  truth <- as.data.frame(prep_sim_analysis(sim)$truth)

  expect_true(all(
    c(
      "y1_lag_y1_within", "y1_lag_y2_within",
      "y2_lag_y1_within", "y2_lag_y2_within",
      "y1_treatmenttreatment:lag_y1_within",
      "y1_treatmenttreatment:lag_y2_within",
      "y2_treatmenttreatment:lag_y1_within",
      "y2_treatmenttreatment:lag_y2_within"
    ) %in% truth$parameter[truth$type == "fixed"]
  ))
  expect_identical(
    truth$truth[truth$parameter == "y1_treatmenttreatment:lag_y2_within"],
    0.1
  )
})

test_that("non-gaussian families use their dpar names and never rescor", {
  poisson_sim <- simulate_data(
    n = 6,
    seed = 404,
    generators = list(
      outcome = gen_outcome(
        y ~ 1,
        family = "poisson",
        params = list(location = list(beta = matrix(
          log(3), nrow = 1, dimnames = list("(Intercept)", "y")
        )))
      )
    )
  )
  poisson_truth <- as.data.frame(prep_sim_analysis(poisson_sim)$truth)
  expect_identical(poisson_truth$parameter, "Intercept")
  expect_identical(poisson_truth$type, "fixed")

  beta_sim <- simulate_data(
    n = 6,
    seed = 405,
    generators = list(
      outcome = gen_outcome(
        y ~ 1,
        scale = phi ~ 1,
        family = "beta",
        params = list(
          location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
          scale = list(beta = matrix(log(10), nrow = 1, dimnames = list("(Intercept)", "y")))
        )
      )
    )
  )
  beta_truth <- as.data.frame(prep_sim_analysis(beta_sim)$truth)
  expect_setequal(beta_truth$parameter, c("Intercept", "phi_Intercept"))
  expect_false(any(beta_truth$type == "rescor"))

  negbin_sim <- simulate_data(
    n = 6,
    seed = 406,
    generators = list(
      outcome = gen_outcome(
        y ~ 1,
        scale = shape ~ 1,
        family = "negbin",
        params = list(
          location = list(beta = matrix(log(3), nrow = 1, dimnames = list("(Intercept)", "y"))),
          scale = list(beta = matrix(log(2), nrow = 1, dimnames = list("(Intercept)", "y")))
        )
      )
    )
  )
  negbin_truth <- as.data.frame(prep_sim_analysis(negbin_sim)$truth)
  expect_setequal(negbin_truth$parameter, c("Intercept", "shape_Intercept"))

  binomial_sim <- simulate_data(
    n = 6,
    seed = 407,
    generators = list(
      outcome = gen_outcome(
        y ~ 1,
        family = "binomial",
        trials = 8,
        params = list(location = list(beta = matrix(
          0, nrow = 1, dimnames = list("(Intercept)", "y")
        )))
      )
    )
  )
  binomial_truth <- as.data.frame(prep_sim_analysis(binomial_sim)$truth)
  expect_identical(binomial_truth$parameter, "Intercept")
  expect_identical(binomial_truth$type, "fixed")

  gamma_sim <- simulate_data(
    n = 6,
    seed = 408,
    generators = list(
      outcome = gen_outcome(
        y ~ 1,
        scale = shape ~ 1,
        family = "gamma",
        params = list(
          location = list(beta = matrix(log(3), nrow = 1, dimnames = list("(Intercept)", "y"))),
          scale = list(beta = matrix(log(2), nrow = 1, dimnames = list("(Intercept)", "y")))
        )
      )
    )
  )
  gamma_truth <- as.data.frame(prep_sim_analysis(gamma_sim)$truth)
  expect_setequal(gamma_truth$parameter, c("Intercept", "shape_Intercept"))
})

test_that("ambiguous sanitized response names abort instead of mislabeling", {
  # brms itself rejects duplicate sanitized responses at formula build; the
  # prefix guard is defense-in-depth against silent truth mislabeling
  expect_error(
    multilevelcoda:::.mlsim_recovery_prefixes(c(y_1 = "y_1", y1 = "y1")),
    "ambiguous"
  )
  expect_error(
    multilevelcoda:::.mlsim_recovery_prefixes(c(y_1 = "y_1", y1 = "y1")),
    "y_1"
  )
  expect_error(
    multilevelcoda:::.mlsim_recovery_prefixes(
      c(a = "z1_11", b = "z11_1", c = "z2_1")
    ),
    "z1_11"
  )
  # univariate models use an empty prefix
  expect_identical(
    multilevelcoda:::.mlsim_recovery_prefixes(c(y = "y")),
    c(y = "")
  )
})

test_that("special term maps translate ILR predictor terms", {
  map <- c("between(ilr1)" = "bz1_1", "within(ilr1)" = "wz1_1")
  expect_identical(
    multilevelcoda:::.mlsim_recovery_term("between(ilr1)", map),
    "bz1_1"
  )
  expect_identical(
    multilevelcoda:::.mlsim_recovery_term("treatmenttreatment:within(ilr1)", map),
    "treatmenttreatment:wz1_1"
  )
  expect_identical(
    multilevelcoda:::.mlsim_recovery_term("between(stress)", character()),
    "stress_between"
  )
  expect_identical(
    multilevelcoda:::.mlsim_recovery_term("(Intercept)", character()),
    "Intercept"
  )
})

test_that("structured random-effect names parse without regexes", {
  parsed <- multilevelcoda:::.mlsim_recovery_parse_random_name(
    "ar|term=ar1()|to=ilr1|from=ilr2",
    scale_label = "scale"
  )
  expect_identical(parsed$component, "ar")
  expect_identical(parsed$term, "ar1()")
  expect_identical(parsed$to, "ilr1")
  expect_identical(parsed$from, "ilr2")

  parsed <- multilevelcoda:::.mlsim_recovery_parse_random_name(
    "scale|outcome=y|term=(Intercept)",
    scale_label = "scale"
  )
  expect_identical(parsed$component, "scale")
  expect_identical(parsed$outcome, "y")
  expect_identical(parsed$term, "(Intercept)")

  expect_error(
    multilevelcoda:::.mlsim_recovery_parse_random_name("nope|foo=bar", scale_label = "scale"),
    "Cannot parse"
  )
})

fake_summary_matrix <- function(names, estimates) {
  matrix(
    c(estimates, rep(0.05, length(names)), estimates - 0.2, estimates + 0.2),
    nrow = length(names),
    dimnames = list(names, c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  )
}

test_that("estimate extraction and joining recover truth and fail loudly", {
  outcomes <- c("a", "b")
  random_names <- sprintf("location|outcome=%s|term=(Intercept)", outcomes)
  covariance <- matrix(c(0.04, 0.01, 0.01, 0.09), 2,
                       dimnames = list(random_names, random_names))
  params <- list(
    location = list(beta = matrix(
      c(0.5, -0.5), nrow = 1, ncol = 2, dimnames = list("(Intercept)", outcomes)
    )),
    scale = list(
      beta = matrix(log(0.3), nrow = 1, ncol = 2, dimnames = list("(Intercept)", outcomes)),
      correlation = matrix(c(1, 0.25, 0.25, 1), 2, dimnames = list(outcomes, outcomes))
    ),
    random = list(ID = list(covariance = covariance))
  )
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 4,
    group_id = "ID",
    seed = 408,
    generators = list(
      outcome = gen_outcome(
        mvbind(a, b) ~ (1 | ID),
        scale = sigma ~ 1,
        params = params
      )
    )
  )
  analysis <- prep_sim_analysis(sim)
  truth <- analysis$truth
  probs <- c(0.025, 0.975)

  fixef_summary <- fake_summary_matrix(
    c("a_Intercept", "b_Intercept", "sigma_a_Intercept", "sigma_b_Intercept"),
    c(0.55, -0.45, log(0.3), log(0.3))
  )
  cor_names <- c("a_Intercept", "b_Intercept")
  cor_array <- array(
    0,
    dim = c(2, 4, 2),
    dimnames = list(cor_names, c("Estimate", "Est.Error", "Q2.5", "Q97.5"), cor_names)
  )
  cor_array[2, , 1] <- cor_array[1, , 2] <- c(0.2, 0.1, -0.1, 0.5)
  varcorr <- list(ID = list(
    sd = fake_summary_matrix(cor_names, c(0.21, 0.31)),
    cor = cor_array
  ))
  rescor <- list(
    summary = fake_summary_matrix("rescor__a__b", 0.3),
    coef1 = "a",
    coef2 = "b"
  )
  variables <- c("cor_ID__a_Intercept__b_Intercept", "rescor__a__b")

  estimates <- multilevelcoda:::.mlsim_recovery_estimates(
    fixef_summary, varcorr, rescor, variables, probs
  )
  joined <- multilevelcoda:::.mlsim_recovery_join(truth, estimates)

  expect_identical(nrow(joined), nrow(truth))
  a_row <- joined[joined$parameter == "a_Intercept" & joined$type == "fixed", ]
  expect_equal(a_row$estimate, 0.55)
  expect_equal(a_row$bias, 0.05, tolerance = 1e-12)
  expect_true(a_row$covered)
  cor_row <- joined[joined$type == "random_cor", ]
  expect_equal(cor_row$estimate, 0.2)
  expect_equal(cor_row$truth, 0.01 / (0.2 * 0.3), tolerance = 1e-12)
  rescor_row <- joined[joined$type == "rescor", ]
  expect_equal(rescor_row$truth, 0.25)
  expect_equal(rescor_row$estimate, 0.3)

  # truth without a fitted estimate is a loud error
  expect_error(
    multilevelcoda:::.mlsim_recovery_join(
      truth,
      multilevelcoda:::.mlsim_recovery_estimates(
        fixef_summary[-1L, , drop = FALSE], varcorr, rescor, variables, probs
      )
    ),
    "not found in the fitted model"
  )

  # a fitted estimate without truth is a loud error
  extra_fixef <- rbind(fixef_summary, fake_summary_matrix("bogus_term", 0))
  expect_error(
    multilevelcoda:::.mlsim_recovery_join(
      truth,
      multilevelcoda:::.mlsim_recovery_estimates(extra_fixef, varcorr, rescor, variables, probs)
    ),
    "no simulation truth"
  )

  # zero-filled VarCorr correlations never silently match: without the
  # posterior variable the pair is not an estimate, so the truth row errors
  expect_error(
    multilevelcoda:::.mlsim_recovery_join(
      truth,
      multilevelcoda:::.mlsim_recovery_estimates(
        fixef_summary, varcorr, rescor, "rescor__a__b", probs
      )
    ),
    "not found in the fitted model"
  )
})

test_that("non-gaussian dpar truth joins fitted estimates end-to-end", {
  beta_sim <- simulate_data(
    n = 6,
    seed = 412,
    generators = list(
      outcome = gen_outcome(
        y ~ 1,
        scale = phi ~ 1,
        family = "beta",
        params = list(
          location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
          scale = list(beta = matrix(log(10), nrow = 1, dimnames = list("(Intercept)", "y")))
        )
      )
    )
  )
  analysis <- prep_sim_analysis(beta_sim)
  fixef_summary <- fake_summary_matrix(c("Intercept", "phi_Intercept"), c(0.1, log(10)))
  estimates <- multilevelcoda:::.mlsim_recovery_estimates(
    fixef_summary, NULL, NULL, character(), c(0.025, 0.975)
  )
  joined <- multilevelcoda:::.mlsim_recovery_join(analysis$truth, estimates)
  expect_setequal(joined$parameter, c("Intercept", "phi_Intercept"))
  expect_identical(
    joined$bias[joined$parameter == "phi_Intercept"],
    0
  )
})

test_that("lag_center = none truth joins fitted estimates end-to-end", {
  params <- list(
    location = list(beta = matrix(0.1, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y"))),
    ar = list(beta = array(0.3, dim = c(1, 1, 1), dimnames = list("ar1()", "y", "y")))
  )
  sim <- simulate_data(
    n = 6,
    time_id = "time",
    seed = 413,
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
  fixef_summary <- fake_summary_matrix(
    c("Intercept", "lag_y", "sigma_Intercept"),
    c(0.07, 0.3, log(0.2))
  )
  estimates <- multilevelcoda:::.mlsim_recovery_estimates(
    fixef_summary, NULL, NULL, character(), c(0.025, 0.975)
  )
  joined <- multilevelcoda:::.mlsim_recovery_join(analysis$truth, estimates)
  expect_identical(nrow(joined), nrow(analysis$truth))
  expect_setequal(joined$parameter, c("Intercept", "lag_y", "sigma_Intercept"))
  expect_identical(joined$bias[joined$parameter == "lag_y"], 0)
})

test_that("sim_recovery validates inputs before touching the fit", {
  expect_error(sim_recovery(NULL, list()), "mlsim_analysis")

  params <- list(
    location = list(beta = matrix(0, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y")))
  )
  sim <- simulate_data(
    n = 6,
    seed = 409,
    generators = list(
      outcome = gen_outcome(y ~ 1, scale = sigma ~ 1, params = params)
    )
  )
  analysis <- prep_sim_analysis(sim)
  expect_error(sim_recovery(NULL, analysis, probs = c(0.5)), "probs")
  expect_error(sim_recovery(NULL, analysis), "brmsfit")

  stripped <- analysis
  stripped$truth <- NULL
  expect_error(sim_recovery(NULL, stripped), "re-run `prep_sim_analysis\\(\\)`")
})

test_that("multiple random-effect terms disable truth with a recorded reason", {
  template <- simulate_data(
    n_groups = 4,
    n_per_group = 4,
    group_id = "ID",
    seed = 410,
    generators = list(
      x = gen_mvn("x", level = "multilevel", fixed_intercept = 0, random_cov = 0.3, residual_cov = 0.7),
      outcome_template = gen_template(
        y ~ within(x) + (1 | ID) + (1 + within(x) | ID),
        scale = sigma ~ 1
      )
    )
  )
  params <- template$generator_metadata$outcome_template$params
  sim <- simulate_data(
    n_groups = 4,
    n_per_group = 4,
    group_id = "ID",
    seed = 410,
    generators = list(
      x = gen_mvn("x", level = "multilevel", fixed_intercept = 0, random_cov = 0.3, residual_cov = 0.7),
      outcome = gen_outcome(
        y ~ within(x) + (1 | ID) + (1 + within(x) | ID),
        scale = sigma ~ 1,
        params = params,
        burnin = 0
      )
    )
  )
  analysis <- prep_sim_analysis(sim)
  expect_null(analysis$truth)
  expect_match(
    analysis$metadata$truth_unavailable_reason,
    "more than one random-effect term"
  )
  expect_error(sim_recovery(NULL, analysis), "more than one random-effect term")
  expect_output(print(analysis), "truth parameters: unavailable")
})

test_that("print methods surface truth availability and coverage", {
  params <- list(
    location = list(beta = matrix(0.5, nrow = 1, dimnames = list("(Intercept)", "y"))),
    scale = list(beta = matrix(log(0.2), nrow = 1, dimnames = list("(Intercept)", "y")))
  )
  sim <- simulate_data(
    n = 6,
    seed = 411,
    generators = list(
      outcome = gen_outcome(y ~ 1, scale = sigma ~ 1, params = params)
    )
  )
  analysis <- prep_sim_analysis(sim)
  expect_output(print(analysis), "truth parameters: 2")

  fixef_summary <- fake_summary_matrix(c("Intercept", "sigma_Intercept"), c(0.6, log(0.2)))
  estimates <- multilevelcoda:::.mlsim_recovery_estimates(
    fixef_summary, NULL, NULL, character(), c(0.025, 0.975)
  )
  joined <- multilevelcoda:::.mlsim_recovery_join(analysis$truth, estimates)
  data.table::setattr(joined, "probs", c(0.025, 0.975))
  data.table::setattr(joined, "class", c("mlsim_recovery", class(data.table::data.table())))
  expect_output(print(joined), "Coverage of the 95% interval: 2/2\\.")
})
