library(testthat)
library(multilevelcoda)

centering_scalar_params <- function() {
  list(
    location = list(beta = matrix(
      c(0, 0.8, 0.3), ncol = 1,
      dimnames = list(c("(Intercept)", "between(x)", "within(x)"), "y")
    )),
    scale = list(beta = matrix(
      log(0.5), nrow = 1, dimnames = list("(Intercept)", "y")
    )),
    random = list(ID = list(covariance = matrix(
      0.4, 1, 1,
      dimnames = list(
        "location|outcome=y|term=(Intercept)",
        "location|outcome=y|term=(Intercept)"
      )
    )))
  )
}

centering_scalar_sim <- function(n_groups = 40, n_per_group = 6, seed = 2026) {
  simulate_data(
    n_groups = n_groups,
    n_per_group = n_per_group,
    group_id = "ID",
    seed = seed,
    generators = list(
      x = gen_mvn(
        "x",
        level = "multilevel",
        fixed_intercept = 0,
        random_cov = 1,
        residual_cov = 1
      ),
      y = gen_outcome(
        y ~ between(x) + within(x) + (1 | ID),
        scale = sigma ~ 1,
        params = centering_scalar_params()
      )
    )
  )
}

centering_composition_params <- function() {
  list(
    location = list(beta = matrix(
      c(0, 0.5, -0.3, 0.2, 0.1), 5, 1,
      dimnames = list(
        c(
          "(Intercept)", "between(ilr1)", "between(ilr2)",
          "within(ilr1)", "within(ilr2)"
        ),
        "y"
      )
    )),
    scale = list(beta = matrix(
      log(0.3), 1, 1, dimnames = list("(Intercept)", "y")
    )),
    random = list(ID = list(covariance = matrix(
      0.1, 1, 1,
      dimnames = list(
        "location|outcome=y|term=(Intercept)",
        "location|outcome=y|term=(Intercept)"
      )
    )))
  )
}

centering_composition_sim <- function(n_groups = 30, n_per_group = 6, seed = 71) {
  simulate_data(
    n_groups = n_groups,
    n_per_group = n_per_group,
    group_id = "ID",
    seed = seed,
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
        params = centering_composition_params()
      )
    )
  )
}

test_that("prep_sim_analysis defaults to manifest centering", {
  sim <- centering_scalar_sim()
  default <- prep_sim_analysis(sim)
  explicit <- prep_sim_analysis(sim, centering = "manifest")

  expect_identical(default$metadata$centering, "manifest")
  expect_equal(default$data, explicit$data)
  expect_equal(default$truth, explicit$truth)
  # manifest between is the realised person mean of the observed predictor
  expect_equal(
    default$data$x_between,
    as.numeric(ave(sim$data$x, sim$data$ID)),
    tolerance = 1e-12
  )
})

test_that("latent centering takes scalar between/within from the generator", {
  sim <- centering_scalar_sim()
  manifest <- prep_sim_analysis(sim)
  latent <- prep_sim_analysis(sim, centering = "latent")

  expect_identical(latent$metadata$centering, "latent")
  expect_equal(latent$data$x_between, sim$data$x_between, tolerance = 1e-12)
  expect_equal(
    latent$data$x_within,
    latent$data$x - latent$data$x_between,
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(latent$data$x_between, manifest$data$x_between)))
  expect_false(isTRUE(all.equal(latent$data$x_within, manifest$data$x_within)))
})

test_that("latent centering resolves level2 predictors from the role column", {
  params <- list(
    location = list(beta = matrix(
      c(0, 0.5), ncol = 1,
      dimnames = list(c("(Intercept)", "between(g)"), "y")
    )),
    scale = list(beta = matrix(
      log(0.5), nrow = 1, dimnames = list("(Intercept)", "y")
    ))
  )
  sim <- simulate_data(
    n_groups = 20,
    n_per_group = 4,
    group_id = "ID",
    seed = 11,
    generators = list(
      g = gen_mvn("g", level = "level2", fixed_intercept = 2, residual_cov = 1),
      y = gen_outcome(y ~ between(g), scale = sigma ~ 1, params = params)
    )
  )
  latent <- prep_sim_analysis(sim, centering = "latent")
  # a level2 variable has no `g_between` column; the role table points at `g`
  expect_false("g_between" %in% names(sim$data))
  expect_equal(latent$data$g_between, sim$data$g, tolerance = 1e-12)
})

test_that("latent centering rebuilds the complr between and within blocks", {
  sim <- centering_composition_sim()
  manifest <- prep_sim_analysis(sim)
  latent <- prep_sim_analysis(sim, centering = "latent")
  out <- latent$complr$output[[1L]]

  # coordinates equal the latent ILR means the simulator recorded
  expect_equal(latent$data$bz1_1, sim$data$ilr1_between,
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(latent$data$bz2_1, sim$data$ilr2_between,
               tolerance = 1e-10, ignore_attr = TRUE)
  # parts equal the latent between composition
  expect_equal(latent$data$bsleep, sim$data$sleep_between,
               tolerance = 1e-10, ignore_attr = TRUE)

  # internal algebra of the rebuilt blocks
  expect_equal(unclass(out$wZ), unclass(out$Z) - unclass(out$bZ),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(unclass(ilr(out$bX, V = out$psi)), unclass(out$bZ),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(unclass(out$wX), unclass(out$X - out$bX),
               tolerance = 1e-10, ignore_attr = TRUE)

  # names are preserved so brmcoda()/substitution() keep working
  expect_identical(names(manifest$data), names(latent$data))
  expect_identical(colnames(manifest$complr$output[[1L]]$bZ), colnames(out$bZ))
  expect_false(isTRUE(all.equal(manifest$data$bz1_1, latent$data$bz1_1)))
})

test_that("manifest compositional centering leaves a non-zero within-person mean", {
  sim <- centering_composition_sim()
  manifest <- prep_sim_analysis(sim)

  # arithmetic-mean closure is not the ILR barycentre, so the within
  # coordinates do not average to zero within person
  person_means <- tapply(manifest$data$wz1_1, manifest$data$ID, mean)
  expect_gt(max(abs(person_means)), 1e-6)
})

centering_ar_sim <- function(with_random, n_groups = 30, n_per_group = 8, seed = 5) {
  ar <- array(0.4, dim = c(1, 1, 1), dimnames = list("ar1()", "y", "y"))
  params <- list(
    location = list(beta = matrix(
      c(0, 0.8), ncol = 1,
      dimnames = list(c("(Intercept)", "between(x)"), "y")
    )),
    scale = list(beta = matrix(
      log(0.5), nrow = 1, dimnames = list("(Intercept)", "y")
    )),
    ar = list(beta = ar)
  )
  formula <- y ~ between(x) + ar1()
  if (isTRUE(with_random)) {
    params$random <- list(ID = list(covariance = matrix(
      0.4, 1, 1,
      dimnames = list(
        "location|outcome=y|term=(Intercept)",
        "location|outcome=y|term=(Intercept)"
      )
    )))
    formula <- y ~ between(x) + ar1() + (1 | ID)
  }
  simulate_data(
    n_groups = n_groups,
    n_per_group = n_per_group,
    group_id = "ID",
    time_id = "day",
    seed = seed,
    generators = list(
      x = gen_mvn("x", level = "multilevel", fixed_intercept = 0,
                  random_cov = 1, residual_cov = 1),
      y = gen_outcome(formula, scale = sigma ~ 1, params = params, burnin = 50)
    )
  )
}

test_that("manifest ar1 prepares under every lag_center", {
  for (with_random in c(FALSE, TRUE)) {
    sim <- centering_ar_sim(with_random)
    for (lc in c("within", "none")) {
      analysis <- prep_sim_analysis(sim, centering = "manifest", lag_center = lc)
      expect_s3_class(analysis, "mlsim_analysis")
      expect_false(is.null(analysis$truth), info = paste(with_random, lc))
      expect_identical(
        analysis$metadata$lag_columns,
        if (identical(lc, "within")) "lag_y_within" else "lag_y",
        info = paste(with_random, lc)
      )
    }
  }
})

test_that("latent centering builds the lag from the recorded residual state", {
  sim <- centering_ar_sim(FALSE)
  analysis <- prep_sim_analysis(sim, centering = "latent")

  expect_identical(analysis$metadata$lag_columns, "lag_y_latent")
  expect_identical(analysis$metadata$lag_center, "within")
  expect_true(analysis$metadata$lag_latent)

  # equals the true residual shifted one time step within group
  expected <- data.table::copy(sim$data)
  expected[, e := sim$generator_metadata$y$residual[, "y"]]
  data.table::setorderv(expected, c("ID", "day"))
  expected[, want := data.table::shift(e, 1L), by = "ID"]
  expect_equal(analysis$data$lag_y_latent, expected$want, tolerance = 1e-12)

  # the contemporaneous residual must not ship inside the analysis data,
  # because y - e == mu would be a fully determining oracle
  expect_false(any(grepl("_ar_residual", names(analysis$data), fixed = TRUE)))
})

test_that("the latent lag is not person-mean-centered", {
  sim <- centering_ar_sim(FALSE)
  latent <- prep_sim_analysis(sim, centering = "latent")
  manifest <- prep_sim_analysis(sim, centering = "manifest", lag_center = "within")

  latent_means <- tapply(latent$data$lag_y_latent, latent$data$ID,
                         mean, na.rm = TRUE)
  manifest_means <- tapply(manifest$data$lag_y_within, manifest$data$ID,
                           mean, na.rm = TRUE)
  # CMC drives the per-person means of the observed lag to ~0; the latent
  # residual is left uncentered, so its per-person means are not all ~0
  expect_gt(max(abs(latent_means)), 1e-6)
  expect_lt(max(abs(manifest_means)), max(abs(latent_means)))
})

test_that("lag_center and centering are independent under ar1", {
  # `lag_center` says whether the lag is centered, `centering` says with what,
  # so `lag_center = "none"` composes with latent between/within predictors
  sim <- centering_ar_sim(FALSE)
  expect_silent(analysis <- prep_sim_analysis(sim, centering = "latent",
                                              lag_center = "none"))
  expect_silent(prep_sim_analysis(sim, centering = "latent"))
  expect_silent(prep_sim_analysis(centering_scalar_sim(), centering = "latent",
                                  lag_center = "none"))

  expect_identical(analysis$metadata$lag_columns, "lag_y")
  expect_false(analysis$metadata$lag_latent)
  # the residual path is skipped entirely, so nothing is attached to drop
  expect_false(any(grepl("_ar_residual", names(analysis$data), fixed = TRUE)))

  # the lag ignores `centering`
  manifest_nc <- prep_sim_analysis(sim, centering = "manifest",
                                   lag_center = "none")
  expect_identical(analysis$data$lag_y, manifest_nc$data$lag_y)

  # the between/within predictors ignore `lag_center`
  latent_cmc <- prep_sim_analysis(sim, centering = "latent")
  expect_identical(analysis$data$x_between, latent_cmc$data$x_between)
  expect_identical(analysis$data$x_within, latent_cmc$data$x_within)
  expect_false(identical(analysis$data$x_between, manifest_nc$data$x_between))
})

test_that("print reports the latent lag source only when the lag is latent", {
  sim <- centering_ar_sim(FALSE)
  expect_output(
    print(prep_sim_analysis(sim, centering = "latent")),
    "lag centering: within \\(latent residual state\\)"
  )
  expect_output(
    print(prep_sim_analysis(sim, centering = "latent", lag_center = "none")),
    "lag centering: none\n"
  )
})

test_that("latent lag respects group boundaries and time gaps", {
  sim <- centering_ar_sim(FALSE)
  analysis <- prep_sim_analysis(sim, centering = "latent")
  first_rows <- analysis$data[, .I[which.min(day)], by = "ID"]$V1
  expect_true(all(is.na(analysis$data$lag_y_latent[first_rows])))

  gapped <- centering_ar_sim(FALSE)
  gapped$data <- gapped$data[!(ID == gapped$data$ID[[1L]] & day == 4)]
  gapped$generator_metadata$y$residual <-
    gapped$generator_metadata$y$residual[-4L, , drop = FALSE]
  gap_analysis <- suppressWarnings(
    prep_sim_analysis(gapped, centering = "latent")
  )
  expect_gt(gap_analysis$metadata$lag_gaps, 0L)
})

test_that("latent lag errors when the residual state is unusable", {
  sim <- centering_ar_sim(FALSE)
  sim$generator_metadata$y$residual <- NULL
  expect_error(
    prep_sim_analysis(sim, centering = "latent"),
    "requires the latent AR residual state"
  )

  short <- centering_ar_sim(FALSE)
  short$generator_metadata$y$residual <-
    short$generator_metadata$y$residual[-1L, , drop = FALSE]
  expect_error(
    prep_sim_analysis(short, centering = "latent"),
    "Refusing to build a misaligned latent lag"
  )
})

test_that("latent lag maps compositional responses through response_map", {
  outcomes <- c("ilr1", "ilr2")
  ar <- array(0, dim = c(1, 2, 2), dimnames = list("ar1()", outcomes, outcomes))
  ar["ar1()", , ] <- matrix(c(0.4, 0.05, 0.05, 0.3), 2, 2, byrow = TRUE)
  params <- list(
    location = list(beta = matrix(
      c(0.2, -0.1), nrow = 1, dimnames = list("(Intercept)", outcomes)
    )),
    scale = list(
      beta = matrix(
        log(c(0.3, 0.25)), nrow = 1, dimnames = list("(Intercept)", outcomes)
      ),
      correlation = matrix(
        c(1, 0.3, 0.3, 1), 2, dimnames = list(outcomes, outcomes)
      )
    ),
    ar = list(beta = ar)
  )
  sim <- simulate_data(
    n_groups = 20,
    n_per_group = 8,
    group_id = "ID",
    time_id = "day",
    seed = 91,
    generators = list(
      y = gen_outcome(
        mvbind(ilr1, ilr2) ~ ar1(),
        scale = sigma ~ 1,
        params = params,
        burnin = 50,
        composition = list(parts = c("sleep", "active", "sedentary"), total = 24)
      )
    )
  )
  analysis <- prep_sim_analysis(sim, centering = "latent")
  expect_identical(analysis$metadata$lag_columns,
                   c("lag_z1_1_latent", "lag_z2_1_latent"))

  # the residual is recorded in simulator ILR space; complr rebuilds the same
  # coordinates from the parts, so the mapping must be by name not position
  expected <- data.table::copy(sim$data)
  expected[, e1 := sim$generator_metadata$y$residual[, "ilr1"]]
  data.table::setorderv(expected, c("ID", "day"))
  expected[, want1 := data.table::shift(e1, 1L), by = "ID"]
  expect_equal(analysis$data$lag_z1_1_latent, expected$want1, tolerance = 1e-10)
})

test_that("latent centering errors when the latent components are absent", {
  sim <- centering_scalar_sim()
  # drop the latent columns to emulate a generator that never emitted them
  sim$data[, c("x_between", "x_within") := NULL]
  expect_error(
    prep_sim_analysis(sim, centering = "latent"),
    "requires the latent between component"
  )
})

test_that("latent centering errors for a non-multilevel composition", {
  sim <- centering_composition_sim()
  sim$data[, c("sleep_between", "active_between", "sedentary_between") := NULL]
  expect_error(
    prep_sim_analysis(sim, centering = "latent"),
    "requires the latent between composition"
  )
})

test_that("centering is inert without between/within terms", {
  params <- list(
    location = list(beta = matrix(
      0, nrow = 1, dimnames = list("(Intercept)", "y")
    )),
    scale = list(beta = matrix(
      log(0.2), nrow = 1, dimnames = list("(Intercept)", "y")
    ))
  )
  sim <- simulate_data(
    n = 10,
    seed = 1,
    generators = list(
      y = gen_outcome(y ~ 1, scale = sigma ~ 1, params = params)
    )
  )
  # no between()/within() terms, so centering is simply inert here
  expect_silent(prep_sim_analysis(sim, centering = "latent"))
})

test_that("centering is reported by print and stored in metadata", {
  sim <- centering_scalar_sim()
  latent <- prep_sim_analysis(sim, centering = "latent")
  expect_output(print(latent), "centering: latent")
  expect_output(print(prep_sim_analysis(sim)), "centering: manifest")
  expect_error(prep_sim_analysis(sim, centering = "geometric"))
})
