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
library(ggplot2)
library(plotly)
library(shiny)

parts <- c("TST", "WAKE", "MVPA", "LPA", "SB")

data(mcompd)
data(sbp)

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

make_par_data <- function(stat, d, use_substitution = FALSE) {
  n_levels <- if (use_substitution) {
    c(`3` = 6L, `4` = 12L, `5` = 20L)[[as.character(d)]]
  } else {
    c(`3` = 7L, `4` = 9L, `5` = 11L)[[as.character(d)]]
  }

  labels <- paste0(if (use_substitution) "S" else "E", seq_len(n_levels))
  out <- data.table(
    Stat = stat,
    JI = factor(rep(c("J1", "J2"), length.out = n_levels), levels = c("J1", "J2")),
    est = seq_len(n_levels) / 50,
    lower = seq_len(n_levels) / 50 - 0.05,
    upper = seq_len(n_levels) / 50 + 0.05
  )

  if (use_substitution) {
    out[, Substitution := factor(labels, levels = labels)]
  } else {
    out[, Estimand := factor(labels, levels = labels)]
  }

  out
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

fit_single <- fit_brmcoda(
  cilr_wide,
  Stress ~ z1_1 + z2_1 + z3_1 + z4_1 + Female + Age
)

test_that("complr covers wide-transform and validation branches", {
  plain_wide <- wide_data[1:5, .(TST, WAKE, MVPA, LPA, SB)]

  expect_message(
    default_ilr <- complr(data = plain_wide, parts = parts, total = 1440),
    "default sbp"
  )
  expect_true(is.complr(default_ilr))
  expect_null(default_ilr$output[[1]]$bX)
  expect_equal(colnames(default_ilr$output[[1]]$Z), paste0("z", 1:4, "_1"))

  wide_alr <- complr(data = plain_wide, parts = parts, total = 1440, transform = "alr")
  expect_equal(dim(wide_alr$output[[1]]$Z), c(5L, 4L))
  expect_null(wide_alr$output[[1]]$sbp)
  expect_null(wide_alr$output[[1]]$psi)

  wide_clr <- complr(data = plain_wide, parts = parts, total = 1440, transform = "clr")
  expect_equal(dim(wide_clr$output[[1]]$Z), c(5L, 5L))

  long_alr <- complr(
    data = ml_data,
    sbp = sbp,
    parts = parts,
    idvar = "ID",
    total = 1440,
    transform = "alr"
  )
  long_clr <- complr(
    data = ml_data,
    sbp = sbp,
    parts = parts,
    idvar = "ID",
    total = 1440,
    transform = "clr"
  )
  expect_equal(dim(long_alr$output[[1]]$Z), c(nrow(ml_data), 4L))
  expect_equal(dim(long_alr$output[[1]]$bZ), c(nrow(ml_data), 4L))
  expect_equal(dim(long_alr$output[[1]]$wZ), c(nrow(ml_data), 4L))
  expect_equal(dim(long_clr$output[[1]]$Z), c(nrow(ml_data), 5L))
  expect_equal(dim(long_clr$output[[1]]$bZ), c(nrow(ml_data), 5L))
  expect_equal(dim(long_clr$output[[1]]$wZ), c(nrow(ml_data), 5L))

  expect_error(
    complr(data = plain_wide, parts = parts, transform = c("ilr", "alr")),
    "only one type of transforms"
  )
  expect_error(
    complr(data = plain_wide, parts = parts, transform = "bad"),
    "'transform' should be one of"
  )
  expect_error(
    complr(data = plain_wide, parts = list()),
    "empty list"
  )
  expect_error(
    complr(data = plain_wide, parts = list(1:2)),
    "character vector"
  )
  expect_error(
    complr(
      data = ml_data,
      parts = list(c("TST", "WAKE"), c("MVPA", "LPA", "SB")),
      total = list(480),
      idvar = "ID"
    ),
    "same length"
  )
  expect_error(
    complr(
      data = ml_data,
      parts = list(c("TST", "WAKE"), c("MVPA", "LPA", "SB")),
      total = list(480, 960),
      sbp = list(as.matrix(data.table(1, -1))),
      idvar = "ID"
    ),
    "parts and sbp should have the same length"
  )

  na_data <- copy(plain_wide)
  na_data[1, TST := NA_real_]
  expect_error(
    complr(data = na_data, parts = parts),
    "missing data"
  )

  dup_name_data <- copy(plain_wide)
  dup_name_data[, z1_1 := 1]
  expect_error(
    complr(data = dup_name_data, sbp = sbp, parts = parts),
    "'data' cannot have any column names the same as logratio variables"
  )

  sbp_df <- structure(sbp, class = "not_matrix")
  expect_message(
    coerced <- complr(data = ml_data, sbp = sbp_df, parts = parts, idvar = "ID", total = 1440),
    "sbp is a"
  )
  expect_true(is.complr(coerced))
})

test_that("complr supports multi-composition objects and deprecated compilr path", {
  pair_sbp <- as.matrix(data.table(1, -1))
  tri_sbp <- build.sbp(c("MVPA", "LPA", "SB"))

  multi_comp <- complr(
    data = ml_data,
    parts = list(c("TST", "WAKE"), c("MVPA", "LPA", "SB")),
    sbp = list(pair_sbp, tri_sbp),
    total = list(480, 960),
    idvar = "ID"
  )

  expect_true(is.complr(multi_comp))
  expect_length(multi_comp$output, 2L)
  expect_true(all(c("z1_1", "z1_2", "z2_2") %in% names(multi_comp$dataout)))
  expect_equal(mean(multi_comp, parts = c("MVPA", "LPA", "SB"))$X |> length(), 3L)

  expect_warning(
    compile_result <- tryCatch(
      compilr(data = ml_data, parts = parts, sbp = sbp, idvar = "ID", total = 1440),
      error = identity
    ),
    "deprecated"
  )
  expect_s3_class(compile_result, "simpleError")
})

test_that("complr methods cover weighting and part selection paths", {
  mean_prop <- mean(cilr_ml, weight = "proportional")
  var_prop <- multilevelcoda:::var.complr(cilr_ml, weight = "proportional")
  df_ml <- as.data.frame(cilr_ml)
  mat_ml <- as.matrix(cilr_ml)

  expect_named(mean_prop, c("X", "bX", "wX", "Z", "bZ", "wZ"))
  expect_named(var_prop, c("X", "bX", "wX"))
  expect_s3_class(df_ml, "data.frame")
  expect_true(is.matrix(mat_ml))
  expect_equal(nrow(df_ml), nrow(cilr_ml$dataout))
})

test_that("build and get helpers cover exported and internal paths", {
  base_one <- build.base(parts)
  base_all <- build.base(parts, type = "one-to-all")
  sbp_default <- build.sbp(parts)
  vars_complr <- get_variables(cilr_ml)
  vars_ml <- get_variables(fit_ml)
  vars_single <- get_variables(fit_single)

  expect_s3_class(base_one, "data.table")
  expect_equal(dim(base_one), c(20L, 5L))
  expect_true(all(rowSums(as.matrix(base_one)) == 0))

  expect_s3_class(base_all, "data.table")
  expect_equal(dim(base_all), c(10L, 5L))
  expect_true(all(abs(rowSums(as.matrix(base_all))) < 1e-10))

  expect_equal(dim(sbp_default), c(4L, 5L))
  expect_equal(unname(sbp_default[1, ]), c(1, -1, -1, -1, -1))

  expect_named(vars_complr, "composition_1")
  expect_equal(sort(vars_ml$fixef_type), c("between", "within"))
  expect_identical(vars_ml$ranef_type, "multilevel")
  expect_identical(vars_single$fixef_type, "aggregate")
  expect_identical(vars_single$ranef_type, "single")

  expect_equal(unname(as.matrix(get_sbp(cilr_ml)[[1]])), unname(as.matrix(sbp)))
  expect_error(get_sbp(list()), "complr")

  expect_identical(multilevelcoda:::.get_parts(cilr_ml, 1), parts)
  expect_identical(multilevelcoda:::.get_parts(cilr_ml, parts), parts)
  expect_error(multilevelcoda:::.get_parts(cilr_ml, c(1, 2)), "single numeric value")
  expect_error(multilevelcoda:::.get_parts(cilr_ml, 99), "between 1 and")
  expect_error(multilevelcoda:::.get_parts(cilr_ml, 1L:2L), "single numeric value")
  expect_error(multilevelcoda:::.get_parts(cilr_ml, list("TST")), "character vector")
  expect_error(suppressWarnings(multilevelcoda:::.get_parts(cilr_ml, c("bad", "parts"))), "not found")
})

test_that("build.rg covers multilevel, single-level, and user-reference branches", {
  rg_between <- build.rg(
    object = fit_ml,
    ref = "grandmean",
    at = NULL,
    parts = parts,
    level = "between",
    weight = "equal"
  )
  rg_single <- build.rg(
    object = fit_single,
    ref = "grandmean",
    at = NULL,
    parts = parts,
    level = "aggregate",
    weight = "equal"
  )
  rg_between_num <- build.rg(
    object = fit_ml,
    ref = "grandmean",
    at = NULL,
    parts = 1,
    level = "between",
    weight = "equal"
  )
  rg_between_prop <- build.rg(
    object = fit_ml,
    ref = "grandmean",
    at = NULL,
    parts = parts,
    level = "between",
    weight = "proportional"
  )
  rg_agg_prop <- build.rg(
    object = fit_ml,
    ref = "grandmean",
    at = NULL,
    parts = parts,
    level = "aggregate",
    weight = "proportional"
  )
  user_ref <- data.table(TST = 500, WAKE = 200, MVPA = 100, LPA = 300, SB = 340, Female = 1, Age = mean(ml_data$Age))
  expect_true(all(c("bz1_1", "wz1_1", "bTST", "wTST", "ID", "Female", "Age", ".wgt.") %in% names(rg_between)))
  expect_true(all(c("z1_1", "tTST", "Female", "Age", ".wgt.") %in% names(rg_single)))
  expect_identical(names(rg_between_num), names(rg_between))
  expect_true(all(c("bz1_1", "wz1_1", "bTST", "wTST") %in% names(rg_between_prop)))
  expect_true(all(c("z1_1", "tTST", "ID") %in% names(rg_agg_prop)))

  expect_error(
    build.rg(
      object = fit_ml,
      ref = user_ref,
      at = NULL,
      parts = parts,
      level = "between",
      weight = "equal",
      fill = FALSE
    ),
    "not meaningful"
  )
  expect_error(
    build.rg(
      object = fit_ml,
      ref = data.table(TST = 500, WAKE = 200, MVPA = 100, LPA = 300, SB = 340, Female = 1),
      at = list(Female = c(0, 1)),
      parts = parts,
      level = "aggregate",
      weight = "equal",
      fill = TRUE
    ),
    "not meaningful"
  )
  expect_error(
    build.rg(
      object = fit_ml,
      ref = data.table(TST = 500, WAKE = 200, MVPA = 100, LPA = 300, Female = 1, Age = mean(ml_data$Age)),
      at = NULL,
      parts = parts,
      level = "between",
      weight = "equal",
      fill = FALSE
    ),
    "column not found"
  )
  expect_error(
    build.rg(
      object = fit_ml,
      ref = data.table(TST = 500, WAKE = 200, MVPA = 100, LPA = 300, SB = 340, Female = 1, Age = mean(ml_data$Age), extra = 1),
      at = NULL,
      parts = parts,
      level = "between",
      weight = "equal",
      fill = FALSE
    ),
    "not meaningful"
  )
})

test_that("internal utilities and plotting helpers exercise untested branches", {
  df_grid <- multilevelcoda:::expand.grid.df(data.frame(a = 1:2), data.frame(b = c("x", "y")))
  dt_grid <- multilevelcoda:::expand.grid.dt(data.frame(a = 1:2), data.frame(b = c("x", "y")))

  expect_s3_class(df_grid, "data.frame")
  expect_s3_class(dt_grid, "data.table")
  expect_equal(nrow(df_grid), 4L)
  expect_equal(nrow(dt_grid), 4L)
  expect_true(multilevelcoda:::is.sequential(1:4))
  expect_false(multilevelcoda:::is.sequential(c(1, 3, 4)))

  expect_s3_class(multilevelcoda:::.par_plot(make_par_data("bias", d = 3), d = 3), "ggplot")
  expect_s3_class(multilevelcoda:::.par_plot(make_par_data("bias", d = 3, use_substitution = TRUE), d = 3), "ggplot")
  expect_s3_class(multilevelcoda:::.par_plot(make_par_data("bias", d = 4), d = 4), "ggplot")
  expect_s3_class(multilevelcoda:::.par_plot(make_par_data("cover", d = 3, use_substitution = TRUE), d = 3), "ggplot")
  expect_s3_class(multilevelcoda:::.par_plot(make_par_data("mse", d = 4), d = 4), "ggplot")
  expect_s3_class(multilevelcoda:::.par_plot(make_par_data("becover", d = 5, use_substitution = TRUE), d = 5), "ggplot")
  expect_s3_class(multilevelcoda:::.par_plot(make_par_data("becover", d = 5), d = 5), "ggplot")
  expect_s3_class(multilevelcoda:::.par_plot(make_par_data("mse", d = 4, use_substitution = TRUE), d = 4), "ggplot")
  expect_s3_class(multilevelcoda:::.par_plot(make_par_data("empse", d = 4, use_substitution = TRUE), d = 4), "ggplot")
  expect_true(inherits(multilevelcoda:::.par_plot(make_par_data("empse", d = 4), shiny = TRUE, d = 4), "plotly"))
})

test_that("shiny entry points build successfully and server outputs are evaluable", {
  expect_error(launch_shinystan(structure(list(model = NULL), class = "brmcoda")))

  app <- suppressWarnings(multilevelcoda_sim())
  expect_s3_class(app, "shiny.appobj")

  server <- app$serverFuncSource()
  expect_true(is.function(server))

  testServer(server, {
    base_inputs <- list(
      par_brmcoda = "All",
      rint_sd_brmcoda = "medium",
      res_sd1_brmcoda = "medium",
      res_sd2_brmcoda = "large",
      res_sd3_brmcoda = "small",
      N_brmcoda = "All",
      K_brmcoda = "All",
      D_brmcoda = 4,
      stat_brmcoda = "bias",
      rint_sd_sub = "medium",
      res_sd1_sub = "medium",
      res_sd2_sub = "large",
      res_sd3_sub = "small",
      N_sub = "All",
      K_sub = "All",
      D_sub = 4,
      delta3_sub = "All",
      delta4_sub = "All",
      delta5_sub = "All",
      level_sub = "All",
      stat_sub = "bias",
      rint_sd_brmcoda_plot = "medium",
      res_sd1_brmcoda_plot = "medium",
      res_sd2_brmcoda_plot = "large",
      res_sd3_brmcoda_plot = "small",
      D_brmcoda_plot = 4,
      stat_brmcoda_plot = "bias",
      rint_sd_sub_plot = "medium",
      res_sd1_sub_plot = "medium",
      res_sd2_sub_plot = "large",
      res_sd3_sub_plot = "small",
      D_sub_plot = 4,
      stat_sub_plot = "bias"
    )

    do.call(session$setInputs, base_inputs)
    expect_true("json" %in% class(output$simsum_brmcoda_table))
    expect_true("json" %in% class(output$simsum_sub_table))
    expect_true("json" %in% class(output$simsum_brmcoda_plot))
    expect_true("json" %in% class(output$simsum_sub_plot))
    expect_equal(names(output$coda), c("src", "heigh", "width"))
    expect_equal(output$coda$width, "50px")

    for (par_choice in c(
      "b_Intercept",
      "b_bilr1",
      "b_bilr2",
      "b_bilr3",
      "b_bilr4",
      "b_wilr1",
      "b_wilr2",
      "b_wilr3",
      "b_wilr4",
      "sd_ID_Intercept",
      "sigma"
    )) {
      do.call(session$setInputs, modifyList(base_inputs, list(par_brmcoda = par_choice)))
      expect_true("json" %in% class(output$simsum_brmcoda_table))
    }

    do.call(session$setInputs, modifyList(base_inputs, list(
      rint_sd_brmcoda = "small",
      D_brmcoda = 3,
      stat_brmcoda = "cover",
      rint_sd_sub = "large",
      D_sub = 3,
      stat_sub = "cover",
      rint_sd_brmcoda_plot = "small",
      D_brmcoda_plot = 3,
      stat_brmcoda_plot = "cover",
      rint_sd_sub_plot = "large",
      D_sub_plot = 3,
      stat_sub_plot = "cover"
    )))
    expect_true("json" %in% class(output$simsum_brmcoda_table))
    expect_true("json" %in% class(output$simsum_sub_plot))

    do.call(session$setInputs, modifyList(base_inputs, list(
      par_brmcoda = "sigma",
      rint_sd_brmcoda = "medium",
      res_sd1_brmcoda = "small",
      N_brmcoda = "50",
      K_brmcoda = "3",
      D_brmcoda = 3,
      stat_brmcoda = "cover",
      rint_sd_sub = "medium",
      res_sd1_sub = "small",
      N_sub = "50",
      K_sub = "3",
      D_sub = 3,
      delta3_sub = "Sleep",
      level_sub = "Within",
      stat_sub = "cover",
      rint_sd_brmcoda_plot = "medium",
      res_sd1_brmcoda_plot = "small",
      D_brmcoda_plot = 3,
      stat_brmcoda_plot = "cover",
      rint_sd_sub_plot = "medium",
      res_sd1_sub_plot = "small",
      D_sub_plot = 3,
      stat_sub_plot = "cover"
    )))
    expect_true("json" %in% class(output$simsum_brmcoda_table))
    expect_true("json" %in% class(output$simsum_sub_table))
    expect_true("json" %in% class(output$simsum_brmcoda_plot))
    expect_true("json" %in% class(output$simsum_sub_plot))

    do.call(session$setInputs, modifyList(base_inputs, list(
      rint_sd_brmcoda = "large",
      D_brmcoda = 5,
      stat_brmcoda = "becover",
      rint_sd_sub = "medium",
      res_sd1_sub = "small",
      D_sub = 5,
      stat_sub = "becover",
      rint_sd_brmcoda_plot = "large",
      D_brmcoda_plot = 5,
      stat_brmcoda_plot = "becover",
      rint_sd_sub_plot = "medium",
      res_sd1_sub_plot = "small",
      D_sub_plot = 5,
      stat_sub_plot = "becover"
    )))
    expect_true("json" %in% class(output$simsum_brmcoda_plot))
    expect_true("json" %in% class(output$simsum_sub_table))

    do.call(session$setInputs, modifyList(base_inputs, list(
      par_brmcoda = "b_bilr1",
      rint_sd_brmcoda = "medium",
      res_sd1_brmcoda = "large",
      N_brmcoda = "30",
      K_brmcoda = "3",
      D_brmcoda = 5,
      stat_brmcoda = "bias",
      rint_sd_sub = "medium",
      res_sd1_sub = "large",
      N_sub = "30",
      K_sub = "3",
      D_sub = 5,
      delta5_sub = "TST",
      level_sub = "Between",
      stat_sub = "bias",
      rint_sd_brmcoda_plot = "medium",
      res_sd1_brmcoda_plot = "large",
      D_brmcoda_plot = 5,
      stat_brmcoda_plot = "bias",
      rint_sd_sub_plot = "medium",
      res_sd1_sub_plot = "large",
      D_sub_plot = 5,
      stat_sub_plot = "bias"
    )))
    expect_true("json" %in% class(output$simsum_brmcoda_table))
    expect_true("json" %in% class(output$simsum_sub_table))
    expect_true("json" %in% class(output$simsum_brmcoda_plot))
    expect_true("json" %in% class(output$simsum_sub_plot))

    do.call(session$setInputs, modifyList(base_inputs, list(
      D_sub = 4,
      delta4_sub = "Sleep",
      stat_sub = "empse"
    )))
    expect_true("json" %in% class(output$simsum_sub_table))
  })
})
