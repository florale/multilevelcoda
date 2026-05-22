#' Internal Generator Specification Helpers
#'
#' Create and normalize generator specification objects used by
#' [simulate_data()]. These helpers are internal implementation details for the
#' public `gen_*()` API.
#'
#' @param type Character scalar naming the generator type.
#' @param vars Character vector of generated variable names.
#' @param level Simulation level: `"single"`, `"level2"`, or `"multilevel"`.
#' @param simulate Function called by [simulate_data()] with a simulation
#'   context.
#' @param ... Additional metadata stored on the generator specification.
#'
#' @return Internal helper outputs used by the simulation engine.
#'
#' @examples
#' spec <- multilevelcoda:::.mlsim_spec(
#'   "custom", "x", "single",
#'   function(context) {
#'     multilevelcoda:::.mlsim_result(rep(0, context$n_rows), "x", list())
#'   }
#' )
#' inherits(spec, "mlsim_generator_spec")
#'
#' @keywords internal
#' @name multilevelcoda-internal-generator-specs
NULL

#' @rdname multilevelcoda-internal-generator-specs
.mlsim_spec <- function(type, vars, level, simulate, ...) {
  structure(
    c(
      list(
        type = type,
        vars = vars,
        level = .mlsim_match_level(level),
        simulate = simulate
      ),
      list(...)
    ),
    class = c(paste0("mlsim_", type, "_spec"), "mlsim_generator_spec")
  )
}

#' Generate Normal Variables
#'
#' Create a generator specification for univariate Gaussian variables.
#'
#' @param vars Character scalar naming the generated variable.
#' @param level Simulation level. `"single"` generates one row-level value per
#'   observation, `"level2"` generates one group-level value and expands it to
#'   each row in the group, and `"multilevel"` generates row-level values from a
#'   random-intercept model.
#' @param mean,sd Direct mean and standard deviation for `"single"` and
#'   `"level2"` generators. Values are recycled to the number of generated
#'   units.
#' @param fixed_intercept Location intercept for multilevel generators. For
#'   normal variables this is on the response scale.
#' @param random_cov Group-level random-intercept covariance. A scalar is
#'   accepted for univariate location-only generators; a matrix is required when
#'   location and scale random effects are modeled jointly.
#' @param residual_var Residual variance for multilevel normal generators when
#'   no scale model is used.
#' @param scale_fixed_intercept Optional log residual standard-deviation
#'   intercept for multilevel generators. When supplied, residual standard
#'   deviations are `exp(scale_fixed_intercept + group_random_effect)`.
#'
#' @return An `mlsim_generator_spec` object for use in [simulate_data()].
#'
#' @details
#' For non-multilevel normal variables, `mean` and `sd` define independent
#' draws directly. For multilevel variables, `fixed_intercept`, `random_cov`,
#' and either `residual_var` or `scale_fixed_intercept` define a
#' random-intercept data-generating model.
#'
#' @examples
#' sim <- simulate_data(
#'   n_groups = 2,
#'   n_per_group = 3,
#'   seed = 1,
#'   generators = list(
#'     x = gen_normal("x", mean = 0, sd = 1),
#'     u = gen_normal(
#'       "u",
#'       level = "multilevel",
#'       fixed_intercept = 0,
#'       random_cov = 0.2,
#'       residual_var = 1
#'     )
#'   )
#' )
#' sim$data
#'
#' @family predictor generators
#' @export
gen_normal <- function(vars, level = c("single", "level2", "multilevel"),
                       mean = 0, sd = NULL,
                       fixed_intercept = NULL, random_cov = NULL,
                       residual_var = NULL, scale_fixed_intercept = NULL) {
  mean_supplied <- !missing(mean)
  sd_supplied <- !missing(sd)
  vars <- .mlsim_check_vars(vars, 1L)
  level <- .mlsim_match_level(level)
  .mlsim_validate_scale_multilevel_args(level, scale_fixed_intercept, random_cov)
  .mlsim_check_multilevel_unused(
    level,
    c(mean = mean_supplied, sd = sd_supplied),
    "normal",
    "`fixed_intercept`, `random_cov`, `residual_var`, and `scale_fixed_intercept`"
  )
  .mlsim_check_multilevel_required(
    level,
    list(fixed_intercept = fixed_intercept),
    "normal"
  )
  if (identical(level, "multilevel")) {
    .mlsim_check_scale_or_arg(residual_var, "residual_var", scale_fixed_intercept, "normal")
  }
  if (!identical(level, "multilevel")) {
    if (!is.null(fixed_intercept)) {
      .mlsim_warn_ignored_args("normal", "mean", "fixed_intercept")
    }
    if (sd_supplied && !is.null(residual_var)) {
      .mlsim_warn_ignored_args("normal", "sd", "residual_var")
    }
  }

  simulate <- function(context) {
    if (identical(level, "multilevel")) {
      fixed <- fixed_intercept
      if (is.null(fixed)) {
        .mlsim_stop("Multilevel normal generator requires `fixed_intercept`.")
      }
      fixed <- .mlsim_check_scalar_number(fixed, "fixed_intercept")
      scale_fixed <- .mlsim_check_scale_fixed_intercept(scale_fixed_intercept, 1L)
      random <- .mlsim_joint_random_intercepts(vars, random_cov, context, scale = !is.null(scale_fixed))
      location_eta <- fixed + random$random_intercepts[context$group_index, 1L]
      if (is.null(scale_fixed)) {
        residual_var <- .mlsim_check_nonnegative(.mlsim_recycle(residual_var, 1L, "residual_var"), "residual_var")
        residual_sd <- sqrt(residual_var)
        residual <- stats::rnorm(context$n_rows, 0, residual_sd)
        row_parameters <- list(residual_var = residual_var, residual_sd = residual_sd)
        scale_eta <- NULL
      } else {
        scale_eta <- scale_fixed + random$scale_random_intercepts[context$group_index, 1L]
        residual_sd <- exp(scale_eta)
        residual <- stats::rnorm(context$n_rows, 0, residual_sd)
        row_parameters <- list(residual_sd = residual_sd, residual_var = residual_sd^2)
      }
      values <- location_eta + residual
      metadata <- c(
        list(
          distribution = "normal",
          level = level,
          vars = vars
        ),
        .mlsim_simple_random_metadata(
          fixed,
          random,
          scale_fixed_intercept = scale_fixed,
          scale_linear_predictor = scale_eta
        ),
        list(
          row_parameters = row_parameters,
          residual_var = row_parameters$residual_var,
          residual_sd = row_parameters$residual_sd
        )
      )
      return(.mlsim_result(values, vars, metadata))
    }

    n <- .mlsim_n_for_level(level, context)
    mean_n <- .mlsim_check_finite_numeric(.mlsim_recycle(mean %||% 0, n, "mean"), "mean")
    sd_n <- .mlsim_check_positive(.mlsim_recycle(sd %||% sqrt(residual_var %||% 1), n, "sd"), "sd")
    values <- stats::rnorm(n, mean = mean_n, sd = sd_n)
    values <- .mlsim_expand_level2_vector(values, level, context)
    .mlsim_result(
      values,
      vars,
      c(
        list(distribution = "normal", level = level, vars = vars),
        .mlsim_parameter_metadata(level, context, mean = mean_n, sd = sd_n)
      )
    )
  }

  .mlsim_spec(
    "normal", vars, level, simulate,
    mean = if (identical(level, "multilevel")) NULL else mean,
    sd = if (identical(level, "multilevel")) NULL else sd,
    fixed_intercept = fixed_intercept,
    random_cov = random_cov,
    residual_var = residual_var,
    scale_fixed_intercept = scale_fixed_intercept
  )
}

#' Generate Multivariate Normal and Compositional Variables
#'
#' Create a generator specification for multivariate Gaussian variables,
#' optionally transformed from ILR coordinates to compositional parts.
#'
#' @param vars Character vector naming the generated variables. For
#'   compositional generators these are ILR coordinate names.
#' @param level Simulation level. `"single"` generates one row-level vector per
#'   observation, `"level2"` generates one group-level vector and expands it to
#'   each row in the group, and `"multilevel"` uses a random-intercept model.
#' @param mean,cov Direct mean vector and covariance matrix for `"single"` and
#'   `"level2"` generators. Defaults are zero means and identity covariance.
#' @param fixed_intercept Location intercept vector for multilevel generators.
#' @param random_cov Group-level random-intercept covariance matrix. When
#'   `scale_fixed_intercept` is supplied this covariance spans location and
#'   scale intercepts jointly.
#' @param residual_cov Residual covariance matrix for multilevel MVN generators
#'   without a row-specific scale model.
#' @param scale_fixed_intercept Optional log residual standard-deviation
#'   intercepts for multilevel generators.
#' @param residual_cor Residual correlation matrix used with
#'   `scale_fixed_intercept`.
#' @param compositional Logical; if `TRUE`, treat `vars` as ILR coordinates and
#'   back-transform to composition parts.
#' @param parts Character vector naming the composition parts. Must have
#'   `length(vars) + 1` entries unless supplied through `sbp` column names.
#' @param total Positive scalar total for closed compositions.
#' @param keep_ilr Logical; if `TRUE`, emit both ILR coordinates and parts. If
#'   `FALSE`, emit only parts.
#' @param sbp Optional sequential binary partition matrix. Columns name parts;
#'   rows define ILR balances using `-1`, `0`, and `1`.
#'
#' @return An `mlsim_generator_spec` object for use in [simulate_data()].
#'
#' @details
#' Compositional MVN generators simulate ILR coordinates first, then use the
#' SBP basis to transform the coordinates into positive composition parts that
#' sum to `total`.
#'
#' @examples
#' sim <- simulate_data(
#'   n = 4,
#'   seed = 2,
#'   generators = list(
#'     z = gen_mvn(
#'       c("z1", "z2"),
#'       mean = c(0, 0),
#'       cov = diag(2),
#'       compositional = TRUE,
#'       parts = c("sleep", "activity", "sedentary")
#'     )
#'   )
#' )
#' sim$data
#' rowSums(as.matrix(sim$data[, c("sleep", "activity", "sedentary"), with = FALSE]))
#'
#' @family predictor generators
#' @export
gen_mvn <- function(vars, level = c("single", "level2", "multilevel"),
                    mean = NULL, cov = NULL,
                    fixed_intercept = NULL, random_cov = NULL, residual_cov = NULL,
                    scale_fixed_intercept = NULL, residual_cor = NULL,
                    compositional = FALSE, parts = NULL, total = 1,
                    keep_ilr = TRUE, sbp = NULL) {
  mean_supplied <- !missing(mean)
  cov_supplied <- !missing(cov)
  vars <- .mlsim_check_vars(vars)
  compositional <- .mlsim_check_scalar_logical(compositional, "compositional")
  keep_ilr <- .mlsim_check_scalar_logical(keep_ilr, "keep_ilr")
  level <- .mlsim_match_level(level)
  d <- length(vars)
  .mlsim_validate_scale_multilevel_args(level, scale_fixed_intercept, random_cov, residual_cor)
  .mlsim_check_multilevel_unused(
    level,
    c(mean = mean_supplied, cov = cov_supplied),
    "MVN",
    "`fixed_intercept`, `random_cov`, `residual_cov`, `scale_fixed_intercept`, and `residual_cor`"
  )
  .mlsim_check_multilevel_required(
    level,
    list(fixed_intercept = fixed_intercept),
    "MVN"
  )
  if (identical(level, "multilevel")) {
    .mlsim_check_scale_or_arg(residual_cov, "residual_cov", scale_fixed_intercept, "MVN")
    if (is.null(scale_fixed_intercept) && !is.null(residual_cor)) {
      .mlsim_stop("`residual_cor` requires `scale_fixed_intercept`.")
    }
  }

  if (compositional) {
    composition <- .mlsim_resolve_composition(parts, sbp, d)
    parts <- composition$parts
    sbp <- composition$sbp
    .mlsim_check_composition_output_names(vars, parts, keep_ilr)
  } else if (!is.null(parts) || !is.null(sbp)) {
    .mlsim_stop("`parts` and `sbp` require `compositional = TRUE`.")
  }
  output_vars <- if (compositional) {
    if (keep_ilr) c(vars, parts) else parts
  } else {
    vars
  }
  if (!identical(level, "multilevel") && !is.null(fixed_intercept)) {
    .mlsim_warn_ignored_args("MVN", c("mean", "cov"), "fixed_intercept")
  }

  simulate <- function(context) {
    if (identical(level, "multilevel")) {
      fixed <- fixed_intercept
      if (is.null(fixed)) {
        .mlsim_stop("Multilevel MVN generator requires `fixed_intercept`.")
      }
      fixed <- .mlsim_check_finite_numeric(
        .mlsim_align_named_vector(fixed, vars, "fixed_intercept"),
        "fixed_intercept"
      )
      scale_fixed <- .mlsim_check_scale_fixed_intercept(scale_fixed_intercept, d, vars)
      random <- .mlsim_joint_random_intercepts(vars, random_cov, context, scale = !is.null(scale_fixed))
      random_intercepts <- random$random_intercepts
      if (is.null(scale_fixed)) {
        residual_cov <- .mlsim_check_position_cov(residual_cov, vars, "residual_cov")
        residual <- .mlsim_rmvnorm(context$n_rows, rep(0, d), residual_cov)
        row_parameters <- list(residual_cov = residual_cov)
        metadata_residual <- list(residual_cov = residual_cov, residual_cor = NULL)
        scale_eta <- NULL
      } else {
        residual_cor <- .mlsim_check_simple_residual_cor(residual_cor, vars)
        scale_eta <- sweep(
          random$scale_random_intercepts[context$group_index, , drop = FALSE],
          2L,
          scale_fixed,
          "+"
        )
        colnames(scale_eta) <- vars
        residual_sd <- exp(scale_eta)
        colnames(residual_sd) <- vars
        residual <- .mlsim_mvn_row_residual(residual_sd, residual_cor)
        row_parameters <- list(
          residual_sd = residual_sd,
          residual_cor = residual_cor,
          residual_cov = .mlsim_mvn_row_covariance_array(residual_sd, residual_cor)
        )
        metadata_residual <- list(
          residual_cov = NULL,
          residual_cor = residual_cor,
          residual_sd = residual_sd,
          scale_fixed_intercept = scale_fixed,
          scale_random_intercepts = random$scale_random_intercepts,
          scale_linear_predictor = scale_eta
        )
      }
      values <- sweep(random_intercepts[context$group_index, , drop = FALSE], 2L, fixed, "+") + residual
      metadata <- c(
        list(distribution = "mvn", level = level, vars = vars),
        .mlsim_simple_random_metadata(
          fixed,
          random,
          scale_fixed_intercept = scale_fixed,
          scale_linear_predictor = scale_eta
        ),
        list(row_parameters = row_parameters),
        metadata_residual
      )
    } else {
      n <- .mlsim_n_for_level(level, context)
      mean <- .mlsim_check_finite_numeric(
        .mlsim_align_named_vector(mean %||% rep(0, d), vars, "mean"),
        "mean"
      )
      cov <- .mlsim_check_named_square_cov(cov, vars, "cov")
      values <- .mlsim_rmvnorm(n, mean, cov)
      values <- .mlsim_expand_level2(values, level, context)
      metadata <- c(
        list(distribution = "mvn", level = level, vars = vars),
        .mlsim_parameter_metadata(level, context, mean = mean, cov = cov)
      )
    }

    colnames(values) <- vars
    out_values <- values
    if (compositional) {
      comp <- .mlsim_ilr_inverse(values, parts, total = total, coordinate_names = vars, sbp = sbp)
      out_values <- if (keep_ilr) cbind(values, comp$values) else comp$values
      metadata$compositional <- TRUE
      metadata$parts <- parts
      metadata$sbp <- comp$sbp
      metadata$basis <- comp$basis
      metadata$ilr_coordinate_map <- comp$coordinate_map
      metadata$total <- comp$total
      metadata$keep_ilr <- keep_ilr
      metadata$ilr_vars <- vars
      metadata$vars <- colnames(out_values)
    } else {
      metadata$compositional <- FALSE
    }

    .mlsim_result(out_values, colnames(out_values), metadata)
  }

  .mlsim_spec(
    "mvn", output_vars, level, simulate,
    mean = if (identical(level, "multilevel")) NULL else mean,
    cov = if (identical(level, "multilevel")) NULL else cov,
    fixed_intercept = fixed_intercept,
    random_cov = random_cov,
    residual_cov = residual_cov,
    scale_fixed_intercept = scale_fixed_intercept,
    residual_cor = residual_cor,
    compositional = compositional,
    ilr_vars = if (compositional) vars else NULL,
    parts = if (compositional) parts else NULL,
    sbp = if (compositional) sbp else NULL,
    total = total,
    keep_ilr = if (compositional) keep_ilr else FALSE
  )
}

#' Generate Categorical Variables
#'
#' Create a generator specification for binary or multicategory categorical
#' variables.
#'
#' @param vars Character scalar naming the generated variable.
#' @param level Simulation level. `"single"` generates row-level categories,
#'   `"level2"` generates group-level categories, and `"multilevel"` uses a
#'   baseline-category logit random-intercept model.
#' @param categories Vector of category values. Defaults to `c(0L, 1L)`.
#' @param prob Category probabilities for `"single"` and `"level2"`
#'   generators. For binary variables a scalar is interpreted as the success
#'   probability for the second category. Multicategory probabilities may be a
#'   vector or row-wise matrix.
#' @param fixed_intercept Baseline-category logits. For `k` categories this has
#'   length `k - 1` and is named or ordered by non-reference category.
#' @param random_cov Multilevel covariance matrix for baseline-category random
#'   intercepts.
#' @param reference Reference category for the baseline-category logit model.
#'   Defaults to the first category.
#' @param output Output type: `"factor"`, `"character"`, or `"integer"`.
#'   Integer output uses zero-based category codes.
#' @param ordered Logical; when `TRUE`, return an ordered factor. Requires
#'   `output = "factor"`.
#'
#' @return An `mlsim_generator_spec` object for use in [simulate_data()].
#'
#' @examples
#' sim <- simulate_data(
#'   n = 6,
#'   seed = 3,
#'   generators = list(
#'     arm = gen_categorical(
#'       "arm",
#'       categories = c("control", "treatment"),
#'       prob = 0.5
#'     )
#'   )
#' )
#' sim$data
#'
#' @family predictor generators
#' @export
gen_categorical <- function(vars, level = c("single", "level2", "multilevel"),
                            categories = NULL, prob = NULL,
                            fixed_intercept = NULL, random_cov = NULL,
                            reference = NULL,
                            output = c("factor", "character", "integer"),
                            ordered = FALSE) {
  prob_supplied <- !missing(prob)
  vars <- .mlsim_check_vars(vars, 1L)
  level <- .mlsim_match_level(level)
  categories <- .mlsim_check_categories(categories)
  category_values <- categories$values
  category_labels <- categories$labels
  k <- length(category_labels)
  reference_index <- .mlsim_reference_index(reference, category_labels)
  nonreference <- category_labels[-reference_index]
  output <- match.arg(output)
  ordered <- .mlsim_check_categorical_ordered(ordered, output)
  fixed_intercept <- .mlsim_check_categorical_logits(fixed_intercept, nonreference, "fixed_intercept")
  if (!identical(level, "multilevel") && !is.null(random_cov)) {
    .mlsim_stop("`random_cov` requires `level = \"multilevel\"`.")
  }
  .mlsim_check_multilevel_unused(
    level,
    c(prob = prob_supplied),
    "categorical",
    "`fixed_intercept` and `random_cov`"
  )
  .mlsim_check_multilevel_required(
    level,
    list(fixed_intercept = fixed_intercept),
    "categorical"
  )
  if (!identical(level, "multilevel") && !is.null(prob) && !is.null(fixed_intercept)) {
    .mlsim_warn_ignored_args("categorical", "prob", "fixed_intercept")
  }

  simulate <- function(context) {
    if (identical(level, "multilevel")) {
      random_cov <- .mlsim_check_categorical_random_cov(random_cov, nonreference)
      random_intercepts <- .mlsim_rmvnorm(context$n_groups, rep(0, length(nonreference)), random_cov)
      colnames(random_intercepts) <- nonreference
      eta <- sweep(random_intercepts[context$group_index, , drop = FALSE], 2L, fixed_intercept, "+")
      prob <- .mlsim_baseline_category_prob(eta, reference_index, category_labels)
      codes <- .mlsim_sample_categories(prob)
      values <- .mlsim_encode_categories(codes, category_labels, output, ordered)
      metadata <- c(
        list(
          distribution = "categorical",
          level = level,
          vars = vars,
          link = "baseline-logit",
          categories = category_values,
          category_labels = category_labels,
          reference = category_values[[reference_index]],
          reference_label = category_labels[[reference_index]],
          output = output,
          ordered = ordered,
          code_map = .mlsim_category_code_map(category_labels),
          fixed_parameters = list(intercept = fixed_intercept),
          group_parameters = list(random_cov = random_cov, random_intercepts = random_intercepts),
          fixed_intercept = fixed_intercept,
          random_cov = random_cov,
          random_intercepts = random_intercepts
        ),
        .mlsim_parameter_metadata(level, context, prob = prob)
      )
      return(.mlsim_result(values, vars, metadata))
    }

    n <- .mlsim_n_for_level(level, context)
    if (is.null(prob)) {
      if (is.null(fixed_intercept)) {
        prob_n <- .mlsim_uniform_category_prob(n, category_labels)
      } else {
        eta <- matrix(rep(fixed_intercept, n), nrow = n, byrow = TRUE)
        colnames(eta) <- nonreference
        prob_n <- .mlsim_baseline_category_prob(eta, reference_index, category_labels)
      }
    } else {
      prob_n <- .mlsim_resolve_category_prob(prob, n, category_labels, "prob")
    }
    codes <- .mlsim_sample_categories(prob_n)
    codes <- .mlsim_expand_level2_vector(codes, level, context)
    values <- .mlsim_encode_categories(codes, category_labels, output, ordered)
    .mlsim_result(
      values,
      vars,
      c(
        list(
          distribution = "categorical",
          level = level,
          vars = vars,
          link = "baseline-logit",
          categories = category_values,
          category_labels = category_labels,
          reference = category_values[[reference_index]],
          reference_label = category_labels[[reference_index]],
          output = output,
          ordered = ordered,
          code_map = .mlsim_category_code_map(category_labels),
          fixed_parameters = if (is.null(fixed_intercept)) NULL else list(intercept = fixed_intercept),
          fixed_intercept = fixed_intercept
        ),
        .mlsim_parameter_metadata(level, context, prob = prob_n)
      )
    )
  }

  .mlsim_spec(
    "categorical", vars, level, simulate,
    categories = category_values,
    reference = category_values[[reference_index]],
    output = output,
    ordered = ordered,
    prob = if (identical(level, "multilevel")) NULL else prob,
    fixed_intercept = fixed_intercept,
    random_cov = random_cov
  )
}

#' Generate Count Variables
#'
#' Create generator specifications for binomial, Poisson, and negative
#' binomial count variables.
#'
#' @param vars Character scalar naming the generated variable.
#' @param level Simulation level. `"single"` generates row-level values,
#'   `"level2"` generates group-level values, and `"multilevel"` uses a
#'   random-intercept model.
#' @param size Binomial trial size or negative-binomial size/dispersion
#'   parameter. May be a scalar, vector, function of `n`, or count-distribution
#'   list where supported.
#' @param prob Binomial success probability for `"single"` and `"level2"`
#'   generators.
#' @param fixed_intercept Link-scale intercept for direct-parameter defaults
#'   and multilevel models. Binomial uses the logit link; Poisson and negative
#'   binomial use the log link.
#' @param random_cov Group-level random-intercept covariance for multilevel
#'   generators.
#' @param lambda Poisson mean for `"single"` and `"level2"` generators.
#' @param mu Negative-binomial mean for `"single"` and `"level2"` generators.
#' @param scale_fixed_intercept Optional log size intercept for multilevel
#'   negative-binomial generators.
#'
#' @return An `mlsim_generator_spec` object for use in [simulate_data()].
#'
#' @details
#' For multilevel count generators, `fixed_intercept` and `random_cov` define
#' the link-scale mean. Non-multilevel generators can use direct distribution
#' parameters or a link-scale `fixed_intercept`.
#'
#' @examples
#' sim <- simulate_data(
#'   n = 5,
#'   seed = 4,
#'   generators = list(
#'     events = gen_poisson("events", lambda = 2),
#'     successes = gen_binomial("successes", size = 4, prob = 0.5),
#'     overdispersed = gen_negbin("overdispersed", size = 3, mu = 2)
#'   )
#' )
#' sim$data
#'
#' @family predictor generators
#' @name count-generators
NULL

#' @rdname count-generators
#' @export
gen_binomial <- function(vars, level = c("single", "level2", "multilevel"),
                         size = 1, prob = NULL, fixed_intercept = NULL,
                         random_cov = NULL) {
  prob_supplied <- !missing(prob)
  vars <- .mlsim_check_vars(vars, 1L)
  level <- .mlsim_match_level(level)
  .mlsim_check_multilevel_only_args(level, c(random_cov = !is.null(random_cov)))
  .mlsim_check_multilevel_unused(
    level,
    c(prob = prob_supplied),
    "binomial",
    "`fixed_intercept`, `random_cov`, and `size`"
  )
	  .mlsim_check_multilevel_required(
	    level,
	    list(fixed_intercept = fixed_intercept),
	    "binomial"
	  )
	  if (!identical(level, "multilevel") && is.null(prob) && is.null(fixed_intercept)) {
	    .mlsim_stop("Binomial generator requires `prob` or `fixed_intercept`.")
	  }
	  if (!identical(level, "multilevel") && !is.null(prob) && !is.null(fixed_intercept)) {
	    .mlsim_warn_ignored_args("binomial", "prob", "fixed_intercept")
	  }

  simulate <- function(context) {
    if (identical(level, "multilevel")) {
      size_n <- .mlsim_resolve_size(size, context$n_rows, "size")
      lp <- .mlsim_random_eta("binomial", fixed_intercept, random_cov, context, vars)
      prob <- .mlsim_mean_from_eta("binomial", lp$eta)
      values <- stats::rbinom(context$n_rows, size = size_n, prob = prob)
      metadata <- c(
        lp$metadata,
        list(distribution = "binomial", level = level, vars = vars, link = "logit"),
        .mlsim_parameter_metadata(level, context, size = size_n, prob = prob)
      )
      return(.mlsim_result(values, vars, metadata))
    }

    n <- .mlsim_n_for_level(level, context)
    size_n <- .mlsim_resolve_size(size, n, "size")
    prob_n <- prob
    if (is.null(prob_n)) {
      prob_n <- .mlsim_mean_from_eta("binomial", .mlsim_recycle(fixed_intercept, 1L, "fixed_intercept"))
    }
    prob_n <- .mlsim_check_prob(.mlsim_recycle(prob_n, n, "prob"))
    values <- stats::rbinom(n, size = size_n, prob = prob_n)
    values <- .mlsim_expand_level2_vector(values, level, context)
    .mlsim_result(
      values,
      vars,
      c(
        list(distribution = "binomial", level = level, vars = vars, link = "logit"),
        .mlsim_parameter_metadata(level, context, size = size_n, prob = prob_n)
      )
    )
  }

  .mlsim_spec(
    "binomial", vars, level, simulate,
    size = size,
    prob = if (identical(level, "multilevel")) NULL else prob,
    fixed_intercept = fixed_intercept,
    random_cov = random_cov
  )
}

#' @rdname count-generators
#' @export
gen_poisson <- function(vars, level = c("single", "level2", "multilevel"),
                        lambda = NULL, fixed_intercept = NULL,
                        random_cov = NULL) {
  lambda_supplied <- !missing(lambda)
  vars <- .mlsim_check_vars(vars, 1L)
  level <- .mlsim_match_level(level)
  .mlsim_check_multilevel_only_args(level, c(random_cov = !is.null(random_cov)))
  .mlsim_check_multilevel_unused(
    level,
    c(lambda = lambda_supplied),
    "poisson",
    "`fixed_intercept` and `random_cov`"
  )
	  .mlsim_check_multilevel_required(
	    level,
	    list(fixed_intercept = fixed_intercept),
	    "poisson"
	  )
	  if (!identical(level, "multilevel") && is.null(lambda) && is.null(fixed_intercept)) {
	    .mlsim_stop("Poisson generator requires `lambda` or `fixed_intercept`.")
	  }
	  if (!identical(level, "multilevel") && !is.null(lambda) && !is.null(fixed_intercept)) {
	    .mlsim_warn_ignored_args("poisson", "lambda", "fixed_intercept")
	  }

  simulate <- function(context) {
    if (identical(level, "multilevel")) {
      lp <- .mlsim_random_eta("poisson", fixed_intercept, random_cov, context, vars)
      lambda <- .mlsim_mean_from_eta("poisson", lp$eta)
      values <- stats::rpois(context$n_rows, lambda = lambda)
      metadata <- c(
        lp$metadata,
        list(distribution = "poisson", level = level, vars = vars, link = "log"),
        .mlsim_parameter_metadata(level, context, lambda = lambda)
      )
      return(.mlsim_result(values, vars, metadata))
    }

    n <- .mlsim_n_for_level(level, context)
    lambda_n <- lambda
    if (is.null(lambda_n)) {
      lambda_n <- .mlsim_mean_from_eta("poisson", .mlsim_recycle(fixed_intercept, 1L, "fixed_intercept"))
    }
    lambda_n <- .mlsim_check_nonnegative(.mlsim_recycle(lambda_n, n, "lambda"), "lambda")
    values <- stats::rpois(n, lambda = lambda_n)
    values <- .mlsim_expand_level2_vector(values, level, context)
    .mlsim_result(
      values,
      vars,
      c(
        list(distribution = "poisson", level = level, vars = vars, link = "log"),
        .mlsim_parameter_metadata(level, context, lambda = lambda_n)
      )
    )
  }

  .mlsim_spec(
    "poisson", vars, level, simulate,
    lambda = if (identical(level, "multilevel")) NULL else lambda,
    fixed_intercept = fixed_intercept,
    random_cov = random_cov
  )
}

#' @rdname count-generators
#' @export
gen_negbin <- function(vars, level = c("single", "level2", "multilevel"),
                       size = NULL, mu = NULL, fixed_intercept = NULL,
                       random_cov = NULL, scale_fixed_intercept = NULL) {
  mu_supplied <- !missing(mu)
  vars <- .mlsim_check_vars(vars, 1L)
  level <- .mlsim_match_level(level)
  .mlsim_validate_scale_multilevel_args(level, scale_fixed_intercept, random_cov)
  .mlsim_check_multilevel_unused(
    level,
    c(mu = mu_supplied),
    "negative binomial",
    "`fixed_intercept`, `random_cov`, `size`, and `scale_fixed_intercept`"
  )
  .mlsim_check_multilevel_required(
    level,
    list(fixed_intercept = fixed_intercept),
    "negative binomial"
  )
	  if (identical(level, "multilevel")) {
	    .mlsim_check_scale_or_arg(size, "size", scale_fixed_intercept, "negative binomial")
	  } else if (is.null(size)) {
	    .mlsim_stop("Negative binomial generator requires `size`.")
	  }
	  if (!identical(level, "multilevel") && is.null(mu) && is.null(fixed_intercept)) {
	    .mlsim_stop("Negative binomial generator requires `mu` or `fixed_intercept`.")
	  }
	  size <- .mlsim_check_positive_or_null(size, "size")
	  if (!identical(level, "multilevel") && !is.null(mu) && !is.null(fixed_intercept)) {
	    .mlsim_warn_ignored_args("negative binomial", "mu", "fixed_intercept")
  }

  simulate <- function(context) {
    if (identical(level, "multilevel")) {
      lp <- .mlsim_random_eta(
        "negbin",
        fixed_intercept,
        random_cov,
        context,
        vars,
        scale_fixed_intercept = scale_fixed_intercept
      )
      size_n <- if (is.null(scale_fixed_intercept)) {
        .mlsim_recycle(size, context$n_rows, "size")
      } else {
        exp(lp$scale_eta)
      }
      mu <- .mlsim_mean_from_eta("negbin", lp$eta)
      values <- stats::rnbinom(context$n_rows, mu = mu, size = size_n)
      metadata <- c(
        lp$metadata,
        list(distribution = "negative_binomial", level = level, vars = vars, link = "log"),
        .mlsim_parameter_metadata(level, context, size = size_n, mu = mu),
        list(size = size_n)
      )
      return(.mlsim_result(values, vars, metadata))
    }

    n <- .mlsim_n_for_level(level, context)
    size_n <- .mlsim_recycle(size, n, "size")
    mu_n <- mu
    if (is.null(mu_n)) {
      eta <- .mlsim_recycle(fixed_intercept, 1L, "fixed_intercept")
      mu_n <- .mlsim_mean_from_eta("negbin", eta)
    }
    mu_n <- .mlsim_check_positive(.mlsim_recycle(mu_n, n, "mu"), "mu")
    values <- stats::rnbinom(n, mu = mu_n, size = size_n)
    values <- .mlsim_expand_level2_vector(values, level, context)
    .mlsim_result(
      values,
      vars,
      c(
        list(distribution = "negative_binomial", level = level, vars = vars, link = if (is.null(mu)) "log" else NA_character_),
        .mlsim_parameter_metadata(level, context, size = size_n, mu = mu_n)
      )
    )
  }

  .mlsim_spec(
    "negbin", vars, level, simulate,
    size = size,
    mu = if (identical(level, "multilevel")) NULL else mu,
    fixed_intercept = fixed_intercept,
    random_cov = random_cov,
    scale_fixed_intercept = scale_fixed_intercept
  )
}

#' Generate Positive and Bounded Continuous Variables
#'
#' Create generator specifications for gamma and beta variables.
#'
#' @param vars Character scalar naming the generated variable.
#' @param level Simulation level. `"single"` generates row-level values,
#'   `"level2"` generates group-level values, and `"multilevel"` uses a
#'   random-intercept model.
#' @param shape Gamma shape parameter.
#' @param mean Mean parameter. For gamma variables this must be positive; for
#'   beta variables it must be in `(0, 1)`.
#' @param rate,scale Alternative gamma parameterizations for `"single"` and
#'   `"level2"` generators. Supply at most one of `rate` or `scale`.
#' @param fixed_intercept Link-scale intercept. Gamma uses the log link; beta
#'   uses the logit link.
#' @param random_cov Group-level random-intercept covariance for multilevel
#'   generators.
#' @param scale_fixed_intercept Optional log shape or log precision intercept
#'   for multilevel generators.
#' @param precision Beta precision parameter.
#'
#' @return An `mlsim_generator_spec` object for use in [simulate_data()].
#'
#' @examples
#' sim <- simulate_data(
#'   n = 5,
#'   seed = 5,
#'   generators = list(
#'     positive = gen_gamma("positive", shape = 2, mean = 1.5),
#'     proportion = gen_beta("proportion", mean = 0.4, precision = 10)
#'   )
#' )
#' sim$data
#'
#' @family predictor generators
#' @name continuous-generators
NULL

#' @rdname continuous-generators
#' @export
gen_gamma <- function(vars, level = c("single", "level2", "multilevel"),
                      shape = NULL, mean = NULL, rate = NULL, scale = NULL,
                      fixed_intercept = NULL, random_cov = NULL,
                      scale_fixed_intercept = NULL) {
  mean_supplied <- !missing(mean)
  rate_supplied <- !missing(rate)
  scale_supplied <- !missing(scale)
  vars <- .mlsim_check_vars(vars, 1L)
  level <- .mlsim_match_level(level)
  .mlsim_validate_scale_multilevel_args(level, scale_fixed_intercept, random_cov)
  .mlsim_check_multilevel_unused(
    level,
    c(mean = mean_supplied, rate = rate_supplied, scale = scale_supplied),
    "gamma",
    "`fixed_intercept`, `random_cov`, `shape`, and `scale_fixed_intercept`"
  )
  .mlsim_check_multilevel_required(
    level,
    list(fixed_intercept = fixed_intercept),
    "gamma"
  )
	  if (identical(level, "multilevel")) {
	    .mlsim_check_scale_or_arg(shape, "shape", scale_fixed_intercept, "gamma")
	  } else if (is.null(shape)) {
	    .mlsim_stop("Gamma generator requires `shape`.")
	  }
	  if (!identical(level, "multilevel") && is.null(mean) && is.null(rate) &&
	      is.null(scale) && is.null(fixed_intercept)) {
	    .mlsim_stop("Gamma generator requires `mean`, `rate`, `scale`, or `fixed_intercept`.")
	  }
	  shape <- .mlsim_check_positive_or_null(shape, "shape")
	  if (!identical(level, "multilevel")) {
	    if (!is.null(rate) || !is.null(scale)) {
      .mlsim_warn_ignored_args(
        "gamma",
        c(if (!is.null(rate)) "rate", if (!is.null(scale)) "scale"),
        c(if (!is.null(mean)) "mean", if (!is.null(fixed_intercept)) "fixed_intercept")
      )
    } else if (!is.null(mean) && !is.null(fixed_intercept)) {
      .mlsim_warn_ignored_args("gamma", "mean", "fixed_intercept")
    }
  }

  simulate <- function(context) {
    if (identical(level, "multilevel")) {
      lp <- .mlsim_random_eta(
        "gamma",
        fixed_intercept,
        random_cov,
        context,
        vars,
        scale_fixed_intercept = scale_fixed_intercept
      )
      mean <- .mlsim_mean_from_eta("gamma", lp$eta)
      shape_n <- if (is.null(scale_fixed_intercept)) {
        .mlsim_recycle(shape, context$n_rows, "shape")
      } else {
        exp(lp$scale_eta)
      }
      values <- stats::rgamma(context$n_rows, shape = shape_n, scale = mean / shape_n)
      metadata <- c(
        lp$metadata,
        list(distribution = "gamma", level = level, vars = vars, link = "log"),
        .mlsim_parameter_metadata(level, context, shape = shape_n, mean = mean),
        list(shape = shape_n)
      )
      return(.mlsim_result(values, vars, metadata))
    }

    n <- .mlsim_n_for_level(level, context)
    shape_n <- .mlsim_recycle(shape, n, "shape")
    if (!is.null(rate) || !is.null(scale)) {
      if (!is.null(rate) && !is.null(scale)) {
        .mlsim_stop("Gamma generator accepts only one of `rate` or `scale`.")
      }
      if (!is.null(rate)) {
        rate_n <- .mlsim_check_positive(.mlsim_recycle(rate, n, "rate"), "rate")
        values <- stats::rgamma(n, shape = shape_n, rate = rate_n)
        metadata_params <- list(rate = rate_n)
      } else {
        scale_n <- .mlsim_check_positive(.mlsim_recycle(scale, n, "scale"), "scale")
        values <- stats::rgamma(n, shape = shape_n, scale = scale_n)
        metadata_params <- list(scale = scale_n)
      }
    } else {
      mean_n <- mean
      if (is.null(mean_n)) {
        mean_n <- .mlsim_mean_from_eta("gamma", .mlsim_recycle(fixed_intercept, 1L, "fixed_intercept"))
      }
      mean_n <- .mlsim_check_positive(.mlsim_recycle(mean_n, n, "mean"), "mean")
      values <- stats::rgamma(n, shape = shape_n, scale = mean_n / shape_n)
      metadata_params <- list(mean = mean_n)
    }
    values <- .mlsim_expand_level2_vector(values, level, context)
    .mlsim_result(
      values,
      vars,
      c(
        list(distribution = "gamma", level = level, vars = vars, link = if (is.null(mean) && is.null(rate) && is.null(scale)) "log" else NA_character_),
        do.call(
          .mlsim_parameter_metadata,
          c(list(level = level, context = context, shape = shape_n), metadata_params)
        )
      )
    )
  }

  .mlsim_spec(
    "gamma", vars, level, simulate,
    shape = shape,
    mean = if (identical(level, "multilevel")) NULL else mean,
    rate = if (identical(level, "multilevel")) NULL else rate,
    scale = if (identical(level, "multilevel")) NULL else scale,
    fixed_intercept = fixed_intercept,
    random_cov = random_cov,
    scale_fixed_intercept = scale_fixed_intercept
  )
}

#' @rdname continuous-generators
#' @export
gen_beta <- function(vars, level = c("single", "level2", "multilevel"),
                     mean = NULL, precision = NULL,
                     fixed_intercept = NULL, random_cov = NULL,
                     scale_fixed_intercept = NULL) {
  mean_supplied <- !missing(mean)
  vars <- .mlsim_check_vars(vars, 1L)
  level <- .mlsim_match_level(level)
  .mlsim_validate_scale_multilevel_args(level, scale_fixed_intercept, random_cov)
  .mlsim_check_multilevel_unused(
    level,
    c(mean = mean_supplied),
    "beta",
    "`fixed_intercept`, `random_cov`, `precision`, and `scale_fixed_intercept`"
  )
  .mlsim_check_multilevel_required(
    level,
    list(fixed_intercept = fixed_intercept),
    "beta"
  )
	  if (identical(level, "multilevel")) {
	    .mlsim_check_scale_or_arg(precision, "precision", scale_fixed_intercept, "beta")
	  } else if (is.null(precision)) {
	    .mlsim_stop("Beta generator requires `precision`.")
	  }
	  if (!identical(level, "multilevel") && is.null(mean) && is.null(fixed_intercept)) {
	    .mlsim_stop("Beta generator requires `mean` or `fixed_intercept`.")
	  }
	  precision <- .mlsim_check_positive_or_null(precision, "precision")
	  if (!identical(level, "multilevel") && !is.null(mean) && !is.null(fixed_intercept)) {
	    .mlsim_warn_ignored_args("beta", "mean", "fixed_intercept")
  }

  simulate <- function(context) {
    if (identical(level, "multilevel")) {
      lp <- .mlsim_random_eta(
        "beta",
        fixed_intercept,
        random_cov,
        context,
        vars,
        scale_fixed_intercept = scale_fixed_intercept
      )
      mean <- .mlsim_check_open_prob(.mlsim_mean_from_eta("beta", lp$eta), "mean")
      precision_n <- if (is.null(scale_fixed_intercept)) {
        .mlsim_recycle(precision, context$n_rows, "precision")
      } else {
        exp(lp$scale_eta)
      }
      values <- stats::rbeta(context$n_rows, shape1 = mean * precision_n, shape2 = (1 - mean) * precision_n)
      metadata <- c(
        lp$metadata,
        list(distribution = "beta", level = level, vars = vars, link = "logit"),
        .mlsim_parameter_metadata(level, context, mean = mean, precision = precision_n),
        list(precision = precision_n)
      )
      return(.mlsim_result(values, vars, metadata))
    }

    n <- .mlsim_n_for_level(level, context)
    mean_n <- mean
    if (is.null(mean_n)) {
      mean_n <- .mlsim_mean_from_eta("beta", .mlsim_recycle(fixed_intercept, 1L, "fixed_intercept"))
    }
    mean_n <- .mlsim_check_open_prob(.mlsim_recycle(mean_n, n, "mean"), "mean")
    precision_n <- .mlsim_recycle(precision, n, "precision")
    values <- stats::rbeta(n, shape1 = mean_n * precision_n, shape2 = (1 - mean_n) * precision_n)
    values <- .mlsim_expand_level2_vector(values, level, context)
    .mlsim_result(
      values,
      vars,
      c(
        list(distribution = "beta", level = level, vars = vars, link = "logit"),
        .mlsim_parameter_metadata(level, context, mean = mean_n, precision = precision_n)
      )
    )
  }

  .mlsim_spec(
    "beta", vars, level, simulate,
    mean = if (identical(level, "multilevel")) NULL else mean,
    precision = precision,
    fixed_intercept = fixed_intercept,
    random_cov = random_cov,
    scale_fixed_intercept = scale_fixed_intercept
  )
}

#' Generate Variables with a User-Supplied Function
#'
#' Wrap a custom generator function so it can participate in [simulate_data()].
#'
#' @param vars Character vector naming the generated variables.
#' @param generator Function called as
#'   `generator(context = context, vars = vars, level = level, ...)`. It must
#'   return a list with `data` and `names` elements, and may include
#'   `metadata`.
#' @param level Simulation level. A level-2 custom generator may return either
#'   one row per group or one row per simulated observation.
#' @param ... Additional arguments passed to `generator`.
#' @param metadata List of wrapper metadata stored on the generator
#'   specification.
#'
#' @return An `mlsim_generator_spec` object for use in [simulate_data()].
#'
#' @examples
#' constant_generator <- function(context, vars, level, value = 1) {
#'   list(data = data.frame(value = rep(value, context$n_rows)), names = vars)
#' }
#'
#' sim <- simulate_data(
#'   n = 3,
#'   generators = list(x = gen_custom("x", constant_generator, value = 7))
#' )
#' sim$data
#'
#' @family predictor generators
#' @export
gen_custom <- function(vars, generator, level = c("single", "level2", "multilevel"),
                       ..., metadata = list()) {
  vars <- .mlsim_check_vars(vars)
  level <- .mlsim_match_level(level)
  if (!is.function(generator)) {
    .mlsim_stop("`generator` must be a function.")
  }
  if (!is.list(metadata)) {
    .mlsim_stop("`metadata` must be a list.")
  }
  user_args <- list(...)

  simulate <- function(context) {
    result <- do.call(
      generator,
      c(
        list(context = context, vars = vars, level = level),
        user_args
      )
    )
    .mlsim_custom_result(result, vars, level, context, metadata)
  }

  .mlsim_spec(
    "custom", vars, level, simulate,
    metadata = metadata,
    argument_names = names(user_args)
  )
}
