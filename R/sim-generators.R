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

#' @rdname multilevelcoda-internal-generator-specs
.mlsim_check_unused_generator_dots <- function(...) {
  dots <- list(...)
  if (length(dots) == 0L) {
    return(invisible(NULL))
  }
  args <- names(dots)
  args[is.na(args) | args == ""] <- "<unnamed>"
  .mlsim_stop(
    "unused argument%s (%s)",
    if (length(args) == 1L) "" else "s",
    paste(args, collapse = ", ")
  )
}

#' @rdname multilevelcoda-internal-generator-specs
.mlsim_gen_distribution <- function(distribution, vars, level,
                                    fixed_intercept = NULL,
                                    random_cov = NULL,
                                    residual_cov = NULL,
                                    scale_fixed_intercept = NULL,
                                    residual_cor = NULL,
                                    compositional = FALSE, parts = NULL,
                                    total = 1, keep_ilr = TRUE, sbp = NULL,
                                    categories = NULL, reference = NULL,
                                    output = c("factor", "character", "integer"),
                                    ordered = FALSE,
                                    size = NULL,
                                    generator = NULL, user_args = list(),
                                    metadata = list()) {
  distribution <- match.arg(
    distribution,
    c("mvn", "categorical", "binomial", "poisson", "negbin", "gamma", "beta", "custom")
  )
  level <- .mlsim_match_level(level)

  require_parameter <- function(generator_label, name, value) {
    if (is.null(value)) {
      .mlsim_stop("%s generator requires `%s`.", generator_label, name)
    }
    invisible(NULL)
  }

  check_random_contract <- function(generator_label) {
    if (identical(level, "multilevel")) {
      require_parameter(paste("Multilevel", generator_label), "random_cov", random_cov)
    } else if (!is.null(random_cov)) {
      .mlsim_stop("`random_cov` requires `level = \"multilevel\"`.")
    }
    invisible(NULL)
  }

  check_fixed_scalar <- function(generator_label) {
    require_parameter(generator_label, "fixed_intercept", fixed_intercept)
    .mlsim_check_scalar_number(fixed_intercept, "fixed_intercept")
  }

  check_scale_scalar <- function(generator_label) {
    require_parameter(generator_label, "scale_fixed_intercept", scale_fixed_intercept)
    .mlsim_check_scale_fixed_intercept(scale_fixed_intercept, 1L)
  }

  linear_predictor <- function(context, fixed, vars, scale_fixed = NULL) {
    fixed <- as.numeric(fixed)
    names(fixed) <- vars

    random <- NULL
    if (identical(level, "multilevel")) {
      random <- .mlsim_joint_random_intercepts(
        vars,
        random_cov,
        context,
        scale = !is.null(scale_fixed)
      )
      eta <- sweep(
        random$random_intercepts[context$group_index, , drop = FALSE],
        2L,
        fixed,
        "+"
      )
    } else {
      n <- .mlsim_n_for_level(level, context)
      eta <- matrix(rep(fixed, each = n), nrow = n, ncol = length(vars))
      colnames(eta) <- vars
    }

    scale_eta <- NULL
    if (!is.null(scale_fixed)) {
      scale_fixed <- as.numeric(scale_fixed)
      names(scale_fixed) <- vars
      if (identical(level, "multilevel") && !is.null(random$scale_random_intercepts)) {
        scale_eta <- sweep(
          random$scale_random_intercepts[context$group_index, , drop = FALSE],
          2L,
          scale_fixed,
          "+"
        )
      } else {
        n <- if (identical(level, "multilevel")) context$n_rows else .mlsim_n_for_level(level, context)
        scale_eta <- matrix(rep(scale_fixed, each = n), nrow = n, ncol = length(vars))
      }
      colnames(scale_eta) <- vars
    }

    list(
      eta = eta,
      scale_eta = scale_eta,
      random = random,
      internal = if (is.null(random)) NULL else .mlsim_multilevel_random_columns(random, context, vars)
    )
  }

  scalar_linear_predictor <- function(context, type, fixed, scale_fixed = NULL) {
    lp <- linear_predictor(context, fixed, vars, scale_fixed = scale_fixed)
    metadata <- list(
      fixed_parameters = list(intercept = fixed),
      fixed_intercept = fixed
    )
    if (!is.null(scale_fixed)) {
      metadata$fixed_parameters$scale_intercept <- scale_fixed
      metadata$scale_fixed_intercept <- scale_fixed
    }
    if (!is.null(lp$random)) {
      metadata$group_parameters <- list(random_cov = lp$random$random_cov)
      metadata$random_cov <- lp$random$random_cov
    }
    c(
      list(
        eta = as.vector(lp$eta),
        scale_eta = if (is.null(lp$scale_eta)) NULL else as.vector(lp$scale_eta),
        random = lp$random,
        internal = lp$internal
      ),
      list(metadata = metadata)
    )
  }

  switch(
    distribution,
    mvn = {
      vars <- .mlsim_check_vars(vars, NULL)
      compositional <- .mlsim_check_scalar_logical(compositional, "compositional")
      keep_ilr <- .mlsim_check_scalar_logical(keep_ilr, "keep_ilr")
      d <- length(vars)
      require_parameter("MVN", "fixed_intercept", fixed_intercept)
      fixed <- .mlsim_check_finite_numeric(
        .mlsim_align_named_vector(fixed_intercept, vars, "fixed_intercept"),
        "fixed_intercept"
      )
      scale_fixed <- .mlsim_check_scale_fixed_intercept(scale_fixed_intercept, d, vars)
      check_random_contract("MVN")

      if (!is.null(residual_cov) && !is.null(scale_fixed)) {
        .mlsim_stop("`residual_cov` cannot be supplied with `scale_fixed_intercept`.")
      }
      if (is.null(residual_cov) && is.null(scale_fixed)) {
        .mlsim_stop("MVN generator requires `residual_cov` or `scale_fixed_intercept`.")
      }
      if (is.null(scale_fixed) && !is.null(residual_cor)) {
        .mlsim_stop("`residual_cor` requires `scale_fixed_intercept`.")
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
        base <- if (keep_ilr) c(vars, parts) else parts
        if (identical(level, "multilevel")) {
          base <- c(
            base,
            paste0(vars, "_between"),
            paste0(vars, "_within"),
            paste0(parts, "_between")
          )
        }
        .mlsim_check_vars(base, NULL, "compositional output columns")
      } else {
        vars
      }

      simulate <- function(context) {
        lp <- linear_predictor(context, fixed, vars, scale_fixed = scale_fixed)
        n <- nrow(lp$eta)

        if (is.null(scale_fixed)) {
          residual_cov_checked <- .mlsim_check_position_cov(residual_cov, vars, "residual_cov")
          residual <- .mlsim_rmvnorm(n, rep(0, d), residual_cov_checked)
          row_parameters <- list(mean = lp$eta, residual_cov = residual_cov_checked)
          metadata_residual <- list(residual_cov = residual_cov_checked, residual_cor = NULL)
        } else {
          residual_cor_checked <- .mlsim_check_simple_residual_cor(residual_cor, vars)
          residual_sd <- exp(lp$scale_eta)
          colnames(residual_sd) <- vars
          residual <- .mlsim_mvn_row_residual(residual_sd, residual_cor_checked)
          row_parameters <- list(
            mean = lp$eta,
            residual_sd = residual_sd,
            residual_cor = residual_cor_checked,
            residual_cov = .mlsim_mvn_row_covariance_array(residual_sd, residual_cor_checked)
          )
          metadata_residual <- list(
            residual_cov = NULL,
            residual_cor = residual_cor_checked,
            residual_sd = residual_sd,
            scale_fixed_intercept = scale_fixed
          )
        }

        values <- lp$eta + residual
        between_values <- NULL
        within_values <- NULL
        if (identical(level, "multilevel")) {
          between_values <- lp$eta
          within_values <- residual
          colnames(between_values) <- paste0(vars, "_between")
          colnames(within_values) <- paste0(vars, "_within")
        }
        if (identical(level, "level2")) {
          values <- values[context$group_index, , drop = FALSE]
        }
        colnames(values) <- vars

        metadata <- c(
          list(distribution = "mvn", level = level, vars = vars, link = "identity"),
          list(
            fixed_parameters = list(intercept = fixed),
            fixed_intercept = fixed
          ),
          if (!is.null(lp$random)) {
            list(
              group_parameters = list(random_cov = lp$random$random_cov),
              random_cov = lp$random$random_cov
            )
          } else {
            list()
          },
          do.call(
            .mlsim_parameter_metadata,
            c(list(level = level, context = context), row_parameters)
          ),
          metadata_residual
        )

        if (!is.null(scale_fixed)) {
          metadata$fixed_parameters$scale_intercept <- scale_fixed
          metadata$scale_fixed_intercept <- scale_fixed
        }

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
          if (!is.null(between_values)) {
            # the latent between composition is the closed back-transform of the
            # group-level ILR means, not the mean of the observed parts
            between_parts <- .mlsim_ilr_inverse(
              lp$eta,
              parts,
              total = total,
              coordinate_names = vars,
              sbp = sbp
            )$values
            colnames(between_parts) <- paste0(parts, "_between")
            out_values <- cbind(out_values, between_values, within_values, between_parts)
            metadata$between_parts_vars <- colnames(between_parts)
            metadata$column_roles <- .mlsim_column_roles(
              column = c(
                if (keep_ilr) vars,
                colnames(between_values),
                colnames(within_values)
              ),
              variable = c(if (keep_ilr) vars, vars, vars),
              component = c(
                if (keep_ilr) rep("observed", length(vars)),
                rep("between", length(vars)),
                rep("within", length(vars))
              ),
              level = c(
                if (keep_ilr) rep("row", length(vars)),
                rep("group", length(vars)),
                rep("row", length(vars))
              )
            )
          } else if (keep_ilr) {
            component <- if (identical(level, "level2")) "between" else "observed"
            role_level <- if (identical(level, "level2")) "group" else "row"
            metadata$column_roles <- .mlsim_column_roles(
              column = vars,
              variable = vars,
              component = rep(component, length(vars)),
              level = rep(role_level, length(vars))
            )
          }
          metadata$vars <- colnames(out_values)
        } else {
          metadata$compositional <- FALSE
          if (!is.null(between_values)) {
            out_values <- cbind(out_values, between_values, within_values)
            metadata$column_roles <- .mlsim_column_roles(
              column = c(vars, colnames(between_values), colnames(within_values)),
              variable = rep(vars, times = 3L),
              component = rep(c("observed", "between", "within"), each = length(vars)),
              level = rep(c("row", "group", "row"), each = length(vars))
            )
          } else if (identical(level, "level2")) {
            metadata$column_roles <- .mlsim_column_roles(
              column = vars,
              variable = vars,
              component = rep("between", length(vars)),
              level = rep("group", length(vars))
            )
          } else {
            metadata$column_roles <- .mlsim_column_roles(
              column = vars,
              variable = vars,
              component = rep("observed", length(vars)),
              level = rep("row", length(vars))
            )
          }
        }

        .mlsim_result(out_values, colnames(out_values), metadata, internal = lp$internal)
      }

      .mlsim_spec(
        "mvn", output_vars, level, simulate,
        fixed_intercept = fixed,
        random_cov = random_cov,
        residual_cov = residual_cov,
        scale_fixed_intercept = scale_fixed,
        residual_cor = residual_cor,
        compositional = compositional,
        ilr_vars = if (compositional) vars else NULL,
        parts = if (compositional) parts else NULL,
        sbp = if (compositional) sbp else NULL,
        total = total,
        keep_ilr = if (compositional) keep_ilr else FALSE
      )
    },
    categorical = {
      vars <- .mlsim_check_vars(vars, 1L)
      categories <- .mlsim_check_categories(categories)
      category_values <- categories$values
      category_labels <- categories$labels
      reference_index <- .mlsim_reference_index(reference, category_labels)
      nonreference <- category_labels[-reference_index]
      output <- match.arg(output)
      ordered <- .mlsim_check_categorical_ordered(ordered, output)
      require_parameter("Categorical", "fixed_intercept", fixed_intercept)
      fixed <- .mlsim_check_categorical_logits(fixed_intercept, nonreference, "fixed_intercept")
      check_random_contract("categorical")

      simulate <- function(context) {
        internal <- NULL
        random_cov_checked <- NULL
        if (identical(level, "multilevel")) {
          random_cov_checked <- .mlsim_check_categorical_random_cov(random_cov, nonreference)
          random_intercepts <- .mlsim_rmvnorm(
            context$n_groups,
            rep(0, length(nonreference)),
            random_cov_checked
          )
          colnames(random_intercepts) <- nonreference
          eta <- sweep(random_intercepts[context$group_index, , drop = FALSE], 2L, fixed, "+")
          internal <- .mlsim_categorical_random_columns(random_intercepts, context, vars)
        } else {
          n <- .mlsim_n_for_level(level, context)
          eta <- matrix(rep(fixed, each = n), nrow = n, ncol = length(nonreference))
          colnames(eta) <- nonreference
        }

        prob <- .mlsim_baseline_category_prob(eta, reference_index, category_labels)
        codes <- .mlsim_sample_categories(prob)
        if (identical(level, "level2")) {
          codes <- codes[context$group_index]
        }
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
            fixed_parameters = list(intercept = fixed),
            fixed_intercept = fixed
          ),
          if (!is.null(random_cov_checked)) {
            list(
              group_parameters = list(random_cov = random_cov_checked),
              random_cov = random_cov_checked
            )
          } else {
            list()
          },
          .mlsim_parameter_metadata(level, context, prob = prob)
        )

        .mlsim_result(values, vars, metadata, internal = internal)
      }

      .mlsim_spec(
        "categorical", vars, level, simulate,
        categories = category_values,
        reference = category_values[[reference_index]],
        output = output,
        ordered = ordered,
        fixed_intercept = fixed,
        random_cov = random_cov
      )
    },
    binomial = {
      vars <- .mlsim_check_vars(vars, 1L)
      require_parameter("Binomial", "size", size)
      fixed <- check_fixed_scalar("Binomial")
      check_random_contract("binomial")

      simulate <- function(context) {
        lp <- scalar_linear_predictor(context, "binomial", fixed)
        n <- length(lp$eta)
        size_n <- .mlsim_resolve_size(size, n, "size")
        prob <- .mlsim_mean_from_eta("binomial", lp$eta)
        values <- stats::rbinom(n, size = size_n, prob = prob)
        if (identical(level, "level2")) {
          values <- values[context$group_index]
        }
        metadata <- c(
          lp$metadata,
          list(distribution = "binomial", level = level, vars = vars, link = "logit"),
          .mlsim_parameter_metadata(level, context, size = size_n, prob = prob)
        )
        .mlsim_result(values, vars, metadata, internal = lp$internal)
      }

      .mlsim_spec(
        "binomial", vars, level, simulate,
        size = size,
        fixed_intercept = fixed,
        random_cov = random_cov
      )
    },
    poisson = {
      vars <- .mlsim_check_vars(vars, 1L)
      fixed <- check_fixed_scalar("Poisson")
      check_random_contract("poisson")

      simulate <- function(context) {
        lp <- scalar_linear_predictor(context, "poisson", fixed)
        lambda <- .mlsim_check_nonnegative(.mlsim_mean_from_eta("poisson", lp$eta), "lambda")
        values <- stats::rpois(length(lambda), lambda = lambda)
        if (identical(level, "level2")) {
          values <- values[context$group_index]
        }
        metadata <- c(
          lp$metadata,
          list(distribution = "poisson", level = level, vars = vars, link = "log"),
          .mlsim_parameter_metadata(level, context, lambda = lambda)
        )
        .mlsim_result(values, vars, metadata, internal = lp$internal)
      }

      .mlsim_spec(
        "poisson", vars, level, simulate,
        fixed_intercept = fixed,
        random_cov = random_cov
      )
    },
    negbin = {
      vars <- .mlsim_check_vars(vars, 1L)
      fixed <- check_fixed_scalar("Negative binomial")
      scale_fixed <- check_scale_scalar("Negative binomial")
      check_random_contract("negative binomial")

      simulate <- function(context) {
        lp <- scalar_linear_predictor(context, "negbin", fixed, scale_fixed = scale_fixed)
        mu <- .mlsim_check_positive(.mlsim_mean_from_eta("negbin", lp$eta), "mu")
        size_n <- .mlsim_check_positive(exp(lp$scale_eta), "size")
        values <- stats::rnbinom(length(mu), mu = mu, size = size_n)
        if (identical(level, "level2")) {
          values <- values[context$group_index]
        }
        metadata <- c(
          lp$metadata,
          list(distribution = "negative_binomial", level = level, vars = vars, link = "log"),
          .mlsim_parameter_metadata(level, context, size = size_n, mu = mu),
          list(size = size_n)
        )
        .mlsim_result(values, vars, metadata, internal = lp$internal)
      }

      .mlsim_spec(
        "negbin", vars, level, simulate,
        fixed_intercept = fixed,
        random_cov = random_cov,
        scale_fixed_intercept = scale_fixed
      )
    },
    gamma = {
      vars <- .mlsim_check_vars(vars, 1L)
      fixed <- check_fixed_scalar("Gamma")
      scale_fixed <- check_scale_scalar("Gamma")
      check_random_contract("gamma")

      simulate <- function(context) {
        lp <- scalar_linear_predictor(context, "gamma", fixed, scale_fixed = scale_fixed)
        mean <- .mlsim_check_positive(.mlsim_mean_from_eta("gamma", lp$eta), "mean")
        shape <- .mlsim_check_positive(exp(lp$scale_eta), "shape")
        values <- stats::rgamma(length(mean), shape = shape, scale = mean / shape)
        if (identical(level, "level2")) {
          values <- values[context$group_index]
        }
        metadata <- c(
          lp$metadata,
          list(distribution = "gamma", level = level, vars = vars, link = "log"),
          .mlsim_parameter_metadata(level, context, shape = shape, mean = mean),
          list(shape = shape)
        )
        .mlsim_result(values, vars, metadata, internal = lp$internal)
      }

      .mlsim_spec(
        "gamma", vars, level, simulate,
        fixed_intercept = fixed,
        random_cov = random_cov,
        scale_fixed_intercept = scale_fixed
      )
    },
    beta = {
      vars <- .mlsim_check_vars(vars, 1L)
      fixed <- check_fixed_scalar("Beta")
      scale_fixed <- check_scale_scalar("Beta")
      check_random_contract("beta")

      simulate <- function(context) {
        lp <- scalar_linear_predictor(context, "beta", fixed, scale_fixed = scale_fixed)
        mean <- .mlsim_check_open_prob(.mlsim_mean_from_eta("beta", lp$eta), "mean")
        precision <- .mlsim_check_positive(exp(lp$scale_eta), "precision")
        values <- stats::rbeta(
          length(mean),
          shape1 = mean * precision,
          shape2 = (1 - mean) * precision
        )
        if (identical(level, "level2")) {
          values <- values[context$group_index]
        }
        metadata <- c(
          lp$metadata,
          list(distribution = "beta", level = level, vars = vars, link = "logit"),
          .mlsim_parameter_metadata(level, context, mean = mean, precision = precision),
          list(precision = precision)
        )
        .mlsim_result(values, vars, metadata, internal = lp$internal)
      }

      .mlsim_spec(
        "beta", vars, level, simulate,
        fixed_intercept = fixed,
        random_cov = random_cov,
        scale_fixed_intercept = scale_fixed
      )
    },
    custom = {
      vars <- .mlsim_check_vars(vars, NULL)
      if (!is.function(generator)) {
        .mlsim_stop("`generator` must be a function.")
      }
      if (!is.list(metadata)) {
        .mlsim_stop("`metadata` must be a list.")
      }
      if (!is.list(user_args)) {
        .mlsim_stop("`user_args` must be a list.")
      }

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
  )
}

#' Generate Normal, Multivariate Normal, and Compositional Variables
#'
#' Create a generator specification for univariate or multivariate Gaussian
#' variables, optionally transformed from ILR coordinates to compositional parts.
#'
#' @param vars Character vector naming the generated variables. For
#'   compositional generators these are ILR coordinate names.
#' @param level Simulation level. `"single"` generates one row-level vector per
#'   observation, `"level2"` generates one group-level vector and expands it to
#'   each row in the group, and `"multilevel"` uses a random-intercept model.
#' @param fixed_intercept Identity-scale location intercept vector.
#' @param residual_cov Residual covariance matrix for Gaussian generators
#'   without a row-specific scale model. For univariate generators, it may be a
#'   scalar residual variance.
#' @param random_cov Group-level random-intercept covariance matrix for
#'   multilevel generators. When `scale_fixed_intercept` is supplied this may be
#'   either a location-only covariance or a joint location-scale covariance.
#' @param scale_fixed_intercept Optional log residual standard-deviation
#'   intercepts.
#' @param residual_cor Residual correlation matrix used with
#'   `scale_fixed_intercept`. Defaults to an identity correlation matrix.
#' @param ... Removed direct distribution parameters are rejected.
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
#' Multilevel generators additionally emit the latent between- and
#' within-group components of each generated variable as visible columns.
#' For non-compositional variables these are `<var>_between` (the group-level
#' mean, fixed plus random intercept) and `<var>_within` (the row-level
#' residual), with `<var> = <var>_between + <var>_within` exactly. For
#' compositional generators the same decomposition is emitted on the ILR
#' scale (`<ilr>_between`, `<ilr>_within`), together with the true between
#' composition `<part>_between`, the closed back-transform of the group-level
#' ILR means. All generated variables are labelled with column roles so that
#' `between()` and `within()` terms in [gen_outcome()] formulas resolve to the
#' latent components; roles (and hence `between()`/`within()`) apply to the
#' ILR coordinates of a compositional generator, not to its parts. When
#' `keep_ilr = FALSE` the observed ILR columns are dropped, but the
#' between/within columns and their roles are still emitted for multilevel
#' generators.
#'
#' Multilevel generators also record their group-level random-effect draws
#' as per-row truth columns: `.mlsim_<var>_random_intercept` for the
#' location intercept and, when `random_cov` is a joint location-scale
#' covariance, `.mlsim_<var>_scale_random_intercept` for the scale
#' intercept. The other multilevel scalar generators ([gen_negbin()],
#' [gen_gamma()], [gen_beta()], [gen_poisson()], [gen_binomial()]) share
#' this machinery and emit the same columns. See the Value section of
#' [simulate_data()] for the full naming contract.
#'
#' @examples
#' sim <- simulate_data(
#'   n = 4,
#'   seed = 2,
#'   generators = list(
#'     x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
#'     z = gen_mvn(
#'       c("z1", "z2"),
#'       fixed_intercept = c(0, 0),
#'       residual_cov = diag(2),
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
                    fixed_intercept = NULL, residual_cov = NULL,
                    random_cov = NULL, ..., scale_fixed_intercept = NULL,
                    residual_cor = NULL,
                    compositional = FALSE, parts = NULL, total = 1,
                    keep_ilr = TRUE, sbp = NULL) {
  .mlsim_check_unused_generator_dots(...)
  .mlsim_gen_distribution(
    distribution = "mvn",
    vars = vars,
    level = level,
    fixed_intercept = fixed_intercept,
    residual_cov = residual_cov,
    random_cov = random_cov,
    scale_fixed_intercept = scale_fixed_intercept,
    residual_cor = residual_cor,
    compositional = compositional,
    parts = parts,
    total = total,
    keep_ilr = keep_ilr,
    sbp = sbp
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
#' @details
#' Multilevel categorical generators record their group-level random-effect
#' draws as per-row truth columns named
#' `.mlsim_<var>_random_intercept_<category>`, one per non-reference
#' category, with `<var>` and `<category>` sanitized by [make.names()]. See
#' the Value section of [simulate_data()] for the full naming contract.
#'
#' @examples
#' sim <- simulate_data(
#'   n = 6,
#'   seed = 3,
#'   generators = list(
#'     arm = gen_categorical(
#'       "arm",
#'       categories = c("control", "treatment"),
#'       fixed_intercept = stats::qlogis(0.5)
#'     )
#'   )
#' )
#' sim$data
#'
#' @family predictor generators
#' @export
gen_categorical <- function(vars, level = c("single", "level2", "multilevel"),
                            categories = NULL, fixed_intercept = NULL,
                            random_cov = NULL, reference = NULL,
                            output = c("factor", "character", "integer"),
                            ordered = FALSE) {
  .mlsim_gen_distribution(
    distribution = "categorical",
    vars = vars,
    level = level,
    categories = categories,
    fixed_intercept = fixed_intercept,
    random_cov = random_cov,
    reference = reference,
    output = output,
    ordered = ordered
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
#' @param size Binomial trial size. May be a scalar, vector, function of `n`, or
#'   count-distribution list. Count-distribution `minimum`/`maximum` bounds
#'   clamp draws to the bounds (censoring, not truncation).
#' @param fixed_intercept Link-scale intercept. Binomial uses the logit link;
#'   Poisson and negative binomial use the log link.
#' @param random_cov Group-level random-intercept covariance for multilevel
#'   generators.
#' @param scale_fixed_intercept Log negative-binomial size intercept.
#' @param ... Removed direct distribution parameters are rejected.
#'
#' @return An `mlsim_generator_spec` object for use in [simulate_data()].
#'
#' @details
#' Count generators use link-scale intercepts at every level. Negative-binomial
#' size is `exp(scale_fixed_intercept)`.
#'
#' @examples
#' sim <- simulate_data(
#'   n = 5,
#'   seed = 4,
#'   generators = list(
#'     events = gen_poisson("events", fixed_intercept = log(2)),
#'     successes = gen_binomial(
#'       "successes",
#'       size = 4,
#'       fixed_intercept = stats::qlogis(0.5)
#'     ),
#'     overdispersed = gen_negbin(
#'       "overdispersed",
#'       fixed_intercept = log(2),
#'       scale_fixed_intercept = log(3)
#'     )
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
                         size = NULL, fixed_intercept = NULL,
                         random_cov = NULL) {
  .mlsim_gen_distribution(
    distribution = "binomial",
    vars = vars,
    level = level,
    size = size,
    fixed_intercept = fixed_intercept,
    random_cov = random_cov
  )
}

#' @rdname count-generators
#' @export
gen_poisson <- function(vars, level = c("single", "level2", "multilevel"),
                        fixed_intercept = NULL, random_cov = NULL) {
  .mlsim_gen_distribution(
    distribution = "poisson",
    vars = vars,
    level = level,
    fixed_intercept = fixed_intercept,
    random_cov = random_cov
  )
}

#' @rdname count-generators
#' @export
gen_negbin <- function(vars, level = c("single", "level2", "multilevel"),
                       fixed_intercept = NULL, ..., scale_fixed_intercept = NULL,
                       random_cov = NULL) {
  .mlsim_check_unused_generator_dots(...)
  .mlsim_gen_distribution(
    distribution = "negbin",
    vars = vars,
    level = level,
    fixed_intercept = fixed_intercept,
    scale_fixed_intercept = scale_fixed_intercept,
    random_cov = random_cov
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
#' @param fixed_intercept Link-scale intercept. Gamma uses the log link and
#'   beta uses the logit link.
#' @param scale_fixed_intercept Log shape intercept for gamma variables or log
#'   precision intercept for beta variables.
#' @param random_cov Group-level random-intercept covariance for multilevel
#'   generators.
#' @param ... Removed direct distribution parameters are rejected.
#'
#' @return An `mlsim_generator_spec` object for use in [simulate_data()].
#'
#' @examples
#' sim <- simulate_data(
#'   n = 5,
#'   seed = 5,
#'   generators = list(
#'     positive = gen_gamma(
#'       "positive",
#'       fixed_intercept = log(1.5),
#'       scale_fixed_intercept = log(2)
#'     ),
#'     proportion = gen_beta(
#'       "proportion",
#'       fixed_intercept = stats::qlogis(0.4),
#'       scale_fixed_intercept = log(10)
#'     )
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
                      fixed_intercept = NULL, ..., scale_fixed_intercept = NULL,
                      random_cov = NULL) {
  .mlsim_check_unused_generator_dots(...)
  .mlsim_gen_distribution(
    distribution = "gamma",
    vars = vars,
    level = level,
    fixed_intercept = fixed_intercept,
    scale_fixed_intercept = scale_fixed_intercept,
    random_cov = random_cov
  )
}

#' @rdname continuous-generators
#' @export
gen_beta <- function(vars, level = c("single", "level2", "multilevel"),
                     fixed_intercept = NULL, ..., scale_fixed_intercept = NULL,
                     random_cov = NULL) {
  .mlsim_check_unused_generator_dots(...)
  .mlsim_gen_distribution(
    distribution = "beta",
    vars = vars,
    level = level,
    fixed_intercept = fixed_intercept,
    scale_fixed_intercept = scale_fixed_intercept,
    random_cov = random_cov
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
  .mlsim_gen_distribution(
    distribution = "custom",
    vars = vars,
    level = level,
    generator = generator,
    user_args = list(...),
    metadata = metadata
  )
}
