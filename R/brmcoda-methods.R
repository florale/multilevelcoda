#' Checks if argument is a \code{brmcoda} object
#'
#' @param x An object of class \code{brmcoda}.
#'
#' @export
is.brmcoda <- function(x) {
  inherits(x, "brmcoda")
}


#' Extract Number of Observations from \code{brmcoda} object
#'
#' @param object A \code{brmcoda} object.
#' @param ... Further arguments to be passed to methods.
#'
#' @importFrom stats nobs
#' @method nobs brmcoda
#' @export
nobs.brmcoda <- function(object, ...) {
  nobs(object$model, ...)
}

#' Extracting the Model Frame from a Formula or Fit from \code{brmcoda} object
#'
#' @param formula A \code{brmcoda} object.
#' @param ... Further arguments to be passed to methods.
#'
#' @importFrom stats model.frame
#' @method model.frame brmcoda
#' @export
model.frame.brmcoda <- function(formula, ...) {
  model.frame(formula$model, ...)
}

#' Index \code{brmcoda} objects
#'
#' @aliases variables nvariables niterations nchains ndraws
#'
#' @param x An object of class \code{brmcoda}.
#' @param ... Arguments passed to individual methods.
#'
#' @name draws-index-brmcoda
#'
NULL

#' @rdname draws-index-brmcoda
#' @importFrom brms variables
#' @method variables brmcoda
#'
#' @seealso \code{\link[brms:variables]{variables.brmsfit}}
#'
#' @export
variables.brmcoda <- function(x, ...) {
  variables(x$model, ...)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms nvariables
#' @method nvariables brmcoda
#'
#' @seealso \code{\link[brms:nvariables.brmsfit]{nvariables.brmsfit}}
#'
#' @export
nvariables.brmcoda <- function(x, ...) {
  nvariables(x$model, ...)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms niterations
#' @method niterations brmcoda
#'
#' @seealso \code{\link[brms:niterations.brmsfit]{niterations.brmsfit}}
#'
#' @export
niterations.brmcoda <- function(x) {
  niterations(x$model)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms nchains
#' @method nchains brmcoda
#'
#' @seealso \code{\link[brms:nchains.brmsfit]{nchains.brmsfit}}
#'
#' @export
nchains.brmcoda <- function(x) {
  nchains(x$model)
}

#' @rdname draws-index-brmcoda
#' @importFrom brms ndraws
#' @method ndraws brmcoda
#'
#' @seealso \code{\link[brms:ndraws.brmsfit]{ndraws.brmsfit}}
#'
#' @export
ndraws.brmcoda <- function(x) {
  ndraws(x$model)
}

# nwarmup.brmcoda <- function(x) {
#   if (inherits(x$model, "brmcoda")) {
#     nwarmup(x$model)
#   }
# }
#
# nthin.brmcoda <- function(x) {
#   if (inherits(x$model, "brmcoda")) {
#     nthin(x$model)
#   }
# }

#' Extract Diagnostic Quantities from \code{brmsfit} Models in \code{brmcoda}
#'
#' @name diagnostic-quantities-brmcoda
#' @aliases log_posterior nuts_params rhat neff_ratio
#'
#' @param x,object A \code{brmcoda} object or another \R object for which
#' the methods are defined.
#' @param ... Arguments passed to individual methods (if applicable).
#'
#' @return The exact form of the output depends on the method.
#'
NULL

#' @rdname diagnostic-quantities-brmcoda
#' @importFrom brms log_posterior
#' @method log_posterior brmcoda
#'
#' @seealso \code{\link[brms:log_posterior.brmsfit]{log_posterior.brmsfit}}
#'
#' @export
log_posterior.brmcoda <- function(object, ...) {
  log_posterior(object$model, ...)
}

#' @rdname diagnostic-quantities-brmcoda
#' @importFrom brms nuts_params
#' @method nuts_params brmcoda
#'
#' @seealso \code{\link[brms:nuts_params.brmsfit]{nuts_params.brmsfit}}
#'
#' @export
nuts_params.brmcoda <- function(object, ...) {
  nuts_params(object$model, ...)
}

#' @rdname diagnostic-quantities-brmcoda
#' @importFrom brms rhat
#' @method rhat brmcoda
#'
#' @seealso \code{\link[brms:rhat.brmsfit]{rhat.brmsfit}}
#'
#' @export
rhat.brmcoda <- function(x, ...) {
  rhat(x$model, ...)
}

#' @rdname diagnostic-quantities-brmcoda
#' @importFrom brms neff_ratio
#' @method neff_ratio brmcoda
#'
#' @seealso \code{\link[brms:neff_ratio.brmsfit]{neff_ratio.brmsfit}}
#'
#' @export
neff_ratio.brmcoda <- function(object, ...) {
  neff_ratio(object$model, ...)
}

#' Bayes Factors from Marginal Likelihoods
#'
#' Compute Bayes factors from marginal likelihoods
#'
#' @param x1 A \code{brmcoda} object.
#' @param x2 Another \code{brmcoda} object based on the same responses.
#' @param ... Other arguments passed to \code{\link[brms:bayes_factor.brmsfit]{bayes_factor.brmsfit}}.
#'
#' @importFrom brms bayes_factor
#' @method bayes_factor brmcoda
#'
#' @seealso \code{\link[brms:bayes_factor.brmsfit]{bayes_factor.brmsfit}}
#'
#' @export
bayes_factor.brmcoda <- function(x1, x2, ...) {
  out <- invisible(bayes_factor(x1 = x1$model, x2 = x2$model, ...))
  
  m1 <- deparse(substitute(x1))
  m2 <- deparse(substitute(x2))
  
  attr(out, "model_names") <- c(m1, m2)
  out
}

#' Extract Priors of a \code{brmsfit} from a \code{brmcoda} object
#'
#' Compute Bayes factors from marginal likelihoods
#'
#' @param object An object of class \code{brmcoda}.
#' @inheritParams brms::prior_summary.brmsfit
#'
#' @importFrom brms prior_summary
#' @method prior_summary brmcoda
#'
#' @seealso \code{\link[brms:prior_summary.brmsfit]{prior_summary.brmsfit}}
#'
#' @export
prior_summary.brmcoda <- function(object, ...) {
  prior_summary(object$model, ...)
}

#' Posterior Predictive Checks for \code{brmcoda} Objects
#'
#' Perform posterior predictive checks with the help of the \pkg{bayesplot} package.
#'
#' @aliases pp_check
#'
#' @param object An object of class \code{brmcoda}.
#' @inheritParams brms::pp_check.brmsfit
#' @inheritParams predict.brmcoda
#' @param parts Optional names or indices of compositional response parts to
#' include when \code{scale = "response"}. If \code{NULL}, all parts are used.
#'
#' @importFrom bayesplot pp_check
#' @method pp_check brmcoda
#'
#' @seealso \code{\link[brms:pp_check.brmsfit]{pp_check.brmsfit}}
#'
#' @export
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   ## fit a model
#'   x <- complr(data = mcompd, sbp = sbp,
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                  idvar = "ID", total = 1440)
#'
#'   m1 <- brmcoda(complr = x,
#'                 formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                    wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'
#'   ## posterior predictive checks
#'   bayesplot::pp_check(m1, ndraws = 5)
#'
#'   ## fit a model with compositional outcome
#'   m2 <- brmcoda(complr = x,
#'                 formula = mvbind(z1_1, z2_1, z3_1, z4_1) ~ 
#'                           bz1_1 + bz2_1 + bz3_1 + bz4_1 + Female + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'
#'   ## posterior predictive checks for compositional outcome -- linear scale
#'   bayesplot::pp_check(m2, resp = "z11", ndraws = 10)
#'   ## posterior predictive checks for compositional outcome -- original response scale
#'   bayesplot::pp_check(m2, parts = "WAKE", scale = "response", ndraws = 10)
#' }}
pp_check.brmcoda <- function(object,
                             type = "dens_overlay",
                             ndraws = NULL,
                             prefix = c("ppc", "ppd"),
                             newdata = NULL,
                             draw_ids = NULL,
                             nsamples = NULL,
                             subset = NULL,
                             scale = c("linear", "response"),
                             parts = NULL,
                             ...) {
  scale <- match.arg(scale)
  prefix <- match.arg(prefix)

  is_mv <- inherits(object[["model"]][["formula"]], "mvbrmsformula")
  brmcoda_vars <- if (is_mv) get_variables(object) else NULL
  is_comp_y <- is_mv && grepl("^multivariate-", brmcoda_vars$resp_type) &&
    !identical(brmcoda_vars$resp_type, "multivariate-non-compositional")

  if (!is_comp_y || identical(scale, "linear")) {
    return(
      pp_check(
        object$model, type = type, ndraws = ndraws, prefix = prefix,
        newdata = newdata, draw_ids = draw_ids, nsamples = nsamples,
        subset = subset, ...
      )
    )
  }

  complr_vars <- get_variables(object$complr)
  idy <- which(
    vapply(complr_vars, function(y) {
      any(sapply(c("Z", "bZ", "wZ"), function(z) {
        identical(sort(y[[z]]), sort(brmcoda_vars$y))
      }))
    }, logical(1))
  )

  if (length(idy) != 1L) {
    stop(
      "The response variables in the brmcoda model do not correspond to any of the compositional parts in the complr object."
    )
  }

  if (identical(brmcoda_vars$y, complr_vars[[idy]]$bZ)) {
    Yn <- complr_vars[[idy]]$bX
  } else if (identical(brmcoda_vars$y, complr_vars[[idy]]$wZ)) {
    Yn <- complr_vars[[idy]]$wX
  } else if (identical(brmcoda_vars$y, complr_vars[[idy]]$Z)) {
    Yn <- complr_vars[[idy]]$X
  } else {
    stop(
      "The response variables in the brmcoda model do not correspond to any of the compositional parts in the complr object."
    )
  }

  data <- if (is.null(newdata)) object$complr$dataout else newdata
  data <- as.data.frame(data)
  missing_y <- setdiff(Yn, names(data))
  if (length(missing_y)) {
    stop(
      "Response-scale posterior predictive checks require the original compositional outcome variables in 'newdata'."
    )
  }

  ndraws <- if (is.null(ndraws)) nsamples else ndraws
  draw_ids <- if (is.null(draw_ids)) subset else draw_ids
  if (is.null(ndraws) && is.null(draw_ids)) {
    ndraws <- 10
  }

  dots <- list(...)
  pred_arg_names <- c(
    "re_formula", "allow_new_levels", "sample_new_levels", "incl_autocor",
    "oos", "resp", "negative_rt", "sort", "ntrys", "cores", "newdata2"
  )
  pred_args <- dots[names(dots) %in% pred_arg_names]
  ppc_args <- dots[!names(dots) %in% pred_arg_names]

  yrep <- do.call(
    predict,
    c(
      list(
        object = object, scale = "response", newdata = newdata,
        ndraws = ndraws, draw_ids = draw_ids, summary = FALSE
      ),
      pred_args
    )
  )

  y <- as.matrix(data[Yn])
  yrep_names <- dimnames(yrep)[[3]]
  part_ids <- seq_along(Yn)
  if (!is.null(parts)) {
    part_ids <- if (is.numeric(parts)) {
      parts
    } else {
      match(parts, Yn)
    }
    if (anyNA(part_ids) && !is.null(yrep_names)) {
      missing <- is.na(part_ids)
      part_ids[missing] <- match(parts[missing], yrep_names)
    }
    if (anyNA(part_ids) || any(part_ids < 1L) || any(part_ids > length(Yn))) {
      stop("'parts' must identify compositional response parts.")
    }
  }

  y <- y[, part_ids, drop = FALSE]
  yrep <- yrep[, , part_ids, drop = FALSE]

  OK <- stats::complete.cases(y)
  if (any(!OK)) {
    warning("NA responses are not shown in 'pp_check'.")
    y <- y[OK, , drop = FALSE]
    yrep <- yrep[, OK, , drop = FALSE]
  }

  y <- as.vector(y)
  yrep <- t(matrix(aperm(yrep, c(2, 3, 1)), ncol = dim(yrep)[1]))

  if (identical(prefix, "ppd")) {
    ppd_fun <- get(paste0("ppd_", type), asNamespace("bayesplot"))
    do.call(ppd_fun, c(list(ypred = yrep), ppc_args))
  } else {
    ppc_fun <- get(paste0("ppc_", type), asNamespace("bayesplot"))
    do.call(pp_check, c(list(object = y, yrep = yrep, fun = ppc_fun), ppc_args))
  }
}

#' #' Posteriors Sampling Diagnostic
#' #'
#' #' Extract diagnostic metrics (Effective Sample Size (`ESS`), `Rhat` and Monte
#' #' Carlo Standard Error `MCSE`).
#' #'
#' #' @param posteriors An object of class \code{brmcoda}.
#' #' @inheritParams bayestestR::diagnostic_posterior
#' #' @param ... Other arguments passed to \code{\link{diagnostic_posterior}}.
#' #'
#' #' @importFrom bayestestR diagnostic_posterior
#' #' @method diagnostic_posterior brmcoda
#' #'
#' #' @seealso \code{\link[bayestestR:diagnostic_posterior]{diagnostic_posterior}}
#' #'
#' #' @export
#' diagnostic_posterior.brmcoda <- function(posteriors, diagnostic = c("ESS", "Rhat"), ...) {
#'   diagnostic_posterior(posterior$model, diagnostic = diagnostic, ...)
#' }
