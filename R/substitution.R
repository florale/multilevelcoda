#' Multilevel Compositional Substitution Analysis
#'
#' Estimate the difference in an outcome
#' when compositional parts are substituted for specific unit(s).
#' The \code{substitution} output encapsulates
#' the substitution results for all compositional parts
#' present in the \code{\link{brmcoda}} object.
#'
#' @param object A fitted \code{\link{brmcoda}} object.
#' @param delta A integer, numeric value or vector indicating the amount of substituted change between compositional parts.
#' @param ref Either a character value or vector or a dataset.
#' Can be \code{"grandmean"} and/or \code{"clustermean"}, or
#' a \code{data.frame} or \code{data.table} of user's specified reference grid consisting
#' of combinations of covariates over which predictions are made.
#' User's specified reference grid is only possible for simple substitution.
#' Single level models are default to \code{"grandmean"}.
#' @param level A character string or vector.
#' Should the estimate of multilevel models focus on the \code{"between"} and/or \code{"within"} or \code{"aggregate"} variance?
#' Single-level models are default to \code{"aggregate"}.
#' @param summary A logical value to obtain summary statistics instead of the raw values. Default is \code{TRUE}.
#' Currently only support outputing raw values for model using grandmean as reference composition.
#' @param at An optional named list of levels for the corresponding variables in the reference grid.
#' @param base An optional base substitution.
#' Can be a \code{data.frame} or \code{data.table} of the base possible substitution of compositional parts,
#' which can be computed using function \code{\link{build.base}}.
#' @param parts A optional character string specifying names of compositional parts that should be considered
#' in the substitution analysis. This should correspond to a single set of names of compositional parts specified
#' in the \code{complr} object. Default to the first composition in the \code{complr} object.
#' @param weight A character value specifying the weight to use in calculation of the reference composition.
#' @param type A character string to indicate the type of substitution.
#' If \code{"one-to-all"}, all possible one-to-remaining reallocations are estimated.
#' If \code{"one-to-one"}, all possible one-to-one reallocations are estimated.
#' If \code{"equal"}, give equal weight to units (e.g., individuals).
#' If \code{"proportional"}, weights in proportion to the frequencies of units being averaged
#' (e.g., observations across individuals).
#' Default to \code{"equal"} for \code{ref = "grandmean"} and \code{"proportional"} for \code{ref = "clustermean"}.
#' @param scale Either \code{"response"} or \code{"linear"}.
#' If \code{"response"}, results are returned on the scale of the response variable.
#' If \code{"linear"}, results are returned on the scale of the linear predictor term,
#' that is without applying the inverse link function or other transformations.
#' @param aorg Internal use. A logical value indicating whether the results should be average across reference grid.
#' @param cores Number of cores to use when executing the chains in parallel,
#' we recommend setting the \code{mc.cores} option
#' to be as many processors as the hardware and RAM allow (up to the number of compositional parts).
#' For non-Windows OS in non-interactive R sessions, forking is used instead of PSOCK clusters.
#' Default to \code{"one-to-one"}.
#' @param ... Further arguments passed to \code{\link[brms:posterior_summary]{posterior_summary}}.
#'
#' @return A list containing the results of multilevel compositional substitution model.
#' The first six lists contain the results of the substitution estimation for a compositional part.
#'   \item{\code{Mean}}{ Posterior means.}
#'   \item{\code{CI_low} and \code{CI_high}}{ 95% credible intervals.}
#'   \item{\code{Delta}}{ Amount substituted across compositional parts.}
#'   \item{\code{From}}{ Compositional part that is substituted from.}
#'   \item{\code{To}}{ Compositional parts that is substituted to.}
#'   \item{\code{Level}}{ Level where changes in composition takes place.}
#'   \item{\code{Reference}}{ Either \code{grandmean}, \code{clustermean}, or \code{users}.}
#'
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo mean.acomp
#'
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   x <- complr(data = mcompd, sbp = sbp,
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                  idvar = "ID", total = 1440)
#'
#'   # model with compositional predictor at between and within-person levels
#'   m1 <- brmcoda(complr = x,
#'                   formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                      wz1_1 + wz2_1 + wz3_1 + wz4_1 +
#'                                      Female + (1 | ID),
#'                   chain = 1, iter = 500, backend = "cmdstanr")
#'
#'   # one to one reallocation at between and within-person levels
#'   sub1 <- substitution(object = m1, delta = 5, level = c("between"))
#'   summary(sub1)
#'
#'   # one to all reallocation at between and within-person levels
#'   sub2 <- substitution(object = m1, delta = 5, level = c("between", "within"),
#'                        type = "one-to-all")
#'   summary(sub2)
#'
#'   # model with compositional predictor at aggregate level
#'   m2 <- brmcoda(complr = x,
#'                   formula = Stress ~ z1_1 + z2_1 + z3_1 + z4_1 + (1 | ID),
#'                   chain = 1, iter = 500, backend = "cmdstanr")
#'   sub3 <- substitution(object = m2, delta = 5, level = c("aggregate"))
#'
#' }}
#' @export
substitution <- function(object,
                         delta,
                         ref = c("grandmean", "clustermean"),
                         level = c("between", "within", "aggregate"),
                         summary = TRUE,
                         at = NULL,
                         parts = 1,
                         base,
                         type,
                         weight = c("equal", "proportional"),
                         scale = c("response", "linear"),
                         aorg = NULL,
                         cores = NULL,
                         ...) {
  if (missing(object)) {
    stop(
      paste(
        "'object' is a required argument and cannot be missing;",
        "  it should be an object of class 'brmcoda'.",
        "  See ?substitution for details.",
        sep = "\n"
      )
    )
  }
  
  if (isFALSE(inherits(object, "brmcoda"))) {
    stop(
      sprintf(
        "Can't handle an object of class (%s)
  It should be a fitted 'brmcoda' object
  See ?substitution for details.",
        class(object)
      )
    )
  }
  
  if (isFALSE(identical(object[["complr"]][["transform"]], "ilr"))) {
    stop(
      sprintf(
        "Can't handle an object of class (%s) in 'substitution',
      'brmcoda' should be fitted with ilr transform to enable substitution analysis.",
        object[["complr"]][["transform"]]
      )
    )
  }
  
  if (isFALSE(missing(delta))) {
    if (isFALSE(is.integer(delta))) {
      if (isFALSE(delta > 0)) {
        stop(" 'delta' should be an positive integer value.")
      }
    }
  } else if (missing(delta)) {
    stop(
      paste(
        "'delta' is a required argument and cannot be missing;",
        "  it should be interger, numeric positive value or vector",
        "  to specify the change in units across compositional parts",
        sep = "\n"
      )
    )
  }
  
  # set default weight to be equal
  if (identical(weight, "proportional")) {
    weight <- "proportional"
  } else {
    weight <- "equal"
  }
  
  # set default ref to be grandmean
  if (identical(ref, "clustermean")) {
    ref <- "clustermean"
  } else {
    ref <- "grandmean"
  }
  
  # default aorg to TRUE when at is null
  if (is.null(at)) {
    aorg <- TRUE
  } else {
    aorg <- FALSE
  }
  
  # get part names
  parts <- .get_parts(object[["complr"]], parts)
  
  # get the index of which index elements of object[["complr"]][["output"]] does the parts correspond to
  idx <- as.integer(which(vapply(lapply(object[["complr"]][["output"]], function(x)
    x$parts), function(p)
      identical(sort(parts), sort(p)), logical(1))))[1]
  
  # get brmcoda variables
  brmcoda_vars <- get_variables(object)
  
  # type
  if (missing(type)) {
    if (missing(base)) {
      base <- build.base(parts = object[["complr"]][["output"]][[idx]][["parts"]])
      type <- "one-to-one"
    } else {
      if (inherits(base, c("data.frame", "data.table", "matrix"))) {
        type <- "one-to-one"
      } else {
        stop("If 'base' is provided, it should be a data frame or data table.")
      }
      # names in base should match parts
      if (isFALSE(identical(sort(parts), sort(colnames(base))))) {
        stop(
          sprintf(
            "'base' should contain the same compositional parts as specified in 'parts' argument: %s",
            paste(parts, collapse = ", ")
          )
        )
      }
    }
  } else {
    if (identical(type, "one-to-one")) {
      base <- build.base(parts = object[["complr"]][["output"]][[idx]][["parts"]])
      type <- "one-to-one"
    } else if (inherits(type, "character") &&
               identical(type, "one-to-all")) {
      base <- build.base(parts = object[["complr"]][["output"]][[idx]][["parts"]], type = "one-to-all")
      type <- "one-to-all"
    }
  }
  
  if (is.null(brmcoda_vars[["fixef_type"]])) {
    stop("No fixed effects of composition in the model to perform substitution analysis.")
  }
  
  ## only grandmean and aggregate for single level model
  if ("aggregate" %in% brmcoda_vars[["fixef_type"]]) {
    if (identical(brmcoda_vars[["ranef_type"]], "single")) {
      if ("clustermean" %in% ref) {
        warning("Can only use grandmean for single level model.")
      }
      level  <- "aggregate"
      ref    <- "grandmean"
      weight <- "equal"
    } else {
      level  <- "aggregate"
      ref    <- ref
    }
  }
  
  ## no between or within for single level model
  if (any(c("between", "within") %in% brmcoda_vars[["fixef_type"]])) {
    if (identical(brmcoda_vars[["ranef_type"]], "single")) {
      stop(
        " between and within substitution analysis cannot be computed on a single level model"
      )
    } else if (identical(brmcoda_vars[["ranef_type"]], "multilevel")) {
      level <- level
      ref <- ref
    }
  }
  
  ## set default to be only between and within if level is not specified
  if (all(c("between", "within", "aggregate") %in% level)) {
    level <- c("between", "within")
  }
  
  ## level args match with coefs in object
  if (any(c("between", "within") %in% level)) {
    if (isFALSE(any(c("between", "within") %in% brmcoda_vars[["fixef_type"]]))) {
      stop(
        sprintf(
          "between and within substitution analysis cannot be computed
  on a model estimated using the (%s) variance of ilrs.
  Please specify the level argument as \"(%s)\" instead or refit 'brmcoda' model.",
          brmcoda_vars[["fixef_type"]],
          brmcoda_vars[["fixef_type"]]
        )
      )
    }
    if ("between" %in% level &&
        isFALSE(all(names(object[["complr"]][["output"]][[idx]][["bZ"]] %in% brmcoda_vars[["x"]])))) {
      stop(
        sprintf(
          "brmcoda model %s should include a complete set of between and within ilr predictors %s",
          "for substitution analysis at between level to be performed",
          brmcoda_vars[["x"]],
          names(object[["complr"]][["output"]][[idx]][["bZ"]])
        )
      )
    }
    
    if ("within" %in% level &&
        isFALSE(all(names(object[["complr"]][["output"]][[idx]][["wZ"]] %in% brmcoda_vars[["x"]])))) {
      stop(
        sprintf(
          "brmcoda model %s should include a complete set of between and within ilr predictors %s",
          "for substitution analysis at within level to be performed",
          brmcoda_vars[["x"]],
          names(object[["complr"]][["output"]][[idx]][["wZ"]])
        )
      )
    }
  }
  if ("aggregate" %in% level) {
    if (isFALSE("aggregate" %in% brmcoda_vars[["fixef_type"]])) {
      stop(
        sprintf(
          "'aggregate' substitution analysis cannot be computed
  on a model estimated using the (%s) variance of ilrs.
  Please specify the level argument as \"(%s)\" instead or refit 'brmcoda' model.",
          brmcoda_vars[["fixef_type"]],
          brmcoda_vars[["fixef_type"]]
        )
      )
    }
    if (isFALSE(all(names(object[["complr"]][["output"]][[idx]][["Z"]]) %in% brmcoda_vars[["x"]]))) {
      stop(
        sprintf(
          "brmcoda model %s should include a complete set of aggregate ilr predictors %s",
          "for substitution analysis at aggregate level to be performed",
          brmcoda_vars[["x"]],
          names(object[["complr"]][["output"]][[idx]][["Z"]])
        )
      )
    }
  }
  
  # deploy to substitution analysis
  bmout <- bout <- NULL
  if ("between" %in% level) {
    if ("grandmean" %in% ref) {
      bout <- bsub(
        object = object,
        delta = delta,
        base = base,
        parts = parts,
        summary = summary,
        aorg = aorg,
        at = at,
        ref = "grandmean",
        level = "between",
        weight = weight,
        scale = scale,
        type = type,
        cores = cores,
        ...
      )
    }
    else if (inherits(ref, c("data.table", "data.frame", "matrix"))) {
      bout <- bsub(
        object = object,
        delta = delta,
        base = base,
        parts = parts,
        summary = summary,
        aorg = aorg,
        at = at,
        ref = ref,
        level = "between",
        weight = weight,
        scale = scale,
        type = type,
        cores = cores,
        ...
      )
    }
    if ("clustermean" %in% ref) {
      bmout <-
        bsubmargin(
          object = object,
          delta = delta,
          base = base,
          parts = parts,
          summary = summary,
          ref = "clustermean",
          level = "between",
          weight = weight,
          scale = scale,
          type = type,
          cores = cores,
          ...
        )
    }
  }
  
  wmout <- wout <- NULL
  if ("within" %in% level) {
    if ("grandmean" %in% ref) {
      wout <- wsub(
        object = object,
        delta = delta,
        base = base,
        parts = parts,
        summary = summary,
        aorg = aorg,
        at = at,
        ref = "grandmean",
        level = "within",
        weight = weight,
        scale = scale,
        type = type,
        cores = cores,
        ...
      )
    }
    else if (inherits(ref, c("data.table", "data.frame", "matrix"))) {
      wout <- wsub(
        object = object,
        delta = delta,
        base = base,
        parts = parts,
        summary = summary,
        aorg = aorg,
        at = at,
        ref = ref,
        level = "within",
        weight = weight,
        scale = scale,
        type = type,
        cores = cores,
        ...
      )
    }
    if ("clustermean" %in% ref) {
      wmout <-
        wsubmargin(
          object = object,
          delta = delta,
          base = base,
          parts = parts,
          summary = summary,
          ref = "clustermean",
          level = "within",
          weight = weight,
          scale = scale,
          type = type,
          cores = cores,
          ...
        )
    }
  }
  
  tmout <- tout <- NULL
  if ("aggregate" %in% level) {
    if ("grandmean" %in% ref) {
      tout <- sub(
        object = object,
        delta = delta,
        base = base,
        parts = parts,
        summary = summary,
        aorg = aorg,
        at = at,
        ref = "grandmean",
        level = "aggregate",
        weight = weight,
        scale = scale,
        type = type,
        cores = cores,
        ...
      )
    }
    else if (inherits(ref, c("data.table", "data.frame", "matrix"))) {
      tout <- sub(
        object = object,
        delta = delta,
        base = base,
        parts = parts,
        summary = summary,
        aorg = aorg,
        at = at,
        ref = ref,
        level = "aggregate",
        weight = weight,
        scale = scale,
        type = type,
        cores = cores,
        ...
      )
    }
    if ("clustermean" %in% ref) {
      tmout <-
        submargin(
          object = object,
          delta = delta,
          base = base,
          parts = parts,
          summary = summary,
          ref = "clustermean",
          level = "aggregate",
          weight = weight,
          scale = scale,
          type = type,
          cores = cores,
          ...
        )
    }
  }
  
  # out
  structure(
    list(
      between_simple_sub = if (exists("bout")) bout  else NULL,
      within_simple_sub = if (exists("wout")) wout else NULL,
      simple_sub = if (exists("tout")) tout else NULL,
      between_avg_sub = if (exists("bmout")) bmout else NULL,
      within_avg_sub = if (exists("wmout")) wmout else NULL,
      avg_sub = if (exists("tmout")) tmout else NULL,
      brmsformula = object$model$formula,
      delta = delta,
      ref = ref,
      level = level,
      weight = weight,
      parts = parts,
      at = at,
      summary = summary,
      type = if (exists("type")) (type) else (NULL)
    ),
    class = "substitution"
  )
  
}