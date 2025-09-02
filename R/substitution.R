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
#' @param aorg A logical value to obtain (a)verage prediction (o)ver the (r)eference (g)rid.
#' Should the estimate at each level of the reference grid (\code{FALSE})
#' or their average (\code{TRUE}) be returned?
#' Default is \code{TRUE}.
#' Only applicable for model with covariates in addition to
#' the isometric log-ratio coordinates (i.e., adjusted model).
#' @param at An optional named list of levels for the corresponding variables in the reference grid.
#' @param type A character string to indicate the type of substitution to be made.
#' If \code{"one-to-all"}, all possible one-to-remaining reallocations are estimated.
#' If \code{"one-to-one"}, all possible one-to-one reallocations are estimated. 
#' @param base An optional base substitution. 
#' Can be a \code{data.frame} or \code{data.table} of the base possible substitution of compositional parts,
#' which can be computed using function \code{\link{build.base}}.
#' @param parts A optional character string specifying names of compositional parts that should be considered
#' in the substitution analysis. This should correspond to a single set of names of compositional parts specified
#' in the \code{complr} object.
#' @param weight A character value specifying the weight to use in calculation of the reference composition.
#' If \code{"equal"}, give equal weight to units (e.g., individuals).
#' If \code{"proportional"}, weights in proportion to the frequencies of units being averaged
#' (e.g., observations across individuals).
#' Default to \code{"equal"} for \code{ref = "grandmean"} and \code{"proportional"} for \code{ref = "clustermean"}.
#' @param scale Either \code{"response"} or \code{"linear"}.
#' If \code{"response"}, results are returned on the scale of the response variable.
#' If \code{"linear"}, results are returned on the scale of the linear predictor term,
#' that is without applying the inverse link function or other transformations.
#' @param cores Number of cores to use when executing the chains in parallel,
#' we recommend setting the \code{mc.cores} option
#' to be as many processors as the hardware and RAM allow (up to the number of compositional parts).
#' For non-Windows OS in non-interactive R sessions, forking is used instead of PSOCK clusters.
#' Default to \code{"one-to-one"}.
#' @param ... currently ignored.
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
#'   fit1 <- brmcoda(complr = x,
#'                   formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                      wz1_1 + wz2_1 + wz3_1 + wz4_1 + Female + (1 | ID),
#'                   chain = 1, iter = 500, backend = "cmdstanr")
#'                   
#'   # one to one reallocation at between and within-person levels
#'   sub1 <- substitution(object = fit1, delta = 5, level = c("between"))
#'   summary(sub1)
#'   
#'   # one to all reallocation at between and within-person levels
#'   sub2 <- substitution(object = fit1, delta = 5, level = c("between", "within"), 
#'                        type = "one-to-all")
#'   summary(sub2) 
#'   
#'   # model with compositional predictor at aggregate level of variance
#'   fit2 <- brmcoda(complr = x,
#'                   formula = Stress ~ z1 + z2 + z3 + z4 + (1 | ID),
#'                   chain = 1, iter = 500, backend = "cmdstanr")
#'   sub3 <- substitution(object = fit2, delta = 5, level = c("aggregate"))
#'
#' }}
#' @export
substitution <- function(object,
                         delta,
                         ref = c("grandmean", "clustermean"),
                         level = c("between", "within", "aggregate"),
                         summary = TRUE,
                         aorg = TRUE,
                         at = NULL,
                         parts = 1,
                         base,
                         weight = c("equal", "proportional"),
                         scale = c("response", "linear"),
                         cores = NULL,
                         type,
                         ...) {
  
  if (missing(object)) {
    stop(paste(
      "'object' is a required argument and cannot be missing;",
      "  it should be an object of class 'brmcoda'.",
      "  See ?bsub for details.",
      sep = "\n"))
  }
  
  if (isFALSE(inherits(object, "brmcoda"))) {
    stop(sprintf(
      "Can't handle an object of class (%s)
  It should be a fitted 'brmcoda' object
  See ?bsub for details.",
      class(object)))
  }
  
  if (isFALSE(identical(object$complr$transform, "ilr"))) {
    stop(sprintf(
      "Can't handle an object of class (%s) in 'substitution',
      'brmcoda' should be fitted with ilr transform to enable substitution analysis.",
      object$complr$transform)
    )
  }
  
  if(isFALSE(missing(delta))) {
    if (isFALSE(is.integer(delta))) {
      if (isFALSE(delta > 0)) {
        stop(" 'delta' should be an positive integer value.")
      }
    }
  } else if (missing(delta)){
    stop(paste(
      "'delta' is a required argument and cannot be missing;",
      "  it should be interger, numeric positive value or vector",
      "  to specify the change in units across compositional parts",
      sep = "\n"))
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
  
  # currently not allowed exporting draws for cluster mean as ref
  if(identical(ref, "clustermean")) {
    if (isFALSE(summary)) {
      stop("Currently can't output raw values when 'clustermean' is used as reference composition.")
    }
  }
  
  # check parts
  if (is.numeric(parts)) {
    if (length(parts) > 1) {
      stop(" 'parts' should be a single numeric value indicating which set of compositional parts to use.")
    }
    if (parts < 1 || parts > length(object$complr$output)) {
      stop(sprintf(
        " 'parts' should be a single numeric value between 1 and %s, corresponding to the number of sets of compositional parts in the 'complr' object.",
        length(object$complr$output)))
    }
    parts <- object$complr$output[[parts]]$parts
    
  } else {
    if (isFALSE(inherits(parts, "character"))) {
      stop(" 'parts' should be a character vector of compositional parts.")
    }
    ## parts should be identical with either one of the parts presented in output of complr
    if (isFALSE((any(vapply(lapply(object$complr$output, function(x) x$parts), function(p) identical(sort(parts), sort(p)), logical(1)))))) {
      stop(sprintf(
        "The specified 'parts' (%s) are not found in the complr object.",
        "  It should corespond to one set of compositional parts, either one of the following:",
        "%s",
        paste(parts, collapse = ", "),
        invisible(lapply(object$complr$output, function(x) cat(paste(x$parts, collapse = ", "), "\n"))),
        sep = "\n"))
    }
  }
  
  ## get the index of which index elements of object$complr$output does the parts correspond to - check w JW
  idx <- which(vapply(lapply(object$complr$output, 
                             function(x) x$parts), 
                      function(p) identical(sort(parts), sort(p)), logical(1)))
  
  # type - check with JW what would be the best way to detect onetoone vs onetoall
  if (missing(type)) {
    if (missing(base)) {
      base <- build.base(parts = object$complr$output[[idx]]$parts)
      names(base) <- object$complr$output[[idx]]$parts
      type <- "one-to-one"
    } else {
      if (inherits(base, c("data.frame", "data.table", "matrix"))) {
        type <- "one-to-one"
      } else {
        stop("If 'base' is provided, it should be a data frame or data table.")
      }
    }
  } else {
    if (identical(type, "one-to-one")) {
      base <- build.base(parts = object$complr$output[[idx]]$parts)
      names(base) <- object$complr$output[[idx]]$parts
      type <- "one-to-one"
    } else if (inherits(type, "character") && identical(type, "one-to-all")) {
      base <- build.base(parts = object$complr$output[[idx]]$parts, type = "one-to-all")
      names(base) <- object$complr$output[[idx]]$parts
      type <- "one-to-all"
    }
  }
  
  # what type of model is being estimated
  model_fixef <- rownames(fixef(object))
  model_ranef <- if(dim(object$model$ranef)[1] > 0) (names(ranef(object))) else (NULL)
  
  model_fixef_level <- model_fixef_coef <- NULL
  
  # grab the correct logratio names
  z_vars  <- get_variables(object$complr)[["logratio", paste0("composition_", idx)]]
  bz_vars <- get_variables(object$complr)[["between_logratio", paste0("composition_", idx)]]
  wz_vars <- get_variables(object$complr)[["within_logratio", paste0("composition_", idx)]]
  
  if (length(grep(paste0(bz_vars, collapse = "|"), model_fixef, value = T)) > 0) {
    model_fixef_level <- append(model_fixef_level, "between")
    model_fixef_coef  <- append(model_fixef_coef,
                                grep(paste0(bz_vars, collapse = "|"), model_fixef, value = T))
  }
  if (length(grep(paste0(wz_vars, collapse = "|"), model_fixef, value = T)) > 0) {
    model_fixef_level <- append(model_fixef_level, "within")
    model_fixef_coef  <- append(model_fixef_coef,
                                grep(paste0(wz_vars, collapse = "|"), model_fixef, value = T))
  }
  if ((length(grep(paste0(z_vars, collapse = "|"), model_fixef, value = T)) > 0) && (length(grep(paste0(c(bz_vars, wz_vars), collapse = "|"), model_fixef, value = T)) == 0)) {
    model_fixef_level <- append(model_fixef_level, "aggregate")
    model_fixef_coef  <- append(model_fixef_coef, setdiff(
      grep(paste0(z_vars, collapse = "|"), model_fixef, value = TRUE),
      grep(paste0(c(bz_vars, wz_vars), collapse = "|"), model_fixef, value = TRUE)
    ))
  }
  z_vars <- bz_vars <- wz_vars <- NULL
  
  # single level or multilevel
  if (length(model_ranef) > 0) {
    model_ranef_level <- "multilevel"
    model_ranef_coef <- model_ranef
  } else {
    model_ranef_level <- "single"
    model_ranef_coef <- NULL
  }
  
  # ensure sub args make sense
  ## only grandmean and aggregate for single level model
  if ("aggregate" %in% model_fixef_level) {
    if (model_ranef_level == "single") {
      if ("clustermean" %in% ref) {
        warning("Can only use grandmean for single level model.")
      }
      level <- "aggregate"
      ref <- "grandmean"
      weight <- "equal"
    } else {
      level <- "aggregate"
      ref <- ref
    }
  }
  
  ## no between or within for single level model
  if (any(c("between", "within") %in% model_fixef_level)) {
    if (model_ranef_level == "single") {
      stop(" 'between' and 'within' substitution analysis cannot be computed on a single level model")
    } else if (model_ranef_level == "multilevel") {
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
    if (isFALSE(any(c("between", "within") %in% model_fixef_level))) {
      stop(sprintf(
        "'between' and 'within' substitution analysis cannot be computed
  on a model estimated using the (%s) variance of ilrs.
  Please specify the level argument as \"(%s)\" instead or refit 'brmcoda' model.",
        model_fixef_level,
        model_fixef_level
      ))
    }
  }
  if ("aggregate" %in% level) {
    if (isFALSE("aggregate" %in% model_fixef_level)) {
      stop(sprintf(
        "'aggregate' substitution analysis cannot be computed
  on a model estimated using the (%s) variance of ilrs.
  Please specify the level argument as \"(%s)\" instead or refit 'brmcoda' model.",
        model_fixef_level,
        model_fixef_level
      ))
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
        aorg = aorg,
        summary = summary,
        ref = "grandmean",
        level = "between",
        weight = weight,
        scale = scale,
        type = type,
        cores = cores,
        ...)
    }
    else if (inherits(ref, c("data.table", "data.frame", "matrix"))) {
      bout <- bsub(
        object = object,
        delta = delta,
        base = base,
        parts = parts,
        aorg = aorg,
        summary = summary,
        ref = ref,
        level = "between",
        weight = weight,
        scale = scale,
        type = type,
        cores = cores,
        ...)
    }
    if ("clustermean" %in% ref) {
      bmout <-
        bsubmargins(
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
          ...)
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
        aorg = aorg,
        summary = summary,
        ref = "grandmean",
        level = "within",
        weight = weight,
        scale = scale,
        type = type,
        cores = cores,
        ...)
    }
    else if (inherits(ref, c("data.table", "data.frame", "matrix"))) {
      wout <- wsub(
        object = object,
        delta = delta,
        base = base,
        parts = parts,
        aorg = aorg,
        summary = summary,
        ref = ref,
        level = "within",
        weight = weight,
        scale = scale,
        type = type,
        cores = cores,
        ...)
    }
    if ("clustermean" %in% ref) {
      wmout <-
        wsubmargins(
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
          ...)
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
        aorg = aorg,
        summary = summary,
        ref = "grandmean",
        level = "aggregate",
        weight = weight,
        scale = scale,
        type = type,
        cores = cores,
        ...)
    }
    else if (inherits(ref, c("data.table", "data.frame", "matrix"))) {
      tout <- sub(
        object = object,
        delta = delta,
        base = base,
        parts = parts,
        aorg = aorg,
        summary = summary,
        ref = ref,
        level = "aggregate",
        weight = weight,
        scale = scale,
        type = type,
        cores = cores,
        ...)
    }
    if ("clustermean" %in% ref) {
      tmout <-
        submargins(
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
          ...)
    }
  }
  
  # out
  structure(
    list(
      between_simple_sub = if(exists("bout")) (bout) else (NULL),
      within_simple_sub = if(exists("wout")) (wout) else (NULL),
      simple_sub = if(exists("tout")) (tout) else (NULL),
      between_avg_sub = if(exists("bmout")) (bmout) else (NULL),
      within_avg_sub = if(exists("wmout")) (wmout) else (NULL),
      avg_sub = if(exists("tmout")) (tmout) else (NULL),
      brmsformula = object$model$formula,
      delta = delta,
      ref = ref,
      level = level,
      weight = weight,
      parts = parts,
      aorg = aorg,
      summary = summary,
      type = if(exists("type")) (type) else (NULL)),
    class = "substitution")
  
}