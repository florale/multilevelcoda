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
#' @param basesub A \code{data.frame} or \code{data.table} of the base possible substitution of compositional parts.
#' This data set can be computed using function \code{\link{basesub}}. 
#' If \code{NULL}, all possible pairwise substitution of compositional parts are used.
#' @param ref Either a character value or vector or a dataset.
#' Can be \code{"grandmean"} and/or \code{"clustermean"}, or
#' a \code{data.frame} or \code{data.table} of user's specified reference grid consisting
#' of combinations of covariates over which predictions are made.
#' User's specified reference grid is only possible for simple substitution.
#' Single level models are default to \code{"grandmean"}.
#' @param summary A logical value. 
#' Should the estimate at each level of the reference grid (\code{FALSE}) 
#' or their average (\code{TRUE}) be returned? 
#' Default is \code{TRUE}.
#' Only applicable for model with covariates in addition to
#' the isometric log-ratio coordinates (i.e., adjusted model).
#' @param level A character string or vector. 
#' Should the estimate of multilevel models focus on the \code{"between"} and/or \code{"within"} or \code{"combined"} 
#' variance?
#' Single-level models are default to \code{"combined"}.
#' @param weight A character value specifying the weight to use in calculation of the reference composition.
#' If \code{"equal"}, give equal weight to units (e.g., individuals).
#' If \code{"proportional"}, weights in proportion to the frequencies of units being averaged 
#' (e.g., observations across individuals).
#' Default to \code{"equal"} for \code{ref = "grandmean"} and \code{"proportional"} for \code{ref = "clustermean"}.
#' @param ... Additional arguments passed to \code{\link{describe_posterior}}.
#' 
#' @return A list containing the results of multilevel compositional substitution model.
#' The first six lists contain the results of the substitution estimation for a compositional part. 
#'   \item{\code{Mean}}{ Posterior means.}
#'   \item{\code{CI_low} and \code{CI_high}}{ 95% credible intervals.}
#'   \item{\code{Delta}}{ Amount substituted across compositional parts.}
#'   \item{\code{From}}{ Compositional part that is substituted from.}
#'   \item{\code{To}}{ Compositional parts that is substituted to.}
#'   \item{\code{Level}}{ Level where changes in composition takes place.
#'   \item{\code{Reference}}{ Either \code{grandmean}, \code{clustermean}, or \code{users}.}
#' 
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo mean.acomp
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   cilr <- complr(data = mcompd, sbp = sbp,
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                  idvar = "ID", total = 1440)
#'   
#'   # model with compositional predictor at between and between-person levels of variance
#'   fit1 <- brmcoda(complr = cilr,
#'                   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                                      wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                   chain = 1, iter = 500, backend = "cmdstanr")
#'   sub1 <- substitution(object = fit1, delta = 5, level = c("between", "within"))
#'   
#'   # model with compositional predictor at combined level of variance
#'   fit2 <- brmcoda(complr = cilr,
#'                   formula = Stress ~ ilr1 + ilr2 + ilr3 + ilr4 + (1 | ID),
#'                   chain = 1, iter = 500, backend = "cmdstanr")
#'   sub2 <- substitution(object = fit2, delta = 5, level = c("combined"))
#'   
#' }}
#' @export
substitution <- function(object,
                         delta,
                         basesub = NULL,
                         summary = TRUE,
                         ref = c("grandmean", "clustermean"),
                         level = c("between", "within", "combined"),
                         weight = c("equal", "proportional"),
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
  
  if (missing(basesub)) {
    count <- length(object$complr$parts)
    n <- count - 2
    
    subvars1 <- c(1, -1)
    subvars2 <- rep(0, n)
    subvars <- c(subvars1, subvars2)
    
    nc <- length(subvars)
    nr <- (nc - 1) * count
    
    basesub <- matrix(0, nrow = nr, ncol = nc, dimnames = list(NULL, object$complr$parts))
    k <- 0
    
    for (i in 1:nc)
      for (j in 1:nc)
        if (i != j) {
          k <- k + 1
          basesub[k, c(i, j)] <- c(1, -1)
        }
    
    basesub <- as.data.table(basesub)
    names(basesub) <- object$complr$parts
    
  } else if(isFALSE(missing(basesub))) {
    if (isFALSE(identical(ncol(basesub), length(object$complr$parts)))) {
      stop(sprintf(
        "The number of columns in 'basesub' (%d) should be the same as the compositional parts in 'parts' (%d).",
        ncol(basesub),
        length(object$complr$parts)
      ))
    }
    if (isFALSE(identical(colnames(basesub), object$complr$parts))) {
      stop(sprintf(
        "The names of compositional parts should be the same in 'basesub' (%s) and 'parts' (%s).",
        colnames(basesub),
        object$complr$parts
      ))
    }
  }
  
  # what type of model is being estimated
  model_fixef <- rownames(fixef(object))
  model_ranef <- if(dim(object$model$ranef)[1] > 0) (names(ranef(object))) else (NULL)
  
  model_fixef_level <- NULL
  model_fixef_coef <- NULL
  
  if (isFALSE(is.null(object$complr$between_logratio)) && length(grep(paste0(names(object$complr$between_logratio), collapse = "|"), model_fixef, value = T)) > 0) {
    model_fixef_level <- append(model_fixef_level, "between")
    model_fixef_coef  <- append(model_fixef_coef, grep(paste0(names(object$complr$between_logratio), collapse = "|"), model_fixef, value = T))
  }
  if (isFALSE(is.null(object$complr$within_logratio)) && length(grep(paste0(names(object$complr$within_logratio), collapse = "|"), model_fixef, value = T)) > 0) {
    model_fixef_level <- append(model_fixef_level, "within")
    model_fixef_coef  <- append(model_fixef_coef, grep(paste0(names(object$complr$within_logratio), collapse = "|"), model_fixef, value = T))
  }
  if (length(grep(paste0(names(object$complr$logratio), collapse = "|"), model_fixef, value = T)) > 0) {
    model_fixef_level <- append(model_fixef_level, "combined")
    model_fixef_coef  <- append(model_fixef_coef, grep(paste0(names(object$complr$logratio), collapse = "|"), model_fixef, value = T))
  }
  
  # single level or multilevel
  if (length(model_ranef) > 0) {
    model_ranef_level <- "multilevel"
    model_ranef_coef <- model_ranef
  } else {
    model_ranef_level <- "single"
    model_ranef_coef <- NULL
  }
  
  # ensure sub args make sense
  ## only grandmean for single level model
  if ("combined" %in% model_fixef_level) {
    if (model_ranef_level == "single") {
      if ("clustermean" %in% ref) {
        warning("Can only use grandmean for single level model.")
      }
      level <- "combined"
      ref <- "grandmean"
      weight <- "equal"
    } else {
      level <- "combined"
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
  if (all(c("between", "within", "combined") %in% level)) {
    level <- c("between", "within")
  }
  
  ## level args match with coefs in object
  if (any(level %in% c("between", "within"))) {
    if (isFALSE(any(model_fixef_level %in% c("between", "within")))) {
      stop(sprintf(
        "'between' and 'within' substitution analysis cannot be computed
  on a model estimated using the (%s) variance of ilrs.
  Please specify the level argument as \"(%s)\" instead or refit 'brmcoda' model.",
  model_fixef_level,
  model_fixef_level
      ))
    }
  }
  if ("combined" %in% level) {
    if (isFALSE("combined" %in% model_fixef_level)) {
      stop(sprintf(
        "'combined' substitution analysis cannot be computed
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
        basesub = basesub,
        summary = summary,
        ref = "grandmean",
        level = "between",
        weight = weight)
    } 
    else if (inherits(ref, c("data.table", "data.frame", "matrix"))) {
      bout <- bsub(
        object = object,
        delta = delta,
        basesub = basesub,
        summary = summary,
        ref = ref,
        level = "between",
        weight = weight)
    }
    if ("clustermean" %in% ref) {
      bmout <-
        bsubmargins(
          object = object,
          delta = delta,
          basesub = basesub,
          ref = "clustermean",
          level = "between",
          weight = weight)
    }
  }
  
  wmout <- wout <- NULL
  if ("within" %in% level) {
    if ("grandmean" %in% ref) {
      wout <- wsub(
        object = object,
        delta = delta,
        basesub = basesub,
        summary = summary,
        ref = "grandmean",
        level = "within",
        weight = weight)
    } 
    else if (inherits(ref, c("data.table", "data.frame", "matrix"))) {
      wout <- wsub(
        object = object,
        delta = delta,
        basesub = basesub,
        summary = summary,
        ref = ref,
        level = "within",
        weight = weight)
    }
    if ("clustermean" %in% ref) {
      wmout <-
        wsubmargins(
          object = object,
          delta = delta,
          basesub = basesub,
          ref = "clustermean",
          level = "within",
          weight = weight)
    }
  }
  
  tmout <- tout <- NULL
  if ("combined" %in% level) {
    if ("grandmean" %in% ref) {
      tout <- sub(
        object = object,
        delta = delta,
        basesub = basesub,
        summary = summary,
        ref = "grandmean",
        level = "combined",
        weight = weight)
    } 
    else if (inherits(ref, c("data.table", "data.frame", "matrix"))) {
      tout <- sub(
        object = object,
        delta = delta,
        basesub = basesub,
        summary = summary,
        ref = ref,
        level = "between",
        weight = weight)
    }
    if ("clustermean" %in% ref) {
      tmout <-
        submargins(
          object = object,
          delta = delta,
          basesub = basesub,
          ref = "clustermean",
          level = "combined",
          weight = weight)
    }
  }
  
  structure(
    list(
      between_simple_sub = if(exists("bout")) (bout) else (NULL),
      within_simple_sub = if(exists("wout")) (wout) else (NULL),
      simple_sub = if(exists("tout")) (tout) else (NULL),
      between_avg_sub = if(exists("bmout")) (bmout) else (NULL),
      within_avg_sub = if(exists("wmout")) (wmout) else (NULL),
      avg_sub = if(exists("tmout")) (tmout) else (NULL),
      delta = delta,
      ref = ref,
      level = level,
      weight = weight,
      parts = object$complr$parts,
      summary = summary),
    class = "substitution")
}