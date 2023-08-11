#' Multilevel Compositional Substitution Model
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
#' @param summary A logical value. 
#' Should the estimate at each level of the reference grid (\code{FALSE}) 
#' or their average (\code{TRUE}) be returned? 
#' Default is \code{TRUE}.
#' Only applicable for model with covariates in addition to
#' the isometric log-ratio coordinates (i.e., adjusted model).
#' @param level A character string or vector. 
#' Should the estimate be at the \code{"between"} and/or \code{"within"} level?
#' @param weight A character value specifying the weight to use in calculation of the reference composition.
#' If \code{"equal"}, give equal weight to units (e.g., individuals).
#' If \code{"proportional"}, weights in proportion to the frequencies of units being averaged 
#' (e.g., observations across individuals)
#' Default is \code{equal}.
#' @param ... Additional arguments passed to \code{\link{describe_posterior}}.
#' 
#' @return A list containing the results of multilevel compositional substitution model.
#' The first four lists contain the results of the substitution estimation for a compositional part. 
#' \itemize{
#'   \item{\code{Mean}}{ Posterior means.}
#'   \item{\code{CI_low} and \code{CI_high}}{ 95% credible intervals.}
#'   \item{\code{Delta}}{ Amount substituted across compositional parts.}
#'   \item{\code{From}}{ Compositional part that is substituted from.}
#'   \item{\code{To}}{ Compositional parts that is substituted to.}
#'   \item{\code{Level}}{ Level where changes in composition takes place. Either \code{between} or \code{within}.}
#'   \item{\code{Reference}}{ Either \code{grandmean}, \code{clustermean}, or \code{users}.}
#' }
#' 
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo mean.acomp
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   cilr <- compilr(data = mcompd, sbp = sbp,
#'                   parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                   idvar = "ID", total = 1440)
#'   
#'   # model with compositional predictor at between and between-person levels
#'   m <- brmcoda(compilr = cilr,
#'                formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                  wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                chain = 1, iter = 500, backend = "cmdstanr")
#'   
#'   subm <- substitution(object = m, delta = 5)
#' }}
#' @export
substitution <- function(object,
                         delta,
                         basesub = NULL,
                         summary = TRUE,
                         ref = c("grandmean", "clustermean"),
                         level = c("between", "within"),
                         weight = c("equal", "proportional"),
                         ...) {
  
  if (isTRUE(missing(object))) {
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
  
  if(isFALSE(missing(delta))) {
    if (isFALSE(is.integer(delta))) {
      if (isFALSE(delta > 0)) {
        stop(" 'delta' must be an positive integer value.")
      }
    }
  } else if (isTRUE(missing(delta))){
    stop(paste(
      "'delta' is a required argument and cannot be missing;",
      "  it should be interger, numeric positive value or vector", 
      "  to specify the change in units across compositional parts", 
      sep = "\n"))
  }
  if (identical(weight, "proportional")) {
    weight <- "proportional"
  } else {
    weight <- "equal"
  }
  
  if (isTRUE(missing(basesub))) {
    count <- length(object$CompILR$parts)
    n <- count - 2
    
    subvars1 <- c(1, -1)
    subvars2 <- rep(0, n)
    subvars <- c(subvars1, subvars2)
    
    nc <- length(subvars)
    nr <- (nc - 1) * count
    
    basesub <- matrix(0, nrow = nr, ncol = nc, dimnames = list(NULL, object$CompILR$parts))
    k <- 0
    
    for (i in 1:nc)
      for (j in 1:nc)
        if (i != j) {
          k <- k + 1
          basesub[k, c(i, j)] <- c(1, -1)
        }
    
    basesub <- as.data.table(basesub)
    names(basesub) <- object$CompILR$parts
    
  } else if(isFALSE(missing(basesub))) {
    if (isFALSE(identical(ncol(basesub), length(object$CompILR$parts)))) {
      stop(sprintf(
        "The number of columns in 'basesub' (%d) must be the same
        as the compositional parts in 'parts' (%d).",
        ncol(basesub),
        length(object$CompILR$parts)))
    }
    if (isFALSE(identical(colnames(basesub), object$CompILR$parts))) {
      stop(sprintf(
        "The names of compositional parts must be the
        same in 'basesub' (%s) and 'parts' (%s).",
        colnames(basesub),
        object$CompILR$parts))
    }
  }
  
  if (isTRUE("between" %in% level)) {
    if (isTRUE("grandmean" %in% ref)) {
      bout <- bsub(
        object = object,
        delta = delta,
        basesub = basesub,
        summary = summary,
        ref = "grandmean",
        level = "between",
        weight = weight)
    } 
    else if (isTRUE(inherits(ref, c("data.table", "data.frame", "matrix")))) {
      bout <- bsub(
        object = object,
        delta = delta,
        basesub = basesub,
        summary = summary,
        ref = ref,
        level = "between",
        weight = weight)
    }
    if (isTRUE("clustermean" %in% ref)) {
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
  
  if (isTRUE("within" %in% level)) {
    if (isTRUE("grandmean" %in% ref)) {
      wout <- wsub(
        object = object,
        delta = delta,
        basesub = basesub,
        summary = summary,
        ref = "grandmean",
        level = "within",
        weight = weight)
    } 
    else if (isTRUE(inherits(ref, c("data.table", "data.frame", "matrix")))) {
      wout <- wsub(
        object = object,
        delta = delta,
        basesub = basesub,
        summary = summary,
        ref = ref,
        level = "within",
        weight = weight)
    }
    if (isTRUE("clustermean" %in% ref)) {
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
  
  structure(
    list(
      BetweenSub = if(exists("bout")) (bout) else (NULL),
      WithinSub = if(exists("wout")) (wout) else (NULL),
      BetweenSubMargins = if(exists("bmout")) (bmout) else (NULL),
      WithinSubMargins = if(exists("wmout")) (wmout) else (NULL),
      delta = delta,
      ref = ref,
      level = level,
      weight = weight,
      parts = object$CompILR$parts,
      summary = summary),
    class = "substitution")
}