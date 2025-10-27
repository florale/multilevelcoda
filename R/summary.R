#' Create a Summary of a \code{complr} object
#'
#' @param object An object of class \code{complr}.
#' @param ... generic argument, not in use.
#'
#' @importFrom compositions summary.acomp summary.rmult clo acomp rmult
#'
#' @method summary complr
#'
#' @examples
#' x1 <- complr(data = mcompd, sbp = sbp,
#'                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                idvar = "ID")
#' summary(x1)
#' x2 <- complr(data = mcompd, sbp = sbp,
#'                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' summary(x2)
#' @export
summary.complr <- function(object, ...) {
  # meta data specific to data
  transform_type <- object$transform
  idvar <- object$idvar
  nobs  <- nrow(object$dataout)
  ngrps <- if (!is.null(object$idvar))
    length(unique(object$dataout[[object$idvar]]))
  else
    nrow(object$dataout)
  
  # meta data specific to composition
  psi = lapply(object$output, `[[`, "psi")
  total = lapply(object$output, `[[`, "total")
  sequential_binary_partition = lapply(object$output, `[[`, "sbp")
  composition_geometry = lapply(object$output, function(x)
    class(x$X))
  logratio_class = lapply(object$output, function(x)
    class(x$Z))
  
  # output
  output_data <- get_variables(object)
  output_data <- lapply(output_data, function(x) Filter(Negate(is.null), x)) # remove NULL
  output_data <- do.call(cbind, output_data)
  meta_data   <- setNames(data.frame(total = unlist(total), row.names = "total"), colnames(output_data))
  output_data <- rbind(output_data, meta_data)
  colnames(output_data) <- paste0("Composition ", 1:ncol(output_data))
  
  # print output
  cat("Summary of Composition and Logratio Variables\n")
  cat("=============================================\n\n")
  cat(sprintf(" Number of observations: %s\n", nobs))
  if (!is.null(ngrps)) {
    cat(sprintf(" Number of groups (%ss): %s\n\n", idvar, ngrps))
  }
  cat(" Transformed variables:\n")
  cat(sprintf(" Composition geometry: %s, logratio ransformation: %s\n\n", 
              unique(unlist(composition_geometry)), 
              transform_type))
  print(output_data, row.names = TRUE)
  
  # save output
  out <- structure(list(
    meta = list(
      nobs = nobs,
      ngrps = ngrps,
      idvar = idvar,
      transform_type = transform_type
    ),
    transform = list(total = total, psi   = psi, sbp   = sequential_binary_partition),
    geometry = list(
      composition_geometry = unique(unlist(composition_geometry)),
      logratio_class = unique(unlist(logratio_class))
    ),
    output = output_data
  ),
  class = "summary.complr")
  invisible(out)
}

#' Print a Summary for a \code{complr} object
#'
#' @param x An object of class \code{complr}.
#' @param ... Other arguments passed to \code{\link{summary.complr}}.
#'
#' @seealso \code{\link{summary.complr}}
#'
#' @examples
#'
#' cilr <- complr(data = mcompd, sbp = sbp,
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                 idvar = "ID")
#' print(cilr)
#' @export
print.complr <- function(x, ...) {
  summary(x, ...)
}

#' Create a Summary of a fitted \code{brmsfit} model in a \code{brmcoda} object
#'
#' @param object An object of class \code{brmcoda}.
#' @param ... Other arguments passed to \code{\link[brms:summary.brmsfit]{summary.brmsfit}}.
#'
#' @method summary brmcoda
#'
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                  idvar = "ID", total = 1440),
#'   formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                      wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'
#'   summary(m)
#' }}
#' @export
summary.brmcoda <- function(object, ...) {
  if (inherits(object$model, "brmsfit")) {
    summary(object$model, ...)
  }
}

#' Print a Summary for a fitted \code{brmsfit} model in a \code{brmcoda} object
#'
#' @param x An object of class \code{brmcoda}.
#' @param ... Other arguments passed to \code{summary.brmcoda}.
#'
#' @seealso \code{\link{summary.brmcoda}}
#'
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID", total = 1440),
#'   formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                      wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'
#'  print(m)
#' }}
#' @export
print.brmcoda <- function(x, ...) {
  print(summary(x$model, ...), ...)
  
}

#' Create a Summary of a fitted \code{brmsfit} model from a \code{pivot_coord} object
#'
#' @param object An object of class \code{pivot_coord}.
#' @param digits A integer value used for number formatting. Default is \code{2}.
#' @param ... currently ignored.
#'
#' @return A data table of results.
#'
#' @method summary pivot_coord
#'
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID", total = 1440),
#'   formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                      wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'
#'   m_pc <- pivot_coord(m, method = "refit")
#'   summary(m_pc)
#' }}
#' @export
summary.pivot_coord <- function(object, digits = 2, ...) {
  
  if (isFALSE(object$summary)) {
    stop("The 'summary' method is currently not available for draws.")
  }
  
  if (isFALSE(digits == "asis")) {
    # object$output[, 1:3] <- round(object$output[, 1:3], digits)
    lapply(object$output, function(X)
      if (is.numeric(X))
        round(X, digits)
      else
        X)
  } else {
    object$output
  }
}

#' Create a Summary of a Substitution Model represented by a \code{substitution} object
#'
#' @param object A \code{substitution} class object.
#' @param delta A integer, numeric value or vector indicating the desired \code{delta}
#' at which substitution results should be summarised.
#' Default to all \code{delta} available in the \code{substitution} object.
#' @param to A character value or vector specifying the names of the compositional parts
#' that were reallocated to in the model.
#' @param from A character value or vector specifying the names of the compositional parts
#' that were reallocated from in the model.
#' @param ref Either a character value or vector ((\code{"grandmean"} and/or \code{"clustermean"} or \code{"users"}),
#' Default to all \code{ref} available in the \code{substitution} object.
#' @param level A character string or vector (\code{"between"} and/or \code{"within"}).
#' Default to all \code{level} available in the \code{substitution} object.
#' @param digits A integer value used for number formatting. Default is \code{2}.
#' @param ... generic argument, not in use.
#'
#' @return A summary of \code{substitution} object.
#'   \item{\code{Mean}}{ Posterior means.}
#'   \item{\code{CI_low} and \code{CI_high}}{ 95% credible intervals.}
#'   \item{\code{Delta}}{ Amount substituted across compositional parts.}
#'   \item{\code{From}}{ Compositional part that is substituted from.}
#'   \item{\code{To}}{ Compositional parts that is substituted to.}
#'   \item{\code{Level}}{ Level where changes in composition takes place. Either \code{between} or \code{within}.}
#'   \item{\code{Reference}}{ Either \code{grandmean}, \code{clustermean}, or \code{users}.}
#'
#' @method summary substitution
#'
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   ## fit a model with compositional predictor at between and between-person levels
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID", total = 1440),
#'   formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                      wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'
#'   subm <- substitution(object = m, delta = 5)
#'   summary(subm)
#' }}
#' @export
summary.substitution <- function(object,
                                 delta,
                                 to,
                                 from,
                                 ref,
                                 level,
                                 digits = 2,
                                 ...) {
  if (isFALSE(object$summary)) {
    stop("The 'summary' method is currently not available for substitution draws.")
  }
  
  if (missing(delta)) {
    delta <- object[["delta"]]
  } else {
    delta <- as.integer(delta)
  }
  
  if (missing(ref)) {
    ref <- object[["ref"]]
  } else {
    if (isFALSE(any(c("grandmean", "clustermean", "users") %in% ref))) {
      stop(
        "'ref' should be either one of the following: \"grandmean\", \"clustermean\", or \"users\"."
      )
    }
    ref <- as.character(ref)
  }
  
  if (missing(level)) {
    level <- object[["level"]]
  } else {
    if (isFALSE(any(c("between", "within", "aggregate") %in% level))) {
      stop(
        "'level' should be either one of the following: \"between\", \"within\", \"aggregate\"."
      )
    }
    level <- as.character(level)
  }
  
  if (missing(to)) {
    to <- object[["parts"]]
  } else {
    if (isFALSE(any(object[["parts"]] %in% to))) {
      stop(
        "'to' should be names of one or more compositional parts present in the substitution model."
      )
    }
    to <- as.character(to)
  }
  
  if (missing(from)) {
    from <- object[["parts"]]
  } else {
    if (isFALSE(any(object[["parts"]] %in% from))) {
      stop(
        "'from' should be names of one or more compositional parts present in the substitution model."
      )
    }
    from <- as.character(from)
  }

  # combine all results
  out <- rbindlist(lapply(
    Filter(Negate(is.null), object[c(
      "between_simple_sub",
      "within_simple_sub",
      "simple_sub",
      "between_avg_sub",
      "within_avg_sub",
      "avg_sub"
    )]),
    function(outi) {
      rbindlist(lapply(outi, function(outii) {
        rbindlist(Map(function(x, nm) cbind(x, Resp = nm), outii, names(outii)))
      }))
    }))

  if (object[["type"]] == "one-to-all") {
    out <- out[abs(Delta) %in% delta & Level %in% level & Reference %in% ref & (To %in% to | From %in% to)]
    out[, Delta := abs(Delta)]
  } else {
    out <- out[Delta %in% delta & Level %in% level & Reference %in% ref & To %in% to & From %in% from]
  }
  
  if (isTRUE(dim(out)[1] == 0)) {
    stop(
      "An empty data.table returned. Please check that the arguments match with your substitution object."
    )
  }
  if (isFALSE(digits == "asis")) {
    # out[, 1:3] <- round(out[, 1:3], digits)
    out[] <- lapply(out, function(X)
      if (is.numeric(X))
        round(X, digits)
      else
        X)
  }
  out
}

#' Print a Summary for a \code{substitution} object
#'
#' @param x A \code{substitution} object.
#' @param ... Additional arguments to be passed to to method \code{summary} of \code{substitution}.
#'
#' @seealso \code{\link{summary.substitution}}
#'
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   ## fit a model with compositional predictor at between and between-person levels
#'   m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
#'                                parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                idvar = "ID", total = 1440),
#'   formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                      wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'
#'   subm <- substitution(object = m, delta = 5)
#'   print(subm)
#' }}
#' @export
print.substitution <- function(x, ...) {
  print(summary(x, ...), ...)
}