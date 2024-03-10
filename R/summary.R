#' Create a Summary of a \code{compilr} object
#' 
#' @param object An object of class \code{compilr}.
#' @param class Optional. Can be \code{"composition"} and/or \code{"logratio"} to
#' specify the geometry of the composition.
#' @param level Optional. Can be \code{"between"}, \code{"within"}, and/or \code{"total"}
#' indicating the level of the geometry.
#' @param weight A character value specifying the weight to use in calculation of the reference composition.
#' If \code{"equal"}, give equal weight to units (e.g., individuals).
#' If \code{"proportional"}, weights in proportion to the frequencies of units being averaged 
#' (e.g., observations across individuals).
#' Default is \code{equal}.
#' @param digits A integer value used for number formatting. Default is \code{3}.
#' @param ... generic argument, not in use.
#' 
#' @importFrom compositions summary.acomp summary.rmult clo acomp rmult
#' @importFrom utils head tail
#' 
#' @method summary compilr
#' 
#' @examples
#' cilr <- compilr(data = mcompd, sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), 
#'                 idvar = "ID")
#' summary(cilr)
#' @export
summary.compilr <- function(object,
                            class = c("composition", "logratio"),
                            level = c("between", "within", "total"),
                            weight = c("equal", "proportional"),
                            digits = 3,
                            ...) {
  
  ## Assemble
  output <- .get.compilr(object)
  
  varn <- c("Compositional Mean", #keep
            "Geometric Mean of the Pairwise Ratios", 
            "Variation Matrix", #keep
            "One-sigma Factor of Pairwise Ratios",
            "Inverse of One-sigma Factor of Pairwise Ratios",
            "Min of Pairwise Ratios",
            "Q1 of Pairwise Ratios",
            "Median of Pairwise Ratios",
            "Q3 of Pairwise Ratios",
            "Max of Pairwise Ratios",
            "Missingness" #keep
  )
  
  output[1:3] <- lapply(output[1:3], function(X) {
    x1 <- list(clo(X[[1]]))
    x2 <- X[-c(1, length(X))]
    x3 <- tail(X, 1)
    
    names(x1) <- head(varn, 1)
    names(x2) <- varn[-c(1, length(varn))]
    names(x3) <- tail(varn, 1)
    
    unlist(list(x1, x2, x3), recursive = FALSE)
  })
  output[1:3] <- lapply(output[1:3], "[", c("Compositional Mean", 
                                            "Variation Matrix"
                                            # , "Missingness"
                                            ))
  
  ## Summary
  # General info
  cat("  Compositional components are: ")
  cat(paste(object$parts, collapse = ", "), "\n")
  cat("  Composition is closed to    : ")
  cat(object$total, "\n")
  cat("  Geometry                    : ")
  cat("relative composition ('acomp') and", "\n")
  cat("                                isometric log-ratios ('real multivariate')", "\n")
  cat("  Number of observations      : ")
  cat(nrow(object$data), "\n")
  cat("  Number of levels            : ")
  cat(length(unique(object$data[[object$idvar]])), "\n")
  
  # cat("\n", "———— Arithmetic Statistics ————", "\n")
  # print(JWileymisc::egltable(BetweenComp[, -1], ...))
  
  if ("composition" %in% class) {
    if ("total" %in% level) {
      cat("\n", " Raw Composition Statistics: ", "\n")
      for (i in seq_along(output$TotalComp)) {
        cat(paste0("\n", names(output$TotalComp)[i], ":", "\n"))
        print(output$TotalComp[[i]], digits = digits)
      }
    }
    if ("between" %in% level) {
      cat("\n", " Between-level Composition: ", "\n")
      for (i in seq_along(output$BetweenComp)) {
        cat(paste0("\n", names(output$BetweenComp)[i], ":", "\n"))
        print(output$BetweenComp[[i]], digits = digits)
      }
    }
    if ("within" %in% level) {
      cat("\n", " Within-level Composition: ", "\n")
      for (i in seq_along(output$WithinComp)) {
        cat(paste0("\n", names(output$WithinComp)[i], ":", "\n"))
        print(output$WithinComp[[i]], digits = digits)
      }
    }
  }
  
  if ("logratio" %in% class) {
    if ("total" %in% level) {
      cat("\n", " Raw Isometric Log-ratios: ", "\n")
      print(output$TotalILR, digits = digits)
    }
    if ("between" %in% level) {
      cat("\n", " Between-level Isometric Log-ratios: ", "\n")
      print(output$BetweenILR, digits = digits)
    }
    if ("within" %in% level) {
      cat("\n", " Within-level Isometric Log-ratios: ", "\n")
      print(output$WithinILR, digits = digits)
    }
  }
  
  ### Return output invisibly
  output <- lapply(output, function(X) {
    row.names(X) <- NULL
    return(X)
  })
  return(invisible(output))
}

#' Print a Summary for a \code{compilr} object
#' 
#' @param x An object of class \code{compilr}.
#' @param ... Other arguments passed to \code{\link{summary.compilr}}.
#' 
#' @seealso \code{\link{summary.compilr}}
#' 
#' @examples
#' 
#' cilr <- compilr(data = mcompd, sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), 
#'                 idvar = "ID")
#' print(cilr)
#' @export
print.compilr <- function(x, ...) {
  summary(x, ...)
}

#' Create a Summary of a fitted \code{brmsfit} model in a \code{brmcoda} object
#' 
#' @param object An object of class \code{brmcoda}.
#' @param ... Other arguments passed to \code{\link{summary.brmsfit}}.
#' 
#' @method summary brmcoda
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   m <- brmcoda(compilr = compilr(data = mcompd, sbp = sbp,
#'                                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                  idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'   
#'   summary(m)
#' }}
#' @export
summary.brmcoda <- function(object, ...) {
  if (inherits(object$Model, "brmsfit")) {
    summary(object$Model, ...)
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
#'   m <- brmcoda(compilr = compilr(data = mcompd, sbp = sbp,
#'                                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                  idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'   
#'  print(m)
#' }}
#' @export
print.brmcoda <- function(x, ...) {
  print(summary(x$Model, ...), ...)
  
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
#'   m <- brmcoda(compilr = compilr(data = mcompd, sbp = sbp,
#'                                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                  idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'   chain = 1, iter = 500,
#'   backend = "cmdstanr")
#'   
#'   subm <- substitution(object = m, delta = 5)
#'   summary(subm)
#' }}
#' @export
summary.substitution <- function(object, delta, to, from,
                                 ref, level,
                                 digits = 2, ...) { 
  
  if (isTRUE(missing(delta))) {
    delta <- object$delta
  } else {
    delta <- as.integer(delta)
  }
  
  if (isTRUE(missing(ref))) {
    ref <- object$ref
  } else {
    if (isFALSE(any(c("grandmean", "clustermean", "users") %in% ref))) {
      stop("'ref' should be grandmean and/or clustermean or users.")
    }
    ref <- as.character(ref)
  }
  if (isTRUE(missing(level))) {
    level <- object$level
  } else {
    if (isFALSE(any(c("between", "within") %in% level))) {
      stop("'level' should be between and/or within.")
    }
    level <- as.character(level)
  }
  
  if (isTRUE(missing(to))) {
    to <- object$parts
  } else {
    if (isFALSE(any(object$parts %in% to))) {
      stop("'to' should be names of one or more compositional parts present in the substitution model.")
    }
    to <- as.character(to)
  }
  if (isTRUE(missing(from))) {
    from <- object$parts
  } else {
    if (isFALSE(any(object$parts %in% from))) {
      stop("'from' should be names of one or more compositional parts present in the substitution model.")
    }
    from <- as.character(from)
  }
  
  out <- lapply(object[1:4], rbindlist)
  out <- rbindlist(out, use.names = TRUE)
  out <- out[Delta %in% delta & Level %in% level & Reference %in% ref & To %in% to & From %in% from]
  
  if(isTRUE(dim(out)[1] == 0)) {
    stop("An empty data.table returned. Are you sure your arguments are correct?")
  }
    
  if(isFALSE(digits == "asis")) {
    # out[, 1:3] <- round(out[, 1:3], digits)
    out[] <- lapply(out, function(X) if(is.numeric(X)) round(X, digits) else X)
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
#'   m <- brmcoda(compilr = compilr(data = mcompd, sbp = sbp,
#'                                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                                  idvar = "ID", total = 1440),
#'   formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'     wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
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