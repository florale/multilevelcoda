#' Create a summary of a \code{compilr} object
#' 
#' @param object An object of class \code{compilr}.
#' @param x  Optional. Can be \code{"composition"} and/or \code{"logratio"} to
#' specify the geometry of the composition.
#' @param level  Optional. Can be \code{"between"}, \code{"within"}, and/or \code{total}
#' indicating the level of the geometry.
#' @param digits A integer value used for number formatting. Default is \coda{3}.
#' 
#' @importFrom compositions summary.acomp summary.rmult
#' @method summary compilr
#' @export
summary.compilr <- function(object,
                            x = c("composition", "logratio"),
                            level = c("between", "within", "total"),
                            digits = 3,
                            ...) {
  # General info
  cat("  Compositional components are: ")
  cat(paste(object$parts, collapse = ", "), "\n")
  cat("  Composition is closed to    : ")
  cat(object$total, "\n")
  cat("  Number of observations      : ")
  cat(nrow(object$data), "\n")
  cat("  Number of levels            : ")
  cat(length(unique(object$data[[object$idvar]])), "\n")

  # comp
  sumc <- list(
    summary(object$BetweenComp, robust = TRUE),
    summary(object$WithinComp, robust = TRUE),
    summary(object$TotalComp, robust = TRUE)
  )
  names(sumc) <- c("bcomp", "wcomp", "tcomp")
  
  varn <- c("Compositional Mean",
            "Geometric Mean of the Pairwise Ratios",
            "Variation Matrix",
            "One-sigma Factor of Pairwise Ratios",
            "Inverse of One-sigma Factor of Pairwise Ratios",
            "Min of Pairwise Ratios",
            "Q1 of Pairwise Ratios",
            "Median of Pairwise Ratios",
            "Q3 of Pairwise Ratios",
            "Max of Pairwise Ratios",
            "Missingness")
  
  sumc <- lapply(sumc, function(X) {
    x1 <- list(as.matrix(t(X[[1]])))
    x2 <- X[-c(1, length(X))]
    x3 <- tail(X, 1)
    
    names(x1) <- head(varn, 1)
    names(x2) <- varn[-c(1, length(varn))]
    names(x3) <- tail(varn, 1)
    
    unlist(list(x1,
                x2,
                x3), recursive = FALSE)
  })
  
  if ("composition" %in% x) {
    if ("between" %in% level) {
      cat("\n", "——— Composition at between-level ———", "\n")
      for (i in seq_along(sumc$bcomp)) {
        cat(paste0("\n", names(sumc$bcomp)[i], ":\n"))
        print(sumc$bcomp[[i]], digits = digits)
      }
    }
    if ("within" %in% level) {
      cat("\n", "——— Composition at within-level ———", "\n")
      for (i in seq_along(sumc$wcomp)) {
        cat(paste0("\n", names(sumc$wcomp)[i], ":\n"))
        print(sumc$wcomp[[i]], digits = digits)
      }
    }
    if ("total" %in% level) {
      cat("\n", "——— Composition with total variance ———", "\n")
      for (i in seq_along(sumc$tcomp)) {
        cat(paste0("\n", names(sumc$tcomp)[i], ":\n"))
        print(sumc$tcomp[[i]], digits = digits)
      }
    }
  }
  
  # log ratio
  sumlr <- list(
    data.frame(summary(object$BetweenILR)),
    data.frame(summary(object$WithinILR)),
    data.frame(summary(object$TotalILR))
  )
  
  names(sumlr) <- c("bilr", "wilr", "tilr")
  
  if ("logratio" %in% x) {
    if ("between" %in% level) {
      cat("\n", "——— Log-ratio at between-level ———", "\n")
      print(sumlr$bilr, digits = digits)
    }
    if ("within" %in% level) {
      cat("\n", "——— Log-ratio at within-level ———", "\n")
      print(sumlr$wilr, digits = digits)
    }
    if ("total" %in% level) {
      cat("\n", "——— Log-ratio with total variance ———", "\n")
      print(sumlr$tilr, digits = digits)
    }
  }
  
  ### Return output invisibly
  output <- lapply(X = list(sumc, sumlr), FUN = function(x) {
    row.names(x) <- NULL
    return(x)
  })
  return(invisible(output))

}

#' Print a summary for a \code{compilr} object
#' 
#' @inheritParams summary.compilr
#' 
#' @seealso \code{\link{summary.compilr}}
#' 
#' @export
print.compilr <- function(x, ...) {
  summary(object <- x, ...)
  
}
#' Create a summary of a fitted \code{brmsfit} model in a \code{brmcoda} object
#' 
#' @param object An object of class \code{brmcoda}.
#' @param ... Other arguments passed to \code{summary.brmsfit}
#' 
#' @method summary brmcoda
#' @export
summary.brmcoda <- function(object, ...) {
  if (inherits(object$Model, "brmsfit")) {
    summary(object$Model, ...)
  }
}

#' Print a summary for a fitted \code{brmsfit} model in a \code{brmcoda} object
#' 
#' @inheritParams summary.brmcoda
#' 
#' @seealso \code{\link{summary.brmcoda}}
#' 
#' @export
print.brmcoda <- function(x, ...) {
  print(summary(x$Model, ...), ...)
  
}

#' Create a summary of a substitution model represented by a \code{\link{substitution}} object
#' 
#' @param object A \code{\link{substitution}} class object.
#' @param delta A integer, numeric value or vector indicating the desired \code{delta} 
#' at which substitution results should be summarised.
#' Default to all \code{delta} available in the \code{\link{substitution}} object.
#' @param to A character value or vector specifying the names of the compositional parts
#' that were reallocated to in the model.
#' @param from A character value or vector specifying the names of the compositional parts
#' that were reallocated from in the model.
#' @param ref Either a character value or vector ((\code{"grandmean"} and/or \code{"clustermean"} or \code{"users"}),
#' Default to all \code{ref} available in the \code{\link{substitution}} object .
#' @param level A character string or vector (\code{"between"} and/or \code{"within"}).
#' Default to all \code{level} available in the \code{\link{substitution}} object.
#' @param digits A integer value used for number formatting. Default is \coda{2}.
#' @param ... generic argument, not in use.
#' 
#' @return A summary of \code{\link{substitution}} object.
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
#' @method summary substitution
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
    out[] <- lapply(out, function(x) if(is.numeric(x)) round(x, digits) else x)
  }
  
  out
  
}

#' Print a summary for a \code{substitution} object
#' 
#' @param x A \code{substitution} object.
#' @param ... Additional arguments to be passed to to method \code{summary} of \code{substitution}.
#' 
#' @seealso \code{\link{summary.substitution}}
#' 
#' @export
print.substitution <- function(x, ...) {
  print(summary(object, ...), ...)
}