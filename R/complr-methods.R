#' Checks if argument is a \code{complr} object
#'
#' @param x An object of class \code{complr}.
#'
#' @export
is.complr <- function(x) {
  inherits(x, "complr")
}

#' Mean amounts and mean compositions presented in a \code{complr} object.
#' 
#' @param x An object of class \code{complr}.
#' @param class  Optional. Can be \code{"composition"} and/or \code{"logratio"} to
#' specify the geometry of the composition.
#' @param level  Optional. Can be \code{"between"}, \code{"within"}, and/or \code{"combined"}
#' indicating the level of the geometry.
#' @param weight A character value specifying the weight to use in calculation of the reference composition.
#' If \code{"equal"}, give equal weight to units (e.g., individuals).
#' If \code{"proportional"}, weights in proportion to the frequencies of units being averaged 
#' (e.g., observations across individuals)
#' Default is \code{equal}.
#' @param digits A integer value used for number formatting. Default is \code{3}.
#' @param ... generic argument, not in use.
#' 
#' @method mean complr
#' 
#' @examples
#' cilr <- complr(data = mcompd, sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), 
#'                 idvar = "ID")
#' mean(cilr)
#' @export
mean.complr <- function(x, ...,
                         class = c("comp", "logratio"),
                         level = c("between", "within", "combined"),
                         weight = c("equal", "proportional"),
                         digits = 3
                         ) {
  
  # Assemble
  output <- .get.complr(object = x)
  
  ## Out
  if ("comp" %in% class) {
    if ("combined" %in% level) {
      cat("\n", "Raw compositionositional Mean:", "\n")
      print(output$comp$mean, digits = digits)
    }
    if ("between" %in% level) {
      cat("\n", "Between-level compositionositional Mean:", "\n")
      print((output$between_comp$mean), digits = digits)
    }
    if ("within" %in% level) {
      cat("\n", "Within-level compositionositional Mean:", "\n")
      print((output$Withincomp$mean), digits = digits)
    }
  }
  if ("logratio" %in% class) {
    if ("combined" %in% level) {
      cat("\n", "Isometric Log-ratio (Real) Mean:", "\n")
      print(as.data.table(output$logratio)[4,], row.names = FALSE, digits = digits)
    }
    if ("between" %in% level) {
      cat("\n", "Between-level Isometric Log-ratio (Real) Mean:", "\n")
      print(as.data.table(output$between_logratio)[4,], row.names = FALSE, digits = digits)
    }
    if ("within" %in% level) {
      cat("\n", "Within-level Isometric Log-ratio (Real) Mean:", "\n")
      print(as.data.table(output$within_logratio)[4,], row.names = FALSE, digits = digits)
    }
  }
  
  ### Return output invisibly
  output <- lapply(output, function(X) {
    row.names(X) <- NULL
    return(X)
  })
  
  return(invisible(output))
}

#' Extract compositionositional Data from \code{complr} object.
#'
#' Extract amounts and compositions in conventional formats
#' as data.frames, matrices, or arrays.
#'
#' @inheritParams mean.complr
#' @param row.names,optional Unused and only added for consistency with
#' the \code{\link[base:as.data.frame]{as.data.frame}} generic.
#'
#' @export
as.data.frame.complr <- function(x, row.names = NULL, optional = TRUE,
                                  class = c("comp", "logratio"),
                                  level = c("between", "within", "combined"),
                                  ...) {
  allout <- lapply(x[1:6], as.data.frame)
  output <- data.frame()
  
  ## Out
  if ("comp" %in% class) {
    if ("combined" %in% level) {
      output <- allout[[1]]
    }
    if ("between" %in% level) {
      output <- cbind(output, allout[[2]])
    }
    if ("within" %in% level) {
      output <- cbind(output, allout[[3]])
    }
  }
  if ("logratio" %in% class) {
    if ("combined" %in% level) {
      output <- cbind(output, allout[[4]])
    }
    if ("between" %in% level) {
      output <- cbind(output, allout[[5]])
    }
    if ("within" %in% level) {
      output <- cbind(output, allout[[6]])
    }
  }
  output
}

#' @rdname as.data.frame.complr
#' @export
as.matrix.complr <- function(x, 
                              class = c("comp", "logratio"),
                              level = c("between", "within", "combined"),
                              ...) {
  as.matrix(as.data.frame(x, class = class, level = level, ...))
}

# ----------------- Extract compositionositional Data -----------------
.get.complr <- function(object, 
                         class = c("comp", "logratio"),
                         level = c("between", "within", "combined"),
                         weight = c("equal", "proportional"),
                         digits = 3,
                         ...) {
  ## Assemble
  if (identical(weight, "proportional")) {
    weight <- "proportional"
  } else {
    weight <- "equal"
  }
  
  if (weight == "equal") {
    bcomp <- cbind(object$data[, object$idvar, with = FALSE], object$between_comp)[!duplicated(get(object$idvar))]
    wcomp <- cbind(object$data[, object$idvar, with = FALSE], object$within_comp)[!duplicated(get(object$idvar))]
    tcomp <- cbind(object$data[, object$idvar, with = FALSE], object$comp)[!duplicated(get(object$idvar))]
    
    bilr <- cbind(object$data[, object$idvar, with = FALSE], object$between_logratio)[!duplicated(get(object$idvar))]
    wilr <- cbind(object$data[, object$idvar, with = FALSE], object$within_logratio)[!duplicated(get(object$idvar))]
    tilr <- cbind(object$data[, object$idvar, with = FALSE], object$logratio)[!duplicated(get(object$idvar))]
    
    output <- list(
      summary(acomp(bcomp[, -1]), robust = TRUE),
      summary(acomp(wcomp[, -1]), robust = TRUE),
      summary(acomp(tcomp[, -1]), robust = TRUE),
      data.frame(summary(rmult(bilr[, -1]))),
      data.frame(summary(rmult(wilr[, -1]))),
      data.frame(summary(rmult(tilr[, -1])))
    )
    
  } else if (weight == "proportional") {
    output <- list(
      summary(object$between_comp, robust = TRUE),
      summary(object$within_comp, robust = TRUE),
      summary(object$comp, robust = TRUE),
      data.frame(summary(object$between_logratio)),
      data.frame(summary(object$within_logratio)),
      data.frame(summary(object$logratio))
    )
  }
  
  names(output) <- c("between_comp", "within_comp", "comp", "between_logratio", "within_logratio", "logratio")
  
  output
}