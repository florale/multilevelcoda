#' Checks if argument is a \code{compilr} object
#'
#' @param x An object of class \code{compilr}.
#'
#' @export
is.compilr <- function(x) {
  inherits(x, "compilr")
}

#' Mean amounts and mean compositions presented in a \code{compilr} object.
#' 
#' @param x An object of class \code{compilr}.
#' @param class  Optional. Can be \code{"composition"} and/or \code{"logratio"} to
#' specify the geometry of the composition.
#' @param level  Optional. Can be \code{"between"}, \code{"within"}, and/or \code{total}
#' indicating the level of the geometry.
#' @param weight A character value specifying the weight to use in calculation of the reference composition.
#' If \code{"equal"}, give equal weight to units (e.g., individuals).
#' If \code{"proportional"}, weights in proportion to the frequencies of units being averaged 
#' (e.g., observations across individuals)
#' Default is \code{equal}.
#' @param digits A integer value used for number formatting. Default is \code{3}.
#' @param ... generic argument, not in use.
#' 
#' @method mean compilr
#' 
#' @examples
#' cilr <- compilr(data = mcompd, sbp = sbp, 
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), 
#'                 idvar = "ID")
#' mean(cilr)
#' @export
mean.compilr <- function(x, ...,
                         class = c("composition", "logratio"),
                         level = c("between", "within", "total"),
                         weight = c("equal", "proportional"),
                         digits = 3
                         ) {
  
  # Assemble
  output <- .get.compilr(object = x)
  
  ## Out
  if ("composition" %in% class) {
    if ("total" %in% level) {
      cat("\n", "Raw Compositional Mean:", "\n")
      print(output$TotalComp$mean, digits = digits)
    }
    if ("between" %in% level) {
      cat("\n", "Between-level Compositional Mean:", "\n")
      print((output$BetweenComp$mean), digits = digits)
    }
    if ("within" %in% level) {
      cat("\n", "Within-level Compositional Mean:", "\n")
      print((output$WithinComp$mean), digits = digits)
    }
  }
  if ("logratio" %in% class) {
    if ("total" %in% level) {
      cat("\n", "Isometric Log-ratio (Real) Mean:", "\n")
      print(as.data.table(output$TotalILR)[4,], row.names = FALSE, digits = digits)
    }
    if ("between" %in% level) {
      cat("\n", "Between-level Isometric Log-ratio (Real) Mean:", "\n")
      print(as.data.table(output$BetweenILR)[4,], row.names = FALSE, digits = digits)
    }
    if ("within" %in% level) {
      cat("\n", "Within-level Isometric Log-ratio (Real) Mean:", "\n")
      print(as.data.table(output$WithinILR)[4,], row.names = FALSE, digits = digits)
    }
  }
  
  ### Return output invisibly
  output <- lapply(output, function(X) {
    row.names(X) <- NULL
    return(X)
  })
  
  return(invisible(output))
}

#' Extract Compositional Data from \code{compilr} object.
#'
#' Extract amounts and compositions in conventional formats
#' as data.frames, matrices, or arrays.
#'
#' @inheritParams mean.compilr
#' @param row.names,optional Unused and only added for consistency with
#' the \code{\link[base:as.data.frame]{as.data.frame}} generic.
#'
#' @export
as.data.frame.compilr <- function(x, row.names = NULL, optional = TRUE,
                                  class = c("composition", "logratio"),
                                  level = c("between", "within", "total"),
                                  ...) {
  allout <- lapply(x[1:6], as.data.frame)
  output <- data.frame()
  
  ## Out
  if ("composition" %in% class) {
    if ("total" %in% level) {
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
    if ("total" %in% level) {
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

#' @rdname as.data.frame.compilr
#' @export
as.matrix.compilr <- function(x, 
                              class = c("composition", "logratio"),
                              level = c("between", "within", "total"),
                              ...) {
  as.matrix(as.data.frame(x, class = class, level = level, ...))
}

# ----------------- Extract Compositional Data -----------------
.get.compilr <- function(object, 
                         class = c("composition", "logratio"),
                         level = c("between", "within", "total"),
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
    bcomp <- cbind(object$data[, object$idvar, with = FALSE], object$BetweenComp)[!duplicated(get(object$idvar))]
    wcomp <- cbind(object$data[, object$idvar, with = FALSE], object$WithinComp)[!duplicated(get(object$idvar))]
    tcomp <- cbind(object$data[, object$idvar, with = FALSE], object$TotalComp)[!duplicated(get(object$idvar))]
    
    bilr <- cbind(object$data[, object$idvar, with = FALSE], object$BetweenILR)[!duplicated(get(object$idvar))]
    wilr <- cbind(object$data[, object$idvar, with = FALSE], object$WithinILR)[!duplicated(get(object$idvar))]
    tilr <- cbind(object$data[, object$idvar, with = FALSE], object$TotalILR)[!duplicated(get(object$idvar))]
    
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
      summary(object$BetweenComp, robust = TRUE),
      summary(object$WithinComp, robust = TRUE),
      summary(object$TotalComp, robust = TRUE),
      data.frame(summary(object$BetweenILR)),
      data.frame(summary(object$WithinILR)),
      data.frame(summary(object$TotalILR))
    )
  }
  
  names(output) <- c("BetweenComp", "WithinComp", "TotalComp", "BetweenILR", "WithinILR", "TotalILR")
  
  output
}