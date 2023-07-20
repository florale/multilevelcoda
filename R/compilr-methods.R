#' Several methods for \code{compilr} objects.
#' 
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
#' @param digits A integer value used for number formatting. Default is \code{3}.
#' 
#' @method mean compilr
#' @export
mean.compilr <- function(x, 
                         class = c("composition", "logratio"),
                         level = c("between", "within", "total"),
                         weight = c("equal", "proportional"),
                         digits = 3,
                         ...) {
  
  # Assemble
  output <- .get.compilr(x)
  
  ## Out
  if ("composition" %in% class) {
    if ("total" %in% level) {
      cat("\n", "Raw Compositional Mean:", "\n")
      print(clo(output$TotalComp), digits = digits)
    }
    if ("between" %in% level) {
      cat("\n", "Between-level Compositional Mean:", "\n")
      print(clo(output$BetweenComp), digits = digits)
    }
    if ("within" %in% level) {
      cat("\n", "Within-level Compositional Mean:", "\n")
      print(clo(output$WithinComp), digits = digits)
    }
  }
  if ("logratio" %in% class) {
    if ("total" %in% level) {
      cat("\n", "Isometric Log-ratio (Real) Mean:", "\n")
      print(as.data.table(t(output$TotalILR)), row.names = FALSE, digits = digits)
    }
    if ("between" %in% level) {
      cat("\n", "Between-level Isometric Log-ratio (Real) Mean:", "\n")
      print(as.data.table(t(output$BetweenILR)), row.names = FALSE, digits = digits)
    }
    if ("within" %in% level) {
      cat("\n", "Within-level Isometric Log-ratio (Real) Mean:", "\n")
      print(as.data.table(t(output$WithinILR)), row.names = FALSE, digits = digits)
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
#'
#' @export
as.data.frame.compilr <- function(x, 
                                  class = c("composition", "logratio"),
                                  level = c("between", "within", "total"),
                                  digits = 3,
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
                              digits = 3,
                              ...) {
  as.matrix(as.data.frame(x, class = class, level = level, digits = digits, ...))
}

# ----------------- Extract Compositional Data -----------------
.get.compilr <- function(x, 
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
  
  names(output) <- c("bcomp", "wcomp", "tcomp", "bilr", "wilr", "tilr")
  
  output
}