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

#' Mean amounts and mean compositions presented in an \code{compilr} object.
#' 
#' @param x An object of class \code{compilr}.
#' @param class  Optional. Can be \code{"composition"} and/or \code{"logratio"} to
#' specify the geometry of the composition.
#' @param level  Optional. Can be \code{"between"}, \code{"within"}, and/or \code{total}
#' indicating the level of the geometry.
#' @param digits A integer value used for number formatting. Default is \coda{3}.
#' 
#' @method mean compilr
#' @export
mean.compilr <- function(x, 
                         class = c("composition", "logratio"),
                         level = c("between", "within", "total"),
                         digits = 3,
                         ...) {
  
  # Assemble
  output <- sapply(object[c("BetweenComp", "WithinComp", "TotalComp",
                            "BetweenILR", "WithinILR", "TotalILR")], function(y) mean(y, robust = TRUE))
  
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

