#' Example Function Title
#'
#' This provides a few sentence description about the example function.
#'
#' @param x A vector.
#' @return The vector coerced to numeric vector in a data table.
#' @importFrom data.table data.table
#' @export
#' @examples
#'
#' examplefun(c(1L, 2L))
#'
#' examplefun(TRUE)
#'
#' examplefun(c(5.1, 2.2))
examplefun <- function(x) {
  if (isFALSE(is.numeric(x)) && isFALSE(is.integer(x)) && isFALSE(is.logical(x))) {
    stop("error x must be a numeric, integer or logical vector")
  }

  data.table(x = as.numeric(x))
}
