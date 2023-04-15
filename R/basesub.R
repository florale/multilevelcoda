#' Base Pairwise Substitution 
#'
#' @description
#' Makes a data set of all possible pairwise substitution of a composition which can be used as 
#' the base for substitution models.
#' @param parts A character vector specifying the names of compositional variables to be used.
#' 
#' @return A data table of all possible pairwise substitution.
#' @importFrom data.table as.data.table 
#' @export
#' @examples
#' 
#' data(mcompd)
#' 
#' ps <- basesub(parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' ps2 <- basesub(c("WAKE", "MVPA", "LPA", "SB"))
#' 
#' print(ps2)
#' 
#' ## cleanup
#' rm(mcompd, ps1, ps2)
basesub <- function(parts) {
  
  count <- length(parts)
  n <- count - 2

  subvars1 <- c(1, -1)
  subvars2 <- rep(0, n)
  subvars <- c(subvars1, subvars2)
  
  nc <- length(subvars)
  nr <- (nc - 1) * count
  
  basesub <- matrix(0, nrow = nr, ncol = nc, dimnames = list(NULL, parts))
  k <- 0
  
  for (i in 1:nc)
    for (j in 1:nc)
      if (i != j) {
        k <- k + 1
        basesub[k, c(i, j)] <- c(1, -1)
      }
  
  basesub <- as.data.table(basesub)
  basesub
  }

