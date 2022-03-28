#' Possible Pairwise Substitution 
#'
#' This function generate a dataset of all possible pairwise substitution of a composition.
#'
#' @param data A \code{data.frame} or \code{data.table}
#' containing the composition used in the analysis. Required.
#' @param composition A character vector specifying the names of compositional variables. Required.
#' 
#' @return A data table of all possible pairwise substitution.
#' @importFrom data.table as.data.table copy
#' @export
#' @examples
#' 
#' data(mcompd)
#' 
#' ps1 <- possub(data = mcompd, composition = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' ps2 <- possub(data = mcompd, composition = c("WAKE", "MVPA", "LPA", "SB"))
#' 
#' print(s2)
#' 
#' ## cleanup
#' rm(mcompd, ps1, ps2)
possub <- function(data, composition) {
  
  count <- length(composition)
  n <- count - 2

  tmp <- copy(data)

  subvars1 <- c(1, -1)
  subvars2 <- rep(0, n)
  subvars <- c(subvars1, subvars2)
  
  nc <- length(subvars)
  nr <- (nc - 1) * count
  
  possub <- matrix(0, nrow = nr, ncol = nc, dimnames = list(NULL, composition))
  k <- 0
  
  for(i in 1:nc)
    for(j in 1:nc)
      if(i != j) {
        k <- k + 1
        possub[k, c(i, j)] <- c(1, -1)
      }
  possub <- as.data.table(possub)
  
  return(possub)
}

