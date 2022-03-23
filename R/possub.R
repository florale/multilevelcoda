#' Possible Pairwise Substitution 
#'
#' This function generate a dataset of all possible pairwise substitution of a composition.
#'
#' @param count A numeric value indicating the number of variables within the composition.
#' @param data A dataset of composition and IDs. Ideally ID is the last column. Required.
#' @param idvar A string character indicating the name of the ID variable. Default is ID.
#' 
#' @return A data table of all possible pairwise substitution.
#' @importFrom data.table as.data.table copy :=
#' @importFrom extraoperators %snin%
#' @export
#' @examples
#' 
#' data(mcompd)
#' 
#' posubtest <- possub(count = 5, data = mcompd[, 1:6], idvar = "ID")
#' 
#' str(posubtest)
#' 
#' ## cleanup
#' rm(mcompd, posubtest)
possub <- function(count, data, idvar) {
  
  # add check for 'count' and number of cols in 'data'
  count <- as.numeric(count)
  n <- count - 2

  d <- copy(data)
  vn <- colnames(d) %snin% idvar
  
  subvars1 <- c(1, -1)
  subvars2 <- rep(0, n)
  subvars <- c(subvars1, subvars2)
  
  nc <- length(subvars)
  nr <- (nc - 1) * count
  
  possub <- matrix(0, nrow = nr, ncol = nc, dimnames = list(NULL, vn))
  k <- 0
  
  for(i in 1:nc)
    for(j in 1:nc)
      if(i != j) {
        k <- k + 1
        possub[k, c(i, j)] <- c(1, -1)
      }
  possub <- as.data.table(possub)
}
