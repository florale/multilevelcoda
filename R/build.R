#' Build Base Pairwise Substitution 
#'
#' @description
#' Make a data set of all possible pairwise substitution of a composition which can be used as 
#' the base for substitution models.
#' @param parts A character vector specifying the names of compositional variables to be used.
#' @param comparison Either \code{"one-to-one"} or \code{"one-to-all"}. Default is \code{"one-to-one"}.
#' 
#' @return A data table of all possible pairwise substitution.
#' @importFrom data.table as.data.table 
#' @export
#' @examples
#' ps1 <- build.basesub(parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' print(ps1)
#' 
#' ps2 <- build.basesub(c("WAKE", "MVPA", "LPA", "SB"), comparison = "one-to-all")
#' print(ps2)
build.basesub <- function(parts, comparison = NULL) {
  
  if (is.null(comparison)) {
    comparison <- "one-to-one"
  }
  
  d <- length(parts)
  n <- d - 2
  
  subvar1 <- c(1, -1)
  subvar2 <- rep(0, n)
  subvar <- c(subvar1, subvar2)
  
  nc <- length(subvar)
  nr <- (nc - 1) * d
  k <- 0
  
  if (comparison == "one-to-all") {
    base_sub_to <- matrix(0, nrow = d, ncol = d, dimnames = list(NULL, parts))
    for (i in 1:nc)
      for (j in 1:nc)
        if (i == j) {
          base_sub_to[i, j] <- 1
          base_sub_to[i, -j] <- -(1/(d-1))
        }
    
    base_sub_from <- matrix(0, nrow = d, ncol = d, dimnames = list(NULL, parts))
    for (i in 1:nc)
      for (j in 1:nc)
        if (i == j) {
          base_sub_from[i, j] <- -1
          base_sub_from[i, -j] <- (1/(d-1))
        }
    
    base_sub <- rbind(base_sub_to, base_sub_from)
    
  } else {
    base_sub <- matrix(0, nrow = nr, ncol = nc, dimnames = list(NULL, parts))
    
    for (i in 1:nc)
      for (j in 1:nc)
        if (i != j) {
          k <- k + 1
          base_sub[k, c(i, j)] <- c(1, -1)
        }
  }
  
  base_sub <- as.data.table(base_sub)
  base_sub
}

#' Build Sequential Binary Partition
#'
#' @description
#' Build a default sequential binary partition for \code{complr} object.
#' The default sequential binary partition is a pivot balance that allows 
#' the effect of this first balance coordinate to be interpreted as the change 
#' in the prediction for the dependent variable 
#' when that given part increases while all remaining parts decrease by a common proportion.
#' @param parts A character vector specifying the names of compositional variables to be used.
#' 
#' @return A matrix sequential binary partition.
#' @importFrom data.table as.data.table 
#' @export
#' @examples
#' sbp1 <- build.sbp(parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' print(sbp1)
#' 
#' sbp2 <- build.sbp(c("WAKE", "MVPA", "LPA", "SB"))
#' print(sbp2)
build.sbp <- function(parts) {
  
  d <- length(parts)
  k <- 0
  
  nc <- d
  nr <- d - 1
  
  base_sbp <- matrix(NA, nrow = nr, ncol = nc, dimnames = list(NULL, parts))
  sbp <- base_sbp
  
  for (i in 1:nr) {
    base_sbp[i,  i] <- 1
    base_sbp[i, -i] <- -1
    base_sbp[i, 0:(i-1)] <- 0
  }
  base_sbp
}

