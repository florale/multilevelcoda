comp <- function(count, names) {  # perhaps default names to be sleep-wake behaviours
   n <- count - 2

   subvars1 <- c(1, -1)
   subvars2 <- rep(0, n)
   subvars <- c(subvars1, subvars2)
   
   nc <- length(subvars)
   K <- (nc-1)*count
  
  comp <- matrix(0, nrow = K, ncol = nc, , dimnames = list(NULL, names))
  k <- 0
  
  for(i in 1:nc) 
  for(j in 1:nc) 
  if(i!=j) {
     k <- k + 1
     comp[k, c(i,j)] <- c(1,-1)
     }
     return(comp)
}
