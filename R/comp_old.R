## # possible within 24-hour real mean total composition 
## poss_vals <- c(1, -1, 0, 0, 0)
## nc <- length(poss_vals)
## K <- factorial(nc)/(nc-length(unique(poss_vals))+1)
## pcomp <- matrix(0,nrow = K, ncol = nc, dimnames = list(NULL,names(t)))
## k <- 0

## for(i in 1:nc) for(j in 1:nc) if(i!=j) {
##   k <- k + 1
##   pcomp[k, c(i,j)] <- c(1,-1)
## }
## pcomp <- as.data.table(pcomp[1:20,])

## save(pcomp, file = "/Users/florale/Documents/GitHub/multilevelcoda/data/pcomp.RData")
