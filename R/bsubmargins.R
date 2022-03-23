#' Margianl Effects of Between-person Substitution
#'
#' Estimate the difference in outcomes
#' when compositional variables are substituted for a specific time period
#' at between-person level.
#'
#' @param data A fitted \code{brms} model data. Required.
#' @param substitute A data frame or data table indicating the possible substitution of variables. This dataset can be computed using \code{possub}. Required.
#' @param minute A integer or numeric value indicating the minute that compositional variable are substituted to/from. Default is 60L.
#' @return A list
#' @importFrom data.table as.data.table copy := 
#' @importFrom compositions acomp ilr clo
#' @importFrom extraoperators %snin% %sin%
#' @export
#' @examples
#' bsubmarginstest <- bsubmargins(data = brmcodatest, substitute = posubtest, minute = 10)
bsubmargins <- function (data, substitute, minute = 60) {

  if(isFALSE(missing(minute))) {
    if (isFALSE(is.integer(minute))) {
      if (isFALSE(is.numeric(minute))) {
        stop("'minute' must be an integer or a numeric value > 0.")
      }
    }
    
  } else if (isTRUE(missing(minute))) {
    minute <- 60L
  }
  
  # Compute between-person composition
  b <- data$CompIlr$BetweenComp
  b <- clo(b, total = 1440)
  b <- as.data.table(b)
  b <- unique(b)
  names(b) <- paste0("B", names(b))
  
  ID <- 1 # why re.form = NA but still needs this?
  min <- as.integer(minute)
  psi <- data$CompIlr$psi
  
  # generate possible substitution
  vn <- colnames(substitute) 
  
  # list to store final output
  allout <- list()
  
  for(i in vn) {
    posub <- copy(substitute)
    posub <- as.data.table(posub)
    posub <- posub[(get(i) != 0)]
    posub <- posub[order(-rank(get(i)))]
    
    # add substitution variable name
    subvar <- colnames(posub) %snin% eval(i)
    iv <- i
    
    # TODO
    result <- NULL
    newcomp <- vector('list')
    subd <- vector("list")
    marginsub <- NULL
    marginsuball <- NULL
    out <- NULL
    comp <- vector('list', length = nrow(b))
    sub <- NULL
    
    for (j in 1:min) {
      sub <- posub * j
      for (k in 1:nrow(sub)) {
        for (l in 1:nrow(b)) {
          newcomp <- b[l, ] + sub[k, ]
          names(newcomp) <- paste0(names(substitute))
          misub <- sub[k, get(i)]
          comp[[l]] <- cbind(b[l, ], newcomp, misub)
          # fix
        }
        subd <- as.data.table(do.call(rbind, comp))

        # names
        colnames(subd)[ncol(subd)] <- "MinSubstituted"
        subd[, Substitute := rep(subvar, length.out = nrow(subd))[k]]
        subd$Predictor <- iv
        subd$MinSubstituted <- as.numeric(subd$MinSubstituted)
        
        ## remove impossible reallocation that result in negative values - TODO
        cols <- colnames(subd) %snin% c("MinSubstituted", "Substitute", "Predictor")
        
        noneg <- function(x){
          res <- ifelse(x < 0, NA, x)
          return(res)
        }
        
        subd[, (cols) := lapply(.SD, noneg), .SDcols = cols]
        subd <- subd[complete.cases(subd), ]
        
        ## add comp and ilr
        bn <- colnames(subd) %sin% names(b)
        tn <- colnames(subd) %sin% names(newcomp)
        
        bcomp <- acomp(subd[, bn, with = FALSE]) 
        tcomp <- acomp(subd[, tn, with = FALSE])
        
        bilr <- ilr(bcomp, V = psi) 
        tilr <- ilr(tcomp, V = psi) 
        wilr <- matrix(0, nrow = nrow(subd), ncol = ncol(bilr))
        wilr <- as.data.table(wilr)
        
        names(bilr) <- c(paste0("bilr", 1:ncol(bilr)))
        names(tilr) <- c(paste0("bilr", 1:ncol(tilr)))
        names(wilr) <- c(paste0("wilr", 1:ncol(wilr)))
        
        ## substitution dataset
        subd <- cbind(subd, tilr, wilr)
        subd$ID <- ID
        
        ## no change dataset
        samed <- cbind(bilr, wilr)
        samed$ID <- ID
        
        # prediction
        ## substitution
        predsub <- as.data.table(fitted(data$Results, newdata = subd, re.form = NA, summary = FALSE))
        
        ## no change
        predsame <- as.data.table(fitted(data$Results, newdata = samed, re.form = NA, summary = FALSE))

        # calculate difference between substitution and no change
        preddif <- predsub - predsame
        preddif <- rowMeans(preddif)
        preddif <- as.data.table(describe_posterior(preddif, centrality = "mean", ci = 0.95, ci_method = "eti"))
        preddif <- preddif[, .(Mean, CI_low, CI_high)]
        preddif$MinSubstituted <- mean(subd[, MinSubstituted])
        marginsub[[k]] <- preddif
      }
      # save results
      marginsuball[[j]] <- do.call(rbind, marginsub)
    }
    
    out <- as.data.table(do.call(rbind, marginsuball))
    out[, Substitute := rep(subvar, length.out = nrow(out))]
    out$Predictor <- iv
    names(out) <- c("Mean", "CI_low", "CI_high", "MinSubstituted", "Substitute", "Predictor")
    
    allout[[i]] <- out
  }
  
  return(allout)
}