#' Within-person Substitution (from sample's compositional mean)
#'
#' Estimate the difference in outcomes
#' when compositional variables are substituted for a specific time period
#' at within-person level.
#'
#' @param data A list resulted from \code{brmcoda} that contains (1) results of a fitted  \code{brm} model and 
#' @param data (2) relevant inputs for substitution model.
#' @param substitute A data frame or data table indicating the possible substitution of variables. This dataset can be computed using \code{possub}. Required.
#' @param minute A integer or numeric value indicating the minute that compositional variable are substituted to/from. Default is 60L.
#' 
#' @return ## TODO
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo
#' @importFrom extraoperators %snin% %sin%
#' @export
#' @examples
#' ## TODO
#' wsubtest <- wsub(data = brmcodatest, substitute = posubtest, minute = 10)
wsub <- function(data, substitute, minute = 60) { 
  
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
  
  # Compute compositional mean
  mcomp <- mean(b)
  mcomp <- clo(mcomp, total = 1440)
  mcomp <- as.data.table(t(mcomp))
  names(mcomp) <- paste0("B", names(mcomp))
  
  # generate input for substitution model
  ID <- 1
  vn <- colnames(substitute) 
  min <- as.integer(paste0(minute))
  psi <- data$CompIlr$psi
  
  # List to store final output
  allout <- list()
  for(i in vn) {
    posub <- copy(substitute)
    posub <- as.data.table(posub)
    posub <- posub[(get(i) != 0)]
    posub <- posub[order(-rank(get(i)))]
    
    # Get substitution variable name for substitution model
    subvar <- colnames(posub) %snin% eval(i)
    iv <- i
    
    # lists to store results - TODO
    result <- NULL
    newcomp <- vector('list')
    subd <- vector("list")
    
    # substitution dataset
    for (j in 1:min) {
      sub <- posub * j
      for (k in 1:nrow(posub)) {
        newcomp <- mcomp + sub[k, ]
        names(newcomp) <- paste0(names(substitute))
        subd[[k]] <- cbind(mcomp, newcomp, sub[k, ][[i]])
      }
      result[[j]] <- do.call(rbind, subd)
    }
    subd <- as.data.table(do.call(rbind, result))
    
    # add names
    colnames(subd)[ncol(subd)] <- "MinSubstituted"
    subd[, Substitute := rep(subvar, length.out = nrow(subd))]
    subd$Predictor <- iv
    
    # ## remove impossible reallocation that result in negative values - TODO
    cols <- colnames(subd) %snin% c("MinSubstituted", "Substitute", "Predictor")
    
    noneg <- function(x){
      res <- ifelse(x < 0, NA, x)
      return(res)
    }
    
    subd[, (cols) := lapply(.SD, noneg), .SDcols = cols]
    subd <- subd[complete.cases(subd), ]
    
    ## add comp and ilr
    
    bn <- colnames(subd) %sin% names(mcomp)
    tn <- colnames(subd) %sin% names(newcomp)
    
    bcomp <- acomp(subd[, bn, with = FALSE]) 
    tcomp <- acomp(subd[, tn, with = FALSE])
    wcomp <- tcomp - bcomp
    
    bilr <- ilr(bcomp, V = psi) 
    wilr <- ilr(wcomp, V = psi) 
    
    names(bilr) <- c(paste0("bilr", 1:ncol(bilr)))
    names(wilr) <- c(paste0("wilr", 1:ncol(wilr)))
    
    subd <- cbind(subd, bilr, wilr)
    subd$ID <- ID
    
    ## dataset for no change
    woilr <- matrix(0, nrow = nrow(subd), ncol = ncol(wilr))
    woilr <- as.data.table(woilr)
    names(woilr) <- c(paste0("wilr", 1:ncol(wilr)))
    
    samed <- cbind(bilr, woilr)
    samed$ID <- ID
    
    # prediction
    ## substitution
    predsub <- as.data.table(fitted(data$Results, newdata = subd, re.form = NA, summary = FALSE))
    
    ## no change
    predsame <- as.data.table(fitted(data$Results, newdata = samed, re.form = NA, summary = FALSE))
    
    # difference between substitution and no change
    preddif <- predsub - predsame
    preddif <- as.data.table(describe_posterior(preddif, centrality = "mean",
                                                ci = 0.95, ci_method = "eti"))
    preddif <- preddif[, .(Mean, CI_low, CI_high)]
    
    # save results
    out <- do.call(cbind, preddif)
    out <- cbind(out, subd[, c("MinSubstituted", "Substitute", "Predictor")])
    out <- as.data.table(out)
    names(out) <- c("Mean", "CI_low", "CI_high", "MinSubstituted", "Substitute", "Predictor")
    
    allout[[i]] <- out
  }  
  return(allout)
}