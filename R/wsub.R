#' Within-person Substitution (from sample's compositional mean)
#'
#' Estimate the difference in outcomes
#' when compositional variables are substituted for a specific time period
#' at within-person level.
#'
#' @param object A fitted \code{brms} model object. Required.
#' @param data A dataset of composition plus a variable containing IDs that was used for the brms model object. Required.
#' @param stutitute A data frame or data table indicating the possible substitution of variables. This dataset can be computed using \code{possub}. Required.
#' @param sbp A signary matrix indicating sequential binary partition. Required.
#' @param minute A integer or numeric value indicating the minute that compositional variable are substituted to/from. Default is 60L.
#' @param idvar A character string indicating the name of the variable containing IDs.
#' @param cov 
#' 
#' @return TODO
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo gsi.buildilrBase
#' @importFrom extraoperators %snin% %sin%
#' @export
#' @examples
#' ## TODO
#' wsubtest <- wsub(object = m, data = mcompd[, 1:6], substitute = posubtest, minute = 10, sbp = sbp)
wsub <- function(object, data, substitute, sbp, minute = 60, idvar = "ID") { # if use compilr do we need to add all argument for it here
  
  if(isFALSE(missing(minute))) {
    if (isFALSE(is.integer(minute))) {
      if (isFALSE(is.numeric(minute))) {
        stop("'minute' must be an integer or a numeric value > 0.")
      }
    }
    
  } else if (isTRUE(missing(minute))) {
    minute <- 60L
  }
  
  # compute between-person composition
  b <- compilr(data = data, sbp = sbp)[[1]]
  psi <- gsi.buildilrBase(t(sbp))
  
  #compute mean composition
  mcomp <- mean(b)
  mcomp <- clo(mcomp, total = 1440)
  mcomp <- as.data.table(t(mcomp))
  names(mcomp) <- paste0("B", names(mcomp))
  
  ID <- 1
  min <- as.integer(paste0(minute))
  
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
    
    # lists to store results - TODO
    result <- NULL
    newcomp <- vector('list')
    subd <- vector("list")
    
    for (j in 1:min) {
      sub <- posub * j
      for (k in 1:nrow(posub)) {
        newcomp <- mcomp + sub[k, ]
        names(newcomp) <- paste0(names(substitute))
        subd[[k]] <- cbind(mcomp, newcomp, sub[k, ][[i]])
      }
      result[[j]] <- do.call(rbind, subd)
    }
    subd <- as.data.table(do.call(rbind,result))
    
    # names
    colnames(subd)[ncol(subd)] <- "MinSubstituted"
    subd[, Substitute := rep(subvar, length.out = nrow(subd))]
    subd$Predictor <- iv
    
    ## remove impossible reallocation that result in negative values - TODO
    subd <- subd %>%
      filter(if_all(-c(Substitute, Minsubvard), ~ . > 0))
    
    ## add comp and ilr
    bn <- colnames(subd) %sin% names(mcomp)
    tn <- colnames(subd) %sin% names(newcomp)
    
    bcomp <- acomp(subd[, bn, with = FALSE]) 
    tcomp <- acomp(subd[, tn, with = FALSE])
    wcomp <- tcomp - bcomp
    
    bilr <- ilr(bcomp, V = psi) 
    wilr <- ilr(wcomp, V = psi) 
    
    subd$bilr1 <- bilr[, 1]
    subd$bilr2 <- bilr[, 2]
    subd$bilr3 <- bilr[, 3]
    subd$bilr4 <- bilr[, 4]
    
    subd$wilr1 <- wilr[, 1]
    subd$wilr2 <- wilr[, 2]
    subd$wilr3 <- wilr[, 3]
    subd$wilr4 <- wilr[, 4]
    
    subd$ID <- ID
    
    ## dataset for no change
    samed <- data.table(bilr1 = bilr[, 1],
                        bilr2 = bilr[, 2],
                        bilr3 = bilr[, 3],
                        bilr4 = bilr[, 4],
                        wilr1 = 0,
                        wilr2 = 0,
                        wilr3 = 0,
                        wilr4 = 0)
    samed$ID <- ID
    
    # prediction
    ## substitution
    predsub <- as.data.table(fitted(object, newdata = subd, re.form = NA, summary = FALSE))
    
    ## no change
    predsame <- as.data.table(fitted(object, newdata = samed, re.form = NA, summary = FALSE))
    
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