#' Margianl Effects of Between-person Substitution
#'
#' Estimate the difference in outcomes
#' when compositional variables are substituted for a specific time period
#' at between-person level.
#'
#' @param object A fitted \code{brms} model object. Required.
#' @param data A dataset of composition plus a variable containing IDs that was used for the brms model object. Required.
#' @param stutitute A data frame or data table indicating the possible substitution of variables. This dataset can be computed using \code{possub}. Required.
#' @param sbp A signary matrix indicating sequential binary partition. Required.
#' @param minute A integer or numeric value indicating the minute that compositional variable are substituted to/from. Default is 60L.
#' @param idvar A character string indicating the name of the variable containing IDs.
#' 
#' @return
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo gsi.buildilrBase
#' @importFrom extraoperators %snin% %sin%
#' @export
#' @examples
#' 
#' bsubmargintest <- bsubmargin(object = m, data = mcompd[, 1:6], stutitute = posubtest, minute = 10, sbp = sbp)
bsubmargin <- function (object, data, substitute, sbp, minute = 60, idvar = "ID") {
  
  if (isTRUE(missing(object))) { ## J to add 
    stop(paste("TODO"
    ))
  }

  if (isTRUE(missing(data))) {
    stop(paste(" 'data' is a required argument and cannot be missing;",
               "it should be dataset of composition and IDs",
               sep = "\n"))
  }

  if (isFALSE(missing(at))) {
    if (isFALSE(is.integer(at) && all(at > 0))) {
      stop("'at' must be integer values > 0.")
    }
  } else if (isTRUE(missing(at))) {
    at <- 1:60L
  }

  if (isFALSE(missing(names))) {
    if (isFALSE(is.character(names))) {
      stop(paste("'names' must be a character string representing the names of the following:",
                 "between-person compostion, total composition, at.",
                 sep = "\n"))
    }
  } else if (isTRUE(missing(names))) {
    names <- c("BTST", "BAWAKE", "BMVPA", "BLPA", "BSB", "TST", "AWAKE", "MVPA", "LPA", "SB", "at")
  }

  if (isTRUE(missing(substitute))) {
    substitute <- names[7:10]
  }

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
  b <- clo(b, total = 1440)
  psi <- gsi.buildilrBase(t(sbp))
  
  # 1 between-person composition per participant
  b <- as.data.table(b)
  b <- unique(b) 
  names(b) <- paste0("B", names(b))
  
  ID <- 1 # why re.form = NA but still needs this?
  min <- as.integer(paste0(minute))
  
  # generate possible substitution
  vn <- colnames(stutitute) 
  
  # list to store final output
  allout <- list()
  
  for(i in vn) {
    posub <- copy(stutitute)
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
          names(newcomp) <- paste0(names(stutitute))
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
        subd <- subd %>%
          filter(if_all(-c(Substitute, MinSubstituted), ~ . > 0))
        
        ## add comp and ilr
        bn <- colnames(subd) %sin% names(b)
        tn <- colnames(subd) %sin% names(newcomp)
        
        bcomp <- acomp(subd[, bn, with = FALSE]) 
        tcomp <- acomp(subd[, tn, with = FALSE])
        
        bilr <- ilr(bcomp, V = psi) 
        tilr <- ilr(tcomp, V = psi) 
        
        subd$bilr1 <- tilr[, 1]
        subd$bilr2 <- tilr[, 2]
        subd$bilr3 <- tilr[, 3]
        subd$bilr4 <- tilr[, 4]
        
        subd$wilr1 <- 0
        subd$wilr2 <- 0
        subd$wilr3 <- 0
        subd$wilr4 <- 0
        
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
        predsub <- as.data.table(fitted(m, newdata = subd, re.form = NA, summary = FALSE))
        
        ## no change
        predsame <- as.data.table(fitted(m, newdata = samed, re.form = NA, summary = FALSE))
        
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