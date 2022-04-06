## make Rcmd check happy
utils::globalVariables(c("Mean",  "CI_low", "CI_high", "Substitute", "MinSubstituted"))

#' Between-person Substitution Model (from sample's compositional mean)
#'
#' Estimate the difference in outcomes
#' when compositional variables are substituted for a specific time period.
#' at between-person level.
#'
#' @param data A \code{brmcoda} object that contains results of a fitted \code{brm} model and 
#' relevant inputs of substitution model from \code{compilr}.
#' @param substitute A \code{data.frame} or \code{data.table} of the possible substitution of variables. 
#' This dataset can be computed using function \code{possub}. Required.
#' @param minute A integer or numeric value specifying the minute that 
#' compositional variable are substituted to/from. Default is \code{60L}.
#' 
#' @return A list of results from substitution models for all compositional variables.
#' @importFrom data.table as.data.table copy := 
#' @importFrom compositions acomp ilr clo mean.acomp
#' @importFrom extraoperators %snin% %sin%
#' @importFrom insight find_predictors
#' @export
#' @examples
#' ## TODO
#' 
#' 
#' data(mcompd)
#' data(sbp)
#' ps <- possub(data = mcompd, composition = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' 
#' bsubtest <- bsub(data = mcm, substitute = ps, minute = 10)
#' 
#' ## cleanup
#' rm(bsubtest, mcompd)
bsub <- function(data, substitute, minute = 60L) { 
  
  if(isFALSE(missing(minute))) {
    if (isFALSE(is.integer(minute))) {
      if (isFALSE(is.numeric(minute))) {
        stop(sprintf("'minute' must be an integer or a numeric value > 0."))
      }
    }

  } else if (isTRUE(missing(minute))) {
    minute <- 60L
  }
  
  if (isFALSE(identical(ncol(substitute), length(data$CompIlr$composition)))) {
    stop(sprintf("The number of columns in substitute (%d) must be the 
                 same as the compositional variables in composition (%d).",
                 ncol(substitute),
                 length(data$CompIlr$composition)))
  }
  
  if (isFALSE(identical(colnames(substitute), data$CompIlr$composition))) {
    stop(sprintf("The names of compositional variables must be the same
                 in substitute (%s) and composition (%s).",
                 colnames(substitute),
                 data$CompIlr$composition))
  }
  
  tmp <- copy(data)
  
  # Compute compositional mean
  b <- tmp$CompIlr$BetweenComp
  
  mcomp <- mean(b, robust = TRUE)
  mcomp <- clo(mcomp, total = 1440)
  mcomp <- as.data.table(t(mcomp))
  names(mcomp) <- paste0("B", names(mcomp))

  # generate input for substitution model
  ID <- 1
  min <- as.integer(minute)
  psi <- tmp$CompIlr$psi

  # List to store final output
  allout <- list()
  for(i in colnames(substitute)) {
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
    newd <- vector("list")
    
    # substitution dataset
    for (j in seq_len(min)) {
      sub <- posub * j
      for (k in seq_len(nrow(sub))) {
        newcomp <- mcomp + sub[k, ]
        names(newcomp) <- paste0(names(substitute))
        newd[[k]] <- cbind(mcomp, newcomp, sub[k, ][[i]])
      }
      result[[j]] <- do.call(rbind, newd)
    }
    newd <- as.data.table(do.call(rbind, result))
    
    # add names
    colnames(subd)[ncol(subd)] <- "MinSubstituted"
    subd[, Substitute := rep(subvar, length.out = nrow(subd))]
    subd$Predictor <- iv

    ## remove impossible reallocation that result in negative values
    cols <- colnames(subd) %snin% c("MinSubstituted", "Substitute", "Predictor")
    subd <- subd[rowSums(subd[, ..cols] < 0) == 0]

    ## add comp and ilr
    
    bn <- colnames(newd) %sin% names(mcomp)
    tn <- colnames(newd) %sin% names(newcomp)
    
    bcomp <- acomp(newd[, bn, with = FALSE]) 
    tcomp <- acomp(newd[, tn, with = FALSE])
    
    bilr <- ilr(bcomp, V = psi)
    tilr <- ilr(tcomp, V = psi)
    wilr <- matrix(0, nrow = nrow(newd), ncol = ncol(bilr))
    wilr <- as.data.table(wilr)
    
    colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
    colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
    colnames(tilr) <- paste0("bilr", seq_len(ncol(tilr)))
    
    ## substitution dataset
    subd <- cbind(newd, tilr, wilr, ID)

    ## no change dataset
    samed <- cbind(newd, bilr, wilr, ID)

    # prediction
    ## substitution
    predsub <- as.data.table(fitted(tmp$BrmModel, newdata = subd, re.form = NA, summary = FALSE))
    
    ## no change
    predsame <- as.data.table(fitted(tmp$BrmModel, newdata = samed, re.form = NA, summary = FALSE))
    
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

    ## final results for entire composition
    allout[[i]] <- out
  }
  
  return(allout)
}
