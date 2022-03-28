## make Rcmd check happy
utils::globalVariables(c("Mean",  "CI_low", "CI_high", "Substitute", "MinSubstituted"))
#' @title Between-person Substitution Model (from sample's compositional mean)
#'
#' @description 
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
#' @return A list
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo mean.acomp
#' @importFrom extraoperators %snin% %sin%
#' @importFrom insight find_predictors
#' @importFrom emmeans ref_grid
#' @export
#' @examples
#'
#' data(mcompd)
#' data(sbp)
#' ps <- possub(data = mcompd, composition = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#'
#' bsubctest2 <- bsubc(data = mcmc, substitute = ps, minute = 10)
#'
#' ## cleanup
#' rm(bsubtest, mcompd, sbp, ps)
bsubc <- function(data, substitute, minute = 60L) {

  if (isFALSE(missing(minute))) {
    if (isFALSE(is.integer(minute))) {
      if (isFALSE(is.numeric(minute))) {
        stop(sprintf("'minute' must be an integer or a numeric value > 0."))
      }
    }

  } else {
    minute <- 60L
  }

  if (isFALSE(identical(ncol(substitute), length(data$CompIlr$composition)))) {
    stop(sprintf("The number of columns in 'substitute' (%d) must be the
  same as the compositional variables in 'composition' (%d).",
                 ncol(substitute),
                 length(data$CompIlr$composition)))
  }

  if (isFALSE(identical(colnames(substitute), data$CompIlr$composition))) {
    stop(sprintf("The names of compositional variables must be the same
  in 'substitute' (%s) and 'composition' (%s).",
                 colnames(substitute),
                 data$CompIlr$composition))
  }
  
  tmp <- copy(data)

  # Compute compositional mean
  mcomp <- .get.mcomp(tmp)
  
  # input for substitution model
  ID <- 1
  min <- as.integer(minute)
  psi <- tmp$CompIlr$psi
  
  iout <- vector("list")
  
  # deal with covariates
  ilrn <- c(names(tmp$CompIlr$BetweenILR), names(tmp$CompIlr$WithinILR))
  vn <- do.call(rbind, find_predictors(tmp$BrmModel))
  
  if (isTRUE(identical(length(vn), length(ilrn)))) {
    
    for(i in colnames(substitute)) {
      posub <- copy(substitute)
      posub <- as.data.table(posub)
      posub <- posub[(get(i) != 0)]
      posub <- posub[order(-rank(get(i)))]
      
      subvar <- colnames(posub) %snin% eval(i) # substitute compositional variables
      iv <- i # central compositional variable
      
      # lists to store results - TODO
      kout <- vector("list", length = nrow(posub)) # a list to store all possible pairwise substitution for 1 specific period
      newd <- vector("list", length = min) # # a list to store all possible pairwise substitution for all period as specified by user 
      
      # Generate substitution dataset for predictions
      for (j in seq_len(min)) { # time level
        sub <- posub * j
        for (k in seq_len(nrow(sub))) { # substitution level
          newcomp <- mcomp + sub[k, ]
          names(newcomp) <- colnames(substitute)

          kout[[k]] <- cbind(mcomp, newcomp, sub[k, ][[i]])
        }
        newd[[j]] <- do.call(rbind, kout)
      }
      # Here is a dataset of all possible substitution from sample's mean across all periods
      newd <- as.data.table(do.call(rbind, newd)) 
      
      # Add more useful information to the final output
      colnames(newd)[ncol(newd)] <- "MinSubstituted"
      newd[, Substitute := rep(subvar, length.out = nrow(newd))]
      newd$Predictor <- iv
      
      # ## remove impossible reallocation that result in negative values - TODO
      cols <- colnames(newd) %snin% c("MinSubstituted", "Substitute", "Predictor")
      
      noneg <- function(x){
        res <- ifelse(x < 0, NA, x)
        return(res)
      }
      
      newd[, (cols) := lapply(.SD, noneg), .SDcols = cols]
      newd <- newd[complete.cases(newd), ]
      
      # Add comp and ilr for predictions      
      bcomp <- acomp(newd[, colnames(tmp$CompIlr$BetweenComp), with = FALSE]) 
      tcomp <- acomp(newd[, tmp$CompIlr$composition, with = FALSE])
      
      bilr <- ilr(bcomp, V = psi) 
      tilr <- ilr(tcomp, V = psi) 
      wilr <- matrix(0, nrow = nrow(newd), ncol = ncol(bilr))
      wilr <- as.data.table(wilr)
      
      colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
      colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
      colnames(tilr) <- paste0("bilr", seq_len(ncol(tilr)))
      
      ## Substitution dataset
      subd <- cbind(newd, tilr, wilr, ID)

      ## No change dataset
      samed <- cbind(newd, bilr, wilr, ID)

      # Generate prediction
      ysub <- as.data.table(fitted(tmp$BrmModel, newdata = subd, re.form = NA, summary = FALSE))
      ysame <- as.data.table(fitted(tmp$BrmModel, newdata = samed, re.form = NA, summary = FALSE))
      
      # Difference between substitution and no change
      yd <- ysub - ysame
      yd <- as.data.table(describe_posterior(yd, centrality = "mean",
                                              ci = 0.95, ci_method = "eti"))
      yd <- yd[, .(Mean, CI_low, CI_high)]
      
      # Save results
      result <- do.call(cbind, yd)
      result <- cbind(result, newd[, c("MinSubstituted", "Substitute", "Predictor")])
      result <- as.data.table(result)
      names(result) <- c("Mean", "CI_low", "CI_high", "MinSubstituted", "Substitute", "Predictor")
      
      # Final results for entire composition
      iout[[i]] <- result
    }
    return(iout)
    
    } else {
      
      refg <- as.data.table(ref_grid(tmp$BrmModel) @grid)
      cv <- colnames(refg) %snin% c(ilrn, ".wgt.")
      
      refg <- refg[, cv, with = FALSE] # reference grid for covariates
      
      for (i in colnames(substitute)) {
        posub <- copy(substitute)
        posub <- as.data.table(posub)
        posub <- posub[(get(i) != 0)]
        posub <- posub[order(-rank(get(i)))]
        
        subvar <- colnames(posub) %snin% eval(i)
        iv <- i 
        
        # lists to store results - TODO
        nd <- NULL
        newd <- vector("list")
        
        # Substitution dataset
        for (j in seq_len(min)) {
          sub <- posub * j
          for (k in seq_len(nrow(sub))) {
            newcomp <- mcomp + sub[k, ]
            names(newcomp) <- paste0(names(substitute))
            nd[[k]] <- cbind(mcomp, newcomp, sub[k, ][[i]])
            }
          newd[[j]] <- do.call(rbind, nd)
          }
        newd <- as.data.table(do.call(rbind, newd))
        
        # Add more useful information to the final output
        colnames(newd)[ncol(newd)] <- "MinSubstituted"
        newd[, Substitute := rep(subvar, length.out = nrow(newd))]
        newd$Predictor <- iv
        
        # Remove impossible reallocation that result in negative values - TODO
        cols <- colnames(newd) %snin% c("MinSubstituted", "Substitute", "Predictor")
        
        noneg <- function(x){
          res <- ifelse(x < 0, NA, x)
          return(res)
          }
        
        newd[, (cols) := lapply(.SD, noneg), .SDcols = cols]
        newd <- newd[complete.cases(newd), ]
        
        # Now add composition and ilr for prediction
        bcomp <- acomp(newd[, colnames(tmp$CompIlr$BetweenComp), with = FALSE]) 
        tcomp <- acomp(newd[, tmp$CompIlr$composition, with = FALSE])
        
        bilr <- ilr(bcomp, V = psi)
        tilr <- ilr(tcomp, V = psi)
        wilr <- matrix(0, nrow = nrow(newd), ncol = ncol(bilr))
        wilr <- as.data.table(wilr)
        
        colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
        colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
        colnames(tilr) <- paste0("bilr", seq_len(ncol(tilr)))
        
        ## Loop through covariates to get AME
        mout <- vector("list", length = nrow(refg))
        
        for(m in seq_len(nrow(refg))) {
          ## Substitution dataset
          subd <- cbind(newd, tilr, wilr, ID, refg[m, ])
          
          ## No change dataset
          samed <- cbind(newd, bilr, wilr, ID, refg[m, ])
          
          # Generate prediction
          ysub <- as.data.table(fitted(tmp$BrmModel, newdata = subd, re.form = NA, summary = FALSE))
          ysame <- as.data.table(fitted(tmp$BrmModel, newdata = samed, re.form = NA, summary = FALSE))
          
          # Difference between substitution and no change
          yd <- ysub - ysame
          if (nrow(refg) == 1) {
            myd <- yd
            } else {
              mout[[m]] <- yd
              myd <- Reduce(`+`, mout) / length(mout)
            }
          }
        # AME across covariates and posteriors
        myd <- as.data.table(describe_posterior(myd, centrality = "mean",
                                                ci = 0.95, ci_method = "eti"))
        myd <- myd[, .(Mean, CI_low, CI_high)]
        
        # Save results
        result <- do.call(cbind, myd)
        result <- cbind(result, newd[, c("MinSubstituted", "Substitute", "Predictor")])
        result <- as.data.table(result)
        names(result) <- c("Mean", "CI_low", "CI_high", "MinSubstituted", "Substitute", "Predictor")
        
        # Finally, entire composition
        iout[[i]] <- result
        }
      return(iout)
    }
  }