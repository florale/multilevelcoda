#' # Functions used only internally
#' @keywords internal
#' @importFrom data.table as.data.table
#' @importFrom compositions clo
#' @noRd
#' ## Remove negative values
.noneg <- function(x){
   res <- ifelse(x < 0, NA, x)
   return(res)
}

#' ## Support Substitution Model
#' # Get substitution data
.adj.subm <- function (sub, b, refg, subvar, iv, tmp, psi, substitute, i, j, k, ID = "1") {
  lout <- vector("list", length = nrow(refg))
  
     for(l in seq_len(nrow(refg))) {
              # participant level
            subk <- sub[k, ]
            subk <- subk[rep(seq_len(nrow(subk)), nrow(b)), ]
            newcomp <- b + subk
            names(newcomp) <- colnames(substitute)

            newd <- cbind(b, newcomp, refg[l, ], sub[k, get(i)])
          
          # Here is a dataset 1 possible substitution per 1 min for the entire sample
          # nrow = nrow(b), ncol = 10

          # Add more useful information to the final output
          colnames(newd)[ncol(newd)] <- "MinSubstituted"
          newd[, Substitute := rep(subvar, length.out = nrow(newd))[k]]
          newd$Predictor <- iv
          newd$MinSubstituted <- as.numeric(newd$MinSubstituted)
          
          ## remove impossible reallocation that result in negative values - TODO
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
          colnames(tilr) <- paste0("bilr", seq_len(ncol(tilr)))
          colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
          
          ## Substitution dataset
            subd <- cbind(newd, tilr, wilr, ID)
            
            ## No change dataset
            samed <- cbind(newd, bilr, wilr, ID)
            
            # Generate prediction
            ysub <- fitted(tmp$BrmModel, newdata = subd, re.form = NA, summary = FALSE)
            ysame <- fitted(tmp$BrmModel, newdata = samed, re.form = NA, summary = FALSE)
            
            # Difference between substitution and no change
            yd <- ysub - ysame
            yd <- rowMeans(yd) # ame of 1 substiution at participant level
            # nrow(yd) = posterior draws, ncol = 1
            lout[[l]] <- yd
          }
}

#' # Get possible substitution
.posub.data <- function(substitute, i) {
  posub <- copy(substitute)
  posub <- as.data.table(posub)
  posub <- posub[(get(i) != 0)]
  posub <- posub[order(-rank(get(i)))]

  posub
}

#' # Get compositional mean
.mcomp <- function(data) {
  
  b <- data$CompIlr$BetweenComp
  mcomp <- mean(b, robust = TRUE)
  mcomp <- clo(mcomp, total = 1440)
  mcomp <- as.data.table(t(mcomp))
  colnames(mcomp) <- colnames(b)
  
  mcomp
}

#' # Unadjusted  substitution model
.unadj.subm <- function(substitute, b, tmp, 
                        min, psi, ID = 1) {
  iout <- vector("list")
  for(i in colnames(substitute)) { # compostion level
    posub <- copy(substitute)
    posub <- as.data.table(posub)
    posub <- posub[(get(i) != 0)]
    posub <- posub[order(-rank(get(i)))]
    
  # add substitution variable name
  subvar <- colnames(posub) %snin% eval(i) # substitute compositional variables
  iv <- i # central compositional variable
  
  # lists to store results - TODO
  lout <- vector("list", length = nrow(b)) # a list to store a pairwise substitution for 1 specific period of all participant
  kout <- vector("list", length = nrow(substitute))
  jout <- vector("list", length = min)
  
  for (j in seq_len(min)) { # time level
    sub <- posub * j
    for (k in seq_len(nrow(sub))) { # substitution level
      subk <- sub[k, ]
      subk <- subk[rep(seq_len(nrow(subk)), nrow(b)), ]
      newcomp <- b + subk
      names(newcomp) <- colnames(substitute)
            
      newd <- cbind(b, newcomp, sub[k, get(i)])
      
      # Add more useful information to the final output
      colnames(newd)[ncol(newd)] <- "MinSubstituted"
      newd[, Substitute := rep(subvar, length.out = nrow(newd))[k]]
      newd$Predictor <- iv
      newd$MinSubstituted <- as.numeric(newd$MinSubstituted)
      
      ## remove impossible reallocation that result in negative values - TODO
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
      colnames(tilr) <- paste0("bilr", seq_len(ncol(tilr)))
      colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
      
      ## Substitution dataset
      subd <- cbind(newd, tilr, wilr, ID)
      
      ## No change dataset
      samed <- cbind(newd, bilr, wilr, ID)
      
      # Generate prediction
      ysub <- as.data.table(fitted(tmp$BrmModel, newdata = subd, re.form = NA, summary = FALSE))
      ysame <- as.data.table(fitted(tmp$BrmModel, newdata = samed, re.form = NA, summary = FALSE))
      
      # Difference between substitution and no change
      yd <- ysub - ysame
      yd <- rowMeans(yd) # ame of 1 substiution at participant level
      # nrow(yd) = 1, ncol = posterior draws
      # probably ame of covariates follows here
      myd <- as.data.table(describe_posterior(yd, centrality = "mean", 
                                              ci = 0.95, ci_method = "eti"))
      myd <- myd[, .(Mean, CI_low, CI_high)]
      myd$MinSubstituted <- mean(newd[, MinSubstituted])
      
      kout[[k]] <- myd # ame at substitution level
    }
    # Save results
    jout[[j]] <- do.call(rbind, kout) 
    # nrow(myds) = nrow(posub), ncol() = 4
  }
  jout <- as.data.table(do.call(rbind, jout))
  jout[, Substitute := rep(subvar, length.out = nrow(jout))]
  jout$Predictor <- iv
  names(jout) <- c("Mean", "CI_low", "CI_high", "MinSubstituted", "Substitute", "Predictor")
  
  # Final results for entire composition
  iout[[i]] <- jout
  
  }
  iout
}

