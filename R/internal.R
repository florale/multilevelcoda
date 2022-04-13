#' # Functions used only internally
#' @keywords internal
#' @importFrom data.table as.data.table
#' @importFrom compositions clo
#' @noRd
#' ## Support Substitution Model
#' # Get possible substitution

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
.get.bsubm <- function(substitute, b, tmp, 
                      min, psi, ID = 1,
                      ysame, refg = NULL, h = 0) {
  iout <- foreach(i = colnames(substitute), .combine = c) %dopar% {
    
    posub <- copy(substitute)
    posub <- as.data.table(posub)
    posub <- posub[(get(i) != 0)]
    posub <- posub[order(-rank(get(i)))]
    
    # get substitution variable names
    subvar <- colnames(posub) %snin% eval(i) # substitute compositional variables
    iv <- i # central compositional variable
    
    # lists to store results - TODO
    kout <- vector("list", length = nrow(posub))
    jout <- vector("list", length = min)
    
    for (j in seq_len(min)) { # time level
      sub <- posub * j
      for (k in seq_len(nrow(sub))) {
        subk <- sub[k, ]
        subk <- subk[rep(seq_len(nrow(subk)), nrow(b)), ]
        newcomp <- b + subk
        MinSubstituted <- subk[, get(i)]
        names(newcomp) <- colnames(substitute)
        
        newd <- cbind(b, newcomp, MinSubstituted)

        if(h == 0) {
          newd <- newd
        } else {
          newd <- newd[rep(seq_len(nrow(newd)), h), ]
        }
        # Add more useful information to the final output
        newd[, Substitute := rep(subvar, length.out = nrow(newd))[k]]
        newd$Predictor <- iv
        newd$MinSubstituted <- as.numeric(newd$MinSubstituted)
        newd <- cbind(newd, refg)
        
        ## remove impossible reallocation that result in negative values 
        cols <- colnames(newd) %sin% c(colnames(b), colnames(substitute))
        newd <- newd[rowSums(newd[, ..cols] < 0) == 0]
        
        # Add comp and ilr for predictions
        bcomp <- acomp(newd[, colnames(tmp$CompIlr$BetweenComp), with = FALSE])
        tcomp <- acomp(newd[, tmp$CompIlr$composition, with = FALSE])
        
        bilr <- ilr(bcomp, V = psi)
        tilr <- ilr(tcomp, V = psi)
        wilr <- matrix(0, nrow = nrow(bilr), ncol = ncol(bilr))
        wilr <- as.data.table(wilr)
        
        colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
        colnames(tilr) <- paste0("bilr", seq_len(ncol(tilr)))
        colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
        
        ## Substitution dataset
        dsub <- cbind(newd, tilr, wilr, ID)
        
        # Generate prediction
        ysub <- fitted(tmp$BrmModel, newdata = dsub, re.form = NA, summary = FALSE)
        ysub <- rowMeans(ysub) # emmeans across participants for 1 possible substitution 
        
        # Y difference between substitution and no change
        ydiff <- ysub - ysame
        
        myd <- as.data.table(describe_posterior(ydiff, centrality = "mean",
                                                ci = 0.95, ci_method = "eti"))
        myd <- myd[, .(Mean, CI_low, CI_high)]
        myd$MinSubstituted <- sub[k, get(i)]
        kout[[k]] <- myd # ame at substitution level
        }
      # Save results
      jout[[j]] <- rbindlist(kout)
      }
    jout <- rbindlist(jout)
    jout[, Substitute := rep(subvar, length.out = nrow(jout))]
    jout$Predictor <- iv
    names(jout) <- c("Mean", "CI_low", "CI_high", "MinSubstituted", "Substitute", "Predictor")
    
    # Final results for entire composition
    jout <- list(jout)
    names(jout) <- i
    jout
  }
  iout
  }

.get.bsubm2 <- function(substitute, b, tmp, 
                      min, psi, ID = 1,
                      ysame, refg = NULL, h = 0) {
  # list to store final output
  iout <- vector("list")
  for(i in colnames(substitute)) { # compostion level
    posub <- copy(substitute)
    posub <- as.data.table(posub)
    posub <- posub[(get(i) != 0)]
    posub <- posub[order(-rank(get(i)))]

    # get substitution variable name
    subvar <- colnames(posub) %snin% eval(i) # substitute compositional variables
    iv <- i # central compositional variable

    # lists to store results
    kout <- vector("list", length = nrow(posub))
    jout <- vector("list", length = min)

    for (j in seq_len(min)) { # time level
      sub <- posub * j
      for (k in seq_len(nrow(sub))) {
        subk <- sub[k, ]
        subk <- subk[rep(seq_len(nrow(subk)), nrow(b)), ]
        newcomp <- b + subk
        MinSubstituted <- subk[, get(i)]
        names(newcomp) <- colnames(substitute)

        newd <- cbind(b, newcomp, MinSubstituted)

        if(h == 0) {
          newd <- newd
        } else {
          newd <- newd[rep(seq_len(nrow(newd)), h), ]
        }
        # Add more useful information to the final output
        newd[, Substitute := rep(subvar, length.out = nrow(newd))[k]]
        newd$Predictor <- iv
        newd$MinSubstituted <- as.numeric(newd$MinSubstituted)
        newd <- cbind(newd, refg)

        kout[[k]] <- newd
      }
      newd <- as.data.table(do.call(rbind, kout))
        ## remove impossible reallocation that result in negative values - TODO
        cols <- colnames(newd) %snin% c("MinSubstituted", "Substitute", "Predictor")
        newd <- newd[rowSums(newd[, ..cols] < 0) == 0]

        # Add comp and ilr for predictions
        tcomp <- acomp(newd[, tmp$CompIlr$composition, with = FALSE])

        tilr <- ilr(tcomp, V = psi)
        wilr <- matrix(0, nrow = nrow(newd), ncol = ncol(tilr))
        wilr <- as.data.table(wilr)

        colnames(tilr) <- paste0("bilr", seq_len(ncol(tilr)))
        colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))

        ## Substitution dataset
        dsub <- cbind(newd, tilr, wilr, ID)
        # Generate prediction
        ysub <- fitted(tmp$BrmModel, newdata = dsub, re.form = NA, summary = FALSE)
        ysub <- cbind(newd[, .(MinSubstituted, Substitute, Predictor)], as.data.table(t(ysub)))
        ysub <- ysub[, lapply(.SD, mean), by = c("MinSubstituted", "Substitute", "Predictor")]
        suppl <- ysub[, .(MinSubstituted, Substitute, Predictor)]
        ysub <- ysub[, -c(1:3)]
        # Difference between substitution and no change
        ydiff <- t(ysub) - ysame

        myd <- apply(ydiff, 2, function(x) {describe_posterior(x, centrality = "mean",
                                                               ci = 0.95, ci_method = "eti")})
        myd <- rbindlist(myd)
        myd <- myd[, .(Mean, CI_low, CI_high)]
        myd <- cbind(myd, suppl)
      # Save results
      jout[[j]] <- myd
    }
    jout <- as.data.table(do.call(rbind, jout))
    names(jout) <- c("Mean", "CI_low", "CI_high", "MinSubstituted", "Substitute", "Predictor")

    # Final results for entire composition
    iout[[i]] <- jout
    }
    iout
    }