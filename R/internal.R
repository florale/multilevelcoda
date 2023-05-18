## make Rcmd check happy
utils::globalVariables(c("i",  "..cols", ".", "To", ".SD", "b", "t",
                         "Mean",  "CI_low", "CI_high", "From", "Delta",
                         "spread", "value", "variable", "ID", "EffectType", "Level",
                         "update"))
#' Functions used only internally
#' @keywords internal
#' @importFrom data.table as.data.table copy := setDT rbindlist
#' @importFrom compositions acomp ilr clo mean.acomp
#' @importFrom bayestestR describe_posterior
#' @importFrom extraoperators %snin% %sin%
#' @importFrom foreach foreach %dopar%
#' @importFrom stats fitted
#' @noRd
# Basic Between-person Substitution model
get.bsub <- function(object, basesub, recomp, 
                     yref, dref,
                     delta, summary,
                     level, type,
                     ID = 1, cv = NULL, regrid = NULL, ...) {
  
  iout <- foreach(i = colnames(basesub), .combine = c) %dopar% {
    
    # possible susbstituion of 1 compositional variable
    posub <- as.data.table(basesub)
    posub <- posub[(get(i) != 0)]
    posub <- posub[order(-rank(get(i)))]
    
    # substitution variable names
    subvar <- colnames(posub) %snin% eval(i)
    iv <- i
    
    kout <- vector("list", length = nrow(posub))
    jout <- vector("list", length = length(delta))
    
    for (j in seq_along(delta)) {
      # delta level
      sub <- posub * delta[j]
      for (k in seq_len(nrow(sub))) {
        newcomp <- recomp + sub[k,]
        names(newcomp) <- object$CompILR$parts
        Delta <- sub[k, get(i)]
        kout[[k]] <- cbind(recomp, newcomp, Delta)
      }
      jout[[j]] <- do.call(rbind, kout)
    }
    dnew <- setDT(do.call(rbind, jout))

    # useful information for the final results
    dnew[, From := rep(subvar, length.out = nrow(dnew))]
    dnew$To <- iv
    dnew$Delta <- as.numeric(dnew$Delta)
    dnew$Level <- level
    dnew$EffectType <- type
    
    # remove impossible reallocation that result in negative values 
    cols <- colnames(dnew) %snin% c("Delta", "From", "To", "Level", "EffectType")
    dnew <- dnew[rowSums(dnew[, ..cols] < 0) == 0]
    
    # compositions and ilrs for predictions
    bcomp <- acomp(dnew[, colnames(object$CompILR$BetweenComp), with = FALSE], total = object$CompILR$total)
    tcomp <- acomp(dnew[, object$CompILR$parts, with = FALSE], total = object$CompILR$total)
    
    bilr <- ilr(bcomp, V = object$CompILR$psi)
    tilr <- ilr(tcomp, V = object$CompILR$psi)
    
    wilr <- dref[, colnames(object$CompILR$WithinILR), with = FALSE]
    
    colnames(tilr) <- colnames(object$CompILR$BetweenILR)
    colnames(wilr) <- colnames(object$CompILR$WithinILR)
    
    # prediction
    if(is.null(regrid)) { # unadjusted
      dsub <- cbind(dnew, tilr, wilr, ID)
      ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
      
      ydiff <- apply(ysub, 2, function(y) {y - yref})
      suppressWarnings(ymean <- apply(ydiff, 2, function(x) {describe_posterior(x, centrality = "mean", ...)}))
      ymean <- rbindlist(ymean)
      ymean <- cbind(ymean[, .(Mean, CI_low, CI_high)], 
                     dsub[, .(Delta, To, From, Level, EffectType)])
      
    } else { # adjusted
      hout <- vector("list", length = nrow(regrid))
      if(isTRUE(summary)) { # averaging over reference grid
        for (h in seq_len(nrow(regrid))) {
          dsub <- cbind(dnew, tilr, wilr, ID, regrid[h, ])
          ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
          ydiff <- ysub - yref[, h]
          hout[[h]] <- ydiff
        }
        
        ymean <- Reduce(`+`, hout) / length(hout)
        suppressWarnings(ymean <- apply(ydiff, 2, function(x) {describe_posterior(x, centrality = "mean", ...)}))
        ymean <- rbindlist(ymean)
        ymean <- cbind(ymean[, .(Mean, CI_low, CI_high)], 
                       dsub[, .(Delta, To, From, Level, EffectType)])
        
      } else { # keeping prediction at each level of reference grid
        for (h in seq_len(nrow(regrid))) {
          dsub <- cbind(dnew, tilr, wilr, ID, regrid[h, ])
          ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
          ydiff <- ysub - yref[, h]
          suppressWarnings(ymean <- apply(ydiff, 2, function(x) {describe_posterior(x, centrality = "mean", ...)}))
          ymean <- rbindlist(ymean)
          ymean <- cbind(ymean[, .(Mean, CI_low, CI_high)],
                         dsub[, c("Delta", "From", "To", "Level", "EffectType", cv), 
                              with = FALSE])
          
          hout[[h]] <- ymean
        }
        ymean <- rbindlist(hout)
      }
    }
    # final results for entire composition
    out <- list(ymean)
    names(out) <- i
    out
  }
  iout
}

# Basic Within-person Substitution model
get.wsub <- function(object, basesub, recomp,
                     yref, dref,
                     delta, summary, 
                     level, type,
                     ID = 1, cv = NULL, regrid = NULL, ...) {
  
  iout <- foreach(i = colnames(basesub), .combine = c) %dopar% {
    
    # possible susbstituion of 1 compositional variable
    posub <- as.data.table(basesub)
    posub <- posub[(get(i) != 0)]
    posub <- posub[order(-rank(get(i)))]
    
    # substitution variable names
    subvar <- colnames(posub) %snin% eval(i)
    iv <- i
    
    kout <- vector("list", length = nrow(posub))
    jout <- vector("list", length = length(delta))
    
    for (j in seq_along(delta)) { # delta level
      sub <- posub * delta[j]
      for (k in seq_len(nrow(sub))) {
        newcomp <- recomp + sub[k, ]
        names(newcomp) <- object$CompILR$parts
        Delta <- sub[k, get(i)]
        kout[[k]] <- cbind(recomp, newcomp, Delta)
      }
      jout[[j]] <- do.call(rbind, kout)
    }
    dnew <- setDT(do.call(rbind, jout))
    
    # useful information for the final results
    dnew[, From := rep(subvar, length.out = nrow(dnew))]
    dnew$To <- iv
    dnew$Delta <- as.numeric(dnew$Delta)
    dnew$Level <- level
    dnew$EffectType <- type
    
    # remove impossible reallocation that result in negative values 
    cols <- colnames(dnew) %snin% c("Delta", "From", "To", "Level", "EffectType")
    dnew <- dnew[rowSums(dnew[, ..cols] < 0) == 0]
    
    # compositions and ilrs for predictions
    bcomp <- acomp(dnew[, colnames(object$CompILR$BetweenComp), with = FALSE], total = object$CompILR$total)
    tcomp <- acomp(dnew[, object$CompILR$parts, with = FALSE], total = object$CompILR$total)
    
    bilr <- ilr(bcomp, V = object$CompILR$psi)
    tilr <- ilr(tcomp, V = object$CompILR$psi)
    wilr <- tilr - bilr
    
    colnames(bilr) <- colnames(object$CompILR$BetweenILR)
    colnames(wilr) <- colnames(object$CompILR$WithinILR)
    
    # prediction
    if(is.null(regrid)) { # unadjusted
      dsub <- cbind(dnew, bilr, wilr, ID)
      ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
      
      ydiff <- apply(ysub, 2, function(y) {y - yref})
      suppressWarnings(ymean <- apply(ydiff, 2, function(x) {describe_posterior(x, centrality = "mean", ...)}))
      ymean <- rbindlist(ymean)
      ymean <- cbind(ymean[, .(Mean, CI_low, CI_high)], 
                     dsub[, .(Delta, To, From, Level, EffectType)])
      
    } else { # adjusted
      hout <- vector("list", length = nrow(regrid))
      if(isTRUE(summary)) { # averaging over reference grid
        for (h in seq_len(nrow(regrid))) {
          dsub <- cbind(dnew, bilr, wilr, ID, regrid[h, ])
          ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
          ydiff <- ysub - yref[, h]
          hout[[h]] <- ydiff
        }
        
        ymean <- Reduce(`+`, hout) / length(hout)
        suppressWarnings(ymean <- apply(ydiff, 2, function(x) {describe_posterior(x, centrality = "mean", ...)}))
        ymean <- rbindlist(ymean)
        ymean <- cbind(ymean[, .(Mean, CI_low, CI_high)], 
                       dsub[, .(Delta, To, From, Level, EffectType)])
        
      } else { # keeping prediction at each level of reference grid
        for (h in seq_len(nrow(regrid))) {
          dsub <- cbind(dnew, bilr, wilr, ID, regrid[h, ])
          ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
          ydiff <- ysub - yref[, h]
          suppressWarnings(ymean <- apply(ydiff, 2, function(x) {describe_posterior(x, centrality = "mean", ...)}))
          ymean <- rbindlist(ymean)
          ymean <- cbind(ymean[, .(Mean, CI_low, CI_high)],
                         dsub[, c("Delta", "From", "To", "Level", "EffectType", cv),
                              with = FALSE])
          
          hout[[h]] <- ymean
        }
        ymean <- rbindlist(hout)
      }
    }
    # final results for entire composition
    out <- list(ymean)
    names(out) <- i
    out
  }
  iout
}

# Between-person Marginal Substitution Model.
.get.bsubmargins <- function(object, basesub, b, 
                             yref, delta, 
                             level, type,
                             ...) {
  
  iout <- foreach(i = colnames(basesub), .combine = c) %dopar% {
    
    # possible susbstituion of 1 compositional variable
    posub <- as.data.table(basesub)
    posub <- posub[(get(i) != 0)]
    posub <- posub[order(-rank(get(i)))]
    
    # substitution variable names
    subvar <- colnames(posub) %snin% eval(i)
    iv <- i
    
    kout <- vector("list", length = nrow(posub))
    jout <- vector("list", length = length(delta))
    
    for (j in seq_along(delta)) { # delta level
      sub <- posub * delta[j]
      for (k in seq_len(nrow(sub))) {
        subk <- sub[k, ]
        subk <- subk[rep(seq_len(nrow(subk)), nrow(b)), ]
        newcomp <- b + subk
        Delta <- subk[, get(i)]
        names(newcomp) <- colnames(basesub)
        
        dnew <- cbind(b, newcomp, object$CompILR$data, Delta)
        
        # useful information for the final results
        dnew[, From := rep(subvar, length.out = nrow(dnew))[k]]
        dnew$To <- iv
        dnew$Delta <- as.numeric(dnew$Delta)
        
        # remove impossible reallocation that result in negative values 
        cols <- colnames(dnew) %sin% c(colnames(b), colnames(basesub))
        dnew <- dnew[rowSums(dnew[, ..cols] < 0) == 0]
        
        # compositions and ilrs for predictions
        bcomp <- acomp(dnew[, colnames(object$CompILR$BetweenComp), with = FALSE], total = object$CompILR$total)
        tcomp <- acomp(dnew[, object$CompILR$parts, with = FALSE], total = object$CompILR$total)
        
        bilr <- ilr(bcomp, V = object$CompILR$psi)
        tilr <- ilr(tcomp, V = object$CompILR$psi)
        
        wilr <- as.data.table(matrix(0, nrow = nrow(tilr), ncol = ncol(tilr)))
        
        colnames(tilr) <- colnames(object$CompILR$BetweenILR)
        colnames(wilr) <- colnames(object$CompILR$WithinILR)
        
        # prediction
        dsub <- cbind(dnew, tilr, wilr)
        ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
        ysub <- rowMeans(ysub)
        
        # difference in outcomes between substitution and no change
        ydiff <- ysub - yref
        
        # posterior means and intervals
        suppressWarnings(ymean <- setDT(describe_posterior(ydiff, centrality = "mean", ...)))
        ymean <- ymean[, .(Mean, CI_low, CI_high)]
        ymean$Delta <- sub[k, get(i)]
        kout[[k]] <- ymean
      }
      jout[[j]] <- rbindlist(kout)
    }
    jout <- rbindlist(jout)
    jout$To <- iv
    jout[, From := rep(subvar, length.out = nrow(jout))]
    jout$Level <- level
    jout$EffectType <- type
    
    names(jout) <- c("Mean", "CI_low", "CI_high", "Delta", "To", "From", 
                     "Level", "EffectType")
    
    # store final results for entire composition
    jout <- list(jout)
    names(jout) <- i
    jout
  }
  iout
}

# Within-person Marginal Substitution Model.
.get.wsubmargins <- function(object, basesub, b,
                             yref, delta, 
                             level, type,
                             ...) {
  
  iout <- foreach(i = colnames(basesub), .combine = c) %dopar% {
    
    posub <- as.data.table(basesub)
    posub <- posub[(get(i) != 0)]
    posub <- posub[order(-rank(get(i)))]
    
    # substitution variable names
    subvar <- colnames(posub) %snin% eval(i)
    iv <- i
    
    kout <- vector("list", length = nrow(posub))
    jout <- vector("list", length = length(delta))
    
    for (j in seq_along(delta)) { # delta level
      sub <- posub * delta[j]
      for (k in seq_len(nrow(sub))) {
        subk <- sub[k, ]
        subk <- subk[rep(seq_len(nrow(subk)), nrow(b)), ]
        newcomp <- b + subk
        Delta <- subk[, get(i)]
        names(newcomp) <- colnames(basesub)
        
        dnew <- cbind(b, newcomp, object$CompILR$data, Delta)
        
        # useful information for the final output
        dnew[, From := rep(subvar, length.out = nrow(dnew))[k]]
        dnew$To <- iv
        dnew$Delta <- as.numeric(dnew$Delta)
        
        # remove impossible reallocation that result in negative values 
        cols <- colnames(dnew) %sin% c(colnames(b), colnames(basesub))
        dnew <- dnew[rowSums(dnew[, ..cols] < 0) == 0]
        
        # compositions and ilr for predictions
        bcomp <- acomp(dnew[, colnames(object$CompILR$BetweenComp), with = FALSE], total = object$CompILR$total)
        tcomp <- acomp(dnew[, object$CompILR$parts, with = FALSE], total = object$CompILR$total)
        
        bilr <- ilr(bcomp, V = object$CompILR$psi)
        tilr <- ilr(tcomp, V = object$CompILR$psi)
        
        wilr <- tilr - bilr 
        
        colnames(bilr) <- colnames(object$CompILR$BetweenILR)
        colnames(wilr) <- colnames(object$CompILR$WithinILR)
        
        # substitution data
        dsub <- cbind(dnew, bilr, wilr)
        
        # prediction
        ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
        ysub <- rowMeans(ysub) 
        
        # difference between substitution and no change
        ydiff <- ysub - yref
        
        # posterior means and intervals
        suppressWarnings(ymean <- setDT(describe_posterior(ydiff, centrality = "mean", ...)))
        ymean <- ymean[, .(Mean, CI_low, CI_high)]
        ymean$Delta <- sub[k, get(i)]
        kout[[k]] <- ymean
      }
      # results
      jout[[j]] <- rbindlist(kout)
    }
    jout <- rbindlist(jout)
    jout$To <- iv
    jout[, From := rep(subvar, length.out = nrow(jout))]
    jout$Level <- level
    jout$EffectType <- type
    
    names(jout) <- c("Mean", "CI_low", "CI_high", "Delta", "To", "From", 
                     "Level", "EffectType")
    
    # final results for entire composition
    jout <- list(jout)
    names(jout) <- i
    jout
  }
  iout
}

# Marginal Substitution Model.
.get.submargins <- function(object, basesub, t,
                            yref, delta,
                            level, type,
                            ...) {
  
  iout <- foreach(i = colnames(basesub), .combine = c) %dopar% {
    
    # possible susbstituion of 1 compositional variable
    posub <- as.data.table(basesub)
    posub <- posub[(get(i) != 0)]
    posub <- posub[order(-rank(get(i)))]
    
    # substitution variable names
    subvar <- colnames(posub) %snin% eval(i)
    iv <- i
    
    kout <- vector("list", length = nrow(posub))
    jout <- vector("list", length = length(delta))
    
    for (j in seq_along(delta)) { # delta level
      sub <- posub * delta[j]
      for (k in seq_len(nrow(sub))) {
        subk <- sub[k, ]
        subk <- subk[rep(seq_len(nrow(subk)), nrow(t)), ]
        newcomp <- t + subk
        Delta <- subk[, get(i)]
        names(newcomp) <- colnames(basesub)
        
        dnew <- cbind(newcomp, object$CompILR$data, Delta)
        
        # useful information for the final results
        dnew[, From := rep(subvar, length.out = nrow(dnew))[k]]
        dnew$To <- iv
        dnew$Delta <- as.numeric(dnew$Delta)
        
        # remove impossible reallocation that result in negative values 
        cols <- colnames(dnew) %sin% c(colnames(t), colnames(basesub))
        dnew <- dnew[rowSums(dnew[, ..cols] < 0) == 0]
        
        # compositions and ilrs for predictions
        tcomp <- acomp(dnew[, object$CompILR$parts, with = FALSE], total = object$CompILR$total)
        tilr <- ilr(tcomp, V = object$CompILR$psi)
        
        colnames(tilr) <- colnames(object$CompILR$TotalILR)
        
        # prediction
        dsub <- cbind(dnew, tilr)
        ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
        ysub <- rowMeans(ysub)
        
        # difference in outcomes between substitution and no change
        ydiff <- ysub - yref
        
        # posterior means and intervals
        suppressWarnings(ymean <- setDT(describe_posterior(ydiff, centrality = "mean", ...)))
        ymean <- ymean[, .(Mean, CI_low, CI_high)]
        ymean$Delta <- sub[k, get(i)]
        kout[[k]] <- ymean
      }
      jout[[j]] <- rbindlist(kout)
    }
    jout <- rbindlist(jout)
    jout$To <- iv
    jout[, From := rep(subvar, length.out = nrow(jout))]
    jout$Level <- level
    jout$EffectType <- type
    
    names(jout) <- c("Mean", "CI_low", "CI_high", "Delta", "To", "From", 
                     "Level", "EffectType")
    
    # store final results for entire composition
    jout <- list(jout)
    names(jout) <- i
    jout
  }
  iout
}