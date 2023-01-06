## make Rcmd check happy
utils::globalVariables(c("i",  "..cols", ".", "To", ".SD",
                         "Mean",  "CI_low", "CI_high", "From", "Delta",
                         "spread", "value", "variable", "ID", "EffectType", "Level"))
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
get.bsub <- function(object, basesub, mcomp, 
                     ysame, delta, summary = summary,
                     level, type,
                     ID = 1, cv = NULL, refg = NULL, ...) {
  
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
        newcomp <- mcomp + sub[k, ]
        names(newcomp) <- object$CompIlr$parts
        Delta <- sub[k, get(i)]
        kout[[k]] <- cbind(mcomp, newcomp, Delta)
        }
      jout[[j]] <- do.call(rbind, kout)
    }
    newd <- setDT(do.call(rbind, jout))

    # useful information for the final results
    newd[, From := rep(subvar, length.out = nrow(newd))]
    newd$To <- iv
    newd$Delta <- as.numeric(newd$Delta)
    newd$Level <- level
    newd$EffectType <- type

    # remove impossible reallocation that result in negative values 
    cols <- colnames(newd) %snin% c("Delta", "From", "To", "Level", "EffectType")
    newd <- newd[rowSums(newd[, ..cols] < 0) == 0]

    # compositions and ilrs for predictions
    bcomp <- acomp(newd[, colnames(object$CompIlr$BetweenComp), with = FALSE])
    tcomp <- acomp(newd[, object$CompIlr$parts, with = FALSE])
    bilr <- ilr(bcomp, V = object$CompIlr$psi)
    tilr <- ilr(tcomp, V = object$CompIlr$psi)
    
    wilr <- matrix(0, nrow = nrow(tilr), ncol = ncol(tilr))
    wilr <- as.data.table(wilr)
    
    colnames(tilr) <- paste0("bilr", seq_len(ncol(tilr)))
    colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
    
    # prediction
    if(is.null(refg)) { # unadjusted
      subd <- cbind(newd, tilr, wilr, ID)
      ysub <- fitted(object$Model, newdata = subd, re_formula = NA, summary = FALSE)
      
      ydiff <- apply(ysub, 2, function(y) {y - ysame})
      ymean <- apply(ydiff, 2, function(x) {describe_posterior(x, centrality = "mean", ...)})
      ymean <- rbindlist(ymean)
      ymean <- cbind(ymean[, .(Mean, CI_low, CI_high)], 
                     subd[, .(Delta, To, From, Level, EffectType)])
      
      } else { # adjusted
        hout <- vector("list", length = nrow(refg))
          if(isTRUE(summary)) { # averaging over reference grid
            for (h in seq_len(nrow(refg))) {
              subd <- cbind(newd, tilr, wilr, ID, refg[h, ])
              ysub <- fitted(object$Model, newdata = subd, re_formula = NA, summary = FALSE)
              ydiff <- ysub - ysame[, h]
              hout[[h]] <- ydiff
            }
            
            ymean <- Reduce(`+`, hout) / length(hout)
            ymean <- apply(ydiff, 2, function(x) {describe_posterior(x, centrality = "mean", ...)})
            ymean <- rbindlist(ymean)
            ymean <- cbind(ymean[, .(Mean, CI_low, CI_high)], 
                           subd[, .(Delta, To, From, Level, EffectType)])
            
            } else { # keeping prediction at each level of reference grid
              for (h in seq_len(nrow(refg))) {
                subd <- cbind(newd, tilr, wilr, ID, refg[h, ])
                ysub <- fitted(object$Model, newdata = subd, re_formula = NA, summary = FALSE)
                ydiff <- ysub - ysame[, h]
                ymean <- apply(ydiff, 2, 
                               function(x) {describe_posterior(x, centrality = "mean", ...)})
                ymean <- rbindlist(ymean)
                ymean <- cbind(ymean[, .(Mean, CI_low, CI_high)],
                               subd[, c("Delta", "From", "To", "Level", "EffectType", cv), 
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
get.wsub <- function(object, basesub, mcomp,
                     ysame, delta, summary = summary, 
                     level, type,
                     ID = 1, cv = NULL, refg = NULL, ...) {
  
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
        newcomp <- mcomp + sub[k, ]
        names(newcomp) <- object$CompIlr$parts
        Delta <- sub[k, get(i)]
        kout[[k]] <- cbind(mcomp, newcomp, Delta)
      }
      jout[[j]] <- do.call(rbind, kout)
    }
    newd <- setDT(do.call(rbind, jout))

    # useful information for the final results
    newd[, From := rep(subvar, length.out = nrow(newd))]
    newd$To <- iv
    newd$Delta <- as.numeric(newd$Delta)
    newd$Level <- level
    newd$EffectType <- type
    
    # remove impossible reallocation that result in negative values 
    cols <- colnames(newd) %snin% c("Delta", "From", "To", "Level", "EffectType")
    newd <- newd[rowSums(newd[, ..cols] < 0) == 0]
    
    # compositions and ilrs for predictions
    bcomp <- acomp(newd[, colnames(object$CompIlr$BetweenComp), with = FALSE])
    tcomp <- acomp(newd[, object$CompIlr$parts, with = FALSE])
    wcomp <- tcomp - bcomp 
    
    bilr <- ilr(bcomp, V = object$CompIlr$psi)
    tilr <- ilr(tcomp, V = object$CompIlr$psi)
    wilr <- ilr(wcomp, V = object$CompIlr$psi)
    
    colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
    colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
    
    # prediction
    if(is.null(refg)) { # unadjusted
      subd <- cbind(newd, bilr, wilr, ID)
      ysub <- fitted(object$Model, newdata = subd, re_formula = NA, summary = FALSE)

      ydiff <- apply(ysub, 2, function(y) {y - ysame})
      ymean <- apply(ydiff, 2, function(x) {describe_posterior(x, centrality = "mean", ...)})
      ymean <- rbindlist(ymean)
      ymean <- cbind(ymean[, .(Mean, CI_low, CI_high)], 
                     subd[, .(Delta, To, From, Level, EffectType)])
      
      } else { # adjusted
        hout <- vector("list", length = nrow(refg))
          if(isTRUE(summary)) { # averaging over reference grid
            for (h in seq_len(nrow(refg))) {
              subd <- cbind(newd, bilr, wilr, ID, refg[h, ])
              ysub <- fitted(object$Model, newdata = subd, re_formula = NA, summary = FALSE)
              ydiff <- ysub - ysame[, h]
              hout[[h]] <- ydiff
            }
            
            ymean <- Reduce(`+`, hout) / length(hout)
            ymean <- apply(ydiff, 2, function(x) {describe_posterior(x, centrality = "mean", ...)})
            ymean <- rbindlist(ymean)
            ymean <- cbind(ymean[, .(Mean, CI_low, CI_high)], 
                           subd[, .(Delta, To, From, Level, EffectType)])
            
            } else { # keeping prediction at each level of reference grid
              for (h in seq_len(nrow(refg))) {
                subd <- cbind(newd, bilr, wilr, ID, refg[h, ])
                ysub <- fitted(object$Model, newdata = subd, re_formula = NA, summary = FALSE)
                ydiff <- ysub - ysame[, h]
                ymean <- apply(ydiff, 2, 
                               function(x) {describe_posterior(x, centrality = "mean", ...)})
                ymean <- rbindlist(ymean)
                ymean <- cbind(ymean[, .(Mean, CI_low, CI_high)],
                               subd[, c("Delta", "From", "To", "Level", "EffectType", cv),
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
                             ysame, delta, 
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
        
        newd <- cbind(b, newcomp, object$CompIlr$data, Delta)
        
        # useful information for the final results
        newd[, From := rep(subvar, length.out = nrow(newd))[k]]
        newd$To <- iv
        newd$Delta <- as.numeric(newd$Delta)
        
        # remove impossible reallocation that result in negative values 
        cols <- colnames(newd) %sin% c(colnames(b), colnames(basesub))
        newd <- newd[rowSums(newd[, ..cols] < 0) == 0]
        
        # compositions and ilrs for predictions
        bcomp <- acomp(newd[, colnames(object$CompIlr$BetweenComp), with = FALSE])
        tcomp <- acomp(newd[, object$CompIlr$parts, with = FALSE])
        bilr <- ilr(bcomp, V = object$CompIlr$psi)
        tilr <- ilr(tcomp, V = object$CompIlr$psi)
        
        wilr <- as.data.table(matrix(0, nrow = nrow(tilr), ncol = ncol(tilr)))
        
        colnames(tilr) <- paste0("bilr", seq_len(ncol(tilr)))
        colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
        
        # prediction
        subd <- cbind(newd, tilr, wilr)
        ysub <- fitted(object$Model, newdata = subd, re_formula = NA, summary = FALSE)
        ysub <- rowMeans(ysub)
        
        # difference in outcomes between substitution and no change
        ydiff <- ysub - ysame
        
        # posterior means and intervals
        ymean <- setDT(describe_posterior(ydiff, centrality = "mean", ...))
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
                             ysame, delta, 
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
        
        newd <- cbind(b, newcomp, object$CompIlr$data, Delta)
        
        # useful information for the final output
        newd[, From := rep(subvar, length.out = nrow(newd))[k]]
        newd$To <- iv
        newd$Delta <- as.numeric(newd$Delta)
        
        # remove impossible reallocation that result in negative values 
        cols <- colnames(newd) %sin% c(colnames(b), colnames(basesub))
        newd <- newd[rowSums(newd[, ..cols] < 0) == 0]
        
        # compositions and ilr for predictions
        bcomp <- acomp(newd[, colnames(object$CompIlr$BetweenComp), with = FALSE])
        tcomp <- acomp(newd[, object$CompIlr$parts, with = FALSE])
        bilr <- ilr(bcomp, V = object$CompIlr$psi)
        tilr <- ilr(tcomp, V = object$CompIlr$psi)
        
        wilr <- tilr - bilr 
        
        colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
        colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
        
        # substitution data
        subd <- cbind(newd, bilr, wilr)
        
        # prediction
        ysub <- fitted(object$Model, newdata = subd, re_formula = NA, summary = FALSE)
        ysub <- rowMeans(ysub) 
        
        # difference between substitution and no change
        ydiff <- ysub - ysame
        
        # posterior means and intervals
        ymean <- setDT(describe_posterior(ydiff, centrality = "mean", ...))
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
                            ysame, delta,
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
        
        newd <- cbind(newcomp, object$CompIlr$data, Delta)
        
        # useful information for the final results
        newd[, From := rep(subvar, length.out = nrow(newd))[k]]
        newd$To <- iv
        newd$Delta <- as.numeric(newd$Delta)
        
        # remove impossible reallocation that result in negative values 
        cols <- colnames(newd) %sin% c(colnames(t), colnames(basesub))
        newd <- newd[rowSums(newd[, ..cols] < 0) == 0]
        
        # compositions and ilrs for predictions
        tcomp <- acomp(newd[, object$CompIlr$parts, with = FALSE])
        tilr <- ilr(tcomp, V = object$CompIlr$psi)
        
        colnames(tilr) <- paste0("ilr", seq_len(ncol(tilr)))

        # prediction
        subd <- cbind(newd, tilr)
        ysub <- fitted(object$Model, newdata = subd, re_formula = NA, summary = FALSE)
        ysub <- rowMeans(ysub)
        
        # difference in outcomes between substitution and no change
        ydiff <- ysub - ysame
        
        # posterior means and intervals
        ymean <- setDT(describe_posterior(ydiff, centrality = "mean", ...))
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