## make Rcmd check happy
utils::globalVariables(c("i",  "..cols", ".", "To", ".SD", "t", "head",
                         "Mean",  "CI_low", "CI_high", "From", "Delta",
                         "spread", "value", "variable", "ID", "EffectType", "Level", "Reference",
                         "update"))
#' Functions used only internally
#' @keywords internal
#' @importFrom data.table as.data.table data.table copy := setDT rbindlist
#' @importFrom compositions acomp ilr clo mean.acomp
#' @importFrom bayestestR describe_posterior
#' @importFrom extraoperators %snin% %sin%
#' @importFrom foreach foreach %dopar%
#' @importFrom stats fitted
#' @noRd
# Grandmean Between-person Substitution model
get.bsub <- function(object, delta, basesub,
                     comp0, y0, d0,
                     summary,
                     level, ref,
                     ...) {
  
  iout <- foreach(i = colnames(basesub), .combine = c) %dopar% {
    
    # dnew - reallocation data
    # possible substitution of 1 compositional variable
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
        newcomp <- comp0 + sub[k,]
        names(newcomp) <- object$CompILR$parts
        Delta <- sub[k, get(i)]
        kout[[k]] <- cbind(comp0, newcomp, Delta)
      }
      jout[[j]] <- do.call(rbind, kout)
    }
    dnew <- setDT(do.call(rbind, jout))
    
    # useful information for the final results
    dnew[, From := rep(subvar, length.out = nrow(dnew))]
    dnew$To <- iv
    dnew$Delta <- as.numeric(dnew$Delta)
    dnew$Level <- level
    dnew$Reference <- ref
    
    # remove impossible reallocation that result in negative values 
    cols <- colnames(dnew) %snin% c("Delta", "From", "To", "Level")
    dnew <- dnew[rowSums(dnew[, ..cols] < 0) == 0]
    
    # compositions and ilrs for predictions
    bcomp0 <- acomp(dnew[, colnames(object$CompILR$BetweenComp), with = FALSE], total = object$CompILR$total)
    bcompsub  <- acomp(dnew[, object$CompILR$parts, with = FALSE], total = object$CompILR$total)
    
    bilrsub <- ilr(bcompsub, V = object$CompILR$psi)
    wilr0 <- d0[1, colnames(object$CompILR$WithinILR), with = FALSE]
    
    colnames(bilrsub) <- colnames(object$CompILR$BetweenILR)
    colnames(wilr0) <- colnames(object$CompILR$WithinILR)
    
    # reference grid 
    ## get covariate + idvar names
    covnames <- colnames(d0) %snin% c(colnames(object$CompILR$BetweenILR),
                                      colnames(object$CompILR$WithinILR),
                                      colnames(object$CompILR$BetweenComp),
                                      colnames(object$CompILR$WithinComp)
    )
    refgrid <- d0[, covnames, with = FALSE]
    
    # predictions
    hout <- vector("list", length = nrow(refgrid))
    if(isTRUE(summary)) { # unadj OR adj averaging over reference grid
      for (h in seq_len(nrow(refgrid))) {
        dsub <- cbind(dnew, bilrsub, wilr0, refgrid[h, ])
        ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
        delta_y <- ysub - y0[, h]
        hout[[h]] <- delta_y
      }
      
      PD_delta_y <- Reduce(`+`, hout) / length(hout)
      suppressWarnings(PD_delta_y <- apply(delta_y, 2, function(x) {describe_posterior(x, centrality = "mean", ...)}))
      PD_delta_y <- rbindlist(PD_delta_y)
      PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)], 
                          dsub[, .(Delta, From, To, Level, Reference)])
      
    } else { # adj keeping prediction at each level of reference grid
      for (h in seq_len(nrow(refgrid))) {
        dsub <- cbind(dnew, bilrsub, wilr0, refgrid[h, ])
        ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
        delta_y <- ysub - y0[, h]
        suppressWarnings(PD_delta_y <- apply(delta_y, 2, function(x) {describe_posterior(x, centrality = "mean", ...)}))
        PD_delta_y <- rbindlist(PD_delta_y)
        PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)], 
                            dsub[, .(Delta, From, To, Level, Reference)],
                            dsub[, colnames(refgrid) %snin% object$CompILR$idvar, with = FALSE])
        
        hout[[h]] <- PD_delta_y
      }
      PD_delta_y <- rbindlist(hout)
    }
    # final results for entire composition
    out <- list(PD_delta_y)
    names(out) <- i
    out
  }
  iout
}

# Grandmean Within-person Substitution model
get.wsub <- function(object, delta, basesub,
                     comp0, y0, d0,
                     summary,
                     level, ref,
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
        newcomp <- comp0 + sub[k, ]
        names(newcomp) <- object$CompILR$parts
        Delta <- sub[k, get(i)]
        kout[[k]] <- cbind(comp0, newcomp, Delta)
      }
      jout[[j]] <- do.call(rbind, kout)
    }
    dnew <- setDT(do.call(rbind, jout))
    
    # useful information for the final results
    dnew[, From := rep(subvar, length.out = nrow(dnew))]
    dnew$To <- iv
    dnew$Delta <- as.numeric(dnew$Delta)
    dnew$Level <- level
    dnew$Reference <- ref
    
    # remove impossible reallocation that result in negative values 
    cols <- colnames(dnew) %snin% c("Delta", "From", "To", "Level")
    dnew <- dnew[rowSums(dnew[, ..cols] < 0) == 0]
    
    # compositions and ilrs for predictions
    bcomp0   <- acomp(dnew[, colnames(object$CompILR$BetweenComp), with = FALSE], total = object$CompILR$total)
    bcompsub <- acomp(dnew[, object$CompILR$parts, with = FALSE], total = object$CompILR$total)
    
    bilr0   <- ilr(bcomp0, V = object$CompILR$psi)
    bilrsub <- ilr(bcompsub, V = object$CompILR$psi)
    wilrsub <- bilrsub - bilr0
    
    colnames(bilr0) <- colnames(object$CompILR$BetweenILR)
    colnames(wilrsub) <- colnames(object$CompILR$WithinILR)
    
    # reference grid 
    ## get covariate + idvar names
    covnames <- colnames(d0) %snin% c(colnames(object$CompILR$BetweenILR),
                                      colnames(object$CompILR$WithinILR),
                                      colnames(object$CompILR$BetweenComp),
                                      colnames(object$CompILR$WithinComp)
    )
    refgrid <- d0[, covnames, with = FALSE]
    
    # predictions
    hout <- vector("list", length = nrow(refgrid))
    if(isTRUE(summary)) { # unadj OR adj averaging over reference grid
      for (h in seq_len(nrow(refgrid))) {
        dsub <- cbind(dnew, bilr0, wilrsub, refgrid[h, ])
        ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
        delta_y <- ysub - y0[, h]
        hout[[h]] <- delta_y
      }
      
      PD_delta_y <- Reduce(`+`, hout) / length(hout)
      suppressWarnings(PD_delta_y <- apply(delta_y, 2, function(x) {describe_posterior(x, centrality = "mean", ...)}))
      PD_delta_y <- rbindlist(PD_delta_y)
      PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)], 
                          dsub[, .(Delta, From, To, Level, Reference)])
      
    } else { # adj keeping prediction at each level of reference grid
      for (h in seq_len(nrow(refgrid))) {
        dsub <- cbind(dnew, bilr0, wilrsub, refgrid[h, ])
        ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
        delta_y <- ysub - y0[, h]
        suppressWarnings(PD_delta_y <- apply(delta_y, 2, function(x) {describe_posterior(x, centrality = "mean", ...)}))
        PD_delta_y <- rbindlist(PD_delta_y)
        PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)], 
                            dsub[, .(Delta, From, To, Level, Reference)],
                            dsub[, colnames(refgrid) %snin% object$CompILR$idvar, with = FALSE])
        
        hout[[h]] <- PD_delta_y
      }
      PD_delta_y <- rbindlist(hout)
    }
    # final results for entire composition
    out <- list(PD_delta_y)
    names(out) <- i
    out
  }
  iout
}

# clustermean Between-person Substitution Model.
.get.bsubmargins <- function(object, delta, basesub,
                             comp0, y0, d0,
                             level, ref,
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
      for (k in seq_len(nrow(sub))) { # reallocation level
        subk <- sub[k, ]
        subk <- subk[rep(seq_len(nrow(subk)), nrow(comp0)), ]
        newcomp <- comp0 + subk
        Delta <- subk[, get(i)]
        names(newcomp) <- object$CompILR$parts
        
        dnew <- cbind(comp0, newcomp, 
                      d0[, colnames(d0) %in% colnames(object$CompILR$data[, -object$CompILR$part, with = FALSE]), with = FALSE], 
                      Delta)
        # useful information for the final results
        dnew[, From := rep(subvar, length.out = nrow(dnew))[k]]
        dnew$To <- iv
        dnew$Delta <- as.numeric(dnew$Delta)
        
        # remove impossible reallocation that result in negative values 
        cols <- colnames(dnew) %sin% c(colnames(comp0), colnames(basesub))
        dnew <- dnew[rowSums(dnew[, ..cols] < 0) == 0]
        
        # compositions and ilrs for predictions
        bcomp0    <- acomp(dnew[, colnames(object$CompILR$BetweenComp), with = FALSE], total = object$CompILR$total)
        bcompsub  <- acomp(dnew[, object$CompILR$parts, with = FALSE], total = object$CompILR$total)
        
        bilr0    <- ilr(bcomp0, V = object$CompILR$psi)
        bilrsub  <- ilr(bcompsub, V = object$CompILR$psi)
        
        wilr0 <- as.data.table(matrix(0, nrow = nrow(bilrsub), ncol = ncol(bilrsub)))
        
        colnames(bilrsub) <- colnames(object$CompILR$BetweenILR)
        colnames(wilr0)   <- colnames(object$CompILR$WithinILR)
        
        # prediction
        dsub <- cbind(dnew, bilrsub, wilr0)
        ysub <- fitted(object$Model, newdata = dsub, re_formula = NULL, summary = FALSE)
        ysub <- rowMeans(ysub)
        
        # difference in outcomes between substitution and no change
        delta_y <- ysub - y0
        
        # posterior means and intervals
        suppressWarnings(PD_delta_y <- setDT(describe_posterior(delta_y, centrality = "mean", ...)))
        PD_delta_y <- PD_delta_y[, .(Mean, CI_low, CI_high)]
        PD_delta_y$Delta <- sub[k, get(i)]
        kout[[k]] <- PD_delta_y
      }
      jout[[j]] <- rbindlist(kout)
    }
    jout <- rbindlist(jout)
    jout$To <- iv
    jout[, From := rep(subvar, length.out = nrow(jout))]
    jout$Level <- level
    jout$Reference <- ref
    
    names(jout) <- c("Mean", "CI_low", "CI_high", 
                     "Delta", "To", "From", "Level", "Reference")
    
    # store final results for entire composition
    jout <- list(jout)
    names(jout) <- i
    jout
  }
  iout
}

# clustermean Within-person Substitution Model.
.get.wsubmargins <- function(object, delta, basesub,
                             comp0, y0, d0,
                             level, ref,
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
        subk <- subk[rep(seq_len(nrow(subk)), nrow(comp0)), ]
        newcomp <- comp0 + subk
        Delta <- subk[, get(i)]
        names(newcomp) <- object$CompILR$parts
        
        dnew <- cbind(comp0, newcomp, 
                      d0[, colnames(d0) %in% colnames(object$CompILR$data[, -object$CompILR$part, with = FALSE]), with = FALSE], 
                      Delta)
        
        # useful information for the final output
        dnew[, From := rep(subvar, length.out = nrow(dnew))[k]]
        dnew$To <- iv
        dnew$Delta <- as.numeric(dnew$Delta)
        
        # remove impossible reallocation that result in negative values 
        cols <- colnames(dnew) %sin% c(colnames(comp0), colnames(basesub))
        dnew <- dnew[rowSums(dnew[, ..cols] < 0) == 0]
        
        # compositions and ilr for predictions
        bcomp0   <- acomp(dnew[, colnames(object$CompILR$BetweenComp), with = FALSE], total = object$CompILR$total)
        bcompsub <- acomp(dnew[, object$CompILR$parts, with = FALSE], total = object$CompILR$total)
        
        bilr0   <- ilr(bcomp0, V = object$CompILR$psi)
        bilrsub <- ilr(bcompsub, V = object$CompILR$psi)
        
        wilrsub <- bilrsub - bilr0 
        
        colnames(bilr0)   <- colnames(object$CompILR$BetweenILR)
        colnames(wilrsub) <- colnames(object$CompILR$WithinILR)
        
        # substitution data
        dsub <- cbind(dnew, bilr0, wilrsub)
        
        # prediction
        ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
        ysub <- rowMeans(ysub) 
        
        # difference between substitution and no change
        delta_y <- ysub - y0
        
        # describe PD of delta y
        suppressWarnings(PD_delta_y <- setDT(describe_posterior(delta_y, centrality = "mean", ...)))
        PD_delta_y <- PD_delta_y[, .(Mean, CI_low, CI_high)]
        PD_delta_y$Delta <- sub[k, get(i)]
        kout[[k]] <- PD_delta_y
      }
      # results
      jout[[j]] <- rbindlist(kout)
    }
    jout <- rbindlist(jout)
    jout$To <- iv
    jout[, From := rep(subvar, length.out = nrow(jout))]
    jout$Level <- level
    jout$Reference <- ref
    
    names(jout) <- c("Mean", "CI_low", "CI_high", 
                     "Delta", "To", "From", "Level", "Reference")
    
    # final results for entire composition
    jout <- list(jout)
    names(jout) <- i
    jout
  }
  iout
}

# clustermean Substitution Model.
.get.submargins <- function(object, basesub, t,
                            y0, delta,
                            level, ref,
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
        names(newcomp) <- object$CompILR$parts
        
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
        ysub <- fitted(object$Model, newdata = dsub, re_formula = NULL, summary = FALSE)
        ysub <- rowMeans(ysub)
        
        # difference in outcomes between substitution and no change
        delta_y <- ysub - y0
        
        # describe PD of delta y
        suppressWarnings(PD_delta_y <- setDT(describe_posterior(delta_y, centrality = "mean", ...)))
        PD_delta_y <- PD_delta_y[, .(Mean, CI_low, CI_high)]
        PD_delta_y$Delta <- sub[k, get(i)]
        kout[[k]] <- PD_delta_y
      }
      jout[[j]] <- rbindlist(kout)
    }
    jout <- rbindlist(jout)
    jout$To <- iv
    jout[, From := rep(subvar, length.out = nrow(jout))]
    jout$Level <- level

    names(jout) <- c("Mean", "CI_low", "CI_high", 
                     "Delta", "To", "From", 
                     "Level")
    
    # store final results for entire composition
    jout <- list(jout)
    names(jout) <- i
    jout
  }
  iout
}

# reference dataset
build.rg <- function(object, 
                     ref,
                     weight, 
                     fill = FALSE) {
  
  cov.grid <- NULL
  
  if(isTRUE(ref == "clustermean")) {
    if(isTRUE(weight == "equal")) {
      comp0 <- object$CompILR$BetweenComp
      d0 <- cbind(object$CompILR$BetweenComp, object$CompILR$data)
      d0 <- d0[, head(.SD, 1), by = eval(object$CompILR$idvar)]
      comp0 <- acomp(d0[, colnames(object$CompILR$BetweenComp), with = FALSE], total = object$CompILR$total)
    }
    if (isTRUE(is.null(weight) || weight == "proportional")) {
      comp0 <- object$CompILR$BetweenComp
      d0 <- object$CompILR$data
      
    }
    
    # d0 ---------------------------
    bcomp0 <- comp0
    bilr0 <- ilr(bcomp0, V = object$CompILR$psi)
    bilr0 <- as.data.table(bilr0)
    
    wcomp0 <- as.data.table(matrix(1, nrow = nrow(bcomp0), ncol = ncol(bcomp0)))
    wilr0 <- as.data.table(matrix(0, nrow = nrow(bilr0), ncol = ncol(bilr0)))
    
    colnames(bilr0) <- colnames(object$CompILR$BetweenILR)
    colnames(wilr0) <- colnames(object$CompILR$WithinILR)
    colnames(bcomp0) <- colnames(object$CompILR$BetweenComp)
    colnames(wcomp0) <- colnames(object$CompILR$WithinComp)
    
    d0 <- cbind(bilr0, wilr0, bcomp0, wcomp0, 
                d0[, colnames(d0) %in% colnames(object$CompILR$data), with = FALSE])
  }
  
  if(isTRUE(ref == "grandmean")) {
    if(isTRUE(is.null(weight) || weight == "equal")) {
      comp0 <- cbind(object$CompILR$BetweenComp, 
                     object$CompILR$data[, object$CompILR$idvar, with = FALSE])
      comp0 <- comp0[, head(.SD, 1), by = eval(object$CompILR$idvar)]
      comp0 <- acomp(comp0[, -object$CompILR$idvar, with = FALSE], total = object$CompILR$total)
      comp0 <- mean(comp0, robust = TRUE)
      comp0 <- acomp(comp0, total = object$CompILR$total)
      comp0 <- as.data.table(t(comp0))
      mcomp <- comp0
    }
    if (isTRUE(weight == "proportional")) {
      comp0 <- mean(object$CompILR$BetweenComp, robust = TRUE)
      comp0 <- acomp(comp0, total = object$CompILR$total)
      comp0 <- as.data.table(t(comp0))
      mcomp <- comp0
    }
    
    if(isTRUE(inherits(ref, c("data.table", "data.frame", "matrix")))) {
      # get compositional mean
      mcomp <- cbind(object$CompILR$BetweenComp, 
                     object$CompILR$data[, object$CompILR$idvar, with = FALSE])
      mcomp <- mcomp[, head(.SD, 1), by = eval(object$CompILR$idvar)]
      mcomp <- acomp(mcomp[, -object$CompILR$idvar, with = FALSE], total = object$CompILR$total)
      mcomp <- mean(mcomp, robust = TRUE)
      mcomp <- acomp(mcomp, total = object$CompILR$total)
      mcomp <- as.data.table(t(mcomp))
      
      if (isTRUE(identical(length(object$CompILR$parts), ncol(ref)))) {
        comp0 <- ref
        
      } else {
        if (object$CompILR$parts %nin% colnames(ref)) {
          stop(
            sprintf(
              "The reference grid should include all compositional components but (%s) are missing.",
              paste0(object$CompILR$parts %nin% colnames(ref), collapse = ", ")
            ))
        }
        comp0 <- ref[, object$CompILR$parts, with = FALSE]
        cov.grid <- ref[, -object$CompILR$parts, with = FALSE]
      }
      
      if(isFALSE(sum(comp0) == object$CompILR$total)) {
        stop(sprintf(
          "The total amount of the reference composition (%s) should be the same as the composition (%s).",
          sum(comp0),
          object$CompILR$total
        ))
      }
      if (isTRUE((any(comp0 > lapply(object$CompILR$data[, object$CompILR$parts, with = FALSE], max)) |
                  any(comp0 < lapply(object$CompILR$data[, object$CompILR$parts, with = FALSE], min))))) {
        stop(paste(
          sprintf(
            "comp0 should be numeric or interger values that are between (%s) and (%s)",
            paste0(round(apply(object$CompILR$data[, object$CompILR$parts, with = FALSE], 2, min)), collapse = ", "),
            paste0(round(apply(object$CompILR$data[, object$CompILR$parts, with = FALSE], 2, max)), collapse = ", ")),
          "\n", 
          " for",
          paste0(object$CompILR$parts, collapse = ", "),
          "respectively"
        ))
      }
      comp0  <- acomp(comp0, total = object$CompILR$total)
      comp0  <- as.data.table(t(comp0))
      colnames(comp0) <- colnames(object$CompILR$BetweenComp)
    }
    
    # reference grid -------------------
    # get possible ilr names
    ilrnames <- c(colnames(object$CompILR$BetweenILR), 
                  colnames(object$CompILR$WithinILR),
                  colnames(object$CompILR$TotalILR))
    
    # get all varnames in model
    varnames <- do.call(rbind, find_predictors(object$Model))
    
    if (isTRUE(all(varnames %in% ilrnames))) { # unadj subsitution model
      # get varnames in model in ilr names
      refgrid <- NULL
      
    } else { # adj subsitution model
      # default reference grid
      rg <- as.data.table(ref_grid(object$Model)@grid)
      covnames <- colnames(rg) %snin% c(ilrnames, ".wgt.")
      
      # user's specified reference grid
      if (isTRUE(is.null(cov.grid))) {
        refgrid <- rg[, covnames, with = FALSE]
      } else {
        gridnames <- colnames(cov.grid) %nin% object$CompILR$parts
        
        if(isFALSE(fill)) {
          if(isFALSE(identical(colnames(cov.grid), covnames))) { # ensure all covs are provided
            stop(paste(
              "'cov.grid' should contains information about",
              "  the covariates in 'brmcoda' model to estimate the substitution model.",
              "  It should not include ILR variables nor any column names starting with 'bilr', 'wilr', or 'ilr',",
              "  as these variables will be calculated by substitution model.",
              "  Please provide a different reference grid.",
              sep = "\n"))
          }
          refgrid <- as.data.table(cov.grid)
        } else {
          # grab covariates in user's specified reference grid
          # and use default if not supplied
          refgrid <- expand.grid.df(cov.grid,
                                    rg[, -gridnames, with = FALSE])}
      } 
    }
    
    # d0 ---------------------------
    # bilr0 is between-person ilr of the ref comp (doesn't have to be compositional mean)
    bcomp0 <- comp0
    bilr0 <- ilr(bcomp0, V = object$CompILR$psi)
    bilr0 <- as.data.table(t(bilr0))
    
    # wcomp0 and wilr0 are the difference between the actual compositional mean of the dataset and bilr
    # is 0 if ref comp is compositional mean
    # but is different if not
    wcomp0 <- bcomp0 - mcomp
    wilr0 <- as.data.table(t(ilr(wcomp0, V = object$CompILR$psi)))
    
    id <- data.table::data.table(1) # to make fitted() happy
    
    colnames(bilr0) <- colnames(object$CompILR$BetweenILR)
    colnames(wilr0) <- colnames(object$CompILR$WithinILR)
    colnames(bcomp0) <- colnames(object$CompILR$BetweenComp)
    colnames(wcomp0) <- colnames(object$CompILR$WithinComp)
    colnames(id) <- object$CompILR$idvar
    
    if(is.null(dim(refgrid))) {
      d0 <- cbind(bilr0, wilr0, bcomp0, wcomp0, id)
    } else {
      d0 <- expand.grid.df(bilr0, wilr0, bcomp0, wcomp0, id, refgrid)
    }}
  
  as.data.table(d0)
}

# expand grid data frame
expand.grid.df <- function(...) Reduce(function(...) merge.data.frame(..., by = NULL, all = TRUE), list(...))
