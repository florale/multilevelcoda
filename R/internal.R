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
# Grandmean Between-person Substitution model
get.bsub <- function(object, basesub, refcomp, 
                     y0, d0,
                     delta, summary,
                     level, type,
                     ID = 1, refgrid = refgrid, ...) {
  
  iout <- foreach(i = colnames(basesub), .combine = c) %dopar% {
    
    # dnew - reallocation data ----------------------------------------------
    # possible substitution of 1 compositional variable
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
        newcomp <- refcomp + sub[k,]
        names(newcomp) <- object$CompILR$parts
        Delta <- sub[k, get(i)]
        kout[[k]] <- cbind(refcomp, newcomp, Delta)
      }
      jout[[j]] <- do.call(rbind, kout)
    }
    dnew <- setDT(do.call(rbind, jout))

    # useful information for the final results
    dnew[, From := rep(subvar, length.out = nrow(dnew))]
    dnew$To <- iv
    dnew$Delta <- as.numeric(dnew$Delta)
    dnew$Level <- level

    # remove impossible reallocation that result in negative values 
    cols <- colnames(dnew) %snin% c("Delta", "From", "To", "Level")
    dnew <- dnew[rowSums(dnew[, ..cols] < 0) == 0]
    
    # compositions and ilrs for predictions
    bcomp <- acomp(dnew[, colnames(object$CompILR$BetweenComp), with = FALSE], total = object$CompILR$total)
    tcomp <- acomp(dnew[, object$CompILR$parts, with = FALSE], total = object$CompILR$total)
    
    bilr <- ilr(bcomp, V = object$CompILR$psi)
    tilr <- ilr(tcomp, V = object$CompILR$psi)
    
    wilr <- d0[, colnames(object$CompILR$WithinILR), with = FALSE]
    
    colnames(tilr) <- colnames(object$CompILR$BetweenILR)
    colnames(wilr) <- colnames(object$CompILR$WithinILR)
    
    # prediction of posterior draws ------------------------
    hout <- vector("list", length = nrow(refgrid))
    if(isTRUE(by.grid)) { # keeping prediction at each level of reference grid
      for (h in seq_len(nrow(refgrid))) {
        dsub <- cbind(dnew, tilr, wilr, ID, refgrid[h, ])
        ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
        delta_y <- ysub - y0[, h]
        
        suppressWarnings(PD_delta_y <- apply(delta_y, 2, function(x)
          {describe_posterior(x, centrality = "mean", ...)}))
        PD_delta_y <- rbindlist(PD_delta_y)
        PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)],
                            dsub[, c("Delta", "From", "To", "Level", colnames(refgrid)), with = FALSE])
        hout[[h]] <- list(delta_y = delta_y,
                          PD_delta_y = PD_delta_y)
      }} else {
        dsub <- expand.grid.df(dnew, tilr, wilr, ID, refgrid)
        ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
        delta_y <- ysub - y0
        
        suppressWarnings(PD_delta_y <- apply(delta_y, 2, function(x)
          {describe_posterior(x, centrality = "mean", ...)}))
        PD_delta_y <- rbindlist(PD_delta_y)
        PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)],
                            dsub[, c("Delta", "From", "To", "Level", colnames(refgrid)), with = FALSE])
        hout <- list(delta_y = delta_y,
                     PD_delta_y = PD_delta_y)
      }

    # final results for entire composition
    out <- if(isTRUE(summary)) (hout$PD_delta_y) else (hout$delta_y)
    names(out) <- i
    out
  }
  iout
}

# Grandmean Within-person Substitution model
get.wsub <- function(object, basesub, refcomp,
                     y0, d0,
                     delta, summary, 
                     level, type,
                     ID = 1, refgrid = NULL, ...) {
  
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
        newcomp <- refcomp + sub[k, ]
        names(newcomp) <- object$CompILR$parts
        Delta <- sub[k, get(i)]
        kout[[k]] <- cbind(refcomp, newcomp, Delta)
      }
      jout[[j]] <- do.call(rbind, kout)
    }
    dnew <- setDT(do.call(rbind, jout))
    
    # useful information for the final results
    dnew[, From := rep(subvar, length.out = nrow(dnew))]
    dnew$To <- iv
    dnew$Delta <- as.numeric(dnew$Delta)
    dnew$Level <- level

    # remove impossible reallocation that result in negative values 
    cols <- colnames(dnew) %snin% c("Delta", "From", "To", "Level")
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
    if(is.null(refgrid)) { # unadjusted
      dsub <- cbind(dnew, bilr, wilr, ID)
      ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
      
      delta_y <- apply(ysub, 2, function(y) {y - y0})
      suppressWarnings(PD_delta_y <- apply(delta_y, 2, function(x) {describe_posterior(x, centrality = "mean", ...)}))
      PD_delta_y <- rbindlist(PD_delta_y)
      PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)], 
                     dsub[, .(Delta, To, From, Level)])
      
    } else { # adjusted
      hout <- vector("list", length = nrow(refgrid))
      if(isTRUE(summary)) { # averaging over reference grid
        for (h in seq_len(nrow(refgrid))) {
          dsub <- cbind(dnew, bilr, wilr, ID, refgrid[h, ])
          ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
          delta_y <- ysub - y0[, h]
          hout[[h]] <- delta_y
        }
        
        PD_delta_y <- Reduce(`+`, hout) / length(hout)
        suppressWarnings(PD_delta_y <- apply(delta_y, 2, function(x) {describe_posterior(x, centrality = "mean", ...)}))
        PD_delta_y <- rbindlist(PD_delta_y)
        PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)], 
                       dsub[, .(Delta, To, From, Level)])
        
      } else { # keeping prediction at each level of reference grid
        for (h in seq_len(nrow(refgrid))) {
          dsub <- cbind(dnew, bilr, wilr, ID, refgrid[h, ])
          ysub <- fitted(object$Model, newdata = dsub, re_formula = NA, summary = FALSE)
          delta_y <- ysub - y0[, h]
          suppressWarnings(PD_delta_y <- apply(delta_y, 2, function(x) {describe_posterior(x, centrality = "mean", ...)}))
          PD_delta_y <- rbindlist(PD_delta_y)
          PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)],
                         dsub[, c("Delta", "From", "To", "Level", colnames(refgrid)),
                              with = FALSE])
          
          hout[[h]] <- PD_delta_y
        }
        PD_delta_y <- rbindlist(hout)
      }
    }
    # final results for entire composition
    out <- list(PD_delta_y)
    names(out) <- i
    out
  }
  iout
}

# UnitMean Between-person Substitution Model.
.get.bsubmargins <- function(object, basesub, b, 
                             y0, delta, 
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
      for (k in seq_len(nrow(sub))) { # reallocation level
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

    names(jout) <- c("Mean", "CI_low", "CI_high", 
                     "Delta", "To", "From", "Level")
    
    # store final results for entire composition
    jout <- list(jout)
    names(jout) <- i
    jout
  }
  iout
}

# Unitmean Within-person Substitution Model.
.get.wsubmargins <- function(object, basesub, b,
                             y0, delta, 
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

    names(jout) <- c("Mean", "CI_low", "CI_high", 
                     "Delta", "To", "From", "Level")
    
    # final results for entire composition
    jout <- list(jout)
    names(jout) <- i
    jout
  }
  iout
}

# Unitmean Substitution Model.
.get.submargins <- function(object, basesub, t,
                            y0, delta,
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

# expand grid data frame
expand.grid.df <- function(...) Reduce(function(...) merge(..., by = NULL, all = TRUE), list(...))

# reference dataset
build.rg <- function(object, 
                     ref = c("grandmean", "unitmean"),
                     weight = NULL, build = FALSE) {
  
  if(isTRUE(ref == "unitmean")) {
    refcomp <- object$CompILR$BetweenComp
  }
  
  if(identical(ref, "grandmean")) {
    if(isTRUE(is.null(weight) || weight == "equal")) {
      refcomp <- cbind(object$CompILR$BetweenComp, 
                       object$CompILR$data[, object$CompILR$idvar, with = FALSE])
      refcomp <- refcomp[, .SD[c(1)], by = eval(object$CompILR$idvar)]
      refcomp <- acomp(refcomp[, -object$CompILR$idvar, with = FALSE], total = object$CompILR$total)
      refcomp <- mean(refcomp, robust = TRUE)
      refcomp <- acomp(refcomp, total = object$CompILR$total)
      refcomp <- as.data.table(t(refcomp))
      mcomp <- refcomp
    }
    else if (isTRUE(weight == "proportional")) {
      # compositional mean
      refcomp <- mean(object$CompILR$BetweenComp, robust = TRUE)
      refcomp <- acomp(refcomp, total = object$CompILR$total)
      refcomp <- as.data.table(t(refcomp))
    }
  }
  
  if(isTRUE(inherits(ref, c("data.table", "data.frame", "matrix")))) {
    if (isTRUE(identical(length(object$CompILR$parts), ncol(refgrid)))) {
      refcomp <- refgrid
    } else {
      if (object$CompILR$parts %nin% colnames(refgrid)) {
        stop(
          sprintf(
            "The reference grid should include all compositional components but (%s) are missing.",
            paste0(object$CompILR$parts %nin% colnames(refgrid), collapse = ", ")
          ))
      }
      refcomp <- refgrid[, object$CompILR$parts, with = FALSE]
    }
    
    if(isFALSE(sum(refcomp) == object$CompILR$total)) {
      stop(sprintf(
        "The total amount of the reference composition (%s) should be the same as the composition (%s).",
        sum(refcomp),
        object$CompILR$total
      ))
    }
    if (isTRUE((any(refcomp > lapply(object$CompILR$data[, object$CompILR$parts, with = FALSE], max)) |
                any(refcomp < lapply(object$CompILR$data[, object$CompILR$parts, with = FALSE], min))))) {
      stop(paste(
        sprintf(
          "refcomp should be numeric or interger values that are between (%s) and (%s)",
          paste0(round(apply(object$CompILR$data[, object$CompILR$parts, with = FALSE], 2, min)), collapse = ", "),
          paste0(round(apply(object$CompILR$data[, object$CompILR$parts, with = FALSE], 2, max)), collapse = ", ")),
        "\n", 
        " for",
        paste0(object$CompILR$parts, collapse = ", "),
        "respectively"
      ))
    }
    refcomp  <- compositions::acomp(refcomp, total = object$CompILR$total)
    refcomp  <- as.data.table(t(refcomp))
    colnames(refcomp) <- colnames(object$CompILR$BetweenComp)
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
    if (isFALSE(is.null(refgrid))) {
      gridnames <- colnames(refgrid) %snin% object$CompILR$parts
      if(isFALSE(build)) {
        if(isFALSE(identical(colnames(refgrid), covnames))) { # ensure all covs are provided
          stop(paste(
            "'refgrid' should contains information about",
            "  the covariates in 'brmcoda' model to estimate the substitution model.",
            "  It should not include ILR variables nor any column names starting with 'bilr', 'wilr', or 'ilr',",
            "  as these variables will be calculated by substitution model.",
            "  Please provide a different reference grid.",
            sep = "\n"))
        }
        refgrid <- as.data.table(refgrid)
      } else {
        # grab covariates in user's specified reference grid
        # and use default if not supplied
        refgrid <- expand.grid.df(refgrid,
                                  rg[, -gridnames, with = FALSE])}
    } else {
      refgrid <- rg[, covnames, with = FALSE]
    }
  }
  
  # d0 ---------------------------
  # bilr is between-person ilr of the ref comp (doesn't have to be compositional mean)
  bilr0 <- ilr(refcomp, V = object$CompILR$psi)
  bilr0 <- as.data.frame(t(bilr0))
  
  # wcomp and wilr are the difference between the actual compositional mean of the dataset and bilr
  # is 0 if ref comp is compositional mean
  # but is different if not
  wcomp <- refcomp - mcomp
  wilr0 <- as.data.frame(t(ilr(wcomp, V = object$CompILR$psi)))
  
  colnames(bilr0) <- colnames(object$CompILR$BetweenILR)
  colnames(wilr0) <- colnames(object$CompILR$WithinILR)
  
  ID <- 1 # to make fitted() happy
  
  if(is.null(dim(refgrid))) {
    cbind(bilr0, wilr0, ID)
  } else {
    expand.grid.df(bilr0, wilr0, ID, refgrid)
  }
}