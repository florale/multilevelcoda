#' Checks if argument is a \code{substitution} object
#'
#' @param x An object of class \code{substitution}.
#'
#' @export
is.substitution <- function(x) {
  inherits(x, "substitution")
}

#' Constructor function for \code{substitution} class.
#'
#' @param between_simple_sub A list of results from \code{bsub} or \code{NULL}
#' @param between_avg_sub A list of results from \code{bsubmargins} or \code{NULL}
#' @param within_simple_sub A list of results from \code{wsub} or \code{NULL}
#' @param within_avg_sub A list of results from \code{wsubmargins} or \code{NULL}
#' @param simple_sub A list of results from \code{sub} or \code{NULL}
#' @param avg_sub A list of results from \code{submargins} or \code{NULL}
#' @param delta A numeric vector of the amount of substitution
#' @param ref A character value specifying the reference grid
#' @param level A character value specifying the level of substitution
#' @param weight The weight to use in calculation of the reference composition
#' @param parts The parts of the composition
#' @param summary A logical value specifying whether to summarize the results
#' 
#' @seealso \code{\link{substitution}}
#' 
#' @return An object of class \code{substitution}
#'
create_substitution <-
  function(between_simple_sub, within_simple_sub, simple_sub,
           between_avg_sub, within_avg_sub, avg_sub,
           delta,
           ref,
           level,
           weight,
           parts,
           summary) {
    
    stopifnot(is.list(between_simple_sub) || is.null(between_simple_sub))
    stopifnot(is.list(within_simple_sub) || is.null(within_simple_sub))
    stopifnot(is.list(simple_sub) || is.null(simple_sub))
    stopifnot(is.list(between_avg_sub) || is.null(between_avg_sub))
    stopifnot(is.list(within_avg_sub) || is.null(within_avg_sub))
    stopifnot(is.list(avg_sub) || is.null(avg_sub))
    
    out <- list(
      between_simple_sub = between_simple_sub,
      within_simple_sub = within_simple_sub,
      simple_sub = simple_sub,
      between_avg_sub = between_avg_sub,
      within_avg_sub = within_avg_sub,
      avg_sub = avg_sub,
      delta = delta,
      ref = ref,
      level = level,
      weight = weight,
      parts = parts,
      summary = summary
    )
    
    class(out) <- "substitution"
    
    return(out)
  }

#' Reference Grid for \code{substitution} model.
#' 
#' Build a dataset for \code{fitted.brmcoda} used in \code{substitution} model
#' 
#' @param object A fitted \code{\link{brmcoda}} object.
#' @param fill Logical value only relevant when \code{ref} is an user's specified reference grid
#' in which information about some, but not all covariates is provided
#' (e.g., models including age and sex as covariate but only age was provided in the reference grid).
#' If \code{TRUE}, the unspecified covariates are filled with the default reference grid.
#' If \code{FALSE}, users will be asked to provide a full reference grid.
#' Currently only support the default to \code{FALSE}.
#' @inheritParams substitution
#' 
#' @importFrom utils head
#' @importFrom data.table as.data.table copy :=
#' @importFrom compositions acomp ilr clo mean.acomp
#' @importFrom insight find_predictors
#' @importFrom emmeans ref_grid
#' @importFrom extraoperators %snin% %sin%
#'
#' @return A reference grid consisting of a combination of covariates in \code{brmcoda}
#' 
#' @export
build.rg <- function(object,
                     ref, 
                     level = NULL,
                     weight, 
                     fill = FALSE) {
  
  cov.grid <- NULL
  
  # get what type of model is being estimated
  model_fixef <- rownames(fixef(object))
  model_ranef <- names(ranef(object))
  
  if (isFALSE(length(grep("^[bilr]", model_fixef, value = T)) == 1)) {
    model_fixef_level <- "between"
    model_fixef_coef <- grep("^[bilr]", model_fixef, value = T)
  }
  if (isFALSE(length(grep("^[wilr]", model_fixef, value = T)) == 1)) {
    model_fixef_level <- "within"
    model_fixef_coef <- grep("^[wilr]", model_fixef, value = T)
  }
  if (isFALSE(length(grep("^[ilr]", model_fixef, value = T)) == 1)) {
    model_fixef_level <- "combined"
    model_fixef_coef <- grep("^[ilr]", model_fixef, value = T)
  }
  
  # single level or multilevel
  if (length(model_ranef) > 0) {
    model_ranef_level <- "multilevel"
    model_ranef_coef <- model_ranef
  } else {
    model_ranef_level <- "single"
    model_ranef_coef <- NULL
  }
  
  # d0 and comp0 for single level model
  if (model_ranef_level == "single") {
    comp0 <- object$complr$comp
    d0 <- object$complr$data
  }
  
  # d0 and comp0 for multilevel level model
  if (model_ranef_level == "multilevel") {
    
    # for combined
    if (level == "combined") {
      comp0 <- object$complr$comp
      d0 <- cbind(object$complr$comp, object$complr$data)
      d0 <- d0[, head(.SD, 1), by = eval(object$complr$idvar)]
      comp0 <- acomp(d0[, colnames(object$complr$comp), with = FALSE], total = object$complr$total)
    }
    
    # for between and/or within
    if (level %in% c("between", "within")) {
      if (isTRUE(ref == "clustermean")) {
        if (isTRUE(weight == "equal")) {
          comp0 <- object$complr$between_comp
          d0 <- cbind(object$complr$between_comp, object$complr$data)
          d0 <- d0[, head(.SD, 1), by = eval(object$complr$idvar)]
          comp0 <- acomp(d0[, colnames(object$complr$between_comp), with = FALSE], total = object$complr$total)
        }
        if (isTRUE(is.null(weight))) {
          comp0 <- object$complr$between_comp
          d0 <- object$complr$data
        }
        if (weight == "proportional") {
          comp0 <- object$complr$between_comp
          d0 <- object$complr$data
        }
        
        # d0 
        bcomp0 <- comp0
        bilr0 <- ilr(bcomp0, V = object$complr$psi)
        bilr0 <- as.data.table(bilr0)
        
        wcomp0 <- as.data.table(matrix(1, nrow = nrow(bcomp0), ncol = ncol(bcomp0)))
        wilr0 <- as.data.table(matrix(0, nrow = nrow(bilr0), ncol = ncol(bilr0)))
        
        colnames(bilr0) <- colnames(object$complr$between_logratio)
        colnames(wilr0) <- colnames(object$complr$within_logratio)
        colnames(bcomp0) <- colnames(object$complr$between_comp)
        colnames(wcomp0) <- colnames(object$complr$within_comp)
        
        d0 <- cbind(bilr0, wilr0, bcomp0, wcomp0,
                    d0[, colnames(d0) %in% colnames(object$complr$data), with = FALSE])
      }
      
      if (isTRUE(ref == "grandmean")) {
        if (isTRUE(is.null(weight) || weight == "equal")) {
          comp0 <- cbind(object$complr$between_comp, 
                         object$complr$data[, object$complr$idvar, with = FALSE])
          comp0 <- comp0[, head(.SD, 1), by = eval(object$complr$idvar)]
          comp0 <- acomp(comp0[, -object$complr$idvar, with = FALSE], total = object$complr$total)
          comp0 <- mean.acomp(comp0, robust = TRUE)
          comp0 <- acomp(comp0, total = object$complr$total)
          comp0 <- as.data.table(t(comp0))
          mcomp <- comp0
        }
        if (isTRUE(weight == "proportional")) {
          comp0 <- mean.acomp(object$complr$between_comp, robust = TRUE)
          comp0 <- acomp(comp0, total = object$complr$total)
          comp0 <- as.data.table(t(comp0))
          mcomp <- comp0
        }
        
        if (isTRUE(inherits(ref, c("data.table", "data.frame", "matrix")))) {
          # get compositional mean
          mcomp <- cbind(object$complr$between_comp, 
                         object$complr$data[, object$complr$idvar, with = FALSE])
          mcomp <- mcomp[, head(.SD, 1), by = eval(object$complr$idvar)]
          mcomp <- acomp(mcomp[, -object$complr$idvar, with = FALSE], total = object$complr$total)
          mcomp <- mean.acomp(mcomp, robust = TRUE)
          mcomp <- acomp(mcomp, total = object$complr$total)
          mcomp <- as.data.table(t(mcomp))
          
          if (isTRUE(identical(length(object$complr$parts), ncol(ref)))) {
            comp0 <- ref
            
          } else {
            if (object$complr$parts %nin% colnames(ref)) {
              stop(
                sprintf(
                  "The reference grid should include all compositional components but (%s) are missing.",
                  paste0(object$complr$parts %nin% colnames(ref), collapse = ", ")
                ))
            }
            comp0 <- ref[, object$complr$parts, with = FALSE]
            cov.grid <- ref[, -object$complr$parts, with = FALSE]
          }
          
          if(isFALSE(sum(comp0) == object$complr$total)) {
            stop(sprintf(
              "The total amount of the reference composition (%s) should be the same as the composition (%s).",
              sum(comp0),
              object$complr$total
            ))
          }
          if (isTRUE((any(comp0 > lapply(object$complr$data[, object$complr$parts, with = FALSE], max)) |
                      any(comp0 < lapply(object$complr$data[, object$complr$parts, with = FALSE], min))))) {
            stop(paste(
              sprintf(
                "comp0 should be numeric or interger values that are between (%s) and (%s)",
                paste0(round(apply(object$complr$data[, object$complr$parts, with = FALSE], 2, min)), collapse = ", "),
                paste0(round(apply(object$complr$data[, object$complr$parts, with = FALSE], 2, max)), collapse = ", ")),
              "\n", 
              " for",
              paste0(object$complr$parts, collapse = ", "),
              "respectively"
            ))
          }
          comp0  <- acomp(comp0, total = object$complr$total)
          comp0  <- as.data.table(t(comp0))
          colnames(comp0) <- colnames(object$complr$between_comp)
        }
        
        # reference grid 
        # get possible ilr names
        ilrnames <- c(colnames(object$complr$between_logratio), 
                      colnames(object$complr$within_logratio),
                      colnames(object$complr$ILR))
        
        # get all varnames in model
        varnames <- do.call(rbind, find_predictors(object$model))
        
        if (isTRUE(all(varnames %in% ilrnames))) { # unadj subsitution model
          # get varnames in model in ilr names
          refgrid <- NULL
          
        } else { # adj subsitution model
          # default reference grid
          rg <- as.data.table(ref_grid(object$model)@grid)
          covnames <- colnames(rg) %snin% c(ilrnames, ".wgt.")
          
          # user's specified reference grid
          if (isTRUE(is.null(cov.grid))) {
            refgrid <- rg[, covnames, with = FALSE]
          } else {
            gridnames <- colnames(cov.grid) %nin% object$complr$parts
            
            if (isFALSE(fill)) {
              if (isFALSE(identical(colnames(cov.grid), covnames))) {
                # ensure all covs are provided
                stop(paste(
                  "'cov.grid' should contains information about",
                  "  the covariates in 'brmcoda' model to estimate substitution",
                  "  It should not include ILR variables nor any column names starting with 'bilr', 'wilr', or 'ilr',",
                  "  as these variables will be computed in substitution analysis.",
                  "  Please provide a different reference grid.",
                  sep = "\n"))
              }
              refgrid <- as.data.table(cov.grid)
            } else {
              # grab covariates in user's specified reference grid
              # and use default if not supplied ######## fill rg here  
              refgrid <- expand.grid.df(ref,
                                        rg[, -gridnames, with = FALSE])}
          } 
        }
        
        # d0 
        # bilr0 is between-person ilr of the ref comp (doesn't have to be compositional mean)
        bcomp0 <- comp0
        bilr0 <- ilr(bcomp0, V = object$complr$psi)
        bilr0 <- as.data.table(t(bilr0))
        
        # wcomp0 and wilr0 are the difference between the actual compositional mean of the dataset and bilr
        # is 0 if ref comp is compositional mean
        # but is different if not
        wcomp0 <- bcomp0 - mcomp
        wilr0 <- as.data.table(t(ilr(wcomp0, V = object$complr$psi)))
        
        id <- data.table::data.table(1) # to make fitted() happy
        
        colnames(bilr0)  <- colnames(object$complr$between_logratio)
        colnames(wilr0)  <- colnames(object$complr$within_logratio)
        colnames(bcomp0) <- colnames(object$complr$between_comp)
        colnames(wcomp0) <- colnames(object$complr$within_comp)
        colnames(id)     <- object$complr$idvar
        
        if (is.null(dim(refgrid))) {
          d0 <- cbind(bilr0, wilr0, bcomp0, wcomp0, id)
        } else {
          d0 <- expand.grid.df(bilr0, wilr0, bcomp0, wcomp0, id, refgrid)
        }}
    }
  }
  as.data.table(d0)
}

#' Helper functions used only internally to estimate substitution model
#' @importFrom data.table as.data.table data.table copy := setDT rbindlist
#' @importFrom compositions acomp ilr clo mean.acomp
#' @importFrom bayestestR describe_posterior
#' @importFrom extraoperators %snin% %sin%
#' @importFrom foreach foreach %dopar%
#' @importFrom doFuture %dofuture%
#' 
#' @name get-substitution
NULL

# Grandmean Between-person Substitution model
.get.bsub <- function(object, delta, basesub,
                      comp0, y0, d0,
                      summary,
                      level, ref,
                      ...) {
  
  iout <- foreach(i = colnames(basesub), .combine = c, 
                  .options.future = list(packages = "multilevelcoda")) %dofuture% {
                    
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
                        names(newcomp) <- object$complr$parts
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
                    bcomp0 <- acomp(dnew[, colnames(object$complr$between_comp), with = FALSE], total = object$complr$total)
                    bcompsub  <- acomp(dnew[, object$complr$parts, with = FALSE], total = object$complr$total)
                    
                    bilrsub <- ilr(bcompsub, V = object$complr$psi)
                    wilr0 <- as.data.table(matrix(0, nrow = nrow(bilrsub), ncol = ncol(bilrsub)))
                    
                    colnames(bilrsub) <- colnames(object$complr$between_logratio)
                    colnames(wilr0) <- colnames(object$complr$within_logratio)
                    
                    # reference grid 
                    ## get covariate + idvar names
                    covnames <- colnames(d0) %snin% c(colnames(object$complr$between_logratio),
                                                      colnames(object$complr$within_logratio),
                                                      colnames(object$complr$between_comp),
                                                      colnames(object$complr$within_comp)
                    )
                    refgrid <- d0[, covnames, with = FALSE]
                    
                    # predictions
                    hout <- vector("list", length = nrow(refgrid))
                    if(isTRUE(summary)) { # unadj OR adj averaging over reference grid
                      for (h in seq_len(nrow(refgrid))) {
                        dsub <- cbind(dnew, bilrsub, wilr0, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        hout[[h]] <- delta_y
                      }
                      
                      PD_delta_y <- Reduce(`+`, hout) / length(hout)
                      suppressWarnings(
                        PD_delta_y <- apply(delta_y, 2, function(x) {
                          describe_posterior(x, centrality = "mean", ...)
                        }))
                      PD_delta_y <- rbindlist(PD_delta_y)
                      PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)], 
                                          dsub[, .(Delta, From, To, Level, Reference)])
                      
                    } else { # adj keeping prediction at each level of reference grid
                      for (h in seq_len(nrow(refgrid))) {
                        dsub <- cbind(dnew, bilrsub, wilr0, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        suppressWarnings(
                          PD_delta_y <- apply(delta_y, 2, function(x) {
                            describe_posterior(x, centrality = "mean", ...)
                          }))
                        PD_delta_y <- rbindlist(PD_delta_y)
                        PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)], 
                                            dsub[, .(Delta, From, To, Level, Reference)],
                                            dsub[, colnames(refgrid) %snin% object$complr$idvar, with = FALSE])
                        
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
.get.wsub <- function(object, delta, basesub,
                      comp0, y0, d0,
                      summary,
                      level, ref,
                      ...) {
  
  iout <- foreach(i = colnames(basesub), .combine = c, 
                  .options.future = list(packages = "multilevelcoda")) %dofuture% {
                    
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
                        names(newcomp) <- object$complr$parts
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
                    bcomp0   <- acomp(dnew[, colnames(object$complr$between_comp), with = FALSE], total = object$complr$total)
                    bcompsub <- acomp(dnew[, object$complr$parts, with = FALSE], total = object$complr$total)
                    
                    bilr0   <- ilr(bcomp0, V = object$complr$psi)
                    bilrsub <- ilr(bcompsub, V = object$complr$psi)
                    wilrsub <- bilrsub - bilr0
                    
                    colnames(bilr0) <- colnames(object$complr$between_logratio)
                    colnames(wilrsub) <- colnames(object$complr$within_logratio)
                    
                    # reference grid 
                    ## get covariate + idvar names
                    covnames <- colnames(d0) %snin% c(colnames(object$complr$between_logratio),
                                                      colnames(object$complr$within_logratio),
                                                      colnames(object$complr$between_comp),
                                                      colnames(object$complr$within_comp)
                    )
                    refgrid <- d0[, covnames, with = FALSE]
                    
                    # predictions
                    hout <- vector("list", length = nrow(refgrid))
                    if(isTRUE(summary)) { # unadj OR adj averaging over reference grid
                      for (h in seq_len(nrow(refgrid))) {
                        dsub <- cbind(dnew, bilr0, wilrsub, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        hout[[h]] <- delta_y
                      }
                      
                      PD_delta_y <- Reduce(`+`, hout) / length(hout)
                      suppressWarnings(
                        PD_delta_y <- apply(delta_y, 2, function(x) {
                          describe_posterior(x, centrality = "mean", ...)
                        }))
                      PD_delta_y <- rbindlist(PD_delta_y)
                      PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)],
                                          dsub[, .(Delta, From, To, Level, Reference)])
                      
                    } else { # adj keeping prediction at each level of reference grid
                      for (h in seq_len(nrow(refgrid))) {
                        dsub <- cbind(dnew, bilr0, wilrsub, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        suppressWarnings(
                          PD_delta_y <- apply(delta_y, 2, function(x) {
                            describe_posterior(x, centrality = "mean", ...)
                          }))
                        PD_delta_y <- rbindlist(PD_delta_y)
                        PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)], 
                                            dsub[, .(Delta, From, To, Level, Reference)],
                                            dsub[, colnames(refgrid) %snin% object$complr$idvar, with = FALSE])
                        
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

# Grandmean Simple Substitution
.get.sub <- function(object, delta, basesub,
                     comp0, y0, d0,
                     summary,
                     level, ref,
                     ...) {
  
  iout <- foreach(i = colnames(basesub), .combine = c, 
                  .options.future = list(packages = "multilevelcoda")) %dofuture% {
                    
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
                        names(newcomp) <- object$complr$parts
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
                    comp0 <- acomp(dnew[, colnames(object$complr$comp), with = FALSE], total = object$complr$total)
                    compsub  <- acomp(dnew[, object$complr$parts, with = FALSE], total = object$complr$total)
                    
                    ilrsub <- ilr(compsub, V = object$complr$psi)
                    colnames(ilrsub) <- colnames(object$complr$ILR)
                    
                    # reference grid 
                    ## get covariate + idvar names
                    covnames <- colnames(d0) %snin% c(colnames(object$complr$ILR),
                                                      colnames(object$complr$comp)
                    )
                    refgrid <- d0[, covnames, with = FALSE]
                    
                    # predictions
                    hout <- vector("list", length = nrow(refgrid))
                    if(isTRUE(summary)) { # unadj OR adj averaging over reference grid
                      for (h in seq_len(nrow(refgrid))) {
                        dsub <- cbind(dnew, ilrsub, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        hout[[h]] <- delta_y
                      }
                      
                      PD_delta_y <- Reduce(`+`, hout) / length(hout)
                      suppressWarnings(
                        PD_delta_y <- apply(delta_y, 2, function(x) {
                          describe_posterior(x, centrality = "mean", ...)
                        }))
                      PD_delta_y <- rbindlist(PD_delta_y)
                      PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)], 
                                          dsub[, .(Delta, From, To, Level, Reference)])
                      
                    } else { # adj keeping prediction at each level of reference grid
                      for (h in seq_len(nrow(refgrid))) {
                        dsub <- cbind(dnew, ilrsub, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        suppressWarnings(
                          PD_delta_y <- apply(delta_y, 2, function(x) {
                            describe_posterior(x, centrality = "mean", ...)
                          }))
                        PD_delta_y <- rbindlist(PD_delta_y)
                        PD_delta_y <- cbind(PD_delta_y[, .(Mean, CI_low, CI_high)], 
                                            dsub[, .(Delta, From, To, Level, Reference)],
                                            dsub[, colnames(refgrid) %snin% object$complr$idvar, with = FALSE])
                        
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

# Clustermean Between-person Substitution model
.get.bsubmargins <- function(object, delta, basesub,
                             comp0, y0, d0,
                             level, ref,
                             ...) {
  
  iout <- foreach(i = colnames(basesub), .combine = c, 
                  .options.future = list(packages = "multilevelcoda")) %dofuture% {
                    
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
                        names(newcomp) <- object$complr$parts
                        
                        dnew <- cbind(comp0, newcomp, 
                                      d0[, colnames(d0) %in% colnames(object$complr$data[, -object$complr$part, with = FALSE]), with = FALSE], 
                                      Delta)
                        # useful information for the final results
                        dnew[, From := rep(subvar, length.out = nrow(dnew))[k]]
                        dnew$To <- iv
                        dnew$Delta <- as.numeric(dnew$Delta)
                        
                        # remove impossible reallocation that result in negative values 
                        cols <- colnames(dnew) %sin% c(colnames(comp0), colnames(basesub))
                        dnew <- dnew[rowSums(dnew[, ..cols] < 0) == 0]
                        
                        # compositions and ilrs for predictions
                        bcomp0    <- acomp(dnew[, colnames(object$complr$between_comp), with = FALSE], total = object$complr$total)
                        bcompsub  <- acomp(dnew[, object$complr$parts, with = FALSE], total = object$complr$total)
                        
                        bilr0    <- ilr(bcomp0, V = object$complr$psi)
                        bilrsub  <- ilr(bcompsub, V = object$complr$psi)
                        
                        wilr0 <- as.data.table(matrix(0, nrow = nrow(bilrsub), ncol = ncol(bilrsub)))
                        
                        colnames(bilrsub) <- colnames(object$complr$between_logratio)
                        colnames(wilr0)   <- colnames(object$complr$within_logratio)
                        
                        # prediction
                        dsub <- cbind(dnew, bilrsub, wilr0)
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NULL,
                            summary = FALSE
                          )
                        ysub <- rowMeans(as.data.frame(ysub))
                        
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

# Clustermean Within-person Substitution model
.get.wsubmargins <- function(object, delta, basesub,
                             comp0, y0, d0,
                             level, ref,
                             ...) {
  
  iout <- foreach(i = colnames(basesub), .combine = c, 
                  .options.future = list(packages = "multilevelcoda")) %dofuture% {    
                    
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
                        names(newcomp) <- object$complr$parts
                        
                        dnew <- cbind(comp0, newcomp, 
                                      d0[, colnames(d0) %in% colnames(object$complr$data[, -object$complr$part, with = FALSE]), with = FALSE], 
                                      Delta)
                        
                        # useful information for the final output
                        dnew[, From := rep(subvar, length.out = nrow(dnew))[k]]
                        dnew$To <- iv
                        dnew$Delta <- as.numeric(dnew$Delta)
                        
                        # remove impossible reallocation that result in negative values 
                        cols <- colnames(dnew) %sin% c(colnames(comp0), colnames(basesub))
                        dnew <- dnew[rowSums(dnew[, ..cols] < 0) == 0]
                        
                        # compositions and ilr for predictions
                        bcomp0   <- acomp(dnew[, colnames(object$complr$between_comp), with = FALSE], total = object$complr$total)
                        bcompsub <- acomp(dnew[, object$complr$parts, with = FALSE], total = object$complr$total)
                        
                        bilr0   <- ilr(bcomp0, V = object$complr$psi)
                        bilrsub <- ilr(bcompsub, V = object$complr$psi)
                        
                        wilrsub <- bilrsub - bilr0 
                        
                        colnames(bilr0)   <- colnames(object$complr$between_logratio)
                        colnames(wilrsub) <- colnames(object$complr$within_logratio)
                        
                        # substitution data
                        dsub <- cbind(dnew, bilr0, wilrsub)
                        
                        # prediction
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NULL,
                            summary = FALSE
                          )
                        ysub <- rowMeans(as.data.frame(ysub))
                        
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

# Clustermean Average Substitution
.get.submargins <- function(object, delta, basesub,
                            comp0, y0, d0,
                            level, ref,
                            ...) {
  
  iout <- foreach(i = colnames(basesub), .combine = c, 
                  .options.future = list(packages = "multilevelcoda")) %dofuture% {    
                    
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
                        subk <- subk[rep(seq_len(nrow(subk)), nrow(comp0)), ]
                        newcomp <- comp0 + subk
                        Delta <- subk[, get(i)]
                        names(newcomp) <- object$complr$parts
                        
                        dnew <- cbind(newcomp, object$complr$data, Delta)
                        
                        # useful information for the final results
                        dnew[, From := rep(subvar, length.out = nrow(dnew))[k]]
                        dnew$To <- iv
                        dnew$Delta <- as.numeric(dnew$Delta)
                        
                        # remove impossible reallocation that result in negative values 
                        cols <- colnames(dnew) %sin% c(colnames(comp0), colnames(basesub))
                        dnew <- dnew[rowSums(dnew[, ..cols] < 0) == 0]
                        
                        # compositions and ilrs for predictions
                        tcomp <- acomp(dnew[, object$complr$parts, with = FALSE], total = object$complr$total)
                        tilr <- ilr(tcomp, V = object$complr$psi)
                        
                        colnames(tilr) <- colnames(object$complr$ILR)
                        
                        # substitution data
                        dsub <- cbind(dnew, tilr)
                        
                        # prediction
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NULL,
                            summary = FALSE
                          )
                        ysub <- rowMeans(as.data.frame(ysub))
                        
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