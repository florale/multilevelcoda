#' Checks if argument is a \code{substitution} object
#'
#' @param x An object of class \code{substitution}.
#'
#' @export
is.substitution <- function(x) {
  inherits(x, "substitution")
}

#' Reference grid for \code{substitution} model.
#' 
#' Build a dataset for \code{fitted.brmcoda} used in \code{substitution} model
#' 
#' @param object A fitted \code{\link{brmcoda}} object.
#' @param ref Either a character value or vector or a dataset.
#' Can be \code{"grandmean"} and/or \code{"clustermean"}, or
#' a \code{data.frame} or \code{data.table} of user's specified reference grid consisting
#' of combinations of covariates over which predictions are made.
#' User's specified reference grid only applicable to substitution model
#' using a single reference composition value.
#' @param weight A character value specifying the weight to use in calculation of the reference composition.
#' If \code{"equal"}, give equal weight to units (e.g., individuals).
#' If \code{"proportional"}, weights in proportion to the frequencies of units being averaged 
#' (e.g., observations across individuals)
#' Default is \code{equal}.
#' @param FALSE Logical value only relevant when \code{ref} is an user's specified reference grid
#' in which information about some, but not all covariates is provided
#' (e.g., models including age and sex as covariate but only age was provided in the reference grid).
#' If \code{TRUE}, the unspecified covariates are filled with the default reference grid.
#' If \code{FALSE}, users will be asked to provide a full reference grid.
#' Default is \code{FALSE}.
#' 
#' @return A reference grid consisting of a combination of covariates in \code{brmcoda}.
#' 
build.rg <- function(object, 
                     ref,
                     weight, 
                     fill = FALSE) {
  
  cov.grid <- NULL
  
  if (isTRUE(ref == "clustermean")) {
    if (isTRUE(weight == "equal")) {
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
  
  if (isTRUE(ref == "grandmean")) {
    if (isTRUE(is.null(weight) || weight == "equal")) {
      comp0 <- cbind(object$CompILR$BetweenComp, 
                     object$CompILR$data[, object$CompILR$idvar, with = FALSE])
      comp0 <- comp0[, head(.SD, 1), by = eval(object$CompILR$idvar)]
      comp0 <- acomp(comp0[, -object$CompILR$idvar, with = FALSE], total = object$CompILR$total)
      comp0 <- mean.acomp(comp0, robust = TRUE)
      comp0 <- acomp(comp0, total = object$CompILR$total)
      comp0 <- as.data.table(t(comp0))
      mcomp <- comp0
    }
    if (isTRUE(weight == "proportional")) {
      comp0 <- mean.acomp(object$CompILR$BetweenComp, robust = TRUE)
      comp0 <- acomp(comp0, total = object$CompILR$total)
      comp0 <- as.data.table(t(comp0))
      mcomp <- comp0
    }
    
    if (isTRUE(inherits(ref, c("data.table", "data.frame", "matrix")))) {
      # get compositional mean
      mcomp <- cbind(object$CompILR$BetweenComp, 
                     object$CompILR$data[, object$CompILR$idvar, with = FALSE])
      mcomp <- mcomp[, head(.SD, 1), by = eval(object$CompILR$idvar)]
      mcomp <- acomp(mcomp[, -object$CompILR$idvar, with = FALSE], total = object$CompILR$total)
      mcomp <- mean.acomp(mcomp, robust = TRUE)
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
        
        if (isFALSE(fill)) {
          if (isFALSE(identical(colnames(cov.grid), covnames))) {
            # ensure all covs are provided
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
    
    colnames(bilr0)  <- colnames(object$CompILR$BetweenILR)
    colnames(wilr0)  <- colnames(object$CompILR$WithinILR)
    colnames(bcomp0) <- colnames(object$CompILR$BetweenComp)
    colnames(wcomp0) <- colnames(object$CompILR$WithinComp)
    colnames(id)     <- object$CompILR$idvar
    
    if (is.null(dim(refgrid))) {
      d0 <- cbind(bilr0, wilr0, bcomp0, wcomp0, id)
    } else {
      d0 <- expand.grid.df(bilr0, wilr0, bcomp0, wcomp0, id, refgrid)
    }}
  
  as.data.table(d0)
}

# expand grid data frame
expand.grid.df <- function(...) Reduce(function(...) merge.data.frame(..., by = NULL, all = TRUE), list(...))

# check sequence of number
is.sequential <- function(x) {
  all(length(x) > 2 & all(abs(diff(x)) == 1))
}