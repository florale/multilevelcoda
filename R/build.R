#' Build Base Pairwise Substitution 
#'
#' @description
#' Make a data set of all possible pairwise substitution of a composition which can be used as 
#' the base for substitution models.
#' @param parts A character vector specifying the names of compositional variables to be used.
#' @param type Either \code{"one-to-one"} or \code{"one-to-all"}. Default is \code{"one-to-one"}.
#' 
#' @return A data table of all possible pairwise substitution.
#' @importFrom data.table as.data.table 
#' @export
#' @examples
#' ps1 <- build.base(parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' print(ps1)
#' 
#' ps2 <- build.base(c("WAKE", "MVPA", "LPA", "SB"), type = "one-to-all")
#' print(ps2)
build.base <- function(parts, type = NULL) {
  
  if (is.null(type)) {
    type <- "one-to-one"
  }
  
  d <- length(parts)
  n <- d - 2
  
  subvar1 <- c(1, -1)
  subvar2 <- rep(0, n)
  subvar <- c(subvar1, subvar2)
  
  nc <- length(subvar)
  nr <- (nc - 1) * d
  k <- 0
  
  if (type == "one-to-all") {
    base_sub_to <- matrix(0, nrow = d, ncol = d, dimnames = list(NULL, parts))
    for (i in 1:nc)
      for (j in 1:nc)
        if (i == j) {
          base_sub_to[i, j] <- 1
          base_sub_to[i, -j] <- -(1/(d-1))
        }
    
    base_sub_from <- matrix(0, nrow = d, ncol = d, dimnames = list(NULL, parts))
    for (i in 1:nc)
      for (j in 1:nc)
        if (i == j) {
          base_sub_from[i, j] <- -1
          base_sub_from[i, -j] <- (1/(d-1))
        }
    
    base_sub <- rbind(base_sub_to, base_sub_from)
    
  } else {
    base_sub <- matrix(0, nrow = nr, ncol = nc, dimnames = list(NULL, parts))
    
    for (i in 1:nc)
      for (j in 1:nc)
        if (i != j) {
          k <- k + 1
          base_sub[k, c(i, j)] <- c(1, -1)
        }
  }
  colnames(base_sub) <- parts
  as.data.table(base_sub)
}

#' Build Sequential Binary Partition
#'
#' @description
#' Build a default sequential binary partition for \code{complr} object.
#' The default sequential binary partition is a pivot balance that allows 
#' the effect of this first balance coordinate to be interpreted as the change 
#' in the prediction for the dependent variable 
#' when that given part increases while all remaining parts decrease by a common proportion.
#' @param parts A character vector specifying the names of compositional variables to be used.
#' 
#' @return A matrix sequential binary partition.
#' @importFrom data.table as.data.table 
#' @export
#' @examples
#' sbp1 <- build.sbp(parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' print(sbp1)
#' 
#' sbp2 <- build.sbp(c("WAKE", "MVPA", "LPA", "SB"))
#' print(sbp2)
build.sbp <- function(parts) {
  
  d <- length(parts)
  k <- 0
  
  nc <- d
  nr <- d - 1
  
  base_sbp <- matrix(NA, nrow = nr, ncol = nc, dimnames = list(NULL, parts))
  sbp <- base_sbp
  
  for (i in 1:nr) {
    base_sbp[i,  i] <- 1
    base_sbp[i, -i] <- -1
    base_sbp[i, 0:(i-1)] <- 0
  }
  base_sbp
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
#' @importFrom emmeans ref_grid
#' @importFrom extraoperators %snin% %sin%
#'
#' @return A reference grid consisting of a combination of covariates in \code{brmcoda}
#'
#' @export
build.rg <- function(object,
                     ref,
                     at,
                     parts,
                     level,
                     weight,
                     fill = FALSE) {
  
  covgrid <- NULL
  
  # get part names
  if (is.numeric(parts)) {
    parts <- .get_parts(object[["complr"]], parts)
  }
  
  # get the index of which index elements of object[["complr"]][["output"]] does the parts correspond to
  idx <- as.integer(which(vapply(lapply(object[["complr"]][["output"]], function(x)
    x$parts), function(p)
      identical(sort(parts), sort(p)), logical(1))))[1]
  
  # grab logratio and composition names
  Xz_vars  <- get_variables(object$complr)[[paste0("composition_", idx)]][["Z"]]
  Xbz_vars <- get_variables(object$complr)[[paste0("composition_", idx)]][["bZ"]]
  Xwz_vars <- get_variables(object$complr)[[paste0("composition_", idx)]][["wZ"]]
  
  Xtx_vars  <- get_variables(object$complr)[[paste0("composition_", idx)]][["X"]]
  Xbx_vars <- get_variables(object$complr)[[paste0("composition_", idx)]][["bX"]]
  Xwx_vars <- get_variables(object$complr)[[paste0("composition_", idx)]][["wX"]]
  
  ## NOTES
  ## ignore weight for clustermean
  ## equal weight is default for grandmean
  
  # what type of model is being estimated
  brmcoda_vars <- get_variables(object)
  
  # d0 and x0 for multilevel level model
  if (identical(brmcoda_vars$ranef_type, "multilevel")) {
    
    # for clustermean
    if ("clustermean" %in% ref) {
      weight <- NULL # ignore weight
      
      ## aggregate
      if ("aggregate" %in% level) {
        d0 <- object[["complr"]][["dataout"]][, head(.SD, 1), by = eval(object[["complr"]][["idvar"]])]
        x0 <- acomp(d0[, Xtx_vars, with = FALSE], total = object[["complr"]][["output"]][[idx]][["total"]])
        
        z0 <- ilr(x0, V = object[["complr"]][["output"]][[idx]][["psi"]])
        z0 <- as.data.table(z0)
        
        colnames(z0) <- Xz_vars
        colnames(x0) <- Xtx_vars
        
        d0 <- cbind(z0, x0, d0[, -colnames(x0), with = FALSE])
      }
      
      ## between and within
      if (any(c("between", "within") %in% level)) {
        d0   <- object[["complr"]][["dataout"]][, head(.SD, 1), by = eval(object[["complr"]][["idvar"]])]
        bx0  <- acomp(d0[, Xbx_vars, with = FALSE], total = object[["complr"]][["output"]][[idx]][["total"]])
        
        bz0 <- ilr(bx0, V = object[["complr"]][["output"]][[idx]][["psi"]])
        bz0 <- as.data.table(bz0)
        
        wx0 <- as.data.table(matrix(1, nrow = nrow(bx0), ncol = ncol(bx0)))
        wz0 <- as.data.table(matrix(0, nrow = nrow(bz0), ncol = ncol(bz0)))
        
        colnames(bz0) <- Xbz_vars
        colnames(wz0) <- Xwz_vars
        colnames(bx0) <- Xbx_vars
        colnames(wx0) <- Xwx_vars
        
        d0 <- cbind(bz0, wz0, bx0, wx0, d0[, colnames(d0) %nin% c(Xbz_vars, Xwz_vars, Xbx_vars, Xwx_vars), with = FALSE])
      }
    } else {
      ## assemble reference grid
      ## get var names
      zs <- c(Xbz_vars, Xwz_vars, Xz_vars)
      
      resp  <- brmcoda_vars[["y"]]
      grp   <- object[["complr"]][["idvar"]]
      preds <- brmcoda_vars[["x"]] %sin% c(Xz_vars, Xbz_vars, Xwz_vars)
      covs  <- brmcoda_vars[["x"]] %snin% c(resp, grp, preds)
      
      ## default reference grid
      refgrid <- as.data.table(ref_grid(object[["model"]], at = at)@grid)
      
      ## reference grid (only covariates and outcome)
      refgrid <- refgrid[, colnames(refgrid) %nin% c(zs), with = FALSE]
      
      ## to make fitted() happy
      id <- data.table::data.table(1) # to make fitted() happy
      colnames(id) <- object[["complr"]][["idvar"]]
      
      # grandmean
      if ("grandmean" %in% ref) {
        
        # aggregate
        if ("aggregate" %in% level) {
          if (weight == "proportional") {
            x0 <- mean.acomp(object[["complr"]][["output"]][[idx]][["X"]], robust = TRUE)
            
          } else {
            x0 <- object[["complr"]][["dataout"]][, head(.SD, 1), by = eval(object[["complr"]][["idvar"]])]
            x0 <- acomp(x0[, Xtx_vars, with = FALSE], total = object[["complr"]][["output"]][[idx]][["total"]])
            x0 <- mean.acomp(x0, robust = TRUE)
          }
          
          x0 <- acomp(x0, total = object[["complr"]][["output"]][[idx]][["total"]])
          x0 <- as.data.table(t(x0))
          
          z0 <- ilr(x0, V = object[["complr"]][["output"]][[idx]][["psi"]])
          z0 <- as.data.table(t(z0))
          
          colnames(z0) <- Xz_vars
          colnames(x0) <- Xtx_vars
          
          d0 <- if (all(dim(refgrid) == 0)) (cbind(z0, x0, id)) else (expand.grid.df(z0, x0, id, refgrid))
        }
        
        # between and/or within
        if (any(c("between", "within") %in% level)) {
          if (weight == "proportional") {
            bx0 <- mean.acomp(object[["complr"]][["output"]][[idx]][["bX"]], robust = TRUE)
            
          } else {
            bx0 <- object[["complr"]][["dataout"]][, head(.SD, 1), by = eval(object[["complr"]][["idvar"]])]
            bx0 <- acomp(bx0[, Xbx_vars, with = FALSE], total = object[["complr"]][["output"]][[idx]][["total"]])
            bx0 <- mean.acomp(bx0, robust = TRUE)
          }
          
          bx0 <- acomp(bx0, total = object[["complr"]][["output"]][[idx]][["total"]])
          bx0 <- as.data.table(t(bx0))
          
          bz0 <- ilr(bx0, V = object[["complr"]][["output"]][[idx]][["psi"]])
          bz0 <- as.data.table(t(bz0))
          
          wx0 <- as.data.table(matrix(1, nrow = nrow(bx0), ncol = ncol(bx0)))
          wz0 <- as.data.table(matrix(0, nrow = nrow(bz0), ncol = ncol(bz0)))
          
          colnames(bz0) <- Xbz_vars
          colnames(wz0) <- Xwz_vars
          colnames(bx0) <- Xbx_vars
          colnames(wx0) <- Xwx_vars
          
          d0 <- if (all(dim(refgrid) == 0)) (cbind(bz0, wz0, bx0, wx0, id)) else (expand.grid.df(bz0, wz0, bx0, wx0, id, refgrid))
        }
      }
      
      # user specified
      if (inherits(ref, c("data.table", "data.frame", "matrix"))) {
        weight <- NULL
        
        if (isFALSE(object[["complr"]][["output"]][[idx]][["parts"]] %in% colnames(ref))) {  # get user's composition
          stop(
            sprintf(
              "The reference grid should include all compositional components but (%s) are missing.",
              paste0(object[["complr"]][["output"]][[idx]][["parts"]] %nin% colnames(ref), collapse = ", ")
            ))
        } else {
          xU <- ref[, object[["complr"]][["output"]][[idx]][["parts"]], with = FALSE]
          xU <- acomp(xU, total = object[["complr"]][["output"]][[idx]][["total"]])
          xU <- as.data.table(t(xU))
        }
        
        # sanity checks
        if (nrow(ref) > 1) {
          stop("Only one reference composition is allowed at a time.")
        }
        if(isFALSE(sum(xU) == object[["complr"]][["output"]][[idx]][["total"]])) {
          stop(sprintf(
            "The total amount of the reference composition (%s) should be the same as the composition (%s).",
            sum(xU),
            object[["complr"]][["output"]][[idx]][["total"]]
          ))
        }
        if (isTRUE((any(xU > lapply(object[["complr"]][["dataout"]][, object[["complr"]][["output"]][[idx]][["parts"]], with = FALSE], max)) |
                    any(xU < lapply(object[["complr"]][["dataout"]][, object[["complr"]][["output"]][[idx]][["parts"]], with = FALSE], min))))) {
          stop(paste(
            sprintf(
              "composition should be numeric or interger values that are between (%s) and (%s)",
              paste0(round(apply(object[["complr"]][["dataout"]][, object[["complr"]][["output"]][[idx]][["parts"]], with = FALSE], 2, min)), collapse = ", "),
              paste0(round(apply(object[["complr"]][["dataout"]][, object[["complr"]][["output"]][[idx]][["parts"]], with = FALSE], 2, max)), collapse = ", ")),
            "\n",
            " for",
            paste0(object[["complr"]][["output"]][[idx]][["parts"]], collapse = ", "),
            "respectively"
          ))
        }
        
        # user's specified reference grid - edit to allow for new var names
        ## any covariates left in the ref
        if (ncol(ref) > ncol(xU)) {
          covgrid <- ref[, -object[["complr"]][["output"]][[idx]][["parts"]], with = FALSE]
          
          if (isFALSE(fill)) {
            if (isFALSE(identical(colnames(covgrid), covs))) {
              # ensure all covs are provided
              stop(paste(
                "'ref' should contain information about",
                "  the covariates in 'brmcoda' model to estimate substitution",
                "  except the logratio variables nor any column names starting with 'z', 'bz', or 'wz',",
                "  as these variables will be computed in substitution analysis.",
                "  Please provide a different reference grid.",
                sep = "\n"))
            }
          } else {
            # grab any covariates in user's specified reference grid
            # and fill refgrid if any is missing
            refgrid <- as.data.table(expand.grid.df(covgrid, refgrid[, -colnames(covgrid), with = FALSE]))
          }
        } else {
          refgrid <- refgrid[, covs, with = FALSE]
        }
        
        if (level == "aggregate") {
          x0 <- xU
          
          z0 <- ilr(x0, V = object[["complr"]][["output"]][[idx]][["psi"]])
          z0 <- as.data.table(t(z0))
          
          colnames(z0) <- Xz_vars
          colnames(x0) <- Xtx_vars
          
          d0 <- if (all(dim(refgrid) == 0)) (cbind(z0, x0, id)) else (expand.grid.df(z0, x0, id, refgrid))
          
        }
        if (level %in% c("between", "within")) {
          x0 <- object[["complr"]][["dataout"]][, head(.SD, 1), by = eval(object[["complr"]][["idvar"]])]
          x0 <- acomp(x0[, Xbx_vars, with = FALSE], total = object[["complr"]][["output"]][[idx]][["total"]])
          x0 <- mean.acomp(x0, robust = TRUE)
          
          # assemble d0
          # bz0 is between-person ilr of the ref comp (doesn't have to be compositional mean)
          bx0 <- xU
          bz0 <- ilr(bx0, V = object[["complr"]][["output"]][[idx]][["psi"]])
          bz0 <- as.data.table(t(bz0))
          
          # wx0 and wz0 are the difference between the actual compositional mean of the dataset and bilr
          # is 0 if ref comp is compositional mean
          # but is different if not
          wx0 <- bx0 - x0
          wz0 <- as.data.table(t(ilr(wx0, V = object[["complr"]][["output"]][[idx]][["psi"]])))
          
          id <- data.table::data.table(1) # to make fitted() happy
          
          colnames(bz0) <- Xbz_vars
          colnames(wz0) <- Xwz_vars
          colnames(bx0) <- Xbx_vars
          colnames(wx0) <- Xwx_vars
          colnames(id)  <- object[["complr"]][["idvar"]]
          
          d0 <- if (all(dim(refgrid) == 0)) (cbind(bz0, wz0, bx0, wx0, id)) else (expand.grid.df(bz0, wz0, bx0, wx0, id, refgrid))
        }
      }
    }
  }
  
  ## d0 and x0 for single level model
  if (identical(brmcoda_vars$ranef_type, "single")) {
    
    x0 <- object[["complr"]][["output"]][[idx]][["X"]]
    x0 <- mean.acomp(x0, robust = TRUE)
    x0 <- acomp(x0, total = object[["complr"]][["output"]][[idx]][["total"]])
    x0 <- as.data.table(t(x0))
    
    z0 <- ilr(x0, V = object[["complr"]][["output"]][[idx]][["psi"]])
    z0 <- as.data.table(t(z0))
    
    colnames(z0) <- Xz_vars
    colnames(x0) <- Xtx_vars
    
    # assemble reference grid
    # get var names
    zs <- c(Xz_vars, Xbz_vars, Xwz_vars)
    
    resp  <- brmcoda_vars[["y"]]
    preds <- Xz_vars
    covs  <- brmcoda_vars[["x"]] %snin% c(resp, zs)
    
    refgrid <- as.data.table(ref_grid(object[["model"]], at = at)@grid)
    
    # reference grid (only covariates and outcome)
    refgrid <- refgrid[, colnames(refgrid) %nin% c(zs), with = FALSE]
    
    d0 <- if (all(dim(refgrid) == 0)) (cbind(z0, x0)) else (expand.grid.df(z0, x0, refgrid))
  }
  
  ## drop rep.measure if not null
  if ("rep.meas" %in% colnames(d0)) {
    d0$rep.meas <- NULL
    d0 <- unique(d0)
  }
  
  as.data.table(d0)
}