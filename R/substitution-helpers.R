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
#' @param aorg A logical value specifying whether to summarize the results (average over the reference grid)
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
           aorg) {
    
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
      brmsformula = object$model$formula,
      delta = delta,
      ref = ref,
      level = level,
      weight = weight,
      parts = parts,
      aorg = aorg
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
  
  ## get the index of which index elements of object$complr$output does the parts correspond to
  parts_index <- which(vapply(lapply(object$complr$output, function(x) x$parts), function(p) identical(parts, p), logical(1)))
  
  ## NOTES
  # # ignore weight for clustermean
  # equal weight is default for grandmean
  
  # what type of model is being estimated
  model_fixef <- rownames(fixef(object))
  model_ranef <- if(dim(object$model$ranef)[1] > 0) (names(ranef(object))) else (NULL)
  
  model_fixef_level <- NULL
  model_fixef_coef  <- NULL
  
  # grab the correct logratio names
  z_vars  <- names(object$complr$output[[parts_index]]$logratio)
  bz_vars <- names(object$complr$output[[parts_index]]$between_logratio)
  wz_vars <- names(object$complr$output[[parts_index]]$within_logratio)
  
  if (length(grep(paste0(bz_vars, collapse = "|"), model_fixef, value = T)) > 0) {
    model_fixef_level <- append(model_fixef_level, "between")
    model_fixef_coef  <- append(model_fixef_coef,
                                grep(paste0(bz_vars, collapse = "|"), model_fixef, value = T))
  }
  if (length(grep(paste0(wz_vars, collapse = "|"), model_fixef, value = T)) > 0) {
    model_fixef_level <- append(model_fixef_level, "within")
    model_fixef_coef  <- append(model_fixef_coef,
                                grep(paste0(wz_vars, collapse = "|"), model_fixef, value = T))
  }
  if ((length(grep(paste0(z_vars, collapse = "|"), model_fixef, value = T)) > 0) && (length(grep(paste0(c(bz_vars, wz_vars), collapse = "|"), model_fixef, value = T)) == 0)) {
    model_fixef_level <- append(model_fixef_level, "aggregate")
    model_fixef_coef  <- append(model_fixef_coef, setdiff(
      grep(paste0(z_vars, collapse = "|"), model_fixef, value = TRUE),
      grep(paste0(c(bz_vars, wz_vars), collapse = "|"), model_fixef, value = TRUE)
    ))
  }
  
  # single level or multilevel
  if (length(model_ranef) > 0) {
    model_ranef_level <- "multilevel"
    model_ranef_coef  <- model_ranef
  } else {
    model_ranef_level <- "single"
    model_ranef_coef  <- NULL
  }
  
  # d0 and x0 for multilevel level model
  if (model_ranef_level == "multilevel") {
    
    # for clustermean
    if ("clustermean" %in% ref) {
      weight <- NULL # ignore weight
      
      # aggregate
      if ("aggregate" %in% level) {
        d0 <- cbind(object$complr$output[[parts_index]]$comp, object$complr$output[[parts_index]]$dataout[, -colnames(object$complr$output[[parts_index]]$comp), with = FALSE])
        d0 <- d0[, head(.SD, 1), by = eval(object$complr$idvar)]
        x0 <- acomp(d0[, colnames(object$complr$output[[parts_index]]$comp), with = FALSE], total = object$complr$output[[parts_index]]$total)
        
        z0 <- ilr(x0, V = object$complr$output[[parts_index]]$psi)
        z0 <- as.data.table(z0)
        
        colnames(z0)  <- colnames(object$complr$output[[parts_index]]$logratio)
        colnames(x0) <- colnames(object$complr$output[[parts_index]]$comp)
        
        d0 <- cbind(z0, x0,
                    d0[, -colnames(x0), with = FALSE])
      }
      
      # between and within
      if (any(c("between", "within") %in% level)) {
        d0 <- cbind(object$complr$output[[parts_index]]$between_comp, 
                    object$complr$output[[parts_index]]$dataout)
        d0 <- d0[, head(.SD, 1), by = eval(object$complr$idvar)]
        x0  <- acomp(d0[, colnames(object$complr$output[[parts_index]]$between_comp), with = FALSE], total = object$complr$output[[parts_index]]$total)
        bx0 <- x0
        
        bz0 <- ilr(bx0, V = object$complr$output[[parts_index]]$psi)
        bz0 <- as.data.table(bz0)
        
        wx0 <- as.data.table(matrix(1, nrow = nrow(bx0), ncol = ncol(bx0)))
        wz0    <- as.data.table(matrix(0, nrow = nrow(bz0), ncol = ncol(bz0)))
        
        colnames(bz0)    <- colnames(object$complr$output[[parts_index]]$between_logratio)
        colnames(wz0)    <- colnames(object$complr$output[[parts_index]]$within_logratio)
        colnames(bx0) <- colnames(object$complr$output[[parts_index]]$between_comp)
        colnames(wx0) <- colnames(object$complr$output[[parts_index]]$within_comp)
        
        d0 <- cbind(bz0, wz0, bx0, wx0,
                    d0[, colnames(d0) %in% colnames(object$complr$output[[parts_index]]$dataout), with = FALSE])
      }
    } else {
      
      # assemble reference grid
      # get var names
      ilrnames <- c(colnames(object$complr$output[[parts_index]]$between_logratio),
                    colnames(object$complr$output[[parts_index]]$within_logratio),
                    colnames(object$complr$output[[parts_index]]$logratio))
      
      vars  <- colnames(model.frame(object))
      resp  <- object$model$formula$formula[[2]]
      grp   <- object$model$ranef$group
      preds <- vars %snin% c(resp, grp)
      covs  <- vars %snin% c(resp, grp, ilrnames)
      
      # default reference grid
      drg <- as.data.table(ref_grid(object$model, at = at)@grid)
      
      # reference grid (only covariates and outcome)
      refgrid <- drg[, colnames(drg) %in% c(covs, resp), with = FALSE]
      
      id <- data.table::data.table(1) # to make fitted() happy
      colnames(id) <- object$complr$idvar
      
      # grandmean
      if ("grandmean" %in% ref) {
        
        # aggregate
        if ("aggregate" %in% level) {
          if (weight == "proportional") {
            x0 <- mean.acomp(object$complr$output[[parts_index]]$comp, robust = TRUE)
            
          } else {
            x0 <- cbind(object$complr$output[[parts_index]]$comp,
                        object$complr$output[[parts_index]]$dataout[, object$complr$idvar, with = FALSE])
            x0 <- x0[, head(.SD, 1), by = eval(object$complr$idvar)]
            x0 <- acomp(x0[, colnames(object$complr$output[[parts_index]]$comp), with = FALSE], total = object$complr$output[[parts_index]]$total)
            x0 <- mean.acomp(x0, robust = TRUE)
          }
          
          x0 <- acomp(x0, total = object$complr$output[[parts_index]]$total)
          x0 <- as.data.table(t(x0))
          
          z0 <- ilr(x0, V = object$complr$output[[parts_index]]$psi)
          z0 <- as.data.table(t(z0))
          
          colnames(z0)  <- colnames(object$complr$output[[parts_index]]$logratio)
          colnames(x0) <- colnames(object$complr$output[[parts_index]]$comp)
          
          d0 <- if (all(dim(refgrid) == 0)) (cbind(z0, x0, id)) else (expand.grid.df(z0, x0, id, refgrid))
        }
        
        # between and/or within
        if (any(c("between", "within") %in% level)) {
          if (weight == "proportional") {
            x0 <- mean.acomp(object$complr$output[[parts_index]]$between_comp, robust = TRUE)
            
          } else {
            x0 <- cbind(object$complr$output[[parts_index]]$between_comp,
                        object$complr$output[[parts_index]]$dataout[, object$complr$idvar, with = FALSE])
            x0 <- x0[, head(.SD, 1), by = eval(object$complr$idvar)]
            x0 <- acomp(x0[, colnames(object$complr$output[[parts_index]]$between_comp), with = FALSE], total = object$complr$output[[parts_index]]$total)
            x0 <- mean.acomp(x0, robust = TRUE)
          }
          
          x0  <- acomp(x0, total = object$complr$output[[parts_index]]$total)
          x0  <- as.data.table(t(x0))
          bx0 <- x0
          
          bz0 <- ilr(bx0, V = object$complr$output[[parts_index]]$psi)
          bz0 <- as.data.table(t(bz0))
          
          wx0 <- as.data.table(matrix(1, nrow = nrow(bx0), ncol = ncol(bx0)))
          wz0 <- as.data.table(matrix(0, nrow = nrow(bz0), ncol = ncol(bz0)))
          
          colnames(bz0) <- colnames(object$complr$output[[parts_index]]$between_logratio)
          colnames(wz0) <- colnames(object$complr$output[[parts_index]]$within_logratio)
          colnames(bx0) <- colnames(object$complr$output[[parts_index]]$between_comp)
          colnames(wx0) <- colnames(object$complr$output[[parts_index]]$within_comp)
          
          d0 <- if (all(dim(refgrid) == 0)) (cbind(bz0, wz0, bx0, wx0, id)) else (expand.grid.df(bz0, wz0, bx0, wx0, id, refgrid))
        }
      }
      
      if (inherits(ref, c("data.table", "data.frame", "matrix"))) {
        weight <- NULL
        
        if (isFALSE(object$complr$output[[parts_index]]$parts %in% colnames(ref))) {  # get user's composition
          stop(
            sprintf(
              "The reference grid should include all compositional components but (%s) are missing.",
              paste0(object$complr$output[[parts_index]]$parts %nin% colnames(ref), collapse = ", ")
            ))
        } else {
          comp_user <- ref[, object$complr$output[[parts_index]]$parts, with = FALSE]
          comp_user <- acomp(comp_user, total = object$complr$output[[parts_index]]$total)
          comp_user <- as.data.table(t(comp_user))
        }
        
        # sanity checks
        if (nrow(ref) > 1) {
          stop("Only one reference composition is allowed at a time.")
        }
        if(isFALSE(sum(comp_user) == object$complr$output[[parts_index]]$total)) {
          stop(sprintf(
            "The total amount of the reference composition (%s) should be the same as the composition (%s).",
            sum(comp_user),
            object$complr$output[[parts_index]]$total
          ))
        }
        if (isTRUE((any(comp_user > lapply(object$complr$output[[parts_index]]$dataout[, object$complr$output[[parts_index]]$parts, with = FALSE], max)) |
                    any(comp_user < lapply(object$complr$output[[parts_index]]$dataout[, object$complr$output[[parts_index]]$parts, with = FALSE], min))))) {
          stop(paste(
            sprintf(
              "composition should be numeric or interger values that are between (%s) and (%s)",
              paste0(round(apply(object$complr$output[[parts_index]]$dataout[, object$complr$output[[parts_index]]$parts, with = FALSE], 2, min)), collapse = ", "),
              paste0(round(apply(object$complr$output[[parts_index]]$dataout[, object$complr$output[[parts_index]]$parts, with = FALSE], 2, max)), collapse = ", ")),
            "\n",
            " for",
            paste0(object$complr$output[[parts_index]]$parts, collapse = ", "),
            "respectively"
          ))
        }
        
        # user's specified reference grid
        ## any covariates left in the ref
        if (ncol(ref) > ncol(comp_user)) {
          covgrid <- ref[, -object$complr$output[[parts_index]]$parts, with = FALSE]
          
          if (isFALSE(fill)) {
            if (isFALSE(identical(colnames(covgrid), covs))) {
              # ensure all covs are provided
              stop(paste(
                "'ref' should contain information about",
                "  the covariates in 'brmcoda' model to estimate substitution",
                "  except the ILR variables nor any column names starting with 'bilr', 'wilr', or 'ilr',",
                "  as these variables will be computed in substitution analysis.",
                "  Please provide a different reference grid.",
                sep = "\n"))
            }
          } else {
            # grab any covariates in user's specified reference grid
            # and fill drg if any is missing
            refgrid <- as.data.table(expand.grid.df(covgrid,
                                                    refgrid[, -colnames(covgrid), with = FALSE]))
          }
        } else {
          refgrid <- drg[, covs, with = FALSE]
        }
        
        if (level == "aggregate") {
          x0 <- comp_user
          
          z0 <- ilr(x0, V = object$complr$output[[parts_index]]$psi)
          z0 <- as.data.table(t(z0))
          
          colnames(z0)  <- colnames(object$complr$output[[parts_index]]$logratio)
          colnames(x0) <- colnames(object$complr$output[[parts_index]]$comp)
          
          d0 <- if (all(dim(refgrid) == 0)) (cbind(z0, x0, id)) else (expand.grid.df(z0, x0, id, refgrid))
          
        }
        if (level %in% c("between", "within")) {
          x0 <- cbind(object$complr$output[[parts_index]]$between_comp,
                      object$complr$output[[parts_index]]$dataout[, object$complr$idvar, with = FALSE])
          x0 <- x0[, head(.SD, 1), by = eval(object$complr$idvar)]
          x0 <- acomp(x0[, colnames(object$complr$output[[parts_index]]$between_comp), with = FALSE], total = object$complr$output[[parts_index]]$total)
          x0 <- mean.acomp(x0, robust = TRUE)
          
          # assemble d0
          # bz0 is between-person ilr of the ref comp (doesn't have to be compositional mean)
          bx0 <- comp_user
          bz0 <- ilr(bx0, V = object$complr$output[[parts_index]]$psi)
          bz0 <- as.data.table(t(bz0))
          
          # wx0 and wz0 are the difference between the actual compositional mean of the dataset and bilr
          # is 0 if ref comp is compositional mean
          # but is different if not
          wx0 <- bx0 - x0
          wz0 <- as.data.table(t(ilr(wx0, V = object$complr$output[[parts_index]]$psi)))
          
          id <- data.table::data.table(1) # to make fitted() happy
          
          colnames(bz0)  <- colnames(object$complr$output[[parts_index]]$between_logratio)
          colnames(wz0)  <- colnames(object$complr$output[[parts_index]]$within_logratio)
          colnames(bx0) <- colnames(object$complr$output[[parts_index]]$between_comp)
          colnames(wx0) <- colnames(object$complr$output[[parts_index]]$within_comp)
          colnames(id)     <- object$complr$idvar
          
          d0 <- if (all(dim(refgrid) == 0)) (cbind(bz0, wz0, bx0, wx0, id)) else (expand.grid.df(bz0, wz0, bx0, wx0, id, refgrid))
        }
      }
    }
  }
  
  ## d0 and x0 for single level model
  if (model_ranef_level == "single") {
    
    x0 <- object$complr$output[[parts_index]]$comp
    x0 <- mean.acomp(x0, robust = TRUE)
    x0 <- acomp(x0, total = object$complr$output[[parts_index]]$total)
    x0 <- as.data.table(t(x0))
    bx0 <- x0
    
    d0 <- object$complr$output[[parts_index]]$dataout
    
    z0 <- ilr(x0, V = object$complr$output[[parts_index]]$psi)
    z0 <- as.data.table(t(z0))
    
    colnames(z0)  <- colnames(object$complr$output[[parts_index]]$logratio)
    colnames(x0) <- colnames(object$complr$output[[parts_index]]$comp)
    
    # assemble reference grid
    # get var names
    ilrnames <- c(colnames(object$complr$output[[parts_index]]$between_logratio),
                  colnames(object$complr$output[[parts_index]]$within_logratio),
                  colnames(object$complr$output[[parts_index]]$logratio))
    
    vars  <- colnames(model.frame(object))
    resp  <- object$model$formula$formula[[2]]
    # grp   <- object$model$ranef$group
    preds <- vars %snin% c(resp)
    covs  <- vars %snin% c(resp, ilrnames)
    
    # default reference grid
    # set binary
    # drg <- model.frame(object)
    # drg[] <- as.data.table(lapply(drg, function(j) if(is.numeric(j) && unique(j) %ain% c(0, 1)) as.factor(j) else j))
    # drg <- as.data.table(insight::get_datagrid(drg,
    #                                            by = paste0(resp),
    #                                            factors = factors,
    #                                            length = NA))
    
    drg <- as.data.table(ref_grid(object$model, at = at)@grid)
    
    # reference grid (only covariates and outcome)
    refgrid <- drg[, colnames(drg) %in% c(covs, resp), with = FALSE]
    
    d0 <- if (all(dim(refgrid) == 0)) (cbind(z0, x0)) else (expand.grid.df(z0, x0, refgrid))
  }
  as.data.table(d0)
}

#' Helper functions used only internally to estimate substitution model
#' @importFrom data.table as.data.table data.table copy := setDT rbindlist .SD
#' @importFrom compositions acomp ilr clo mean.acomp
#' @importFrom bayestestR describe_posterior
#' @importFrom extraoperators %snin% %sin%
#' @importFrom foreach foreach %dopar%
#' @importFrom doFuture %dofuture%
#' @importFrom future plan multisession sequential
#'
#' @name get-substitution
NULL

# Grandmean Between-person Substitution model
.get.bsub <- function(object,
                      delta,
                      base,
                      x0,
                      y0,
                      d0,
                      aorg,
                      summary,
                      level,
                      ref,
                      scale,
                      type,
                      cores,
                      
                      ...) {
  
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }
  oopts <- options(future.globals.maxSize = +Inf, future.globals.onReference = NULL)
  on.exit(options(oopts))
  
  iout <- foreach(i = colnames(base), .combine = c,
                  .options.future = list(packages = "multilevelcoda", seed = TRUE)) %dofuture% {
                    
                    # substitution variables
                    if (type == "one-to-all") {
                      
                      # one to remaining
                      base_tmp <- as.data.table(base)
                      base_tmp <- base_tmp[(get(i) %in% c(1, -1))]
                      
                      sub_from_var <- c("remaining", i)
                      sub_to_var  <- c(i, "remaining")
                    }
                    else {
                      # possible pairwise substitution of 1 compositional variable
                      # one to one
                      base_tmp <- as.data.table(base)
                      base_tmp <- base_tmp[(get(i) != 0)]
                      base_tmp <- base_tmp[order(-rank(get(i)))]
                      
                      sub_from_var <- colnames(base_tmp) %snin% eval(i)
                      sub_to_var   <- i
                    }
                    
                    # loop substitution
                    kout <- vector("list", length = nrow(base_tmp))
                    jout <- vector("list", length = length(delta))
                    
                    for (j in seq_along(delta)) { # delta level
                      sub_tmp_j <- base_tmp * delta[j]
                      for (k in seq_len(nrow(sub_tmp_j))) {
                        x1 <- x0 + sub_tmp_j[k,]
                        names(x1) <- object$complr$output[[parts_index]]$parts
                        Delta <- sub_tmp_j[k, get(i)]
                        kout[[k]] <- cbind(x0, x1, Delta)
                      }
                      jout[[j]] <- do.call(rbind, kout)
                    }
                    d1 <- setDT(do.call(rbind, jout))
                    
                    # useful information for the final results
                    d1[, From := rep(sub_from_var, length.out = nrow(d1))]
                    d1[, To := rep(sub_to_var, length.out = nrow(d1))]
                    d1[, Delta := as.numeric(Delta)]
                    d1[, Level := level]
                    d1[, Reference := ref]
                    
                    # remove impossible reallocation that result in negative values
                    cols <- colnames(d1) %snin% c("Delta", "From", "To", "Level")
                    d1 <- d1[rowSums(d1[, ..cols] < 0) == 0]
                    
                    # compositions and ilrs for predictions
                    bx0 <- acomp(d1[, colnames(object$complr$output[[parts_index]]$between_comp), with = FALSE], total = object$complr$output[[parts_index]]$total)
                    bcompsub  <- acomp(d1[, object$complr$output[[parts_index]]$parts, with = FALSE], total = object$complr$output[[parts_index]]$total)
                    
                    bilrsub <- ilr(bcompsub, V = object$complr$output[[parts_index]]$psi)
                    wz0 <- as.data.table(matrix(0, nrow = nrow(bilrsub), ncol = ncol(bilrsub)))
                    
                    colnames(bilrsub) <- colnames(object$complr$output[[parts_index]]$between_logratio)
                    colnames(wz0) <- colnames(object$complr$output[[parts_index]]$within_logratio)
                    
                    # reference grid
                    ## get covariate + idvar names
                    covs <- colnames(d0) %snin% c(colnames(object$complr$output[[parts_index]]$between_logratio),
                                                  colnames(object$complr$output[[parts_index]]$within_logratio),
                                                  colnames(object$complr$output[[parts_index]]$between_comp),
                                                  colnames(object$complr$output[[parts_index]]$within_comp)
                    )
                    refgrid <- d0[, covs, with = FALSE]
                    
                    # predictions
                    hout <- vector("list", length = nrow(d0))
                    if (aorg) { # unadj OR adj averaging over reference grid
                      for (h in seq_len(nrow(d0))) {
                        dsub <- cbind(d1, bilrsub, wz0, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            scale = scale,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        hout[[h]] <- delta_y
                      }
                      delta_y_avg <- Reduce(`+`, hout) / length(hout)
                      
                      # posterior summary or draws?
                      if(summary) {
                        suppressWarnings(posterior_delta_y <- apply(delta_y_avg, 2, function(x) describe_posterior(x, centrality = "mean", ...)))
                        posterior_delta_y <- rbindlist(posterior_delta_y)
                        posterior_delta_y <- cbind(posterior_delta_y[, .(Mean, CI_low, CI_high)],
                                                   dsub[, .(Delta, From, To, Level, Reference)])
                      } else {
                        posterior_delta_y <- t(delta_y_avg)
                        posterior_delta_y <- cbind(dsub[, .(Delta, From, To, Level, Reference)], posterior_delta_y)
                      }
                      
                    } else { # adj keeping prediction at each level of reference grid
                      for (h in seq_len(nrow(d0))) {
                        dsub <- cbind(d1, bilrsub, wz0, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            scale = scale,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        
                        # posterior summary or draws?
                        if(summary) {
                          suppressWarnings(posterior_delta_y <- apply(delta_y, 2, function(x) describe_posterior(x, centrality = "mean", ...)))
                          posterior_delta_y <- rbindlist(posterior_delta_y)
                          posterior_delta_y <- cbind(posterior_delta_y[, .(Mean, CI_low, CI_high)],
                                                     dsub[, .(Delta, From, To, Level, Reference)])
                        } else {
                          posterior_delta_y <- t(delta_y)
                          posterior_delta_y <- cbind(dsub[, .(Delta, From, To, Level, Reference)], posterior_delta_y)
                        }
                        
                        hout[[h]] <- posterior_delta_y
                      }
                      posterior_delta_y <- rbindlist(hout)
                    }
                    
                    # final results for entire composition
                    out <- list(posterior_delta_y)
                    names(out) <- i
                    out
                  }
  iout
}

# Grandmean Within-person Substitution model
.get.wsub <- function(object,
                      delta,
                      base,
                      x0,
                      y0,
                      d0,
                      aorg,
                      summary,
                      level,
                      ref,
                      scale,
                      type,
                      cores,
                      ...) {
  
  
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }
  oopts <- options(future.globals.maxSize = +Inf, future.globals.onReference = NULL)
  on.exit(options(oopts))
  
  iout <- foreach(i = colnames(base), .combine = c,
                  .options.future = list(packages = "multilevelcoda", seed = TRUE)) %dofuture% {
                    
                    # substitution variables
                    if (type == "one-to-all") {
                      
                      # one to remaining
                      base_tmp <- as.data.table(base)
                      base_tmp <- base_tmp[(get(i) %in% c(1, -1))]
                      
                      sub_from_var <- c("remaining", i)
                      sub_to_var <- c(i, "remaining")
                    }
                    else {
                      # possible pairwise substitution of 1 compositional variable
                      # one to one
                      base_tmp <- as.data.table(base)
                      base_tmp <- base_tmp[(get(i) != 0)]
                      base_tmp <- base_tmp[order(-rank(get(i)))]
                      
                      sub_from_var <- colnames(base_tmp) %snin% eval(i)
                      sub_to_var <- i
                    }
                    
                    # loop substitution
                    kout <- vector("list", length = nrow(base_tmp))
                    jout <- vector("list", length = length(delta))
                    
                    for (j in seq_along(delta)) { # delta level
                      sub_tmp_j <- base_tmp * delta[j]
                      for (k in seq_len(nrow(sub_tmp_j))) {
                        x1 <- x0 + sub_tmp_j[k, ]
                        names(x1) <- object$complr$output[[parts_index]]$parts
                        Delta <- sub_tmp_j[k, get(i)]
                        kout[[k]] <- cbind(x0, x1, Delta)
                      }
                      jout[[j]] <- do.call(rbind, kout)
                    }
                    d1 <- setDT(do.call(rbind, jout))
                    
                    # useful information for the final results
                    d1[, From := rep(sub_from_var, length.out = nrow(d1))]
                    d1[, To := rep(sub_to_var, length.out = nrow(d1))]
                    d1[, Delta := as.numeric(Delta)]
                    d1[, Level := level]
                    d1[, Reference := ref]
                    
                    # remove impossible reallocation that result in negative values
                    cols <- colnames(d1) %snin% c("Delta", "From", "To", "Level")
                    d1 <- d1[rowSums(d1[, ..cols] < 0) == 0]
                    
                    # compositions and ilrs for predictions
                    bx0   <- acomp(d1[, colnames(object$complr$output[[parts_index]]$between_comp), with = FALSE], total = object$complr$output[[parts_index]]$total)
                    bcompsub <- acomp(d1[, object$complr$output[[parts_index]]$parts, with = FALSE], total = object$complr$output[[parts_index]]$total)
                    
                    bz0   <- ilr(bx0, V = object$complr$output[[parts_index]]$psi)
                    bilrsub <- ilr(bcompsub, V = object$complr$output[[parts_index]]$psi)
                    wilrsub <- bilrsub - bz0
                    
                    colnames(bz0) <- colnames(object$complr$output[[parts_index]]$between_logratio)
                    colnames(wilrsub) <- colnames(object$complr$output[[parts_index]]$within_logratio)
                    
                    # reference grid
                    ## get covariate + idvar names
                    covs <- colnames(d0) %snin% c(colnames(object$complr$output[[parts_index]]$between_logratio),
                                                  colnames(object$complr$output[[parts_index]]$within_logratio),
                                                  colnames(object$complr$output[[parts_index]]$between_comp),
                                                  colnames(object$complr$output[[parts_index]]$within_comp)
                    )
                    refgrid <- d0[, covs, with = FALSE]
                    
                    # predictions
                    hout <- vector("list", length = nrow(d0))
                    if (aorg) { # unadj OR adj averaging over reference grid
                      for (h in seq_len(nrow(d0))) {
                        dsub <- cbind(d1, bz0, wilrsub, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            scale = scale,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        hout[[h]] <- delta_y
                      }
                      delta_y_avg <- Reduce(`+`, hout) / length(hout)
                      
                      # posterior summary or draws?
                      if(summary) {
                        suppressWarnings(posterior_delta_y <- apply(delta_y_avg, 2, function(x) describe_posterior(x, centrality = "mean", ...)))
                        posterior_delta_y <- rbindlist(posterior_delta_y)
                        posterior_delta_y <- cbind(posterior_delta_y[, .(Mean, CI_low, CI_high)],
                                                   dsub[, .(Delta, From, To, Level, Reference)])
                      } else {
                        posterior_delta_y <- t(delta_y_avg)
                        posterior_delta_y <- cbind(dsub[, .(Delta, From, To, Level, Reference)], posterior_delta_y)
                      }
                      
                    } else { # adj keeping prediction at each level of reference grid
                      for (h in seq_len(nrow(d0))) {
                        dsub <- cbind(d1, bz0, wilrsub, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            scale = scale,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        
                        # posterior summary or draws?
                        if(summary) {
                          suppressWarnings(posterior_delta_y <- apply(delta_y, 2, function(x) describe_posterior(x, centrality = "mean", ...)))
                          posterior_delta_y <- rbindlist(posterior_delta_y)
                          posterior_delta_y <- cbind(posterior_delta_y[, .(Mean, CI_low, CI_high)],
                                                     dsub[, .(Delta, From, To, Level, Reference)])
                        } else {
                          posterior_delta_y <- t(delta_y)
                          posterior_delta_y <- cbind(dsub[, .(Delta, From, To, Level, Reference)], posterior_delta_y)
                        }
                        
                        hout[[h]] <- posterior_delta_y
                      }
                      posterior_delta_y <- rbindlist(hout)
                    }
                    # final results for entire composition
                    out <- list(posterior_delta_y)
                    names(out) <- i
                    out
                  }
  iout
}

# Grandmean Simple Substitution
.get.sub <- function(object,
                     delta,
                     base,
                     x0,
                     y0,
                     d0,
                     aorg,
                     summary,
                     level,
                     ref,
                     scale,
                     type,
                     cores,
                     ...) {
  
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }
  oopts <- options(future.globals.maxSize = +Inf, future.globals.onReference = NULL)
  on.exit(options(oopts))
  
  iout <- foreach(i = colnames(base), .combine = c,
                  .options.future = list(packages = "multilevelcoda", seed = TRUE)) %dofuture% {
                    
                    # substitution variables
                    if (type == "one-to-all") {
                      
                      # one to remaining
                      base_tmp <- as.data.table(base)
                      base_tmp <- base_tmp[(get(i) %in% c(1, -1))]
                      
                      sub_from_var <- c("remaining", i)
                      sub_to_var <- c(i, "remaining")
                    }
                    else {
                      # possible pairwise substitution of 1 compositional variable
                      # one to one
                      base_tmp <- as.data.table(base)
                      base_tmp <- base_tmp[(get(i) != 0)]
                      base_tmp <- base_tmp[order(-rank(get(i)))]
                      
                      sub_from_var <- colnames(base_tmp) %snin% eval(i)
                      sub_to_var <- i
                    }
                    
                    # loop substitution
                    kout <- vector("list", length = nrow(base_tmp))
                    jout <- vector("list", length = length(delta))
                    
                    for (j in seq_along(delta)) { # delta level
                      sub_tmp_j <- base_tmp * delta[j]
                      for (k in seq_len(nrow(sub_tmp_j))) {
                        x1 <- x0 + sub_tmp_j[k,]
                        names(x1) <- object$complr$output[[parts_index]]$parts
                        Delta <- sub_tmp_j[k, get(i)]
                        kout[[k]] <- cbind(x1, Delta)
                      }
                      jout[[j]] <- do.call(rbind, kout)
                    }
                    d1 <- setDT(do.call(rbind, jout))
                    
                    # useful information for the final results
                    d1[, From := rep(sub_from_var, length.out = nrow(d1))]
                    d1[, To := sub_to_var]
                    d1[, Delta := as.numeric(Delta)]
                    d1[, Level := level]
                    d1[, Reference := ref]
                    
                    # remove impossible reallocation that result in negative values
                    cols <- colnames(d1) %snin% c("Delta", "From", "To", "Level")
                    d1 <- d1[rowSums(d1[, ..cols] < 0) == 0]
                    
                    # compositions and ilrs for predictions
                    compsub  <- acomp(d1[, object$complr$output[[parts_index]]$parts, with = FALSE], total = object$complr$output[[parts_index]]$total)
                    
                    ilrsub <- ilr(compsub, V = object$complr$output[[parts_index]]$psi)
                    colnames(ilrsub) <- colnames(object$complr$output[[parts_index]]$logratio)
                    
                    # reference grid
                    ## get covariate + idvar names
                    covs <- colnames(d0) %snin% c(colnames(object$complr$output[[parts_index]]$logratio),
                                                  colnames(object$complr$output[[parts_index]]$comp)
                    )
                    refgrid <- d0[, covs, with = FALSE]
                    
                    # predictions
                    hout <- vector("list", length = nrow(d0))
                    if (aorg) { # unadj OR adj averaging over reference grid
                      for (h in seq_len(nrow(d0))) {
                        dsub <- cbind(d1, ilrsub, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            scale = scale,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        hout[[h]] <- delta_y
                      }
                      delta_y_avg <- Reduce(`+`, hout) / length(hout)
                      
                      # posterior summary or draws?
                      if(summary) {
                        suppressWarnings(posterior_delta_y <- apply(delta_y_avg, 2, function(x) describe_posterior(x, centrality = "mean", ...)))
                        posterior_delta_y <- rbindlist(posterior_delta_y)
                        posterior_delta_y <- cbind(posterior_delta_y[, .(Mean, CI_low, CI_high)],
                                                   dsub[, .(Delta, From, To, Level, Reference)])
                      } else {
                        posterior_delta_y <- t(delta_y_avg)
                        posterior_delta_y <- cbind(dsub[, .(Delta, From, To, Level, Reference)], posterior_delta_y)
                      }
                      
                    } else { # adj keeping prediction at each level of reference grid
                      for (h in seq_len(nrow(d0))) {
                        dsub <- cbind(d1, ilrsub, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            scale = scale,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        
                        # posterior summary or draws?
                        if(summary) {
                          suppressWarnings(posterior_delta_y <- apply(delta_y, 2, function(x) describe_posterior(x, centrality = "mean", ...)))
                          posterior_delta_y <- rbindlist(posterior_delta_y)
                          posterior_delta_y <- cbind(posterior_delta_y[, .(Mean, CI_low, CI_high)],
                                                     dsub[, .(Delta, From, To, Level, Reference)])
                        } else {
                          posterior_delta_y <- t(delta_y)
                          posterior_delta_y <- cbind(dsub[, .(Delta, From, To, Level, Reference)], posterior_delta_y)
                        }
                        
                        hout[[h]] <- posterior_delta_y
                      }
                      posterior_delta_y <- rbindlist(hout)
                    }
                    # final results for entire composition
                    out <- list(posterior_delta_y)
                    names(out) <- i
                    out
                  }
  iout
}

# Clustermean Between-person Substitution model
.get.bsubmargins <- function(object,
                             delta,
                             base,
                             x0,
                             y0,
                             d0,
                             summary,
                             level,
                             ref,
                             scale,
                             type,
                             cores,
                             ...) {
  
  
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }
  oopts <- options(future.globals.maxSize = +Inf, future.globals.onReference = NULL)
  on.exit(options(oopts))
  
  iout <- foreach(i = colnames(base), .combine = c,
                  .options.future = list(packages = "multilevelcoda", seed = TRUE)) %dofuture% {
                    
                    # substitution variables
                    if (type == "one-to-all") {
                      
                      # one to remaining
                      base_tmp <- as.data.table(base)
                      base_tmp <- base_tmp[(get(i) %in% c(1, -1))]
                      
                      sub_from_var <- c("remaining", i)
                      sub_to_var <- c(i, "remaining")
                    }
                    else {
                      # possible pairwise substitution of 1 compositional variable
                      # one to one
                      base_tmp <- as.data.table(base)
                      base_tmp <- base_tmp[(get(i) != 0)]
                      base_tmp <- base_tmp[order(-rank(get(i)))]
                      
                      sub_from_var <- colnames(base_tmp) %snin% eval(i)
                      sub_to_var <- i
                    }
                    
                    # loop substitution
                    kout <- vector("list", length = nrow(base_tmp))
                    jout <- vector("list", length = length(delta))
                    
                    for (j in seq_along(delta)) { # delta level
                      sub_tmp_j <- base_tmp * delta[j]
                      for (k in seq_len(nrow(sub_tmp_j))) { # reallocation level
                        sub_tmp_k <- sub_tmp_j[k, ]
                        sub_tmp_k <- sub_tmp_k[rep(seq_len(nrow(sub_tmp_k)), nrow(x0)), ]
                        x1 <- x0 + sub_tmp_k
                        names(x1) <- object$complr$output[[parts_index]]$parts
                        
                        Delta <- sub_tmp_k[, get(i)]
                        
                        d1 <- cbind(x0, x1,
                                    d0[, colnames(d0) %in% colnames(object$complr$output[[parts_index]]$dataout[, -object$complr$output[[parts_index]]$part, with = FALSE]), with = FALSE],
                                    Delta)
                        
                        # # useful information for the final results
                        # d1[, From := rep(sub_from_var, length.out = nrow(d1))[k]]
                        # d1[, To := rep(sub_to_var, length.out = nrow(d1))[k]]
                        # d1[, Delta := abs(as.numeric(Delta))]
                        
                        # remove impossible reallocation that result in negative values
                        cols <- colnames(d1) %sin% c(colnames(x0), colnames(base))
                        d1 <- d1[rowSums(d1[, ..cols] < 0) == 0]
                        
                        # compositions and ilrs for predictions
                        bx0    <- acomp(d1[, colnames(object$complr$output[[parts_index]]$between_comp), with = FALSE], total = object$complr$output[[parts_index]]$total)
                        bcompsub  <- acomp(d1[, object$complr$output[[parts_index]]$parts, with = FALSE], total = object$complr$output[[parts_index]]$total)
                        
                        bz0    <- ilr(bx0, V = object$complr$output[[parts_index]]$psi)
                        bilrsub  <- ilr(bcompsub, V = object$complr$output[[parts_index]]$psi)
                        
                        wz0 <- as.data.table(matrix(0, nrow = nrow(bilrsub), ncol = ncol(bilrsub)))
                        
                        colnames(bilrsub) <- colnames(object$complr$output[[parts_index]]$between_logratio)
                        colnames(wz0)   <- colnames(object$complr$output[[parts_index]]$within_logratio)
                        
                        # prediction
                        dsub <- cbind(d1, bilrsub, wz0)
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NULL,
                            scale = scale,
                            summary = FALSE
                          )
                        ysub <- rowMeans(as.data.frame(ysub))
                        
                        # difference in outcomes between substitution and no change
                        delta_y <- ysub - y0
                        
                        # posterior means and intervals
                        suppressWarnings(posterior_delta_y <- setDT(describe_posterior(delta_y, centrality = "mean", ...)))
                        posterior_delta_y <- posterior_delta_y[, .(Mean, CI_low, CI_high)]
                        posterior_delta_y$Delta <- sub_tmp_j[k, get(i)]
                        kout[[k]] <- posterior_delta_y
                      }
                      jout[[j]] <- rbindlist(kout)
                    }                        
                    
                    jout <- rbindlist(jout)
                    jout[, Delta := as.numeric(Delta)]
                    jout[, From := rep(sub_from_var, length.out = nrow(jout))]
                    jout[, To := sub_to_var]
                    jout[, Level := level]
                    jout[, Reference := ref]
                    
                    names(jout) <- c("Mean", "CI_low", "CI_high",
                                     "Delta", "From", "To", "Level", "Reference")
                    
                    # store final results for entire composition
                    jout <- list(jout)
                    names(jout) <- i
                    jout
                  }
  iout
}

# Clustermean Within-person Substitution model
.get.wsubmargins <- function(object,
                             delta,
                             base,
                             x0,
                             y0,
                             d0,
                             summary,
                             level,
                             ref,
                             scale,
                             type,
                             cores,
                             ...) {
  
  
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }
  oopts <- options(future.globals.maxSize = +Inf, future.globals.onReference = NULL)
  on.exit(options(oopts))
  
  iout <- foreach(i = colnames(base), .combine = c,
                  .options.future = list(packages = "multilevelcoda", seed = TRUE)) %dofuture% {
                    
                    # substitution variables
                    if (type == "one-to-all") {
                      
                      # one to remaining
                      base_tmp <- as.data.table(base)
                      base_tmp <- base_tmp[(get(i) %in% c(1, -1))]
                      
                      sub_from_var <- c("remaining", i)
                      sub_to_var <- c(i, "remaining")
                    }
                    else {
                      # possible pairwise substitution of 1 compositional variable
                      # one to one
                      base_tmp <- as.data.table(base)
                      base_tmp <- base_tmp[(get(i) != 0)]
                      base_tmp <- base_tmp[order(-rank(get(i)))]
                      
                      sub_from_var <- colnames(base_tmp) %snin% eval(i)
                      sub_to_var <- i
                    }
                    
                    # loop substitution
                    kout <- vector("list", length = nrow(base_tmp))
                    jout <- vector("list", length = length(delta))
                    
                    for (j in seq_along(delta)) { # delta level
                      sub_tmp_j <- base_tmp * delta[j]
                      for (k in seq_len(nrow(sub_tmp_j))) {
                        sub_tmp_k <- sub_tmp_j[k, ]
                        sub_tmp_k <- sub_tmp_k[rep(seq_len(nrow(sub_tmp_k)), nrow(x0)), ]
                        x1 <- x0 + sub_tmp_k
                        names(x1) <- object$complr$output[[parts_index]]$parts
                        
                        Delta <- sub_tmp_k[, get(i)]
                        
                        d1 <- cbind(x0, x1,
                                    d0[, colnames(d0) %in% colnames(object$complr$output[[parts_index]]$dataout[, -object$complr$output[[parts_index]]$part, with = FALSE]), with = FALSE],
                                    Delta)
                        
                        # # useful information for the final output
                        # d1[, From := rep(sub_from_var, length.out = nrow(d1))[k]]
                        # d1[, To := rep(sub_to_var, length.out = nrow(d1))[k]]
                        # d1[, Delta := abs(as.numeric(Delta))]
                        
                        # remove impossible reallocation that result in negative values
                        cols <- colnames(d1) %sin% c(colnames(x0), colnames(base))
                        d1 <- d1[rowSums(d1[, ..cols] < 0) == 0]
                        
                        # compositions and ilr for predictions
                        bx0   <- acomp(d1[, colnames(object$complr$output[[parts_index]]$between_comp), with = FALSE], total = object$complr$output[[parts_index]]$total)
                        bcompsub <- acomp(d1[, object$complr$output[[parts_index]]$parts, with = FALSE], total = object$complr$output[[parts_index]]$total)
                        
                        bz0   <- ilr(bx0, V = object$complr$output[[parts_index]]$psi)
                        bilrsub <- ilr(bcompsub, V = object$complr$output[[parts_index]]$psi)
                        
                        wilrsub <- bilrsub - bz0
                        
                        colnames(bz0)   <- colnames(object$complr$output[[parts_index]]$between_logratio)
                        colnames(wilrsub) <- colnames(object$complr$output[[parts_index]]$within_logratio)
                        
                        # substitution data
                        dsub <- cbind(d1, bz0, wilrsub)
                        
                        # prediction
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NULL,
                            scale = scale,
                            summary = FALSE
                          )
                        ysub <- rowMeans(as.data.frame(ysub))
                        
                        # difference between substitution and no change
                        delta_y <- ysub - y0
                        
                        # describe PD of delta y
                        suppressWarnings(posterior_delta_y <- setDT(describe_posterior(delta_y, centrality = "mean", ...)))
                        posterior_delta_y <- posterior_delta_y[, .(Mean, CI_low, CI_high)]
                        posterior_delta_y$Delta <- sub_tmp_j[k, get(i)]
                        kout[[k]] <- posterior_delta_y
                      }
                      # results
                      jout[[j]] <- rbindlist(kout)
                    }
                    
                    jout <- rbindlist(jout)
                    jout[, Delta := as.numeric(Delta)]
                    jout[, From := rep(sub_from_var, length.out = nrow(jout))]
                    jout[, To := sub_to_var]
                    jout[, Level := level]
                    jout[, Reference := ref]
                    
                    names(jout) <- c("Mean", "CI_low", "CI_high",
                                     "Delta", "From", "To", "Level", "Reference")
                    
                    # final results for entire composition
                    jout <- list(jout)
                    names(jout) <- i
                    jout
                  }
  iout
}

# Clustermean Average Substitution
.get.submargins <- function(object,
                            delta,
                            base,
                            x0,
                            y0,
                            d0,
                            summary,
                            level,
                            ref,
                            scale,
                            type,
                            cores,
                            ...) {
  
  
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }
  oopts <- options(future.globals.maxSize = +Inf, future.globals.onReference = NULL)
  on.exit(options(oopts))
  
  iout <- foreach(i = colnames(base), .combine = c,
                  .options.future = list(packages = "multilevelcoda", seed = TRUE)) %dofuture% {
                    
                    # substitution variables
                    if (type == "one-to-all") {
                      
                      # one to remaining
                      base_tmp <- as.data.table(base)
                      base_tmp <- base_tmp[(get(i) %in% c(1, -1))]
                      
                      sub_from_var <- c("remaining", i)
                      sub_to_var <- c(i, "remaining")
                    }
                    else {
                      # possible pairwise substitution of 1 compositional variable
                      # one to one
                      base_tmp <- as.data.table(base)
                      base_tmp <- base_tmp[(get(i) != 0)]
                      base_tmp <- base_tmp[order(-rank(get(i)))]
                      
                      sub_from_var <- colnames(base_tmp) %snin% eval(i)
                      sub_to_var <- i
                    }
                    
                    # loop substitution
                    kout <- vector("list", length = nrow(base_tmp))
                    jout <- vector("list", length = length(delta))
                    
                    for (j in seq_along(delta)) { # delta level
                      sub_tmp_j <- base_tmp * delta[j]
                      for (k in seq_len(nrow(sub_tmp_j))) {
                        sub_tmp_k <- sub_tmp_j[k, ]
                        sub_tmp_k <- sub_tmp_k[rep(seq_len(nrow(sub_tmp_k)), nrow(x0)), ]
                        x1 <- x0 + sub_tmp_k
                        names(x1) <- object$complr$output[[parts_index]]$parts
                        
                        Delta <- sub_tmp_k[, get(i)]
                        
                        d1 <- cbind(x1,
                                    d0[, colnames(d0) %in% colnames(object$complr$output[[parts_index]]$dataout[, -object$complr$output[[parts_index]]$part, with = FALSE]), with = FALSE],
                                    Delta)
                        
                        # # useful information for the final results
                        # d1[, From := rep(sub_from_var, length.out = nrow(d1))[k]]
                        # d1[, To := rep(sub_to_var, length.out = nrow(d1))[k]]
                        # d1[, Delta := abs(as.numeric(Delta))]
                        
                        # remove impossible reallocation that result in negative values
                        cols <- colnames(d1) %sin% c(colnames(x0), colnames(base))
                        d1 <- d1[rowSums(d1[, ..cols] < 0) == 0]
                        
                        # compositions and ilrs for predictions
                        tcomp <- acomp(d1[, object$complr$output[[parts_index]]$parts, with = FALSE], total = object$complr$output[[parts_index]]$total)
                        tilr <- ilr(tcomp, V = object$complr$output[[parts_index]]$psi)
                        
                        colnames(tilr) <- colnames(object$complr$output[[parts_index]]$logratio)
                        
                        # substitution data
                        dsub <- cbind(d1, tilr)
                        
                        # prediction
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NULL,
                            scale = scale,
                            summary = FALSE
                          )
                        ysub <- rowMeans(as.data.frame(ysub))
                        
                        # difference in outcomes between substitution and no change
                        delta_y <- ysub - y0
                        
                        # describe PD of delta y
                        suppressWarnings(posterior_delta_y <- setDT(describe_posterior(delta_y, centrality = "mean", ...)))
                        posterior_delta_y <- posterior_delta_y[, .(Mean, CI_low, CI_high)]
                        posterior_delta_y$Delta <- sub_tmp_j[k, get(i)]
                        kout[[k]] <- posterior_delta_y
                      }
                      jout[[j]] <- rbindlist(kout)
                    }
                    
                    jout <- rbindlist(jout)
                    jout[, Delta := as.numeric(Delta)]
                    jout[, From := rep(sub_from_var, length.out = nrow(jout))]
                    jout[, To := sub_to_var]
                    jout[, Level := level]
                    jout[, Reference := ref]
                    
                    names(jout) <- c("Mean", "CI_low", "CI_high",
                                     "Delta", "From", "To", "Level", "Reference")
                    
                    # store final results for entire composition
                    jout <- list(jout)
                    names(jout) <- i
                    jout
                  }
  iout
}
