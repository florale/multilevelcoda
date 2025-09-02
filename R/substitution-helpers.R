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
  
  ## get the index of which index elements of object$complr$output do the parts correspond to
  idx <- which(vapply(lapply(object$complr$output, function(x) x$parts), function(p) identical(parts, p), logical(1)))
  
  # grab logratio and composition names
  z_vars  <- names(object$complr)[["logratio", paste0("composition_", idx)]]
  bz_vars <- names(object$complr)[["between_logratio", paste0("composition_", idx)]]
  wz_vars <- names(object$complr)[["within_logratio", paste0("composition_", idx)]]
  
  x_vars  <- names(object$complr)[["composition", paste0("composition_", idx)]]
  bx_vars <- names(object$complr)[["between_composition", paste0("composition_", idx)]]
  wx_vars <- names(object$complr)[["within_composition", paste0("composition_", idx)]]
  
  ## NOTES
  ## ignore weight for clustermean
  ## equal weight is default for grandmean
  
  # what type of model is being estimated
  model_fixef <- rownames(fixef(object))
  model_ranef <- if(dim(object$model$ranef)[1] > 0) (names(ranef(object))) else (NULL)
  
  model_fixef_level <- model_fixef_coef <- NULL
  
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
  if (identical(model_ranef_level, "multilevel")) {
    
    # for clustermean
    if ("clustermean" %in% ref) {
      weight <- NULL # ignore weight
      
      ## aggregate
      if ("aggregate" %in% level) {
        d0 <- object$complr$dataout[, head(.SD, 1), by = eval(object$complr$idvar)]
        x0 <- acomp(d0[, x_vars, with = FALSE], total = object$complr$output[[idx]]$total)
        
        z0 <- ilr(x0, V = object$complr$output[[idx]]$psi)
        z0 <- as.data.table(z0)
        
        colnames(z0) <- z_vars
        colnames(x0) <- x_vars
        
        d0 <- cbind(z0, x0, d0[, -colnames(x0), with = FALSE])
      }
      
      ## between and within
      if (any(c("between", "within") %in% level)) {
        d0  <- object$complr$dataout[, head(.SD, 1), by = eval(object$complr$idvar)]
        bx0  <- acomp(d0[, bx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
        
        bz0 <- ilr(bx0, V = object$complr$output[[idx]]$psi)
        bz0 <- as.data.table(bz0)
        
        wx0 <- as.data.table(matrix(1, nrow = nrow(bx0), ncol = ncol(bx0)))
        wz0 <- as.data.table(matrix(0, nrow = nrow(bz0), ncol = ncol(bz0)))
        
        colnames(bz0) <- bz_vars
        colnames(wz0) <- wz_vars
        colnames(bx0) <- bx_vars
        colnames(wx0) <- wx_vars
        
        d0 <- cbind(bz0, wz0, bx0, wx0, d0[, colnames(d0) %in% colnames(object$complr$dataout), with = FALSE])
      }
    } else {
      ## assemble reference grid
      ## get var names
      zs <- c(bz_vars, wz_vars, z_vars)
      
      vars  <- colnames(model.frame(object))
      resp  <- object$model$formula$formula[[2]]
      grp   <- object$model$ranef$group
      preds <- vars %snin% c(resp, grp)
      covs  <- vars %snin% c(resp, grp, zs)
      
      ## default reference grid
      drg <- as.data.table(ref_grid(object$model, at = at)@grid)
      
      ## reference grid (only covariates and outcome) _ FL comment out to prep weighted summary??
      refgrid <- drg[, colnames(drg) %in% c(covs, resp), with = FALSE]
      
      ## to make fitted() happy
      id <- data.table::data.table(1) # to make fitted() happy
      colnames(id) <- object$complr$idvar
      
      # grandmean
      if ("grandmean" %in% ref) {
        
        # aggregate
        if ("aggregate" %in% level) {
          if (weight == "proportional") {
            x0 <- mean.acomp(object$complr$output[[idx]]$X, robust = TRUE)
            
          } else {
            x0 <- object$complr$dataout[, head(.SD, 1), by = eval(object$complr$idvar)]
            x0 <- acomp(x0[, x_vars, with = FALSE], total = object$complr$output[[idx]]$total)
            x0 <- mean.acomp(x0, robust = TRUE)
          }
          
          x0 <- acomp(x0, total = object$complr$output[[idx]]$total)
          x0 <- as.data.table(t(x0))
          
          z0 <- ilr(x0, V = object$complr$output[[idx]]$psi)
          z0 <- as.data.table(t(z0))
          
          colnames(z0) <- z_vars
          colnames(x0) <- x_vars
          
          d0 <- if (all(dim(refgrid) == 0)) (cbind(z0, x0, id)) else (expand.grid.df(z0, x0, id, refgrid))
        }
        
        # between and/or within
        if (any(c("between", "within") %in% level)) {
          if (weight == "proportional") {
            bx0 <- mean.acomp(object$complr$output[[idx]]$bX, robust = TRUE)
            
          } else {
            bx0 <- object$complr$dataout[, head(.SD, 1), by = eval(object$complr$idvar)]
            bx0 <- acomp(bx0[, bx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
            bx0 <- mean.acomp(bx0, robust = TRUE)
          }
          
          bx0 <- acomp(bx0, total = object$complr$output[[idx]]$total)
          bx0 <- as.data.table(t(bx0))
          
          bz0 <- ilr(bx0, V = object$complr$output[[idx]]$psi)
          bz0 <- as.data.table(t(bz0))
          
          wx0 <- as.data.table(matrix(1, nrow = nrow(bx0), ncol = ncol(bx0)))
          wz0 <- as.data.table(matrix(0, nrow = nrow(bz0), ncol = ncol(bz0)))
          
          colnames(bz0) <- bz_vars
          colnames(wz0) <- wz_vars
          colnames(bx0) <- bx_vars
          colnames(wx0) <- wx_vars
          
          d0 <- if (all(dim(refgrid) == 0)) (cbind(bz0, wz0, bx0, wx0, id)) else (expand.grid.df(bz0, wz0, bx0, wx0, id, refgrid))
        }
      }
      
      # user specified
      if (inherits(ref, c("data.table", "data.frame", "matrix"))) {
        weight <- NULL
        
        if (isFALSE(object$complr$output[[idx]]$parts %in% colnames(ref))) {  # get user's composition
          stop(
            sprintf(
              "The reference grid should include all compositional components but (%s) are missing.",
              paste0(object$complr$output[[idx]]$parts %nin% colnames(ref), collapse = ", ")
            ))
        } else {
          xU <- ref[, object$complr$output[[idx]]$parts, with = FALSE]
          xU <- acomp(xU, total = object$complr$output[[idx]]$total)
          xU <- as.data.table(t(xU))
        }
        
        # sanity checks
        if (nrow(ref) > 1) {
          stop("Only one reference composition is allowed at a time.")
        }
        if(isFALSE(sum(xU) == object$complr$output[[idx]]$total)) {
          stop(sprintf(
            "The total amount of the reference composition (%s) should be the same as the composition (%s).",
            sum(xU),
            object$complr$output[[idx]]$total
          ))
        }
        if (isTRUE((any(xU > lapply(object$complr$dataout[, object$complr$output[[idx]]$parts, with = FALSE], max)) |
                    any(xU < lapply(object$complr$dataout[, object$complr$output[[idx]]$parts, with = FALSE], min))))) {
          stop(paste(
            sprintf(
              "composition should be numeric or interger values that are between (%s) and (%s)",
              paste0(round(apply(object$complr$dataout[, object$complr$output[[idx]]$parts, with = FALSE], 2, min)), collapse = ", "),
              paste0(round(apply(object$complr$dataout[, object$complr$output[[idx]]$parts, with = FALSE], 2, max)), collapse = ", ")),
            "\n",
            " for",
            paste0(object$complr$output[[idx]]$parts, collapse = ", "),
            "respectively"
          ))
        }
        
        # user's specified reference grid - edit to allow for new var names
        ## any covariates left in the ref
        if (ncol(ref) > ncol(xU)) {
          covgrid <- ref[, -object$complr$output[[idx]]$parts, with = FALSE]
          
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
            # and fill drg if any is missing
            refgrid <- as.data.table(expand.grid.df(covgrid,
                                                    refgrid[, -colnames(covgrid), with = FALSE]))
          }
        } else {
          refgrid <- drg[, covs, with = FALSE]
        }
        
        if (level == "aggregate") {
          x0 <- xU
          
          z0 <- ilr(x0, V = object$complr$output[[idx]]$psi)
          z0 <- as.data.table(t(z0))
          
          colnames(z0) <- z_vars
          colnames(x0) <- x_vars
          
          d0 <- if (all(dim(refgrid) == 0)) (cbind(z0, x0, id)) else (expand.grid.df(z0, x0, id, refgrid))
          
        }
        if (level %in% c("between", "within")) {
          x0 <- object$complr$dataout[, head(.SD, 1), by = eval(object$complr$idvar)]
          x0 <- acomp(x0[, bx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
          x0 <- mean.acomp(x0, robust = TRUE)
          
          # assemble d0
          # bz0 is between-person ilr of the ref comp (doesn't have to be compositional mean)
          bx0 <- xU
          bz0 <- ilr(bx0, V = object$complr$output[[idx]]$psi)
          bz0 <- as.data.table(t(bz0))
          
          # wx0 and wz0 are the difference between the actual compositional mean of the dataset and bilr
          # is 0 if ref comp is compositional mean
          # but is different if not
          wx0 <- bx0 - x0
          wz0 <- as.data.table(t(ilr(wx0, V = object$complr$output[[idx]]$psi)))
          
          id <- data.table::data.table(1) # to make fitted() happy
          
          colnames(bz0) <- bz_vars
          colnames(wz0) <- wz_vars
          colnames(bx0) <- bx_vars
          colnames(wx0) <- wx_vars
          colnames(id)  <- object$complr$idvar
          
          d0 <- if (all(dim(refgrid) == 0)) (cbind(bz0, wz0, bx0, wx0, id)) else (expand.grid.df(bz0, wz0, bx0, wx0, id, refgrid))
        }
      }
    }
  }
  
  ## d0 and x0 for single level model
  if (model_ranef_level == "single") {
    
    x0 <- object$complr$output[[idx]]$X
    x0 <- mean.acomp(x0, robust = TRUE)
    x0 <- acomp(x0, total = object$complr$output[[idx]]$total)
    x0 <- as.data.table(t(x0))
    
    z0 <- ilr(x0, V = object$complr$output[[idx]]$psi)
    z0 <- as.data.table(t(z0))
    
    colnames(z0) <- z_vars
    colnames(x0) <- x_vars
    
    # assemble reference grid
    # get var names
    zs <- c(z_vars, bz_vars, wz_vars)
    
    vars  <- colnames(model.frame(object))
    resp  <- object$model$formula$formula[[2]]
    # grp   <- object$model$ranef$group
    preds <- vars %snin% c(resp)
    covs  <- vars %snin% c(resp, zs)
    
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
                      parts,
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
  
  ## get the index of which index elements of object$complr$output do the parts correspond to
  idx <- which(vapply(lapply(object$complr$output, function(x) x$parts), function(p) identical(parts, p), logical(1)))
  
  # grab logratio and composition names
  z_vars  <- names(object$complr)[["logratio", paste0("composition_", idx)]]
  bz_vars <- names(object$complr)[["between_logratio", paste0("composition_", idx)]]
  wz_vars <- names(object$complr)[["within_logratio", paste0("composition_", idx)]]
  
  x_vars  <- names(object$complr)[["composition", paste0("composition_", idx)]]
  bx_vars <- names(object$complr)[["between_composition", paste0("composition_", idx)]]
  wx_vars <- names(object$complr)[["within_composition", paste0("composition_", idx)]]
  
  grid <- d0[, colnames(d0) %nin% c(z_vars, bz_vars, wz_vars, x_vars, bx_vars, wx_vars, object$complr$idvar), with = FALSE]
  
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
                      sub_to_var   <- c(i, "remaining")
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
                        xsub <- x0 + sub_tmp_j[k,]
                        names(xsub) <- x_vars
                        Delta <- sub_tmp_j[k, get(i)]
                        kout[[k]] <- cbind(x0, xsub, Delta)
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
                    d1   <- d1[rowSums(d1[, ..cols] < 0) == 0]
                    
                    # compositions and ilrs for predictions
                    bx0   <- acomp(d1[, bx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
                    bxsub <- acomp(d1[, x_vars, with = FALSE], total = object$complr$output[[idx]]$total)
                    
                    bzsub <- ilr(bxsub, V = object$complr$output[[idx]]$psi)
                    wz0   <- as.data.table(matrix(0, nrow = nrow(bzsub), ncol = ncol(bzsub)))
                    
                    colnames(bzsub) <- bz_vars
                    colnames(wz0)   <- wz_vars
                    
                    # reference grid
                    ## get covariate + idvar names
                    covs <- colnames(d0) %snin% c(bz_vars, wz_vars, bx_vars, wx_vars)
                    refgrid <- d0[, covs, with = FALSE]
                    
                    # predictions
                    hout <- vector("list", length = nrow(d0))
                    for (h in seq_len(nrow(d0))) {
                      dsub <- cbind(d1, bzsub, wz0, refgrid[h, ]) # dif for b w t sub
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
                    
                    if (aorg) { # unadj OR adj averaging over reference grid
                      posterior_delta_y <- list(Reduce(`+`, hout) / length(hout))
                    } else { # adj keeping prediction at each level of reference grid
                      posterior_delta_y <- hout
                    }
                    posterior_delta_y <- lapply(posterior_delta_y, function(x) cbind(dsub[, .(Delta, From, To, Level, Reference)], t(x)))
                    
                    # final results for entire composition
                    list(posterior_delta_y)
                  }
  
  if(summary) { 
    ## sub1 <- substitution(object = fit1, delta = 5, level = c("between"), aorg = FALSE, summary = TRUE)
    ## sub1 <- substitution(object = fit1, delta = 5, level = c("between"), aorg = TRUE, summary = TRUE)
    iout <- lapply(iout, function(y) {
      do.call(rbind, Map(function(x, i) {
        dmeta  <- x[,  c("Delta", "From", "To", "Level", "Reference")]
        result <- x[, -c("Delta", "From", "To", "Level", "Reference")]
        result <- rbindlist(lapply(as.data.table(t(result)), describe_posterior, centrality = "mean", ...))
        if(aorg) (cbind(result, dmeta)) else (cbind(result, dmeta, grid[i, ]))
      }, y, seq_along(y)))
    })
    
  } else { 
    ## sub1 <- substitution(object = fit1, delta = 5, level = c("between"), aorg = FALSE, summary = FALSE)
    ## sub1 <- substitution(object = fit1, delta = 5, level = c("between"), aorg = TRUE, summary = FALSE)
    iout <- lapply(seq_along(iout), function(i) {
      if(aorg) (as.data.table(iout[[i]])) else (list(posterior = iout[[i]], grid = as.data.table(grid)))
    })
  }
  
  names(iout) <- parts
  iout
}

# Grandmean Within-person Substitution model
.get.wsub <- function(object,
                      delta,
                      base,
                      parts,
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
  
  ## get the index of which index elements of object$complr$output do the parts correspond to
  idx <- which(vapply(lapply(object$complr$output, function(x) x$parts), function(p) identical(parts, p), logical(1)))
  
  # grab logratio and composition names
  z_vars  <- names(object$complr)[["logratio", paste0("composition_", idx)]]
  bz_vars <- names(object$complr)[["between_logratio", paste0("composition_", idx)]]
  wz_vars <- names(object$complr)[["within_logratio", paste0("composition_", idx)]]
  
  x_vars  <- names(object$complr)[["composition", paste0("composition_", idx)]]
  bx_vars <- names(object$complr)[["between_composition", paste0("composition_", idx)]]
  wx_vars <- names(object$complr)[["within_composition", paste0("composition_", idx)]]
  
  grid <- d0[, colnames(d0) %nin% c(z_vars, bz_vars, wz_vars, x_vars, bx_vars, wx_vars, object$complr$idvar), with = FALSE]
  
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
                      sub_to_var   <- c(i, "remaining")
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
                        xsub <- x0 + sub_tmp_j[k, ]
                        names(xsub) <- x_vars
                        Delta <- sub_tmp_j[k, get(i)]
                        kout[[k]] <- cbind(x0, xsub, Delta)
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
                    bx0   <- acomp(d1[, bx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
                    bxsub <- acomp(d1[, x_vars, with = FALSE], total = object$complr$output[[idx]]$total)
                    
                    bz0   <- ilr(bx0, V = object$complr$output[[idx]]$psi)
                    bzsub <- ilr(bxsub, V = object$complr$output[[idx]]$psi)
                    wzsub <- bzsub - bz0
                    
                    colnames(bz0)   <- bz_vars
                    colnames(wzsub) <- wz_vars
                    
                    # reference grid
                    ## get covariate + idvar names
                    covs <- colnames(d0) %snin% c(bz_vars, wz_vars, bx_vars, wx_vars)
                    refgrid <- d0[, covs, with = FALSE]
                    
                    # predictions
                    hout <- vector("list", length = nrow(d0))
                    for (h in seq_len(nrow(d0))) {
                      dsub <- cbind(d1, bz0, wzsub, refgrid[h, ])
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
                    
                    if (aorg) { # unadj OR adj averaging over reference grid
                      posterior_delta_y <- list(Reduce(`+`, hout) / length(hout))
                    } else { # adj keeping prediction at each level of reference grid
                      posterior_delta_y <- hout
                    }
                    posterior_delta_y <- lapply(posterior_delta_y, function(x) cbind(dsub[, .(Delta, From, To, Level, Reference)], t(x)))
                    
                    # final results for entire composition
                    list(posterior_delta_y)
                  }
  
  if(summary) { 
    iout <- lapply(iout, function(y) {
      do.call(rbind, Map(function(x, i) {
        dmeta  <- x[,  c("Delta", "From", "To", "Level", "Reference")]
        result <- x[, -c("Delta", "From", "To", "Level", "Reference")]
        result <- rbindlist(lapply(as.data.table(t(result)), describe_posterior, centrality = "mean", ...))
        if(aorg) (cbind(result, dmeta)) else (cbind(result, dmeta, grid[i, ]))
      }, y, seq_along(y)))
    })
    
  } else { 
    iout <- lapply(seq_along(iout), function(i) {
      if(aorg) (as.data.table(iout[[i]])) else (list(posterior = iout[[i]], grid = as.data.table(grid)))
    })
  }
  
  names(iout) <- parts
  iout
}

# Grandmean Simple Substitution
.get.sub <- function(object,
                     delta,
                     base,
                     parts,
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
  
  ## get the index of which index elements of object$complr$output do the parts correspond to
  idx <- which(vapply(lapply(object$complr$output, function(x) x$parts), function(p) identical(parts, p), logical(1)))
  
  # grab logratio and composition names
  z_vars  <- names(object$complr)[["logratio", paste0("composition_", idx)]]
  bz_vars <- names(object$complr)[["between_logratio", paste0("composition_", idx)]]
  wz_vars <- names(object$complr)[["within_logratio", paste0("composition_", idx)]]
  
  x_vars  <- names(object$complr)[["composition", paste0("composition_", idx)]]
  bx_vars <- names(object$complr)[["between_composition", paste0("composition_", idx)]]
  wx_vars <- names(object$complr)[["within_composition", paste0("composition_", idx)]]
  
  grid <- d0[, colnames(d0) %nin% c(z_vars, bz_vars, wz_vars, x_vars, bx_vars, wx_vars, object$complr$idvar), with = FALSE]
  
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
                      sub_to_var   <- c(i, "remaining")
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
                        xsub <- x0 + sub_tmp_j[k,]
                        names(xsub) <- x_vars
                        Delta <- sub_tmp_j[k, get(i)]
                        kout[[k]] <- cbind(xsub, Delta)
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
                    compsub  <- acomp(d1[, x_vars, with = FALSE], total = object$complr$output[[idx]]$total)
                    
                    zsub <- ilr(compsub, V = object$complr$output[[idx]]$psi)
                    colnames(zsub) <- z_vars
                    
                    # reference grid
                    ## get covariate + idvar names
                    covs <- colnames(d0) %snin% c(z_vars, x_vars)
                    refgrid <- d0[, covs, with = FALSE]
                    
                    
                    # predictions
                    hout <- vector("list", length = nrow(d0))
                    for (h in seq_len(nrow(d0))) {
                      dsub <- cbind(d1, zsub, refgrid[h, ])
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
                    
                    if (aorg) { # unadj OR adj averaging over reference grid
                      posterior_delta_y <- list(Reduce(`+`, hout) / length(hout))
                    } else { # adj keeping prediction at each level of reference grid
                      posterior_delta_y <- hout
                    }
                    posterior_delta_y <- lapply(posterior_delta_y, function(x) cbind(dsub[, .(Delta, From, To, Level, Reference)], t(x)))
                    
                    # final results for entire composition
                    list(posterior_delta_y)
                  }
  
  if(summary) { 
    iout <- lapply(iout, function(y) {
      do.call(rbind, Map(function(x, i) {
        dmeta  <- x[,  c("Delta", "From", "To", "Level", "Reference")]
        result <- x[, -c("Delta", "From", "To", "Level", "Reference")]
        result <- rbindlist(lapply(as.data.table(t(result)), describe_posterior, centrality = "mean", ...))
        if(aorg) (cbind(result, dmeta)) else (cbind(result, dmeta, grid[i, ]))
      }, y, seq_along(y)))
    })
    
  } else { 
    iout <- lapply(seq_along(iout), function(i) {
      if(aorg) (as.data.table(iout[[i]])) else (list(posterior = iout[[i]], grid = as.data.table(grid)))
    })
  }
  
  names(iout) <- parts
  iout
}                  

# Clustermean Between-person Substitution model
.get.bsubmargins <- function(object,
                             delta,
                             base,
                             parts,
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
  
  ## get the index of which index elements of object$complr$output do the parts correspond to
  idx <- which(vapply(lapply(object$complr$output, function(x) x$parts), function(p) identical(parts, p), logical(1)))
  
  # grab logratio and composition names
  z_vars  <- names(object$complr)[["logratio", paste0("composition_", idx)]]
  bz_vars <- names(object$complr)[["between_logratio", paste0("composition_", idx)]]
  wz_vars <- names(object$complr)[["within_logratio", paste0("composition_", idx)]]
  
  x_vars  <- names(object$complr)[["composition", paste0("composition_", idx)]]
  bx_vars <- names(object$complr)[["between_composition", paste0("composition_", idx)]]
  wx_vars <- names(object$complr)[["within_composition", paste0("composition_", idx)]]
  
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
                      sub_to_var   <- c(i, "remaining")
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
                        xsub <- x0 + sub_tmp_k
                        names(xsub) <- x_vars
                        
                        Delta <- sub_tmp_k[, get(i)]
                        
                        d1 <- cbind(x0, xsub,
                                    d0[, colnames(d0) %in% colnames(object$complr$dataout[, -x_vars, with = FALSE]), with = FALSE],
                                    Delta)
                        
                        # remove impossible reallocation that result in negative values
                        cols <- colnames(d1) %sin% c(colnames(x0), colnames(base))
                        d1 <- d1[rowSums(d1[, ..cols] < 0) == 0]
                        
                        # compositions and ilrs for predictions
                        bx0   <- acomp(d1[, bx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
                        bxsub <- acomp(d1[, x_vars, with = FALSE], total = object$complr$output[[idx]]$total)
                        
                        bz0   <- ilr(bx0, V = object$complr$output[[idx]]$psi)
                        bzsub <- ilr(bxsub, V = object$complr$output[[idx]]$psi)
                        
                        wz0 <- as.data.table(matrix(0, nrow = nrow(bzsub), ncol = ncol(bzsub)))
                        
                        colnames(bzsub) <- bz_vars
                        colnames(wz0)   <- wz_vars
                        
                        # prediction
                        dsub <- cbind(d1, bzsub, wz0)
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
                             parts,
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
  
  ## get the index of which index elements of object$complr$output do the parts correspond to
  idx <- which(vapply(lapply(object$complr$output, function(x) x$parts), function(p) identical(parts, p), logical(1)))
  
  # grab logratio and composition names
  z_vars  <- names(object$complr)[["logratio", paste0("composition_", idx)]]
  bz_vars <- names(object$complr)[["between_logratio", paste0("composition_", idx)]]
  wz_vars <- names(object$complr)[["within_logratio", paste0("composition_", idx)]]
  
  x_vars  <- names(object$complr)[["composition", paste0("composition_", idx)]]
  bx_vars <- names(object$complr)[["between_composition", paste0("composition_", idx)]]
  wx_vars <- names(object$complr)[["within_composition", paste0("composition_", idx)]]
  
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
                      sub_to_var   <- c(i, "remaining")
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
                        xsub <- x0 + sub_tmp_k
                        names(xsub) <- x_vars
                        
                        Delta <- sub_tmp_k[, get(i)]
                        
                        d1 <- cbind(x0, xsub,
                                    d0[, colnames(d0) %in% colnames(object$complr$dataout[, -x_vars, with = FALSE]), with = FALSE],
                                    Delta)
                        
                        # # useful information for the final output
                        # d1[, From := rep(sub_from_var, length.out = nrow(d1))[k]]
                        # d1[, To := rep(sub_to_var, length.out = nrow(d1))[k]]
                        # d1[, Delta := abs(as.numeric(Delta))]
                        
                        # remove impossible reallocation that result in negative values
                        cols <- colnames(d1) %sin% c(colnames(x0), colnames(base))
                        d1 <- d1[rowSums(d1[, ..cols] < 0) == 0]
                        
                        # compositions and ilr for predictions
                        bx0   <- acomp(d1[, bx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
                        bxsub <- acomp(d1[, x_vars, with = FALSE], total = object$complr$output[[idx]]$total)
                        
                        bz0   <- ilr(bx0, V = object$complr$output[[idx]]$psi)
                        bzsub <- ilr(bxsub, V = object$complr$output[[idx]]$psi)
                        
                        wzsub <- bzsub - bz0
                        
                        colnames(bz0)   <- bz_vars
                        colnames(wzsub) <- wz_vars
                        
                        # substitution data
                        dsub <- cbind(d1, bz0, wzsub)
                        
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
                            parts,
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
  
  ## get the index of which index elements of object$complr$output do the parts correspond to
  idx <- which(vapply(lapply(object$complr$output, function(x) x$parts), function(p) identical(parts, p), logical(1)))
  
  # grab logratio and composition names
  z_vars  <- names(object$complr)[["logratio", paste0("composition_", idx)]]
  bz_vars <- names(object$complr)[["between_logratio", paste0("composition_", idx)]]
  wz_vars <- names(object$complr)[["within_logratio", paste0("composition_", idx)]]
  
  x_vars  <- names(object$complr)[["composition", paste0("composition_", idx)]]
  bx_vars <- names(object$complr)[["between_composition", paste0("composition_", idx)]]
  wx_vars <- names(object$complr)[["within_composition", paste0("composition_", idx)]]
  
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
                      sub_to_var   <- c(i, "remaining")
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
                        sub_tmp_k <- sub_tmp_j[k, ]
                        sub_tmp_k <- sub_tmp_k[rep(seq_len(nrow(sub_tmp_k)), nrow(x0)), ]
                        xsub <- x0 + sub_tmp_k
                        names(xsub) <- x_vars
                        
                        Delta <- sub_tmp_k[, get(i)]
                        
                        d1 <- cbind(xsub,
                                    d0[, colnames(d0) %in% colnames(object$complr$dataout[, -x_vars, with = FALSE]), with = FALSE],
                                    Delta)
                        
                        # # useful information for the final results
                        # d1[, From := rep(sub_from_var, length.out = nrow(d1))[k]]
                        # d1[, To := rep(sub_to_var, length.out = nrow(d1))[k]]
                        # d1[, Delta := abs(as.numeric(Delta))]
                        
                        # remove impossible reallocation that result in negative values
                        cols <- colnames(d1) %sin% c(colnames(x0), colnames(base))
                        d1   <- d1[rowSums(d1[, ..cols] < 0) == 0]
                        
                        # compositions and ilrs for predictions
                        xsub <- acomp(d1[, x_vars, with = FALSE], total = object$complr$output[[idx]]$total)
                        zsub <- ilr(xsub, V = object$complr$output[[idx]]$psi)
                        
                        colnames(zsub) <- z_vars
                        
                        # substitution data
                        dsub <- cbind(d1, zsub)
                        
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

