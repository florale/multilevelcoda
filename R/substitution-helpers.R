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
  function(between_simple_sub,
           within_simple_sub,
           simple_sub,
           between_avg_sub,
           within_avg_sub,
           avg_sub,
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

#' Helper functions used only internally to estimate substitution model
#' @importFrom data.table as.data.table data.table copy := setDT rbindlist .SD
#' @importFrom compositions acomp ilr clo mean.acomp
#' @importFrom brms posterior_summary
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
                      at,
                      aorg,
                      summary,
                      level,
                      ref,
                      scale,
                      type,
                      cores,
                      ...) {
  ## get the index of which index elements of object$complr$output do the parts correspond to
  idx <- which(vapply(lapply(object$complr$output, function(x)
    x$parts), function(p)
      identical(parts, p), logical(1)))
  
  # grab logratio and composition names
  z_vars  <- get_variables(object$complr)[["Z", paste0("composition_", idx)]]
  bz_vars <- get_variables(object$complr)[["bZ", paste0("composition_", idx)]]
  wz_vars <- get_variables(object$complr)[["wZ", paste0("composition_", idx)]]
  
  x_vars  <- get_variables(object$complr)[["X", paste0("composition_", idx)]]
  bx_vars <- get_variables(object$complr)[["bX", paste0("composition_", idx)]]
  wx_vars <- get_variables(object$complr)[["wX", paste0("composition_", idx)]]
  
  sx_vars <- paste0("s", object$complr$output[[idx]]$parts)
  
  grid <- d0[, colnames(d0) %nin% c(z_vars,
                                    bz_vars,
                                    wz_vars,
                                    x_vars,
                                    bx_vars,
                                    wx_vars,
                                    object$complr$idvar), with = FALSE]
  grid[, at := if (!is.null(at))
    names(at)
    else
      NA]
  
  # setup parallel processing
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }
  
  oopts <- options(future.globals.maxSize = +Inf,
                   future.globals.onReference = NULL)
  on.exit(options(oopts))
  
  # substitution loop
  iout <- foreach(
    i = colnames(base),
    .combine = c,
    .options.future = list(packages = "multilevelcoda", seed = TRUE)
  ) %dofuture% {
    # substitution variables
    if (type == "one-to-all") {
      # one to remaining
      base_i <- as.data.table(base)
      base_i <- base_i[(get(i) %in% c(1, -1))]
      
      sub_from_var <- c("remaining", i)
      sub_to_var   <- c(i, "remaining")
    }
    else {
      # possible pairwise substitution of 1 compositional variable
      # one to one
      base_i <- as.data.table(base)
      base_i <- base_i[(get(i) != 0)]
      base_i <- base_i[order(-rank(get(i)))]
      
      sub_from_var <- colnames(base_i) %snin% eval(i)
      sub_to_var   <- i
    }
    
    # loop substitution
    kout <- vector("list", length = nrow(base_i))
    jout <- vector("list", length = length(delta))
    
    for (j in seq_along(delta)) {
      # delta level
      sub_delta_j <- base_i * delta[j]
      for (k in seq_len(nrow(sub_delta_j))) {
        xsub <- x0 + sub_delta_j[k, ]
        x0_xsub_delta_k <- cbind(x0, xsub, sub_delta_j[k, get(i)])
        x0_xsub_delta_k <- setNames(x0_xsub_delta_k, c(bx_vars, sx_vars, "Delta"))
        kout[[k]]     <- x0_xsub_delta_k
      }
      jout[[j]] <- rbindlist(kout)
    }
    d1 <- rbindlist(jout)
    
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
    bxsub <- acomp(d1[, sx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
    
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
    
    if (aorg) {
      # unadj OR adj averaging over reference grid
      weight <- grid$.wgt. / sum(grid$.wgt.)

      weighted_hout <- Map(function(x, w) x * w, hout, weight)
      posterior_delta_y <- list(Reduce(`+`, weighted_hout) / length(weighted_hout))
      
    } else {
      # adj keeping prediction at each level of at
      grid[, wgt_at := sum(.wgt.), by = names(at)]
      at_weight <- grid$wgt_at / sum(grid$wgt_at)
      at_levels <- grid[, names(at), with = FALSE]
      at_id <- at_levels[, idx := .I][, .(idx_list = list(idx)), by = names(at)]$idx_list
      
      at_weighted_hout <- Map(function(x, w) x * w, hout, at_weight)
      posterior_delta_y <- lapply(at_id, function(idx) {
        Reduce(`+`, at_weighted_hout[idx]) / length(idx)
      })
      
      unique_at_levels <- unique(at_levels[, names(at), with = FALSE])
      names(posterior_delta_y) <- apply(unique_at_levels, 1, function(x)
        paste(paste0(colnames(unique_at_levels), x), collapse = "_"))
    }
    
    posterior_delta_y <- lapply(posterior_delta_y, function(x)
      cbind(dsub[, .(Delta, From, To, Level, Reference)], t(x)))
    
    # final results for entire composition
    list(posterior_delta_y)
  }
  
  if (summary) {
    ## sub1 <- substitution(object = fit1, delta = 5, level = c("between"), at = list(Female = c(0,1)), summary = TRUE)
    ## sub1 <- substitution(object = fit1, delta = 5, level = c("between"), aorg = TRUE, summary = TRUE)
    iout <- lapply(iout, function(y) {
      do.call(rbind, Map(function(x, i) {
        dmeta  <- x[, c("Delta", "From", "To", "Level", "Reference")]
        result <- apply(x[, -c("Delta", "From", "To", "Level", "Reference")], 1, posterior_summary, ...)
        row.names(result) <- c("Estimate", "Est.Error", "CI_low", "CI_high")
        if (aorg)
          cbind(t(result), dmeta)
        else
          cbind(t(result), dmeta, grid[i, names(at), with = FALSE])
      }, y, seq_along(y)))
    })
    
  } else {
    ## sub1 <- substitution(object = fit1, delta = 5, level = c("between"), aorg = FALSE, summary = FALSE)
    ## sub1 <- substitution(object = fit1, delta = 5, level = c("between"), aorg = TRUE, summary = FALSE)
    iout <- lapply(seq_along(iout), function(i) {
      if (aorg)
        as.data.table(iout[[i]])
      else
        list(posterior = iout[[i]], grid = as.data.table(grid[i, names(at), with = FALSE]))
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
                      at,
                      aorg,
                      summary,
                      level,
                      ref,
                      scale,
                      type,
                      cores,
                      ...) {
  ## get the index of which index elements of object$complr$output do the parts correspond to
  idx <- which(vapply(lapply(object$complr$output, function(x)
    x$parts), function(p)
      identical(parts, p), logical(1)))
  
  # grab logratio and composition names
  z_vars  <- get_variables(object$complr)[["Z", paste0("composition_", idx)]]
  bz_vars <- get_variables(object$complr)[["bZ", paste0("composition_", idx)]]
  wz_vars <- get_variables(object$complr)[["wZ", paste0("composition_", idx)]]
  
  x_vars  <- get_variables(object$complr)[["X", paste0("composition_", idx)]]
  bx_vars <- get_variables(object$complr)[["bX", paste0("composition_", idx)]]
  wx_vars <- get_variables(object$complr)[["wX", paste0("composition_", idx)]]
  
  sx_vars <- paste0("s", object$complr$output[[idx]]$parts)
  
  grid <- d0[, colnames(d0) %nin% c(z_vars,
                                    bz_vars,
                                    wz_vars,
                                    x_vars,
                                    bx_vars,
                                    wx_vars,
                                    object$complr$idvar), with = FALSE]
  grid[, at := if (!is.null(at))
    names(at)
    else
      NA]
  
  # setup parallel processing
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }
  
  oopts <- options(future.globals.maxSize = +Inf,
                   future.globals.onReference = NULL)
  on.exit(options(oopts))
  
  # substitution loop
  iout <- foreach(
    i = colnames(base),
    .combine = c,
    .options.future = list(packages = "multilevelcoda", seed = TRUE)
  ) %dofuture% {
    # substitution variables
    if (type == "one-to-all") {
      # one to remaining
      base_i <- as.data.table(base)
      base_i <- base_i[(get(i) %in% c(1, -1))]
      
      sub_from_var <- c("remaining", i)
      sub_to_var   <- c(i, "remaining")
    }
    else {
      # possible pairwise substitution of 1 compositional variable
      # one to one
      base_i <- as.data.table(base)
      base_i <- base_i[(get(i) != 0)]
      base_i <- base_i[order(-rank(get(i)))]
      
      sub_from_var <- colnames(base_i) %snin% eval(i)
      sub_to_var <- i
    }
    
    # loop substitution
    kout <- vector("list", length = nrow(base_i))
    jout <- vector("list", length = length(delta))
    
    for (j in seq_along(delta)) {
      # delta level
      sub_delta_j <- base_i * delta[j]
      for (k in seq_len(nrow(sub_delta_j))) {
        xsub <- x0 + sub_delta_j[k, ]
        x0_xsub_delta_k <- cbind(x0, xsub, sub_delta_j[k, get(i)])
        x0_xsub_delta_k <- setNames(x0_xsub_delta_k, c(bx_vars, sx_vars, "Delta"))
        kout[[k]]       <- x0_xsub_delta_k
      }
      jout[[j]] <- rbindlist(kout)
    }
    d1 <- rbindlist(jout)
    
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
    bxsub <- acomp(d1[, sx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
    
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
    
    if (aorg) {
      # unadj OR adj averaging over reference grid
      weight <- grid$.wgt. / sum(grid$.wgt.)
      
      weighted_hout <- Map(function(x, w) x * w, hout, weight)
      posterior_delta_y <- list(Reduce(`+`, weighted_hout) / length(weighted_hout))
      
    } else {
      # adj keeping prediction at each level of at
      grid[, wgt_at := sum(.wgt.), by = names(at)]
      at_weight <- grid$wgt_at / sum(grid$wgt_at)
      at_levels <- grid[, names(at), with = FALSE]
      at_id <- at_levels[, idx := .I][, .(idx_list = list(idx)), by = names(at)]$idx_list
      
      at_weighted_hout <- Map(function(x, w) x * w, hout, at_weight)
      posterior_delta_y <- lapply(at_id, function(idx) {
        Reduce(`+`, at_weighted_hout[idx]) / length(idx)
      })
      
      unique_at_levels <- unique(at_levels[, names(at), with = FALSE])
      names(posterior_delta_y) <- apply(unique_at_levels, 1, function(x)
        paste(paste0(colnames(unique_at_levels), x), collapse = "_"))
    }
    
    posterior_delta_y <- lapply(posterior_delta_y, function(x)
      cbind(dsub[, .(Delta, From, To, Level, Reference)], t(x)))
    
    # final results for entire composition
    list(posterior_delta_y)
  }
  
  if (summary) {
    iout <- lapply(iout, function(y) {
      do.call(rbind, Map(function(x, i) {
        dmeta  <- x[, c("Delta", "From", "To", "Level", "Reference")]
        result <- apply(x[, -c("Delta", "From", "To", "Level", "Reference")], 1, posterior_summary, ...)
        row.names(result) <- c("Estimate", "Est.Error", "CI_low", "CI_high")
        if (aorg)
          cbind(t(result), dmeta)
        else
          cbind(t(result), dmeta, grid[i, names(at), with = FALSE])
      }, y, seq_along(y)))
    })
    
  } else {
    iout <- lapply(seq_along(iout), function(i) {
      if (aorg)
        as.data.table(iout[[i]])
      else
        list(posterior = iout[[i]], grid = as.data.table(grid[i, names(at), with = FALSE]))
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
                     at,
                     aorg,
                     summary,
                     level,
                     ref,
                     scale,
                     type,
                     cores,
                     ...) {
  ## get the index of which index elements of object$complr$output do the parts correspond to
  idx <- which(vapply(lapply(object$complr$output, function(x)
    x$parts), function(p)
      identical(parts, p), logical(1)))
  
  # grab logratio and composition names
  z_vars  <- get_variables(object$complr)[["Z", paste0("composition_", idx)]]
  bz_vars <- get_variables(object$complr)[["bZ", paste0("composition_", idx)]]
  wz_vars <- get_variables(object$complr)[["wZ", paste0("composition_", idx)]]
  
  x_vars  <- get_variables(object$complr)[["X", paste0("composition_", idx)]]
  bx_vars <- get_variables(object$complr)[["bX", paste0("composition_", idx)]]
  wx_vars <- get_variables(object$complr)[["wX", paste0("composition_", idx)]]
  
  sx_vars <- paste0("s", object$complr$output[[idx]]$parts)
  
  grid <- d0[, colnames(d0) %nin% c(z_vars,
                                    bz_vars,
                                    wz_vars,
                                    x_vars,
                                    bx_vars,
                                    wx_vars,
                                    object$complr$idvar), with = FALSE]
  grid[, at := if (!is.null(at))
    names(at)
    else
      NA]
  
  # setup parallel processing
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }
  
  oopts <- options(future.globals.maxSize = +Inf,
                   future.globals.onReference = NULL)
  on.exit(options(oopts))
  
  # substitution loop
  iout <- foreach(
    i = colnames(base),
    .combine = c,
    .options.future = list(packages = "multilevelcoda", seed = TRUE)
  ) %dofuture% {
    # substitution variables
    if (type == "one-to-all") {
      # one to remaining
      base_i <- as.data.table(base)
      base_i <- base_i[(get(i) %in% c(1, -1))]
      
      sub_from_var <- c("remaining", i)
      sub_to_var   <- c(i, "remaining")
    }
    else {
      # possible pairwise substitution of 1 compositional variable
      # one to one
      base_i <- as.data.table(base)
      base_i <- base_i[(get(i) != 0)]
      base_i <- base_i[order(-rank(get(i)))]
      
      sub_from_var <- colnames(base_i) %snin% eval(i)
      sub_to_var <- i
    }
    
    # loop substitution
    kout <- vector("list", length = nrow(base_i))
    jout <- vector("list", length = length(delta))
    
    for (j in seq_along(delta)) {
      # delta level
      sub_delta_j <- base_i * delta[j]
      for (k in seq_len(nrow(sub_delta_j))) {
        xsub <- x0 + sub_delta_j[k, ]
        xsub_delta_k <- cbind(xsub, sub_delta_j[k, get(i)])
        xsub_delta_k <- setNames(xsub_delta_k, c(sx_vars, "Delta"))
        kout[[k]]  <- xsub_delta_k
      }
      jout[[j]] <- rbindlist(kout)
    }
    d1 <- rbindlist(jout)
    
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
    xsub  <- acomp(d1[, sx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
    zsub  <- ilr(xsub, V = object$complr$output[[idx]]$psi)
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
    
    if (aorg) {
      # unadj OR adj averaging over reference grid
      weight <- grid$.wgt. / sum(grid$.wgt.)
      
      weighted_hout <- Map(function(x, w) x * w, hout, weight)
      posterior_delta_y <- list(Reduce(`+`, weighted_hout) / length(weighted_hout))
      
    } else {
      # adj keeping prediction at each level of at
      grid[, wgt_at := sum(.wgt.), by = names(at)]
      at_weight <- grid$wgt_at / sum(grid$wgt_at)
      at_levels <- grid[, names(at), with = FALSE]
      at_id <- at_levels[, idx := .I][, .(idx_list = list(idx)), by = names(at)]$idx_list
      
      at_weighted_hout <- Map(function(x, w) x * w, hout, at_weight)
      posterior_delta_y <- lapply(at_id, function(idx) {
        Reduce(`+`, at_weighted_hout[idx]) / length(idx)
      })
      
      unique_at_levels <- unique(at_levels[, names(at), with = FALSE])
      names(posterior_delta_y) <- apply(unique_at_levels, 1, function(x)
        paste(paste0(colnames(unique_at_levels), x), collapse = "_"))
    }
    
    posterior_delta_y <- lapply(posterior_delta_y, function(x)
      cbind(dsub[, .(Delta, From, To, Level, Reference)], t(x)))
    
    # final results for entire composition
    list(posterior_delta_y)
  }
  
  if (summary) {
    iout <- lapply(iout, function(y) {
      do.call(rbind, Map(function(x, i) {
        dmeta  <- x[, c("Delta", "From", "To", "Level", "Reference")]
        result <- apply(x[, -c("Delta", "From", "To", "Level", "Reference")], 1, posterior_summary, ...)
        row.names(result) <- c("Estimate", "Est.Error", "CI_low", "CI_high")
        if (aorg)
          cbind(t(result), dmeta)
        else
          cbind(t(result), dmeta, grid[i, names(at), with = FALSE])
      }, y, seq_along(y)))
    })
    
  } else {
    iout <- lapply(seq_along(iout), function(i) {
      if (aorg)
        as.data.table(iout[[i]])
      else
        list(posterior = iout[[i]], grid = as.data.table(grid[i, names(at), with = FALSE]))
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
  idx <- which(vapply(lapply(object$complr$output, function(x)
    x$parts), function(p)
      identical(parts, p), logical(1)))
  
  # grab logratio and composition names
  z_vars  <- get_variables(object$complr)[["Z", paste0("composition_", idx)]]
  bz_vars <- get_variables(object$complr)[["bZ", paste0("composition_", idx)]]
  wz_vars <- get_variables(object$complr)[["wZ", paste0("composition_", idx)]]
  
  x_vars  <- get_variables(object$complr)[["X", paste0("composition_", idx)]]
  bx_vars <- get_variables(object$complr)[["bX", paste0("composition_", idx)]]
  wx_vars <- get_variables(object$complr)[["wX", paste0("composition_", idx)]]
  
  vars <- c(z_vars, bz_vars, wz_vars, x_vars, bx_vars, wx_vars)
  
  sx_vars <- paste0("s", object$complr$output[[idx]]$parts)
  
  # setup parallel processing
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }
  
  oopts <- options(future.globals.maxSize = +Inf,
                   future.globals.onReference = NULL)
  on.exit(options(oopts))
  
  # substitution loop
  iout <- foreach(
    i = colnames(base),
    .combine = c,
    .options.future = list(packages = "multilevelcoda", seed = TRUE)
  ) %dofuture% {
    # substitution variables
    if (type == "one-to-all") {
      # one to remaining
      base_i <- as.data.table(base)
      base_i <- base_i[(get(i) %in% c(1, -1))]
      
      sub_from_var <- c("remaining", i)
      sub_to_var   <- c(i, "remaining")
    }
    else {
      # possible pairwise substitution of 1 compositional variable
      # one to one
      base_i <- as.data.table(base)
      base_i <- base_i[(get(i) != 0)]
      base_i <- base_i[order(-rank(get(i)))]
      
      sub_from_var <- colnames(base_i) %snin% eval(i)
      sub_to_var <- i
    }
    
    # loop substitution
    kout <- vector("list", length = nrow(base_i))
    jout <- vector("list", length = length(delta))
    
    for (j in seq_along(delta)) {
      # delta level
      sub_delta_j <- base_i * delta[j]
      for (k in seq_len(nrow(sub_delta_j))) {
        # reallocation level
        sub_delta_k <- sub_delta_j[k, ]
        sub_delta_k <- sub_delta_k[rep(seq_len(nrow(sub_delta_k)), nrow(x0)), ]
        xsub <- x0 + sub_delta_k
        
        x0_xsub_delta_k <- cbind(x0, xsub, sub_delta_k[, get(i)])
        x0_xsub_delta_k <- setNames(x0_xsub_delta_k, c(bx_vars, sx_vars, "Delta"))
        d1 <- cbind(x0_xsub_delta_k, d0[, colnames(d0) %nin% vars, with = FALSE])
        
        # remove impossible reallocation that result in negative values
        cols <- colnames(d1) %sin% c(colnames(x0), colnames(base))
        d1 <- d1[rowSums(d1[, ..cols] < 0) == 0]
        
        # compositions and ilrs for predictions
        bx0   <- acomp(d1[, bx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
        bxsub <- acomp(d1[, sx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
        
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
        kout[[k]] <- cbind(as.data.table(t(delta_y)), dsub[1, .(Delta)]) # take first delta row as all rows are identical
      }
      jout[[j]] <- rbindlist(kout)
    }
    jout <- lapply(jout, function(out) {
      out[, Delta := as.numeric(Delta)]
      out[, From := rep(sub_from_var, length.out = nrow(out))]
      out[, To := sub_to_var]
      out[, Level := level]
      out[, Reference := ref]
    })
    jout
  }
  
  # split by part
  iout <- split(iout, factor(vapply(iout, function(x)
    unique(x$To), character(1)), levels = parts))
  
  if (summary) {
    iout <- lapply(iout, function(y) {
      rbindlist(lapply(y, function(x) {
        dmeta  <- unique(x[, c("Delta", "From", "To", "Level", "Reference")])
        result <- apply(x[, -c("Delta", "From", "To", "Level", "Reference")], 1, posterior_summary, ...)
        row.names(result) <- c("Estimate", "Est.Error", "CI_low", "CI_high")
        cbind(t(result), dmeta)
      }))
    })
  } else {
    iout <- lapply(iout, function(x)
      rbindlist(lapply(x, as.data.table)))
  }
  
  names(iout) <- parts
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
  idx <- which(vapply(lapply(object$complr$output, function(x)
    x$parts), function(p)
      identical(parts, p), logical(1)))
  
  # grab logratio and composition names
  z_vars  <- get_variables(object$complr)[["Z", paste0("composition_", idx)]]
  bz_vars <- get_variables(object$complr)[["bZ", paste0("composition_", idx)]]
  wz_vars <- get_variables(object$complr)[["wZ", paste0("composition_", idx)]]
  
  x_vars  <- get_variables(object$complr)[["X", paste0("composition_", idx)]]
  bx_vars <- get_variables(object$complr)[["bX", paste0("composition_", idx)]]
  wx_vars <- get_variables(object$complr)[["wX", paste0("composition_", idx)]]
  
  sx_vars <- paste0("s", object$complr$output[[idx]]$parts)
  
  # setup parallel processing
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }
  
  oopts <- options(future.globals.maxSize = +Inf,
                   future.globals.onReference = NULL)
  on.exit(options(oopts))
  
  # substitution loop
  iout <- foreach(
    i = colnames(base),
    .combine = c,
    .options.future = list(packages = "multilevelcoda", seed = TRUE)
  ) %dofuture% {
    # substitution variables
    if (type == "one-to-all") {
      # one to remaining
      base_i <- as.data.table(base)
      base_i <- base_i[(get(i) %in% c(1, -1))]
      
      sub_from_var <- c("remaining", i)
      sub_to_var   <- c(i, "remaining")
    }
    else {
      # possible pairwise substitution of 1 compositional variable
      # one to one
      base_i <- as.data.table(base)
      base_i <- base_i[(get(i) != 0)]
      base_i <- base_i[order(-rank(get(i)))]
      
      sub_from_var <- colnames(base_i) %snin% eval(i)
      sub_to_var <- i
    }
    
    # loop substitution
    kout <- vector("list", length = nrow(base_i))
    jout <- vector("list", length = length(delta))
    
    for (j in seq_along(delta)) {
      # delta level
      sub_delta_j <- base_i * delta[j]
      for (k in seq_len(nrow(sub_delta_j))) {
        sub_delta_k <- sub_delta_j[k, ]
        sub_delta_k <- sub_delta_k[rep(seq_len(nrow(sub_delta_k)), nrow(x0)), ]
        xsub <- x0 + sub_delta_k
        
        x0_xsub_delta_k <- cbind(x0, xsub, sub_delta_k[, get(i)])
        x0_xsub_delta_k <- setNames(x0_xsub_delta_k, c(bx_vars, sx_vars, "Delta"))
        d1 <- cbind(x0_xsub_delta_k, d0[, colnames(d0) %nin% colnames(x0_xsub_delta_k), with = FALSE])
        
        # remove impossible reallocation that result in negative values
        cols <- colnames(d1) %sin% c(colnames(x0), colnames(base))
        d1 <- d1[rowSums(d1[, ..cols] < 0) == 0]
        
        # compositions and ilr for predictions
        bx0   <- acomp(d1[, bx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
        bxsub <- acomp(d1[, sx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
        
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
        
        # difference in outcomes between substitution and no change
        delta_y <- ysub - y0
        kout[[k]] <- cbind(as.data.table(t(delta_y)), dsub[1, .(Delta)]) # take first delta row as all rows are identical
      }
      jout[[j]] <- rbindlist(kout)
    }
    jout <- lapply(jout, function(out) {
      out[, Delta := as.numeric(Delta)]
      out[, From := rep(sub_from_var, length.out = nrow(out))]
      out[, To := sub_to_var]
      out[, Level := level]
      out[, Reference := ref]
    })
    jout
  }
  
  # split by part
  iout <- split(iout, factor(vapply(iout, function(x)
    unique(x$To), character(1)), levels = parts))
  
  if (summary) {
    iout <- lapply(iout, function(y) {
      rbindlist(lapply(y, function(x) {
        dmeta  <- unique(x[, c("Delta", "From", "To", "Level", "Reference")])
        result <- apply(x[, -c("Delta", "From", "To", "Level", "Reference")], 1, posterior_summary, ...)
        row.names(result) <- c("Estimate", "Est.Error", "CI_low", "CI_high")
        cbind(t(result), dmeta)
      }))
    })
  } else {
    iout <- lapply(iout, function(x)
      rbindlist(lapply(x, as.data.table)))
  }
  
  names(iout) <- parts
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
  idx <- which(vapply(lapply(object$complr$output, function(x)
    x$parts), function(p)
      identical(parts, p), logical(1)))
  
  # grab logratio and composition names
  z_vars  <- get_variables(object$complr)[["Z", paste0("composition_", idx)]]
  bz_vars <- get_variables(object$complr)[["bZ", paste0("composition_", idx)]]
  wz_vars <- get_variables(object$complr)[["wZ", paste0("composition_", idx)]]
  
  x_vars  <- get_variables(object$complr)[["X", paste0("composition_", idx)]]
  bx_vars <- get_variables(object$complr)[["bX", paste0("composition_", idx)]]
  wx_vars <- get_variables(object$complr)[["wX", paste0("composition_", idx)]]
  
  sx_vars <- paste0("s", object$complr$output[[idx]]$parts)
  
  # setup parallel processing
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }
  
  oopts <- options(future.globals.maxSize = +Inf,
                   future.globals.onReference = NULL)
  on.exit(options(oopts))
  
  # substitution loop
  iout <- foreach(
    i = colnames(base),
    .combine = c,
    .options.future = list(packages = "multilevelcoda", seed = TRUE)
  ) %dofuture% {
    # substitution variables
    if (type == "one-to-all") {
      # one to remaining
      base_i <- as.data.table(base)
      base_i <- base_i[(get(i) %in% c(1, -1))]
      
      sub_from_var <- c("remaining", i)
      sub_to_var   <- c(i, "remaining")
    }
    else {
      # possible pairwise substitution of 1 compositional variable
      # one to one
      base_i <- as.data.table(base)
      base_i <- base_i[(get(i) != 0)]
      base_i <- base_i[order(-rank(get(i)))]
      
      sub_from_var <- colnames(base_i) %snin% eval(i)
      sub_to_var   <- i
    }
    
    # loop substitution
    kout <- vector("list", length = nrow(base_i))
    jout <- vector("list", length = length(delta))
    
    for (j in seq_along(delta)) {
      # delta level
      sub_delta_j <- base_i * delta[j]
      for (k in seq_len(nrow(sub_delta_j))) {
        sub_delta_k <- sub_delta_j[k, ]
        sub_delta_k <- sub_delta_k[rep(seq_len(nrow(sub_delta_k)), nrow(x0)), ]
        xsub <- x0 + sub_delta_k
        
        xsub_delta_k <- cbind(xsub, sub_delta_k[, get(i)])
        xsub_delta_k <- setNames(xsub_delta_k, c(sx_vars, "Delta"))
        d1 <- cbind(xsub_delta_k, d0[, colnames(d0) %nin% colnames(xsub_delta_k), with = FALSE])
        
        # remove impossible reallocation that result in negative values
        cols <- colnames(d1) %sin% c(colnames(x0), colnames(base))
        d1   <- d1[rowSums(d1[, ..cols] < 0) == 0]
        
        # compositions and ilrs for predictions
        xsub <- acomp(d1[, sx_vars, with = FALSE], total = object$complr$output[[idx]]$total)
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
        kout[[k]] <- cbind(as.data.table(t(delta_y)), dsub[1, .(Delta)]) # take first delta row as all rows are identical
      }
      jout[[j]] <- rbindlist(kout)
    }
    jout <- lapply(jout, function(out) {
      out[, Delta := as.numeric(Delta)]
      out[, From := rep(sub_from_var, length.out = nrow(out))]
      out[, To := sub_to_var]
      out[, Level := level]
      out[, Reference := ref]
    })
    jout
  }
  
  # split by part
  iout <- split(iout, factor(vapply(iout, function(x)
    unique(x$To), character(1)), levels = parts))
  
  if (summary) {
    iout <- lapply(iout, function(y) {
      rbindlist(lapply(y, function(x) {
        dmeta  <- unique(x[, c("Delta", "From", "To", "Level", "Reference")])
        result <- apply(x[, -c("Delta", "From", "To", "Level", "Reference")], 1, posterior_summary, ...)
        row.names(result) <- c("Estimate", "Est.Error", "CI_low", "CI_high")
        cbind(t(result), dmeta)
      }))
    })
  } else {
    iout <- lapply(iout, function(x)
      rbindlist(lapply(x, as.data.table)))
  }
  
  names(iout) <- parts
  iout
}
