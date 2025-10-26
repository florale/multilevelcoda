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
#' @param between_avg_sub A list of results from \code{bsubmargin} or \code{NULL}
#' @param within_simple_sub A list of results from \code{wsub} or \code{NULL}
#' @param within_avg_sub A list of results from \code{wsubmargin} or \code{NULL}
#' @param simple_sub A list of results from \code{sub} or \code{NULL}
#' @param avg_sub A list of results from \code{submargin} or \code{NULL}
#' @param brmsformula A \code{brmsformula} object
#' @param delta A numeric vector of the amount of substitution
#' @param ref A character value specifying the reference grid
#' @param level A character value specifying the level of substitution
#' @param parts The parts of the composition
#' @param weight The weight to use in calculation of the reference composition
#' @param at An named list of levels for the corresponding variables in the reference grid or \code{NULL}
#' @param type A character value specifying the type of substitution.
#'
#' @seealso \code{\link{substitution}}
#'
#' @return An object of class \code{substitution}
#'
#' @noRd
create_substitution <- function(between_simple_sub,
                                within_simple_sub,
                                simple_sub,
                                between_avg_sub,
                                within_avg_sub,
                                avg_sub,
                                brmsformula,
                                delta,
                                ref,
                                level,
                                parts,
                                weight,
                                at,
                                type) {
  stopifnot(is.list(between_simple_sub) || is.null(between_simple_sub))
  stopifnot(is.list(within_simple_sub) || is.null(within_simple_sub))
  stopifnot(is.list(simple_sub) || is.null(simple_sub))
  stopifnot(is.list(between_avg_sub) || is.null(between_avg_sub))
  stopifnot(is.list(within_avg_sub) || is.null(within_avg_sub))
  stopifnot(is.list(avg_sub) || is.null(avg_sub))

  structure(
    list(
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
      parts = parts,
      weight = weight,
      at = at,
      type = type
    ),
    class = "substitution"
  )
}


#' Substitution analysis helper functions
#'
#' Functions used only internally to estimate substitution model
#'
#' @importFrom data.table as.data.table data.table copy := setDT rbindlist .SD .I
#' @importFrom stats setNames
#' @importFrom compositions acomp ilr clo mean.acomp
#' @importFrom brms posterior_summary
#' @importFrom extraoperators %snin% %sin% %ain%
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
  # extract variables from complr and brmcoda objects for use in substitution models
  get_vars <- .get.subvars(object = object, parts = parts, scale = scale)
  
  grid <- d0[, colnames(d0) %nin% c(get_vars[["Xxz"]], object[["complr"]][["idvar"]]), with = FALSE]
  # grid[, at := if (!is.null(at)) {
  #   names(at)
  # } else {
  #   NA
  # }]

  # setup parallel processing
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }

  oopts <- options(
    future.globals.maxSize = +Inf,
    future.globals.onReference = NULL
  )
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
      sub_to_var <- c(i, "remaining")
    } else {
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
        x0_xsub_delta_k <- setNames(x0_xsub_delta_k, c(get_vars[["XbX"]], get_vars[["sX"]], "Delta"))
        kout[[k]] <- x0_xsub_delta_k
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
    bx0 <- acomp(d1[, get_vars[["XbX"]], with = FALSE], total = object[["complr"]][["output"]][[get_vars[["idx"]]]][["total"]])
    bxsub <- acomp(d1[, get_vars[["sX"]], with = FALSE], total = object[["complr"]][["output"]][[get_vars[["idx"]]]][["total"]])

    bzsub <- ilr(bxsub, V = object[["complr"]][["output"]][[get_vars[["idx"]]]][["psi"]])
    wz0 <- as.data.table(matrix(0, nrow = nrow(bzsub), ncol = ncol(bzsub)))

    colnames(bzsub) <- get_vars[["XbZ"]]
    colnames(wz0) <- get_vars[["XwZ"]]

    # reference grid
    ## get covariate + idvar names
    covs <- colnames(d0) %snin% c(get_vars[["XbZ"]], get_vars[["XwZ"]], get_vars[["XbX"]], get_vars[["XwX"]])
    refgrid <- d0[, covs, with = FALSE]

    # predictions
    hout <- vector("list", length = nrow(d0))
    for (h in seq_len(nrow(d0))) {
      dsub <- cbind(d1, bzsub, wz0, refgrid[h, ]) # dif for b w t sub
      ysub <- fitted(
        object,
        newdata = dsub,
        re_formula = NA,
        scale = scale,
        summary = FALSE
      )
      if (length(get_vars[["brmcoda"]][["y"]]) > 1) {
        # when multiple outcomes
        ## when outcome is compositional
        ## delta y is acomp change
        if (all(sort(dimnames(ysub)[[3]]) %in% c(
          sort(get_vars$XX),
          sort(get_vars$XwX),
          sort(get_vars$XbX)
        ))) {
          delta_y <- lapply(seq(dim(ysub)[2]), function(j) {
            ## j = delta, h = ref grid
            acomp(ysub[, j, ], total = object[["complr"]][["output"]][[get_vars[["idy"]]]][["total"]])
            - acomp(y0[, h, ], total = object[["complr"]][["output"]][[get_vars[["idy"]]]][["total"]])
          })
        } else { ## when outcome is not compositional
          delta_y <- lapply(seq(dim(ysub)[3]), function(v) {
            ysub[, , v] - y0[, h, v]
          })
        }
      } else { # when single outcome
        delta_y <- list(ysub - y0[, h])
      }
      hout[[h]] <- delta_y
    }

    ## restructure hout so that
    ## first level is the outcome
    ## second level is the reference grid
    hout <- lapply(seq_along(get_vars[["Yn"]]), function(m) {
      Map(`[[`, hout, m)
    })

    if (aorg) {
      # unadj OR adj averaging over reference grid
      weight <- grid$.wgt. / sum(grid$.wgt.)

      posterior_delta_y <- lapply(hout, function(h) {
        weighted_hout <- Map(function(x, w) x * w, h, weight)
        list(Reduce(`+`, weighted_hout))
      })
    } else {
      # adj keeping prediction at each level of at
      grid[, wgt_at := sum(.wgt.), by = names(at)]
      at_weight <- grid$wgt_at / sum(grid$wgt_at)
      at_levels <- grid[, names(at), with = FALSE]
      at_id <- at_levels[, ida := .I][, .(ida_list = list(ida)), by = names(at)]$ida_list
      unique_at_levels <- unique(at_levels[, names(at), with = FALSE])

      # for each outcome, weight the hout by at_weight
      posterior_delta_y <- lapply(hout, function(h) {
        at_weighted_hout <- Map(function(x, w) x * w, h, at_weight)
        pdy <- lapply(at_id, function(ida) {
          Reduce(`+`, at_weighted_hout[ida])
        })
        names(pdy) <- apply(unique_at_levels, 1, function(x) {
          paste(paste0(colnames(
            unique_at_levels
          ), x), collapse = "_")
        })
        pdy
      })
    }
    posterior_delta_y <- lapply(posterior_delta_y, function(x) {
      lapply(x, function(z) {
        cbind(dsub[, .(Delta, From, To, Level, Reference)], t(z))
      })
    })

    # final results for entire composition
    list(posterior_delta_y)
  }

  if (summary) {
    ## sub1 <- substitution(object = m, delta = 5, level = c("between"), at = list(Female = c(0,1)), summary = TRUE)
    ## sub2 <- substitution(object = m2, delta = 5, level = c("between"), summary = TRUE)
    iout <- lapply(iout, function(iouti) {
      iouti <- lapply(iouti, function(ioutii) {
        do.call(rbind, Map(function(d, ida) {
          dmeta <- d[, c("Delta", "From", "To", "Level", "Reference")]
          result <- apply(d[, -c("Delta", "From", "To", "Level", "Reference")], 1, posterior_summary, ...)
          row.names(result) <- c("Estimate", "Est.Error", "CI_low", "CI_high")
          if (aorg) {
            cbind(t(result), dmeta)
          } else {
            cbind(t(result), dmeta, grid[ida, names(at), with = FALSE])
          }
        }, ioutii, seq_along(ioutii)))
      })
      names(iouti) <- get_vars[["Yn"]]
      iouti
    })
  } else {
    ## sub1 <- substitution(object = m, delta = 5, level = c("between"), at = list(Female = c(0,1)), summary = FALSE)
    ## sub1 <- substitution(object = m, delta = 5, level = c("between"), summary = FALSE)
    iout <- lapply(iout, function(iouti) {
      iouti <- lapply(seq_along(iouti), function(i) {
        if (aorg) {
          iouti[[i]]
        } else {
          list(posterior = iouti[[i]], grid = as.data.table(grid[, names(at), with = FALSE]))
        }
      })
      names(iouti) <- get_vars[["Yn"]]
      iouti
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
  # extract variables from complr and brmcoda objects for use in substitution models
  get_vars <- .get.subvars(object = object, parts = parts, scale = scale)

  grid <- d0[, colnames(d0) %nin% c(get_vars[["Xxz"]], object[["complr"]][["idvar"]]), with = FALSE]
  # grid[, at := if (!is.null(at)) {
  #   names(at)
  # } else {
  #   NA
  # }]

  # setup parallel processing
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }

  oopts <- options(
    future.globals.maxSize = +Inf,
    future.globals.onReference = NULL
  )
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
      sub_to_var <- c(i, "remaining")
    } else {
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
        x0_xsub_delta_k <- setNames(x0_xsub_delta_k, c(get_vars[["XbX"]], get_vars[["sX"]], "Delta"))
        kout[[k]] <- x0_xsub_delta_k
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
    bx0 <- acomp(d1[, get_vars[["XbX"]], with = FALSE], total = object[["complr"]][["output"]][[get_vars[["idx"]]]][["total"]])
    bxsub <- acomp(d1[, get_vars[["sX"]], with = FALSE], total = object[["complr"]][["output"]][[get_vars[["idx"]]]][["total"]])

    bz0 <- ilr(bx0, V = object[["complr"]][["output"]][[get_vars[["idx"]]]][["psi"]])
    bzsub <- ilr(bxsub, V = object[["complr"]][["output"]][[get_vars[["idx"]]]][["psi"]])
    wzsub <- bzsub - bz0

    colnames(bz0) <- get_vars[["XbZ"]]
    colnames(wzsub) <- get_vars[["XwZ"]]

    # reference grid
    ## get covariate + idvar names
    covs <- colnames(d0) %snin% c(get_vars[["XbZ"]], get_vars[["XwZ"]], get_vars[["XbX"]], get_vars[["XwX"]])
    refgrid <- d0[, covs, with = FALSE]

    # predictions
    hout <- vector("list", length = nrow(d0))
    for (h in seq_len(nrow(d0))) {
      dsub <- cbind(d1, bz0, wzsub, refgrid[h, ])
      ysub <- fitted(
        object,
        newdata = dsub,
        re_formula = NA,
        scale = scale,
        summary = FALSE
      )
      if (length(get_vars[["brmcoda"]][["y"]]) > 1) {
        # when multiple outcomes
        ## when outcome is compositional
        ## delta y is acomp change
        if (all(sort(dimnames(ysub)[[3]]) %in% c(
          sort(get_vars$XX),
          sort(get_vars$XwX),
          sort(get_vars$XbX)
        ))) {
          delta_y <- lapply(seq(dim(ysub)[2]), function(j) {
            ## j = delta, h = ref grid
            acomp(ysub[, j, ], total = object[["complr"]][["output"]][[get_vars[["idy"]]]][["total"]])
            - acomp(y0[, h, ], total = object[["complr"]][["output"]][[get_vars[["idy"]]]][["total"]])
          })
        } else { ## when outcome is not compositional
          delta_y <- lapply(seq(dim(ysub)[3]), function(v) {
            ysub[, , v] - y0[, h, v]
          })
        }
      } else { # when single outcome
        delta_y <- list(ysub - y0[, h])
      }
      hout[[h]] <- delta_y
    }

    ## restructure hout so that
    ## first level is the outcome
    ## second level is the reference grid
    hout <- lapply(seq_along(get_vars[["Yn"]]), function(m) {
      Map(`[[`, hout, m)
    })

    if (aorg) {
      # unadj OR adj averaging over reference grid
      weight <- grid$.wgt. / sum(grid$.wgt.)
      
      posterior_delta_y <- lapply(hout, function(h) {
        weighted_hout <- Map(function(x, w) x * w, h, weight)
        list(Reduce(`+`, weighted_hout))
      })
    } else {
      # adj keeping prediction at each level of at
      grid[, wgt_at := sum(.wgt.), by = names(at)]
      at_weight <- grid$wgt_at / sum(grid$wgt_at)
      at_levels <- grid[, names(at), with = FALSE]
      at_id <- at_levels[, ida := .I][, .(ida_list = list(ida)), by = names(at)]$ida_list
      unique_at_levels <- unique(at_levels[, names(at), with = FALSE])
      
      # for each outcome, weight the hout by at_weight
      posterior_delta_y <- lapply(hout, function(h) {
        at_weighted_hout <- Map(function(x, w) x * w, h, at_weight)
        pdy <- lapply(at_id, function(ida) {
          Reduce(`+`, at_weighted_hout[ida])
        })
        names(pdy) <- apply(unique_at_levels, 1, function(x) {
          paste(paste0(colnames(
            unique_at_levels
          ), x), collapse = "_")
        })
        pdy
      })
    }
    posterior_delta_y <- lapply(posterior_delta_y, function(x) {
      lapply(x, function(z) {
        cbind(dsub[, .(Delta, From, To, Level, Reference)], t(z))
      })
    })

    # final results for entire composition
    list(posterior_delta_y)
  }

  if (summary) {
    ## sub1 <- substitution(object = m, delta = 5, level = c("between"), at = list(Female = c(0,1)), summary = TRUE)
    ## sub1 <- substitution(object = m, delta = 5, level = c("between"), summary = TRUE)
    iout <- lapply(iout, function(iouti) {
      iouti <- lapply(iouti, function(ioutii) {
        do.call(rbind, Map(function(d, ida) {
          dmeta <- d[, c("Delta", "From", "To", "Level", "Reference")]
          result <- apply(d[, -c("Delta", "From", "To", "Level", "Reference")], 1, posterior_summary, ...)
          row.names(result) <- c("Estimate", "Est.Error", "CI_low", "CI_high")
          if (aorg) {
            cbind(t(result), dmeta)
          } else {
            cbind(t(result), dmeta, grid[ida, names(at), with = FALSE])
          }
        }, ioutii, seq_along(ioutii)))
      })
      names(iouti) <- get_vars[["Yn"]]
      iouti
    })
  } else {
    ## sub1 <- substitution(object = m, delta = 5, level = c("between"), at = list(Female = c(0,1)), summary = FALSE)
    ## sub1 <- substitution(object = m, delta = 5, level = c("between"), summary = FALSE)
    iout <- lapply(iout, function(iouti) {
      iouti <- lapply(seq_along(iouti), function(i) {
        if (aorg) {
          iouti[[i]]
        } else {
          list(posterior = iouti[[i]], grid = as.data.table(grid[, names(at), with = FALSE]))
        }
      })
      names(iouti) <- get_vars[["Yn"]]
      iouti
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
  # extract variables from complr and brmcoda objects for use in substitution models
  get_vars <- .get.subvars(object = object, parts = parts, scale = scale)
  
  grid <- d0[, colnames(d0) %nin% c(get_vars[["Xxz"]], object[["complr"]][["idvar"]]), with = FALSE]
  # grid[, at := if (!is.null(at)) {
  #   names(at)
  # } else {
  #   NA
  # }]

  # setup parallel processing
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }

  oopts <- options(
    future.globals.maxSize = +Inf,
    future.globals.onReference = NULL
  )
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
      sub_to_var <- c(i, "remaining")
    } else {
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
        xsub_delta_k <- setNames(xsub_delta_k, c(get_vars[["sX"]], "Delta"))
        kout[[k]] <- xsub_delta_k
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
    xsub <- acomp(d1[, get_vars[["sX"]], with = FALSE], total = object[["complr"]][["output"]][[get_vars[["idx"]]]][["total"]])
    zsub <- ilr(xsub, V = object[["complr"]][["output"]][[get_vars[["idx"]]]][["psi"]])
    colnames(zsub) <- get_vars[["XZ"]]

    # reference grid
    ## get covariate + idvar names
    covs <- colnames(d0) %snin% c(get_vars[["XZ"]], get_vars$XX)
    refgrid <- d0[, covs, with = FALSE]

    # predictions
    hout <- vector("list", length = nrow(d0))
    for (h in seq_len(nrow(d0))) {
      dsub <- cbind(d1, zsub, refgrid[h, ])
      ysub <- fitted(
        object,
        newdata = dsub,
        re_formula = NA,
        scale = scale,
        summary = FALSE
      )
      if (length(get_vars[["brmcoda"]][["y"]]) > 1) {
        # when multiple outcomes
        ## when outcome is compositional
        ## delta y is acomp change
        if (all(sort(dimnames(ysub)[[3]]) %in% c(
          sort(get_vars$XX),
          sort(get_vars$XwX),
          sort(get_vars$XbX)
        ))) {
          delta_y <- lapply(seq(dim(ysub)[2]), function(j) {
            ## j = delta, h = ref grid
            acomp(ysub[, j, ], total = object[["complr"]][["output"]][[get_vars[["idy"]]]][["total"]])
            - acomp(y0[, h, ], total = object[["complr"]][["output"]][[get_vars[["idy"]]]][["total"]])
          })
        } else { ## when outcome is not compositional
          delta_y <- lapply(seq(dim(ysub)[3]), function(v) {
            ysub[, , v] - y0[, h, v]
          })
        }
      } else { # when single outcome
        delta_y <- list(ysub - y0[, h])
      }
      hout[[h]] <- delta_y
    }

    ## restructure hout so that
    ## first level is the outcome
    ## second level is the reference grid
    hout <- lapply(seq_along(get_vars[["Yn"]]), function(m) {
      Map(`[[`, hout, m)
    })

    if (aorg) {
      # unadj OR adj averaging over reference grid
      weight <- grid$.wgt. / sum(grid$.wgt.)
      
      posterior_delta_y <- lapply(hout, function(h) {
        weighted_hout <- Map(function(x, w) x * w, h, weight)
        list(Reduce(`+`, weighted_hout))
      })
    } else {
      # adj keeping prediction at each level of at
      grid[, wgt_at := sum(.wgt.), by = names(at)]
      at_weight <- grid$wgt_at / sum(grid$wgt_at)
      at_levels <- grid[, names(at), with = FALSE]
      at_id <- at_levels[, ida := .I][, .(ida_list = list(ida)), by = names(at)]$ida_list
      unique_at_levels <- unique(at_levels[, names(at), with = FALSE])
      
      # for each outcome, weight the hout by at_weight
      posterior_delta_y <- lapply(hout, function(h) {
        at_weighted_hout <- Map(function(x, w) x * w, h, at_weight)
        pdy <- lapply(at_id, function(ida) {
          Reduce(`+`, at_weighted_hout[ida])
        })
        names(pdy) <- apply(unique_at_levels, 1, function(x) {
          paste(paste0(colnames(
            unique_at_levels
          ), x), collapse = "_")
        })
        pdy
      })
    }
    posterior_delta_y <- lapply(posterior_delta_y, function(x) {
      lapply(x, function(z) {
        cbind(dsub[, .(Delta, From, To, Level, Reference)], t(z))
      })
    })

    # final results for entire composition
    list(posterior_delta_y)
  }
  if (summary) {
    iout <- lapply(iout, function(iouti) {
      iouti <- lapply(iouti, function(ioutii) {
        do.call(rbind, Map(function(d, ida) {
          dmeta <- d[, c("Delta", "From", "To", "Level", "Reference")]
          result <- apply(d[, -c("Delta", "From", "To", "Level", "Reference")], 1, posterior_summary, ...)
          row.names(result) <- c("Estimate", "Est.Error", "CI_low", "CI_high")
          if (aorg) {
            cbind(t(result), dmeta)
          } else {
            cbind(t(result), dmeta, grid[ida, names(at), with = FALSE])
          }
        }, ioutii, seq_along(ioutii)))
      })
      names(iouti) <- get_vars[["Yn"]]
      iouti
    })
  } else {
    iout <- lapply(iout, function(iouti) {
      iouti <- lapply(seq_along(iouti), function(i) {
        if (aorg) {
          iouti[[i]]
        } else {
          list(posterior = iouti[[i]], grid = as.data.table(grid[, names(at), with = FALSE]))
        }
      })
      names(iouti) <- get_vars[["Yn"]]
      iouti
    })
  }
  names(iout) <- parts
  iout
}

# Clustermean Between-person Substitution model
.get.bsubmargin <- function(object,
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
  # extract variables from complr and brmcoda objects for use in substitution models
  get_vars <- .get.subvars(object = object, parts = parts, scale = scale)

  # setup parallel processing
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }

  oopts <- options(
    future.globals.maxSize = +Inf,
    future.globals.onReference = NULL
  )
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
      sub_to_var <- c(i, "remaining")
    } else {
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
        x0_xsub_delta_k <- setNames(x0_xsub_delta_k, c(get_vars[["XbX"]], get_vars[["sX"]], "Delta"))
        d1 <- cbind(x0_xsub_delta_k, d0[, colnames(d0) %nin% get_vars[["Xxz"]], with = FALSE])

        # useful information for the final results
        d1[, From := rep(sub_from_var, length.out = nrow(sub_delta_j))[k]]
        d1[, To := sub_to_var]
        d1[, Delta := as.numeric(Delta)]
        d1[, Level := level]
        d1[, Reference := ref]

        # remove impossible reallocation that result in negative values
        cols <- colnames(d1) %sin% c(colnames(x0), colnames(base))
        d1 <- d1[rowSums(d1[, ..cols] < 0) == 0]

        # compositions and ilrs for predictions
        bx0 <- acomp(d1[, get_vars[["XbX"]], with = FALSE], total = object[["complr"]][["output"]][[get_vars[["idx"]]]][["total"]])
        bxsub <- acomp(d1[, get_vars[["sX"]], with = FALSE], total = object[["complr"]][["output"]][[get_vars[["idx"]]]][["total"]])

        bz0 <- ilr(bx0, V = object[["complr"]][["output"]][[get_vars[["idx"]]]][["psi"]])
        bzsub <- ilr(bxsub, V = object[["complr"]][["output"]][[get_vars[["idx"]]]][["psi"]])

        wz0 <- as.data.table(matrix(0, nrow = nrow(bzsub), ncol = ncol(bzsub)))

        colnames(bzsub) <- get_vars[["XbZ"]]
        colnames(wz0) <- get_vars[["XwZ"]]

        # prediction
        dsub <- cbind(d1, bzsub, wz0)
        ysub <- fitted(
          object,
          newdata = dsub,
          re_formula = NULL,
          scale = scale,
          summary = FALSE
        )
        if (length(get_vars[["brmcoda"]][["y"]]) > 1) {
          # when multiple outcomes
          ## when outcome is compositional
          ## delta y is acomp change
          if (all(sort(dimnames(ysub)[[3]]) %in% c(
            sort(get_vars$XX),
            sort(get_vars$XwX),
            sort(get_vars$XbX)
          ))) {
            delta_y <- lapply(seq(dim(ysub)[2]), function(j) {
              # cal delta by participants, then avg across participants
              acomp(ysub[, j, ], total = object[["complr"]][["output"]][[get_vars[["idy"]]]][["total"]])
              - acomp(y0[, j, ], total = object[["complr"]][["output"]][[get_vars[["idy"]]]][["total"]])
            })
            delta_y <- array(unlist(delta_y), dim = c(dim(ysub)[1], dim(ysub)[2], dim(ysub)[3]))
            delta_y <- lapply(seq(dim(delta_y)[3]), function(k) {
              rowMeans(delta_y[, , k])
            })
          } else {
            ## when outcome is not compositional
            delta_y <- lapply(seq(dim(ysub)[3]), function(v) {
              # loop outcomes, then avg across participants
              rowMeans(ysub[, , v] - y0[, , v])
            })
          }
        } else {
          delta_y <- list(rowMeans(ysub - y0))
        }
        kout[[k]] <- lapply(delta_y, function(dy) {
          cbind.data.frame(unique(d1[, c("Delta", "From", "To", "Level", "Reference")]), t(dy))
        })
      }

      ## restructure so that
      ## first level is the outcome
      ## second level is delta
      kout <- lapply(seq_along(get_vars[["Yn"]]), function(m) {
        rbindlist(lapply(kout, `[[`, m))
      })
      jout[[j]] <- kout
    }
    jout <- lapply(seq_along(get_vars[["Yn"]]), function(m) {
      rbindlist(lapply(jout, `[[`, m))
    })
    list(jout)
  }

  if (summary) {
    iout <- lapply(iout, function(iouti) {
      iouti <- lapply(iouti, function(ioutii) {
        dmeta <- unique(ioutii[, c("Delta", "From", "To", "Level", "Reference")])
        result <- apply(ioutii[, -c("Delta", "From", "To", "Level", "Reference")], 1, posterior_summary, ...)
        row.names(result) <- c("Estimate", "Est.Error", "CI_low", "CI_high")
        cbind(t(result), dmeta)
      })
      names(iouti) <- get_vars[["Yn"]]
      iouti
    })
  } else {
    iout <- lapply(iout, function(iouti) {
      iouti <- lapply(seq_along(iouti), function(i) {
        iouti[[i]]
      })
      names(iouti) <- get_vars[["Yn"]]
      iouti
    })
  }
  names(iout) <- parts
  iout
}

# Clustermean Within-person Substitution model
.get.wsubmargin <- function(object,
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
  # extract variables from complr and brmcoda objects for use in substitution models
  get_vars <- .get.subvars(object = object, parts = parts, scale = scale)

  # setup parallel processing
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }

  oopts <- options(
    future.globals.maxSize = +Inf,
    future.globals.onReference = NULL
  )
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
      sub_to_var <- c(i, "remaining")
    } else {
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
        x0_xsub_delta_k <- setNames(x0_xsub_delta_k, c(get_vars[["XbX"]], get_vars[["sX"]], "Delta"))
        d1 <- cbind(x0_xsub_delta_k, d0[, colnames(d0) %nin% get_vars[["Xxz"]], with = FALSE])

        # useful information for the final results
        d1[, From := rep(sub_from_var, length.out = nrow(sub_delta_j))[k]]
        d1[, To := sub_to_var]
        d1[, Delta := as.numeric(Delta)]
        d1[, Level := level]
        d1[, Reference := ref]

        # remove impossible reallocation that result in negative values
        cols <- colnames(d1) %sin% c(colnames(x0), colnames(base))
        d1 <- d1[rowSums(d1[, ..cols] < 0) == 0]

        # compositions and ilr for predictions
        bx0 <- acomp(d1[, get_vars[["XbX"]], with = FALSE], total = object[["complr"]][["output"]][[get_vars[["idx"]]]][["total"]])
        bxsub <- acomp(d1[, get_vars[["sX"]], with = FALSE], total = object[["complr"]][["output"]][[get_vars[["idx"]]]][["total"]])

        bz0 <- ilr(bx0, V = object[["complr"]][["output"]][[get_vars[["idx"]]]][["psi"]])
        bzsub <- ilr(bxsub, V = object[["complr"]][["output"]][[get_vars[["idx"]]]][["psi"]])

        wzsub <- bzsub - bz0

        colnames(bz0) <- get_vars[["XbZ"]]
        colnames(wzsub) <- get_vars[["XwZ"]]

        # substitution data
        dsub <- cbind(d1, bz0, wzsub)

        # prediction
        ysub <- fitted(
          object,
          newdata = dsub,
          re_formula = NULL,
          scale = scale,
          summary = FALSE
        )
        if (length(get_vars[["brmcoda"]][["y"]]) > 1) {
          # when multiple outcomes
          ## when outcome is compositional
          ## delta y is acomp change
          if (all(sort(dimnames(ysub)[[3]]) %in% c(
            sort(get_vars$XX),
            sort(get_vars$XwX),
            sort(get_vars$XbX)
          ))) {
            delta_y <- lapply(seq(dim(ysub)[2]), function(j) {
              # cal delta by participants, then avg across participants
              acomp(ysub[, j, ], total = object[["complr"]][["output"]][[get_vars[["idy"]]]][["total"]])
              - acomp(y0[, j, ], total = object[["complr"]][["output"]][[get_vars[["idy"]]]][["total"]])
            })
            delta_y <- array(unlist(delta_y), dim = c(dim(ysub)[1], dim(ysub)[2], dim(ysub)[3]))
            delta_y <- lapply(seq(dim(delta_y)[3]), function(k) {
              rowMeans(delta_y[, , k])
            })
          } else {
            ## when outcome is not compositional
            delta_y <- lapply(seq(dim(ysub)[3]), function(v) {
              # loop outcomes, then avg across participants
              rowMeans(ysub[, , v] - y0[, , v])
            })
          }
        } else {
          delta_y <- list(rowMeans(ysub - y0))
        }
        kout[[k]] <- lapply(delta_y, function(dy) {
          cbind.data.frame(unique(d1[, c("Delta", "From", "To", "Level", "Reference")]), t(dy))
        })
      }

      ## restructure so that
      ## first level is the outcome
      ## second level is delta
      kout <- lapply(seq_along(get_vars[["Yn"]]), function(m) {
        rbindlist(lapply(kout, `[[`, m))
      })
      jout[[j]] <- kout
    }
    jout <- lapply(seq_along(get_vars[["Yn"]]), function(m) {
      rbindlist(lapply(jout, `[[`, m))
    })
    list(jout)
  }

  if (summary) {
    iout <- lapply(iout, function(iouti) {
      iouti <- lapply(iouti, function(ioutii) {
        dmeta <- unique(ioutii[, c("Delta", "From", "To", "Level", "Reference")])
        result <- apply(ioutii[, -c("Delta", "From", "To", "Level", "Reference")], 1, posterior_summary, ...)
        row.names(result) <- c("Estimate", "Est.Error", "CI_low", "CI_high")
        cbind(t(result), dmeta)
      })
      names(iouti) <- get_vars[["Yn"]]
      iouti
    })
  } else {
    iout <- lapply(iout, function(iouti) {
      iouti <- lapply(seq_along(iouti), function(i) {
        iouti[[i]]
      })
      names(iouti) <- get_vars[["Yn"]]
      iouti
    })
  }
  names(iout) <- parts
  iout
}

# Clustermean Average Substitution
.get.submargin <- function(object,
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
  # extract variables from complr and brmcoda objects for use in substitution models
  get_vars <- .get.subvars(object = object, parts = parts, scale = scale)

  # setup parallel processing
  if (isFALSE(is.null(cores))) {
    oplan <- plan(multisession, workers = cores)
    on.exit(plan(oplan))
  } else {
    plan(sequential)
  }

  oopts <- options(
    future.globals.maxSize = +Inf,
    future.globals.onReference = NULL
  )
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
      sub_to_var <- c(i, "remaining")
    } else {
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

        xsub_delta_k <- cbind(xsub, sub_delta_k[, get(i)])
        xsub_delta_k <- setNames(xsub_delta_k, c(get_vars[["sX"]], "Delta"))
        d1 <- cbind(xsub_delta_k, d0[, colnames(d0) %nin% get_vars[["Xxz"]], with = FALSE])


        # useful information for the final results
        d1[, From := rep(sub_from_var, length.out = nrow(sub_delta_j))[k]]
        d1[, To := sub_to_var]
        d1[, Delta := as.numeric(Delta)]
        d1[, Level := level]
        d1[, Reference := ref]

        # remove impossible reallocation that result in negative values
        cols <- colnames(d1) %sin% c(colnames(x0), colnames(base))
        d1 <- d1[rowSums(d1[, ..cols] < 0) == 0]

        # compositions and ilrs for predictions
        xsub <- acomp(d1[, get_vars[["sX"]], with = FALSE], total = object[["complr"]][["output"]][[get_vars[["idx"]]]][["total"]])
        zsub <- ilr(xsub, V = object[["complr"]][["output"]][[get_vars[["idx"]]]][["psi"]])

        colnames(zsub) <- get_vars[["XZ"]]

        # substitution data
        dsub <- cbind(d1, zsub)

        # prediction
        ysub <- fitted(
          object,
          newdata = dsub,
          re_formula = NULL,
          scale = scale,
          summary = FALSE
        )
        if (length(get_vars[["brmcoda"]][["y"]]) > 1) {
          # when multiple outcomes
          ## when outcome is compositional
          ## delta y is acomp change
          if (all(sort(dimnames(ysub)[[3]]) %in% c(
            sort(get_vars$XX),
            sort(get_vars$XwX),
            sort(get_vars$XbX)
          ))) {
            delta_y <- lapply(seq(dim(ysub)[2]), function(j) {
              # cal delta by participants, then avg across participants
              acomp(ysub[, j, ], total = object[["complr"]][["output"]][[get_vars[["idy"]]]][["total"]])
              - acomp(y0[, j, ], total = object[["complr"]][["output"]][[get_vars[["idy"]]]][["total"]])
            })
            delta_y <- array(unlist(delta_y), dim = c(dim(ysub)[1], dim(ysub)[2], dim(ysub)[3]))
            delta_y <- lapply(seq(dim(delta_y)[3]), function(k) {
              rowMeans(delta_y[, , k])
            })
          } else {
            ## when outcome is not compositional
            delta_y <- lapply(seq(dim(ysub)[3]), function(v) {
              # loop outcomes, then avg across participants
              rowMeans(ysub[, , v] - y0[, , v])
            })
          }
        } else {
          delta_y <- list(rowMeans(ysub - y0))
        }
        kout[[k]] <- lapply(delta_y, function(dy) {
          cbind.data.frame(unique(d1[, c("Delta", "From", "To", "Level", "Reference")]), t(dy))
        })
      }

      ## restructure so that
      ## first level is the outcome
      ## second level is delta
      kout <- lapply(seq_along(get_vars[["Yn"]]), function(m) {
        rbindlist(lapply(kout, `[[`, m))
      })
      jout[[j]] <- kout
    }
    jout <- lapply(seq_along(get_vars[["Yn"]]), function(m) {
      rbindlist(lapply(jout, `[[`, m))
    })
    list(jout)
  }

  if (summary) {
    iout <- lapply(iout, function(iouti) {
      iouti <- lapply(iouti, function(ioutii) {
        dmeta <- unique(ioutii[, c("Delta", "From", "To", "Level", "Reference")])
        result <- apply(ioutii[, -c("Delta", "From", "To", "Level", "Reference")], 1, posterior_summary, ...)
        row.names(result) <- c("Estimate", "Est.Error", "CI_low", "CI_high")
        cbind(t(result), dmeta)
      })
      names(iouti) <- get_vars[["Yn"]]
      iouti
    })
  } else {
    iout <- lapply(iout, function(iouti) {
      iouti <- lapply(seq_along(iouti), function(i) {
        iouti[[i]]
      })
      names(iouti) <- get_vars[["Yn"]]
      iouti
    })
  }
  names(iout) <- parts
  iout
}
