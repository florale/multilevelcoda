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
      brmsformula = object$model$formula,
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
#' @importFrom emmeans ref_grid
#' @importFrom extraoperators %snin% %sin%
#'
#' @return A reference grid consisting of a combination of covariates in \code{brmcoda}
#'
#' @export
build.rg <- function(object,
                     ref,
                     level,
                     weight,
                     fill = FALSE) {

  covgrid <- NULL

  ## NOTES
  # # ignore weight for clustermean
  # equal weight is default for grandmean

  # what type of model is being estimated
  model_fixef <- rownames(fixef(object))
  model_ranef <- if(dim(object$model$ranef)[1] > 0) (names(ranef(object))) else (NULL)

  model_fixef_level <- NULL
  model_fixef_coef <- NULL

  if (length(grep("bilr", model_fixef, value = T)) > 0) {
    model_fixef_level <- append(model_fixef_level, "between")
    model_fixef_coef  <- append(model_fixef_coef, grep(".*bilr", model_fixef, value = T))
  }
  if (length(grep("wilr", model_fixef, value = T)) > 0) {
    model_fixef_level <- append(model_fixef_level, "within")
    model_fixef_coef  <- append(model_fixef_coef, grep(".*wilr", model_fixef, value = T))
  }
  if ((length(grep("ilr", model_fixef, value = T)) > 0) && (length(grep("[b|w]ilr", model_fixef, value = T)) == 0)) {
    model_fixef_level <- append(model_fixef_level, "aggregate")
    model_fixef_coef  <- append(model_fixef_coef, grep(paste0(names(object$complr$logratio), collapse = "|"), model_fixef, value = T))
  }

  # single level or multilevel
  if (length(model_ranef) > 0) {
    model_ranef_level <- "multilevel"
    model_ranef_coef <- model_ranef
  } else {
    model_ranef_level <- "single"
    model_ranef_coef <- NULL
  }

  # d0 and comp0 for multilevel level model
  if (model_ranef_level == "multilevel") {

    # for clustermean
    if ("clustermean" %in% ref) {
      weight <- NULL # ignore weight

      # aggregate
      if ("aggregate" %in% level) {
        d0 <- cbind(object$complr$comp, object$complr$data[, -colnames(object$complr$comp), with = FALSE])
        d0 <- d0[, head(.SD, 1), by = eval(object$complr$idvar)]
        comp0 <- acomp(d0[, colnames(object$complr$comp), with = FALSE], total = object$complr$total)

        ilr0 <- ilr(comp0, V = object$complr$psi)
        ilr0 <- as.data.table(ilr0)

        colnames(ilr0)  <- colnames(object$complr$logratio)
        colnames(comp0) <- colnames(object$complr$comp)

        d0 <- cbind(ilr0, comp0,
                    d0[, -colnames(comp0), with = FALSE])
      }

      # between and within
      if (any(c("between", "within") %in% level)) {
        d0 <- cbind(object$complr$between_comp, object$complr$data)
        d0 <- d0[, head(.SD, 1), by = eval(object$complr$idvar)]
        comp0 <- acomp(d0[, colnames(object$complr$between_comp), with = FALSE], total = object$complr$total)
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
    } else {

      # assemble reference grid
      # get var names
      ilrnames <- c(colnames(object$complr$between_logratio),
                    colnames(object$complr$within_logratio),
                    colnames(object$complr$logratio))

      vars  <- colnames(model.frame(object))
      resp  <- object$model$formula$formula[[2]]
      grp   <- object$model$ranef$group
      preds <- vars %snin% c(resp, grp)
      covs  <- vars %snin% c(resp, grp, ilrnames)

      # default reference grid
      # set binary
      # drg <- model.frame(object)
      # drg[] <- as.data.table(lapply(drg, function(j) if(is.numeric(j) && unique(j) %ain% c(0, 1)) as.factor(j) else j))
      # drg <- as.data.table(insight::get_datagrid(drg,
      #                                            by = paste0(resp),
      #                                            factors = factors,
      #                                            length = NA))
      
      drg <- as.data.table(ref_grid(object$model)@grid)

      # reference grid (only covariates and outcome)
      refgrid <- drg[, colnames(drg) %in% c(covs, resp), with = FALSE]

      id <- data.table::data.table(1) # to make fitted() happy
      colnames(id) <- object$complr$idvar

      # grandmean
      if ("grandmean" %in% ref) {

        # aggregate
        if ("aggregate" %in% level) {
          if (weight == "proportional") {
            comp0 <- mean.acomp(object$complr$comp, robust = TRUE)

          } else {
            comp0 <- cbind(object$complr$comp,
                           object$complr$data[, object$complr$idvar, with = FALSE])
            comp0 <- comp0[, head(.SD, 1), by = eval(object$complr$idvar)]
            comp0 <- acomp(comp0[, colnames(object$complr$comp), with = FALSE], total = object$complr$total)
            comp0 <- mean.acomp(comp0, robust = TRUE)
          }

          comp0 <- acomp(comp0, total = object$complr$total)
          comp0 <- as.data.table(t(comp0))

          ilr0 <- ilr(comp0, V = object$complr$psi)
          ilr0 <- as.data.table(t(ilr0))

          colnames(ilr0)  <- colnames(object$complr$logratio)
          colnames(comp0) <- colnames(object$complr$comp)

          d0 <- if (all(dim(refgrid) == 0)) (cbind(ilr0, comp0, id)) else (expand.grid.df(ilr0, comp0, id, refgrid))
        }

        # between and/or within
        if (any(c("between", "within") %in% level)) {
          if (weight == "proportional") {
            comp0 <- mean.acomp(object$complr$between_comp, robust = TRUE)

          } else {
            comp0 <- cbind(object$complr$between_comp,
                           object$complr$data[, object$complr$idvar, with = FALSE])
            comp0 <- comp0[, head(.SD, 1), by = eval(object$complr$idvar)]
            comp0 <- acomp(comp0[, colnames(object$complr$between_comp), with = FALSE], total = object$complr$total)
            comp0 <- mean.acomp(comp0, robust = TRUE)
          }

          comp0 <- acomp(comp0, total = object$complr$total)
          comp0 <- as.data.table(t(comp0))
          bcomp0 <- comp0

          bilr0 <- ilr(bcomp0, V = object$complr$psi)
          bilr0 <- as.data.table(t(bilr0))

          wcomp0 <- as.data.table(matrix(1, nrow = nrow(bcomp0), ncol = ncol(bcomp0)))
          wilr0 <- as.data.table(matrix(0, nrow = nrow(bilr0), ncol = ncol(bilr0)))

          colnames(bilr0)  <- colnames(object$complr$between_logratio)
          colnames(wilr0)  <- colnames(object$complr$within_logratio)
          colnames(bcomp0) <- colnames(object$complr$between_comp)
          colnames(wcomp0) <- colnames(object$complr$within_comp)

          d0 <- if (all(dim(refgrid) == 0)) (cbind(bilr0, wilr0, bcomp0, wcomp0, id)) else (expand.grid.df(bilr0, wilr0, bcomp0, wcomp0, id, refgrid))
        }
      }

      if (inherits(ref, c("data.table", "data.frame", "matrix"))) {
        weight <- NULL

        if (isFALSE(object$complr$parts %in% colnames(ref))) {  # get user's composition
          stop(
            sprintf(
              "The reference grid should include all compositional components but (%s) are missing.",
              paste0(object$complr$parts %nin% colnames(ref), collapse = ", ")
            ))
        } else {
          comp_user <- ref[, object$complr$parts, with = FALSE]
          comp_user <- acomp(comp_user, total = object$complr$total)
          comp_user <- as.data.table(t(comp_user))
        }

        # sanity checks
        if (nrow(ref) > 1) {
          stop("Only one reference composition is allowed at a time.")
        }
        if(isFALSE(sum(comp_user) == object$complr$total)) {
          stop(sprintf(
            "The total amount of the reference composition (%s) should be the same as the composition (%s).",
            sum(comp_user),
            object$complr$total
          ))
        }
        if (isTRUE((any(comp_user > lapply(object$complr$data[, object$complr$parts, with = FALSE], max)) |
                    any(comp_user < lapply(object$complr$data[, object$complr$parts, with = FALSE], min))))) {
          stop(paste(
            sprintf(
              "composition should be numeric or interger values that are between (%s) and (%s)",
              paste0(round(apply(object$complr$data[, object$complr$parts, with = FALSE], 2, min)), collapse = ", "),
              paste0(round(apply(object$complr$data[, object$complr$parts, with = FALSE], 2, max)), collapse = ", ")),
            "\n",
            " for",
            paste0(object$complr$parts, collapse = ", "),
            "respectively"
          ))
        }

        # user's specified reference grid
        ## any covariates left in the ref
        if (ncol(ref) > ncol(comp_user)) {
          covgrid <- ref[, -object$complr$parts, with = FALSE]

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
          comp0 <- comp_user

          ilr0 <- ilr(comp0, V = object$complr$psi)
          ilr0 <- as.data.table(t(ilr0))

          colnames(ilr0)  <- colnames(object$complr$logratio)
          colnames(comp0) <- colnames(object$complr$comp)

          d0 <- if (all(dim(refgrid) == 0)) (cbind(ilr0, comp0, id)) else (expand.grid.df(ilr0, comp0, id, refgrid))

        }
        if (level %in% c("between", "within")) {
          comp0 <- cbind(object$complr$between_comp,
                         object$complr$data[, object$complr$idvar, with = FALSE])
          comp0 <- comp0[, head(.SD, 1), by = eval(object$complr$idvar)]
          comp0 <- acomp(comp0[, colnames(object$complr$between_comp), with = FALSE], total = object$complr$total)
          comp0 <- mean.acomp(comp0, robust = TRUE)

          # assemble d0
          # bilr0 is between-person ilr of the ref comp (doesn't have to be compositional mean)
          bcomp0 <- comp_user
          bilr0 <- ilr(bcomp0, V = object$complr$psi)
          bilr0 <- as.data.table(t(bilr0))

          # wcomp0 and wilr0 are the difference between the actual compositional mean of the dataset and bilr
          # is 0 if ref comp is compositional mean
          # but is different if not
          wcomp0 <- bcomp0 - comp0
          wilr0 <- as.data.table(t(ilr(wcomp0, V = object$complr$psi)))

          id <- data.table::data.table(1) # to make fitted() happy

          colnames(bilr0)  <- colnames(object$complr$between_logratio)
          colnames(wilr0)  <- colnames(object$complr$within_logratio)
          colnames(bcomp0) <- colnames(object$complr$between_comp)
          colnames(wcomp0) <- colnames(object$complr$within_comp)
          colnames(id)     <- object$complr$idvar

          d0 <- if (all(dim(refgrid) == 0)) (cbind(bilr0, wilr0, bcomp0, wcomp0, id)) else (expand.grid.df(bilr0, wilr0, bcomp0, wcomp0, id, refgrid))
        }
      }
    }
  }

  ## d0 and comp0 for single level model
  if (model_ranef_level == "single") {

    comp0 <- object$complr$comp
    comp0 <- mean.acomp(comp0, robust = TRUE)
    comp0 <- acomp(comp0, total = object$complr$total)
    comp0 <- as.data.table(t(comp0))
    bcomp0 <- comp0

    d0 <- object$complr$data

    ilr0 <- ilr(comp0, V = object$complr$psi)
    ilr0 <- as.data.table(t(ilr0))

    colnames(ilr0)  <- colnames(object$complr$logratio)
    colnames(comp0) <- colnames(object$complr$comp)

    # assemble reference grid
    # get var names
    ilrnames <- c(colnames(object$complr$between_logratio),
                  colnames(object$complr$within_logratio),
                  colnames(object$complr$logratio))

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
    
    drg <- as.data.table(ref_grid(object$model)@grid)

    # reference grid (only covariates and outcome)
    refgrid <- drg[, colnames(drg) %in% c(covs, resp), with = FALSE]

    d0 <- if (all(dim(refgrid) == 0)) (cbind(ilr0, comp0)) else (expand.grid.df(ilr0, comp0, refgrid))
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
.get.bsub <- function(object, delta, basesub,
                      comp0, y0, d0,
                      summary,
                      level, ref, scale,
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

  iout <- foreach(i = colnames(basesub), .combine = c,
                  .options.future = list(packages = "multilevelcoda")) %dofuture% {

                    # dnew - reallocation data
                    # possible substitution of 1 compositional variable
                    basesub_tmp <- as.data.table(basesub)
                    basesub_tmp <- basesub_tmp[(get(i) != 0)]
                    basesub_tmp <- basesub_tmp[order(-rank(get(i)))]

                    # substitution variable names
                    subvar <- colnames(basesub_tmp) %snin% eval(i)
                    iv <- i

                    kout <- vector("list", length = nrow(basesub_tmp))
                    jout <- vector("list", length = length(delta))

                    for (j in seq_along(delta)) { # delta level
                      sub_tmp_j <- basesub_tmp * delta[j]
                      for (k in seq_len(nrow(sub_tmp_j))) {
                        comp1 <- comp0 + sub_tmp_j[k,]
                        names(comp1) <- object$complr$parts
                        Delta <- sub_tmp_j[k, get(i)]
                        kout[[k]] <- cbind(comp0, comp1, Delta)
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
                    covs <- colnames(d0) %snin% c(colnames(object$complr$between_logratio),
                                                  colnames(object$complr$within_logratio),
                                                  colnames(object$complr$between_comp),
                                                  colnames(object$complr$within_comp)
                    )
                    refgrid <- d0[, covs, with = FALSE]

                    # predictions
                    hout <- vector("list", length = nrow(d0))
                    if (summary) { # unadj OR adj averaging over reference grid
                      for (h in seq_len(nrow(d0))) {
                        dsub <- cbind(dnew, bilrsub, wilr0, refgrid[h, ])
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

                      suppressWarnings(posterior_delta_y <- apply(delta_y_avg, 2, function(x) describe_posterior(x, centrality = "mean", ...)))
                      posterior_delta_y <- rbindlist(posterior_delta_y)
                      posterior_delta_y <- cbind(posterior_delta_y[, .(Mean, CI_low, CI_high)],
                                                 dsub[, .(Delta, From, To, Level, Reference)])

                    } else { # adj keeping prediction at each level of reference grid
                      for (h in seq_len(nrow(d0))) {
                        dsub <- cbind(dnew, bilrsub, wilr0, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            scale = scale,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        suppressWarnings(
                          posterior_delta_y <- apply(delta_y, 2, function(x) {
                            describe_posterior(x, centrality = "mean", ...)
                          }))
                        posterior_delta_y <- rbindlist(posterior_delta_y)
                        posterior_delta_y <- cbind(posterior_delta_y[, .(Mean, CI_low, CI_high)],
                                                   dsub[, .(Delta, From, To, Level, Reference)],
                                                   dsub[, colnames(refgrid) %snin% object$complr$idvar, with = FALSE])

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
.get.wsub <- function(object, delta, basesub,
                      comp0, y0, d0,
                      summary,
                      level, ref, scale,
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

  iout <- foreach(i = colnames(basesub), .combine = c,
                  .options.future = list(packages = "multilevelcoda")) %dofuture% {

                    # possible susbstituion of 1 compositional variable
                    basesub_tmp <- as.data.table(basesub)
                    basesub_tmp <- basesub_tmp[(get(i) != 0)]
                    basesub_tmp <- basesub_tmp[order(-rank(get(i)))]

                    # substitution variable names
                    subvar <- colnames(basesub_tmp) %snin% eval(i)
                    iv <- i

                    kout <- vector("list", length = nrow(basesub_tmp))
                    jout <- vector("list", length = length(delta))

                    for (j in seq_along(delta)) { # delta level
                      sub_tmp_j <- basesub_tmp * delta[j]
                      for (k in seq_len(nrow(sub_tmp_j))) {
                        comp1 <- comp0 + sub_tmp_j[k, ]
                        names(comp1) <- object$complr$parts
                        Delta <- sub_tmp_j[k, get(i)]
                        kout[[k]] <- cbind(comp0, comp1, Delta)
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
                    covs <- colnames(d0) %snin% c(colnames(object$complr$between_logratio),
                                                  colnames(object$complr$within_logratio),
                                                  colnames(object$complr$between_comp),
                                                  colnames(object$complr$within_comp)
                    )
                    refgrid <- d0[, covs, with = FALSE]

                    # predictions
                    hout <- vector("list", length = nrow(d0))
                    if (summary) { # unadj OR adj averaging over reference grid
                      for (h in seq_len(nrow(d0))) {
                        dsub <- cbind(dnew, bilr0, wilrsub, refgrid[h, ])
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

                      suppressWarnings(posterior_delta_y <- apply(delta_y_avg, 2, function(x) describe_posterior(x, centrality = "mean", ...)))
                      posterior_delta_y <- rbindlist(posterior_delta_y)
                      posterior_delta_y <- cbind(posterior_delta_y[, .(Mean, CI_low, CI_high)],
                                                 dsub[, .(Delta, From, To, Level, Reference)])

                    } else { # adj keeping prediction at each level of reference grid
                      for (h in seq_len(nrow(d0))) {
                        dsub <- cbind(dnew, bilr0, wilrsub, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            scale = scale,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        suppressWarnings(
                          posterior_delta_y <- apply(delta_y, 2, function(x) {
                            describe_posterior(x, centrality = "mean", ...)
                          }))
                        posterior_delta_y <- rbindlist(posterior_delta_y)
                        posterior_delta_y <- cbind(posterior_delta_y[, .(Mean, CI_low, CI_high)],
                                                   dsub[, .(Delta, From, To, Level, Reference)],
                                                   dsub[, colnames(refgrid) %snin% object$complr$idvar, with = FALSE])

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
.get.sub <- function(object, delta, basesub,
                     comp0, y0, d0,
                     summary,
                     level, ref, scale,
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

  iout <- foreach(i = colnames(basesub), .combine = c,
                  .options.future = list(packages = "multilevelcoda")) %dofuture% {

                    # dnew - reallocation data
                    # possible substitution of 1 compositional variable
                    basesub_tmp <- as.data.table(basesub)
                    basesub_tmp <- basesub_tmp[(get(i) != 0)]
                    basesub_tmp <- basesub_tmp[order(-rank(get(i)))]

                    # substitution variable names
                    subvar <- colnames(basesub_tmp) %snin% eval(i)
                    iv <- i

                    kout <- vector("list", length = nrow(basesub_tmp))
                    jout <- vector("list", length = length(delta))

                    for (j in seq_along(delta)) { # delta level
                      sub_tmp_j <- basesub_tmp * delta[j]
                      for (k in seq_len(nrow(sub_tmp_j))) {
                        comp1 <- comp0 + sub_tmp_j[k,]
                        names(comp1) <- object$complr$parts
                        Delta <- sub_tmp_j[k, get(i)]
                        kout[[k]] <- cbind(comp1, Delta)
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
                    compsub  <- acomp(dnew[, object$complr$parts, with = FALSE], total = object$complr$total)

                    ilrsub <- ilr(compsub, V = object$complr$psi)
                    colnames(ilrsub) <- colnames(object$complr$logratio)

                    # reference grid
                    ## get covariate + idvar names
                    covs <- colnames(d0) %snin% c(colnames(object$complr$logratio),
                                                  colnames(object$complr$comp)
                    )
                    refgrid <- d0[, covs, with = FALSE]

                    # predictions
                    hout <- vector("list", length = nrow(d0))
                    if (summary) { # unadj OR adj averaging over reference grid
                      for (h in seq_len(nrow(d0))) {
                        dsub <- cbind(dnew, ilrsub, refgrid[h, ])
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

                      suppressWarnings(posterior_delta_y <- apply(delta_y_avg, 2, function(x) describe_posterior(x, centrality = "mean", ...)))
                      posterior_delta_y <- rbindlist(posterior_delta_y)
                      posterior_delta_y <- cbind(posterior_delta_y[, .(Mean, CI_low, CI_high)],
                                                 dsub[, .(Delta, From, To, Level, Reference)])

                    } else { # adj keeping prediction at each level of reference grid
                      for (h in seq_len(nrow(d0))) {
                        dsub <- cbind(dnew, ilrsub, refgrid[h, ])
                        ysub <-
                          fitted(
                            object,
                            newdata = dsub,
                            re_formula = NA,
                            scale = scale,
                            summary = FALSE
                          )
                        delta_y <- ysub - y0[, h]
                        suppressWarnings(
                          posterior_delta_y <- apply(delta_y, 2, function(x) {
                            describe_posterior(x, centrality = "mean", ...)
                          }))
                        posterior_delta_y <- rbindlist(posterior_delta_y)
                        posterior_delta_y <- cbind(posterior_delta_y[, .(Mean, CI_low, CI_high)],
                                                   dsub[, .(Delta, From, To, Level, Reference)],
                                                   dsub[, colnames(refgrid) %snin% object$complr$idvar, with = FALSE])

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
.get.bsubmargins <- function(object, delta, basesub,
                             comp0, y0, d0,
                             level, ref, scale,
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

  iout <- foreach(i = colnames(basesub), .combine = c,
                  .options.future = list(packages = "multilevelcoda")) %dofuture% {

                    # possible susbstituion of 1 compositional variable
                    basesub_tmp <- as.data.table(basesub)
                    basesub_tmp <- basesub_tmp[(get(i) != 0)]
                    basesub_tmp <- basesub_tmp[order(-rank(get(i)))]

                    # substitution variable names
                    subvar <- colnames(basesub_tmp) %snin% eval(i)
                    iv <- i

                    kout <- vector("list", length = nrow(basesub_tmp))
                    jout <- vector("list", length = length(delta))

                    for (j in seq_along(delta)) { # delta level
                      sub_tmp_j <- basesub_tmp * delta[j]
                      for (k in seq_len(nrow(sub_tmp_j))) { # reallocation level
                        sub_tmp_k <- sub_tmp_j[k, ]
                        sub_tmp_k <- sub_tmp_k[rep(seq_len(nrow(sub_tmp_k)), nrow(comp0)), ]
                        comp1 <- comp0 + sub_tmp_k
                        names(comp1) <- object$complr$parts

                        Delta <- sub_tmp_k[, get(i)]

                        dnew <- cbind(comp0, comp1,
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
                             level, ref, scale,
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

  iout <- foreach(i = colnames(basesub), .combine = c,
                  .options.future = list(packages = "multilevelcoda")) %dofuture% {

                    basesub_tmp <- as.data.table(basesub)
                    basesub_tmp <- basesub_tmp[(get(i) != 0)]
                    basesub_tmp <- basesub_tmp[order(-rank(get(i)))]

                    # substitution variable names
                    subvar <- colnames(basesub_tmp) %snin% eval(i)
                    iv <- i

                    kout <- vector("list", length = nrow(basesub_tmp))
                    jout <- vector("list", length = length(delta))

                    for (j in seq_along(delta)) { # delta level
                      sub_tmp_j <- basesub_tmp * delta[j]
                      for (k in seq_len(nrow(sub_tmp_j))) {
                        sub_tmp_k <- sub_tmp_j[k, ]
                        sub_tmp_k <- sub_tmp_k[rep(seq_len(nrow(sub_tmp_k)), nrow(comp0)), ]
                        comp1 <- comp0 + sub_tmp_k
                        names(comp1) <- object$complr$parts

                        Delta <- sub_tmp_k[, get(i)]

                        dnew <- cbind(comp0, comp1,
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
                            level, ref, scale,
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

  iout <- foreach(i = colnames(basesub), .combine = c,
                  .options.future = list(packages = "multilevelcoda")) %dofuture% {

                    # possible susbstituion of 1 compositional variable
                    basesub_tmp <- as.data.table(basesub)
                    basesub_tmp <- basesub_tmp[(get(i) != 0)]
                    basesub_tmp <- basesub_tmp[order(-rank(get(i)))]

                    # substitution variable names
                    subvar <- colnames(basesub_tmp) %snin% eval(i)
                    iv <- i

                    kout <- vector("list", length = nrow(basesub_tmp))
                    jout <- vector("list", length = length(delta))

                    for (j in seq_along(delta)) { # delta level
                      sub_tmp_j <- basesub_tmp * delta[j]
                      for (k in seq_len(nrow(sub_tmp_j))) {
                        sub_tmp_k <- sub_tmp_j[k, ]
                        sub_tmp_k <- sub_tmp_k[rep(seq_len(nrow(sub_tmp_k)), nrow(comp0)), ]
                        comp1 <- comp0 + sub_tmp_k
                        names(comp1) <- object$complr$parts

                        Delta <- sub_tmp_k[, get(i)]

                        dnew <- cbind(comp1,
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
                        tcomp <- acomp(dnew[, object$complr$parts, with = FALSE], total = object$complr$total)
                        tilr <- ilr(tcomp, V = object$complr$psi)

                        colnames(tilr) <- colnames(object$complr$logratio)

                        # substitution data
                        dsub <- cbind(dnew, tilr)

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
