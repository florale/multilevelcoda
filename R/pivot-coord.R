#' Estimate pivot balance coordinates
#' @param object An object of class \code{brmcoda}.
#' @param method A character string.
#' Should the pivot balance coordinates be estimated by \code{"rotate"} the sequential binary partition 
#' using the same \code{brmcoda} object or \code{"refit"} the \code{brmcoda} object?
#' Default is \code{"rotate"}.
#' @param summary Should summary statistics be returned instead of the raw values? Default is \code{TRUE}.
#' @param ... currently ignored.
#' 
#' @return A list of \code{\link{brmcoda}} for each pivot balance coordinate.
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   cilr <- complr(data = mcompd, sbp = sbp,
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID",
#'                  total = 1440)
#'   
#'   # inspects ILRs before passing to brmcoda
#'   names(cilr$between_logratio)
#'   names(cilr$within_logratio)
#'   names(cilr$logratio)
#'   
#'   # model with compositional predictor at between and within-person levels
#'   m <- brmcoda(complr = cilr,
#'                 formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                                    wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   m_pivot_coord <- pivot_coord(m)
#'   summary(m_pivot_coord)
#'   }}
#' @export
pivot_coord <- function (object, summary = TRUE, 
                         method = c("rotate", "refit"),
                         ...) {
  
  if (all(c("rotate", "refit") %in% method)) {
    method <- "rotate"
  }
  if (method == "rotate") {
    out <- pivot_coord_rotate(object = object,
                              summary = summary,
                              ...)
  }
  if (method == "refit") {
    out <- pivot_coord_refit(object = object,
                             ...)
  }
  structure(out, class = "pivot_coord")
  out
}

#' Estimate pivot balance coordinates by rotating sequential binary partition.
#' 
#' @param object An object of class \code{brmcoda}.
#' @param summary Should summary statistics be returned instead of the raw values? Default is \code{TRUE}.
#' @param ... currently ignored.
#' 
#' @return A list of \code{\link{brmcoda}} for each pivot balance coordinate.
#' 
#' @importFrom posterior summarise_draws as_draws_array
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   cilr <- complr(data = mcompd, sbp = sbp,
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID",
#'                  total = 1440)
#'   
#'   m <- brmcoda(complr = cilr,
#'                 formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                                    wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   m_pivot_coord_rotate <- pivot_coord_rotate(m)
#'   summary(m_pivot_coord_rotate)
#'   
#'   m_pivot_coord_raw <-  pivot_coord_rotate(m, summary = FALSE)
#'   posterior::summarise_draws(posterior::as_draws_array(m_pivot_coord_raw$output))
#'   }}
#' @export
pivot_coord_rotate <- function (object, summary = TRUE, ...) {
  
  out <- vector("list")
  
  b_sbp_0 <- fixef(object,
                   # newdata = model.frame(object),
                   summary = FALSE,
                   ...
  )
  
  # what type of model is being estimated
  model_fixef <- rownames(fixef(object))
  # model_ranef <- if(dim(object$model$ranef)[1] > 0) (names(ranef(object))) else (NULL)
  
  ilr_vars  <- grep("ilr", model_fixef, value = T)
  bilr_vars <- grep(".*bilr", model_fixef, value = T)
  wilr_vars <- grep(".*wilr", model_fixef, value = T)
  
  ilr_sbp_0  <- object$complr$logratio
  bilr_sbp_0 <- object$complr$between_logratio
  wilr_sbp_0 <- object$complr$within_logratio
  
  b_ilr_sbp_0  <- b_sbp_0[, colnames(b_sbp_0) %in% ilr_vars]
  b_bilr_sbp_0 <- b_sbp_0[, colnames(b_sbp_0) %in% bilr_vars]
  b_wilr_sbp_0 <- b_sbp_0[, colnames(b_sbp_0) %in% wilr_vars]
  
  b_sbp_summary_d <- vector("list")
  for (d in object$complr$parts) {
    parts_d <- append(d, grep(d, object$complr$parts, value = T, invert = T))
    sbp_d   <- build.sbp(parts_d)
    sbp_d   <- sbp_d[, object$complr$parts]
    
    clr_d   <- complr(data  = object$complr$data, 
                      sbp   = sbp_d,
                      parts = object$complr$parts,
                      idvar = object$complr$idvar,
                      total = object$complr$total)
    
    R <- crossprod(object$complr$psi, clr_d$psi)

    pars <- c("intercept", "between_logratio", "within_logratio", "logratio")
    b_sbp_target_i <- vector("list", length = length(pars))
    names(b_sbp_target_i) <- pars
    
    for (i in seq_along(b_sbp_target_i)) {
      if (i == 1) {
        b_sbp_target_i[[i]] <- b_sbp_0[, "Intercept"]
      } else {
        
        ## bilr
        if (i == 2) {
          if (length(grep("bilr", model_fixef, value = T)) > 0) {
            b_sbp_target_i[[i]] <- (b_bilr_sbp_0) %*% R
          } 
        }
        
        ## wilr
        if (i == 3) {
          if (length(grep("wilr", model_fixef, value = T)) > 0 ) {
            b_sbp_target_i[[i]] <- (b_wilr_sbp_0) %*% R
          }
        }
        
        ## ilr
        if (i == 4) {
          if ((length(grep("ilr", model_fixef, value = T)) > 0) 
              && (length(grep("[b|w]ilr", model_fixef, value = T)) == 0)) {
            b_sbp_target_i[[i]] <- (b_ilr_sbp_0) %*% R
          }
        }
      }
    }

    # take only non-empty elements (between vs within vs aggregate results)
    b_sbp_target_i <- Filter(Negate(is.null), b_sbp_target_i)
    
    if ("logratio" %in% names(b_sbp_target_i)) {
      level <- "aggregate"
      varnames <- c(bilr_vars, wilr_vars)
    } else {
      level <- c("between", "within")
      varnames <- c(ilr_vars)
    }
    
    # summarise posteriors
    if (isTRUE(summary)) {
      b_sbp_target_summary <- lapply(b_sbp_target_i, posterior_summary)
      b_sbp_target_summary <- do.call(rbind, b_sbp_target_summary)
      
      rownames(b_sbp_target_summary) <- c("Intercept", varnames)
      b_sbp_target_summary <- b_sbp_target_summary[rownames(b_sbp_target_summary) %in% c("bilr1", "wilr1", "ilr1"), ]
      
      # assemble output table
      b_sbp_target_summary <- cbind.data.frame(`Pivot coordinate` = paste0(d, "_vs_remaining"),
                                               Level = level,
                                               b_sbp_target_summary)
    } else {
      b_sbp_target_summary <- matrix(unlist(b_sbp_target_i), ncol = ndraws(object), byrow = TRUE)
      rownames(b_sbp_target_summary) <- c("Intercept", varnames)
      b_sbp_target_summary <- b_sbp_target_summary[rownames(b_sbp_target_summary) %in% c("bilr1", "wilr1", "ilr1"), ]
    }
    b_sbp_summary_d[[d]] <- b_sbp_target_summary
  }
  
  if (isTRUE(summary)) {
    out <- do.call(rbind, b_sbp_summary_d)
    rownames(out) <- NULL
  } else {
    out <- array(unlist(b_sbp_summary_d), 
                 dim = c(length(level),
                         ndraws(object),
                         length(object$complr$parts))
    )
    # out <- brms::do_call(abind::abind, c(out, along = 3))
    out <- aperm(out, c(2, 1, 3)) 
    
    dimnames(out)[[1]] <- 1:ndraws(object)
    dimnames(out)[[2]] <- level
    dimnames(out)[[3]] <- object$complr$parts
  }
  out <- structure(list(output = out, method = "rotate"), class = "pivot_coord")
  out
}

#' Estimate pivot balance coordinates by refitting model.
#' 
#' @param object An object of class \code{brmcoda}.
#' @param ... Further arguments passed to \code{\link[brms:brm]{brm}}.
#' 
#' @return A list of \code{\link{brmcoda}} for each pivot balance coordinate.
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   cilr <- complr(data = mcompd, sbp = sbp,
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID",
#'                  total = 1440)
#'   
#'   m <- brmcoda(complr = cilr,
#'                 formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                                    wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   m_pivot_coord_refit <- pivot_coord_refit(m)
#'   summary(m_pivot_coord_refit)
#'   }}
#' @export
pivot_coord_refit <- function (object, ...) {
  
  out_d <- vector("list")
  
  # loop through parts
  for (d in object$complr$parts) {
    parts_d <- append(d, grep(d, object$complr$parts, value = T, invert = T))
    sbp_d   <- build.sbp(parts_d)
    sbp_d   <- sbp_d[, object$complr$parts]
    
    clr_d <- complr(data  = object$complr$data, 
                    sbp   = sbp_d,
                    parts = object$complr$parts,
                    idvar = object$complr$idvar,
                    total = object$complr$total)
    
    dat_d <-  cbind(clr_d$data,
                    clr_d$between_logratio,
                    clr_d$within_logratio,
                    clr_d$logratio)
    
    fit_d <- update(object$model,
                    newdata = dat_d,
                    ...)
    
    brmcoda_d <- structure(list(complr = clr_d,
                                model  = fit_d),
                           class = "brmcoda")
    
    out_d[[d]] <-   brmcoda_d
  }
  out <- structure(list(output = out_d, method = "refit"), class = "pivot_coord")
  out
}
