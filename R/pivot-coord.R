#' Estimate pivot balance coordinates
#' 
#' This function estimates pivot balance coordinates for each compositional part by
#' either \code{"rotate"} the sequential binary partition using the same \code{brmcoda} object
#' or \code{"refit"} the \code{brmcoda} object.
#' 
#' @param object An object of class \code{brmcoda}.
#' @param method A character string.
#' Should the pivot balance coordinates be estimated by \code{"rotate"} the sequential binary partition 
#' using the same \code{brmcoda} object or \code{"refit"} the \code{brmcoda} object?
#' Default is \code{"rotate"}.
#' @param parts A optional character string specifying names of compositional parts that should be considered
#' in the substitution analysis. This should correspond to a single set of names of compositional parts specified
#' in the \code{complr} object. Default to the first composition in the \code{complr} object.
#' @param summary Should summary statistics be returned instead of the raw values? Default is \code{TRUE}.
#' @param ... Further arguments passed to \code{\link{posterior_summary}}.
#' 
#' @return Estimated pivot balance coordinates representing the
#' effect of increasing one compositional part relative to the remaining compositional parts.
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   x <- complr(data = mcompd, sbp = sbp,
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID",
#'                  total = 1440)
#'   
#'   # inspects ILRs before passing to brmcoda
#'   names(x$between_logratio)
#'   names(x$within_logratio)
#'   names(x$logratio)
#'   
#'   # model with compositional predictor at between and within-person levels
#'   m <- brmcoda(complr = x,
#'                 formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                    wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   m_pivot_coord <- pivot_coord(m)
#'   summary(m_pivot_coord)
#'   }}
#' @export
pivot_coord <- function (object,
                         summary = TRUE,
                         method = c("rotate", "refit"),
                         parts = 1,
                         ...) {
  if (all(c("rotate", "refit") %in% method)) {
    method <- "rotate"
  }
  if (method == "rotate") {
    out <- pivot_coord_rotate(object = object, summary = summary, parts = parts, ...)
  }
  if (method == "refit") {
    out <- pivot_coord_refit(object = object, parts = parts, ...)
  }
  structure(out, class = "pivot_coord")
  out
}

#' Estimate pivot balance coordinates by rotating sequential binary partition.
#' 
#' This function is an alias of \code{\link{pivot_coord}} to estimates the
#' pivot balance coordinates by \code{"rotate"} the sequential binary partition on the
#' same \code{brmcoda} object.
#' 
#' @seealso \code{\link{pivot_coord}}
#' 
#' @inheritParams pivot_coord
#' 
#' @inherit pivot_coord return
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   x <- complr(data = mcompd, sbp = sbp,
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID",
#'                  total = 1440)
#'   
#'   m <- brmcoda(complr = x,
#'                 formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                    wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   m_pivot_coord_rotate <- pivot_coord_rotate(m)
#'   summary(m_pivot_coord_rotate)
#'   
#'   m_pivot_coord_raw <-  pivot_coord_rotate(m, summary = FALSE)
#'   lapply(m_pivot_coord_raw$output, posterior_summary)
#'   }}
#' @export
pivot_coord_rotate <- function (object,
                                summary = TRUE,
                                parts = 1,
                                ...) {
  # get parts
  parts <- get_parts(object$complr, parts)
  
  ## get the index of which index elements of object$complr$output does the parts correspond to
  idx <- which(vapply(lapply(object$complr$output, function(x)
    x$parts), function(p)
      identical(sort(parts), sort(p)), logical(1)))
  
  
  b_sbp_0 <- fixef(object,
                   summary = FALSE,
                   ...
  )
  # grab the correct logratio names
  z_vars  <- get_variables(object$complr)[["Z", paste0("composition_", idx)]]
  bz_vars <- get_variables(object$complr)[["bZ", paste0("composition_", idx)]]
  wz_vars <- get_variables(object$complr)[["wZ", paste0("composition_", idx)]]
  
  b_z_sbp_0  <- b_sbp_0[, colnames(b_sbp_0) %in% z_vars]
  b_bz_sbp_0 <- b_sbp_0[, colnames(b_sbp_0) %in% bz_vars]
  b_wz_sbp_0 <- b_sbp_0[, colnames(b_sbp_0) %in% wz_vars]
  
  out <- vector("list")
  for (target in parts) {
    
    # new complr object with rotated sbp
    parts_target <- append(target, grep(target, parts, value = T, invert = T))
    sbp_target   <- build.sbp(parts_target)
    sbp_rotate   <- sbp_target[, parts]
    
    clr_target   <- complr(
      data  = object$complr$datain,
      sbp   = sbp_rotate,
      parts = parts,
      idvar = if(!is.null(object$complr$idvar)) object$complr$idvar else NULL,
      total = object$complr$output[[idx]]$total
    )
    
    # rotation matrix
    R <- crossprod(object$complr$output[[idx]]$psi, clr_target$output[[1]]$psi)
    
    # multiply posterior samples with rotation matrix
    b_sbp_target <- list(
      b_a_sbp_target  = b_sbp_0[, "Intercept"],
      b_bz_sbp_target = if (!is.null(colnames(b_bz_sbp_0)) && all(colnames(b_bz_sbp_0) %in% colnames(b_sbp_0))) as.matrix(b_bz_sbp_0) %*% R else NULL,
      b_wz_sbp_target = if (!is.null(colnames(b_wz_sbp_0)) && all(colnames(b_wz_sbp_0) %in% colnames(b_sbp_0))) as.matrix(b_wz_sbp_0) %*% R else NULL,
      b_z_sbp_target  = if (!is.null(colnames(b_z_sbp_0))  && all(colnames(b_z_sbp_0)  %in% colnames(b_sbp_0))) as.matrix(b_z_sbp_0) %*% R else NULL
    )
    
    # take only non-empty elements (between vs within vs aggregate results)
    b_sbp_target <- Filter(Negate(is.null), b_sbp_target)
    
    # name new variables the same as in the original model
    b_sbp_target <- lapply(names(b_sbp_target), function(n) {
      d <- as.matrix(b_sbp_target[[n]])
      if (n == "b_a_sbp_target")  colnames(d) <- "Intercept"
      if (n == "b_bz_sbp_target") colnames(d) <- colnames(b_bz_sbp_0)
      if (n == "b_wz_sbp_target") colnames(d) <- colnames(b_wz_sbp_0)
      if (n == "b_z_sbp_target")  colnames(d) <- colnames(b_z_sbp_0)
      d
    })
    b_sbp_target <- do.call(cbind, b_sbp_target)

    # summarise posteriors
    if (summary) {
      b_sbp_target_summary <- apply(b_sbp_target, 2, posterior_summary, ...)
      b_sbp_target_summary <- b_sbp_target_summary[, grep("z1", colnames(b_sbp_target), value = T), drop = FALSE]
      dimnames(b_sbp_target_summary) <- list(c("Estimate", "Est.Error", "CI_low", "CI_high"), grep("z1", colnames(b_sbp_target), value = TRUE))
    } else {
      b_sbp_target_summary <- b_sbp_target[, grep("z1", colnames(b_sbp_target), value = T), drop = FALSE]
    }
    out[[target]] <- b_sbp_target_summary
  }
  names(out) <- parts
  structure(list(output = out, method = "rotate", summary = summary), class = "pivot_coord")
}

#' Estimate pivot balance coordinates by refitting model.
#' 
#' This function is an alias of \code{\link{pivot_coord}} to estimates the
#' pivot balance coordinates by \code{"refit"} the \code{brmcoda} object.
#' 
#' @seealso \code{\link{pivot_coord}}
#' 
#' @inheritParams pivot_coord
#' 
#' @inherit pivot_coord return
#' 
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#'   x <- complr(data = mcompd, sbp = sbp,
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID",
#'                  total = 1440)
#'   
#'   m <- brmcoda(complr = x,
#'                 formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
#'                                    wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   m_pivot_coord_refit <- pivot_coord_refit(m)
#'   summary(m_pivot_coord_refit)
#'   
#'   m_pivot_coord_raw <-  pivot_coord_refit(m, summary = FALSE)
#'   lapply(m_pivot_coord_raw$output, posterior_summary)
#'   
#'   }}
#' @export
pivot_coord_refit <- function (object,
                               summary = TRUE,
                               parts = 1,
                               ...) {
  
  # get parts
  parts <- get_parts(object$complr, parts)
  
  ## get the index of which index elements of object$complr$output does the parts correspond to
  idx <- which(vapply(lapply(object$complr$output, function(x)
    x$parts), function(p)
      identical(sort(parts), sort(p)), logical(1)))
  
  b_sbp_0 <- fixef(object,
                   summary = FALSE,
                   ...
  )
  # grab the correct logratio names
  z_vars  <- get_variables(object$complr)[["Z", paste0("composition_", idx)]]
  bz_vars <- get_variables(object$complr)[["bZ", paste0("composition_", idx)]]
  wz_vars <- get_variables(object$complr)[["wZ", paste0("composition_", idx)]]
  
  b_z_sbp_0  <- b_sbp_0[, colnames(b_sbp_0) %in% z_vars]
  b_bz_sbp_0 <- b_sbp_0[, colnames(b_sbp_0) %in% bz_vars]
  b_wz_sbp_0 <- b_sbp_0[, colnames(b_sbp_0) %in% wz_vars]
  
  out <- vector("list")
  for (target in parts) {
    parts_target <- append(target, grep(target, parts, value = T, invert = T))
    sbp_target   <- build.sbp(parts_target)
    # sbp_target   <- sbp_target[, parts]
    
    clr_target   <- complr(
      data  = object$complr$datain,
      sbp   = sbp_target,
      parts = parts_target,
      idvar = if(!is.null(object$complr$idvar)) (object$complr$idvar) else NULL,
      total = object$complr$output[[idx]]$total
    )
    new_data <- cbind(clr_target$dataout, object$complr$dataout[, colnames(object$complr$dataout) %nin% colnames(clr_target$dataout), with = FALSE])
    
    brmcoda_target <- update(object$model, newdata = new_data, ...)
    
    b_sbp_target <- fixef(brmcoda_target,
                          summary = FALSE,
                          ...)
    
    # summarise posteriors
    if (summary) {
      b_sbp_target_summary <- apply(b_sbp_target, 2, posterior_summary, ...)
      b_sbp_target_summary <- b_sbp_target_summary[, grep("z1", colnames(b_sbp_target), value = T), drop = FALSE]
      dimnames(b_sbp_target_summary) <- list(c("Estimate", "Est.Error", "CI_low", "CI_high"), grep("z1", colnames(b_sbp_target), value = TRUE))
    } else {
      b_sbp_target_summary <- b_sbp_target[, grep("z1", colnames(b_sbp_target), value = T), drop = FALSE]
    }
    out[[target]] <- b_sbp_target_summary
  }
  names(out) <- parts
  structure(list(output = out, method = "refit", summary = summary), class = "pivot_coord")
}

