get_variables = function(object) UseMethod("get_variables")

#' Extract variable names from a \code{complr} object.
#' @param object A \code{complr} object
#' 
#' @method get_variables complr
#' 
#' @rdname get_variables
#' @export
get_variables.complr <- function(object) {
  out <- lapply(object$output, function(x) {
    list(X = names(x$X),
         bX = names(x$bX),
         wX = names(x$wX),
         Z = names(x$Z),
         bZ = names(x$bZ),
         wZ =names(x$wZ))
  })
  names(out) <- paste0("composition_", seq_along(out))
  # print(as.data.frame(do.call(cbind, out)))
  out
}

#' Extract variable names from a \code{brmcoda} object.
#' @param object A \code{brmcoda} object
#'
#' @method get_variables brmcoda
#' @rdname get_variables
#' @export
get_variables.brmcoda <- function(object) {
  if (!is.brmcoda(object)) {
    stop("object must be of class 'brmcoda'")
  }
  
  # grab complr varnames
  complr_vars <- list(
    X  = unlist(lapply(
      get_variables(object$complr), `[[`, "X"
    ), use.names = FALSE),
    bX = unlist(lapply(
      get_variables(object$complr), `[[`, "bX"
    ), use.names = FALSE),
    wX = unlist(lapply(
      get_variables(object$complr), `[[`, "wX"
    ), use.names = FALSE),
    Z  = unlist(lapply(
      get_variables(object$complr), `[[`, "Z"
    ), use.names = FALSE),
    bZ = unlist(lapply(
      get_variables(object$complr), `[[`, "bZ"
    ), use.names = FALSE),
    wZ = unlist(lapply(
      get_variables(object$complr), `[[`, "wZ"
    ), use.names = FALSE)
  )
  
  # what type of model is being estimated
  model_ranef <- if (dim(object$model$ranef)[1] > 0)
    (names(ranef(object)))
  else
    (NULL)
  
  # single level or multilevel
  if (length(model_ranef) > 0) {
    model_ranef_type <- "multilevel"
  } else {
    model_ranef_type <- "single"
  }
  model_fixef_type <- model_resp_type <- NULL
  
  if (inherits(object$model$formula, "mvbrmsformula")) {
    # extract variable names from brms formula
    formulas <- object[["model"]][["formula"]][["forms"]]
    model_vars <- lapply(formulas, function(f) {
      split_formula <- strsplit(as.character(f), "~")[[1]]
      lhs <- trimws(split_formula[1])
      rhs <- trimws(split_formula[2])
      rhs <- trimws(unlist(strsplit(rhs, "\\+|\\*|\\^|\\(|\\)|\\|")))
      rhs <- rhs[rhs != ""]
      rhs <- rhs %sin% colnames(model.frame(object))
      list(lhs = lhs, rhs = rhs)
    })
    
    ## variables in brms formula
    model_resp  <- unique(unlist(lapply(model_vars, function(x)
      x$lhs)))
    model_pred  <- unique(unlist(lapply(model_vars, function(x)
      x$rhs)))
    
    ## response type
    if (length(complr_vars$bZ) > 0 && all(complr_vars$bZ %in% model_resp)) {
      model_resp_type <- append(model_resp_type, "multivariate-between")
    }
    if (length(complr_vars$wZ) > 0 && all(complr_vars$wZ %in% model_resp)) {
      model_resp_type <- append(model_resp_type, "multivariate-within")
    }
    if (length(complr_vars$Z) > 0 && all(complr_vars$Z %in% model_resp)) {
      model_resp_type <- append(model_resp_type, "multivariate-aggregate")
    }

    if (!is.null(model_resp_type) && length(model_resp_type) > 1) {
      stop(
        "Multivariate response models with different levels (between, within, aggregate) are not supported."
      )
    }
    if (is.null(model_resp_type)) {
      model_resp_type <- "multivariate-non-compositional"
    }
    
    ## predictors type
    ### check if the predictors include z variables
    if (any(c(complr_vars$bZ, complr_vars$wZ, complr_vars$Z) %in% model_pred)) {
      model_pred_type <- "compositional"
      
      # if so which ones
      bZ_pred  <- model_pred[which(model_pred %in% c(complr_vars$bZ))]
      wZ_pred  <- model_pred[which(model_pred %in% c(complr_vars$wZ))]
      Z_pred   <- model_pred[which(model_pred %in% c(complr_vars$Z))]
      nZ_pred  <- model_pred[which(model_pred %nin% c(complr_vars$bZ, complr_vars$wZ, complr_vars$Z))]
      
      if (length(bZ_pred) > 0 && all(bZ_pred %in% complr_vars$bZ)) {
        model_fixef_type <- append(model_fixef_type, "between")
      }
      if (length(wZ_pred) > 0 && all(wZ_pred %in% complr_vars$wZ)) {
        model_fixef_type <- append(model_fixef_type, "within")
      }
      if (length(Z_pred) > 0 && all(Z_pred %in% complr_vars$Z)) {
        model_fixef_type <- append(model_fixef_type, "aggregate")
      }
    } else {
      model_pred_type <- "non-compositional"
    }
    
  } else {
    formulas <- as.character(object[["model"]][["formula"]])[[1]]
    split_formula <- strsplit(formulas, "~")[[1]]
    lhs <- trimws(split_formula[1])
    rhs <- trimws(split_formula[2])
    rhs <- trimws(unlist(strsplit(rhs, "\\+|\\*|\\^|\\(|\\)|\\|")))
    rhs <- rhs[rhs != ""]
    rhs <- rhs %sin% colnames(model.frame(object))
    # rhs <- sapply(rhs, function(term) {
    #   # If term is a function call, extract what's inside the parentheses, and anything after a '|'
    #   m <- regmatches(term, regexec("\\(([^)]+)\\)", term))
    #   if (length(m[[1]]) > 1) {
    #     inside <- m[[1]][2]
    #     if (grepl("\\|", inside)) {
    #       trimws(strsplit(inside, "\\|")[[1]][2])
    #     } else {
    #       trimws(inside)
    #     }
    #   } else {
    #     term
    #   }
    # })
    model_vars <- list(lhs = lhs, rhs = rhs)
    
    ## z variables in brms formula
    model_resp  <- model_vars$lhs
    model_pred  <- model_vars$rhs
    
    # response type
    model_resp_type <- "univariate"
    
    ## predictors type
    ### check if the predictors include z variables
    if (any(c(complr_vars$bZ, complr_vars$wZ, complr_vars$Z) %in% model_pred)) {
      model_pred_type <- "compositional"
      
      # if so which ones
      bZ_pred  <- model_pred[which(model_pred %in% complr_vars$bZ)]
      wZ_pred  <- model_pred[which(model_pred %in% complr_vars$wZ)]
      Z_pred   <- model_pred[which(model_pred %in% complr_vars$Z)]
      nZ_pred  <- model_pred[which(model_pred %nin% c(complr_vars$bZ, complr_vars$wZ, complr_vars$Z))]
      
      if (length(bZ_pred) > 0 && all(bZ_pred %in% complr_vars$bZ)) {
        model_fixef_type <- append(model_fixef_type, "between")
      }
      if (length(wZ_pred) > 0 && all(wZ_pred %in% complr_vars$wZ)) {
        model_fixef_type <- append(model_fixef_type, "within")
      }
      if (length(Z_pred) > 0 && all(Z_pred %in% complr_vars$Z)) {
        model_fixef_type <- append(model_fixef_type, "aggregate")
      }
    } else {
      model_pred_type <- "non-compositional"
    }
  }
  list(
    resp_type  = model_resp_type,
    fixef_type = model_fixef_type,
    ranef_type = model_ranef_type,
    y = model_resp,
    x = model_pred
  )
}

#' Extract Sequential Binary Partition from a \code{complr} object.
#' 
#' @param object A \code{complr} object
#' 
#' @export
get_sbp <- function(object) {
  if (isFALSE(inherits(object, "complr"))) {
    stop(sprintf(
      "Can't handle an object of class (%s)
  It should be a 'complr' object
  See ?complr for details.",
      class(object)))
  }
  lapply(object$output, function(x) x$sbp)
}

#' Extract variables from complr and brmcoda objects for use in substitution models
#' 
#' Internal use only
#' @noRd
.get.subvars <- function(object, parts, scale) {
  brmcoda_vars <- get_variables(object)
  complr_vars  <- get_variables(object$complr)
  
  ## get the index of which index elements of object$complr$output do the parts correspond to
  idx <- as.integer(which(vapply(lapply(object[["complr"]][["output"]], function(x)
    x$parts), function(p)
      identical(sort(parts), sort(p)), logical(1))))[1]
  
  idy <- as.integer(which(vapply(complr_vars, function(y) {
    any(sapply(c("Z", "bZ", "wZ"), function(z) {
      identical(sort(y[[z]]), sort(brmcoda_vars[["y"]]))
    }))
  }, logical(1))))[1]
  
  # grab logratio and composition names
  z_vars  <- complr_vars[[paste0("composition_", idx)]][["Z"]]
  bz_vars <- complr_vars[[paste0("composition_", idx)]][["bZ"]]
  wz_vars <- complr_vars[[paste0("composition_", idx)]][["wZ"]]
  x_vars  <- complr_vars[[paste0("composition_", idx)]][["X"]]
  bx_vars <- complr_vars[[paste0("composition_", idx)]][["bX"]]
  wx_vars <- complr_vars[[paste0("composition_", idx)]][["wX"]]
  
  xz_vars <- c(z_vars, bz_vars, wz_vars, x_vars, bx_vars, wx_vars)
  sx_vars <- paste0("s", object[["complr"]][["output"]][[idx]][["parts"]])
  
  if (inherits(object[["model"]][["formula"]], "mvbrmsformula") &&
      (
        identical(brmcoda_vars[["y"]], z_vars)  ||
        identical(brmcoda_vars[["y"]], bz_vars) ||
        identical(brmcoda_vars[["y"]], wz_vars)
      )) {
    if (identical(scale, "response")) {
      y_vars <- complr_vars[[paste0("composition_", idy)]][["X"]]
    } else {
      y_vars <- complr_vars[[paste0("composition_", idy)]][["Z"]]
    }
  } else {
    y_vars <- brmcoda_vars[["y"]]
  }
  
  list(brmcoda_vars = brmcoda_vars,
       complr_vars = complr_vars,
       z_vars = z_vars,
       bz_vars = bz_vars,
       wz_vars = wz_vars,
       x_vars = x_vars,
       bx_vars = bx_vars,
       wx_vars = wx_vars,
       xz_vars = xz_vars,
       sx_vars = sx_vars,
       y_vars = y_vars,
       idx = idx,
       idy = if(!is.na(idy)) idy else 0L
  )
}

#' Extract names of parts from a \code{complr} object.
#' 
#' Internal use only
#' 
#' @param object A \code{complr} object
#' @param parts A optional character string specifying names of compositional parts that should be considered
#' in the substitution analysis. This should correspond to a single set of names of compositional parts specified
#' in the \code{complr} object. Default to the first composition in the \code{complr} object.
#' 
#' @keywords internal
#' @noRd
.get_parts <- function(object, parts = 1) {
  
  if (isFALSE(inherits(object, "complr"))) {
    stop(sprintf(
      "Can't handle an object of class (%s)
  It should be a 'complr' object
  See ?complr for details.",
      class(object)))
  }
  
  if (is.numeric(parts)) {
    if (length(parts) > 1) {
      stop(" 'parts' should be a single numeric value indicating which set of compositional parts to use.")
    }
    if (parts < 1 || parts > length(object$output)) {
      stop(sprintf(
        " 'parts' should be a single numeric value between 1 and %s, corresponding to the number of sets of compositional parts in the 'complr' object.",
        length(object$output)))
    }
    parts <- object$output[[parts]]$parts
    
  } else {
    if (isFALSE(inherits(parts, "character"))) {
      stop(" 'parts' should be a character vector of compositional parts.")
    }
    ## parts should be identical with either one of the parts presented in output of complr
    if (isFALSE((any(vapply(lapply(object$output, function(x) x$parts), function(p) identical(sort(parts), sort(p)), logical(1)))))) {
      stop(sprintf(
        "The specified 'parts' (%s) are not found in the complr object.",
        "  It should corespond to one set of compositional parts, either one of the following:",
        "%s",
        paste(parts, collapse = ", "),
        invisible(lapply(object$output, function(x) cat(paste(x$parts, collapse = ", "), "\n"))),
        sep = "\n"))
    }
  }
  parts
}