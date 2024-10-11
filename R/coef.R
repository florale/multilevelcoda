#' Calculate "coefficients" based on substitutions for each compositional part
#'
#' @param object An object of class \code{brmcoda}.
#' @param level A character string specifying the level of the coefficients to be calculated.
#'   Either \dQuote{between} or \dQuote{within}.
#' @param h A numeric value specifying the step size for the substitution.
#' @return A data table of results.
#' @importFrom data.table as.data.table copy
#' @importFrom stats fitted model.frame
#' @importFrom testthat expect_equal
#' @importFrom brms posterior_summary
#' @importFrom nlme fixef
#' @export
#' @examples
#' \donttest{
#' if(requireNamespace("cmdstanr")){
#' sbp2 <- sbp
#' sbp2[1, ] <- c(-1, 1, -1, -1, -1)
#' sbp2[2, ] <- c( 1, 0, -1, -1, -1)
#' sbp2[3, ] <- c( 0, 0,  1, -1, -1)
#' sbp2[4, ] <- c( 0, 0,  0,  1, -1)
#' 
#' cilr <- complr(data = mcompd, sbp = sbp2, 
#'   parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'   idvar = "ID")
#' 
#' m1 <- brmcoda(complr = cilr,
#'               formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                                  wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'               chain = 4, iter = 1000, cores = 4L,
#'               backend = "cmdstanr")
#' substition_coef(m1, level = "between", h = 10)
#' substition_coef(m1, level = "within", h = 10)
#' rm(sbp2, cilr, m1) ## cleanup
#' }
#' }
substition_coef <- function(object, level = c("between", "within"), h = 10) {
  level <- match.arg(level)
  expect_s3_class(object, "brmcoda")
  
  if (object$model$family$family == "gaussian" && object$model$family$link == "identity") {
    linear <- TRUE
  } else {
    linear <- FALSE
  }
  
  parts <- object$complr$parts
  x <- object$complr$data[, ..parts]
  
  if (isFALSE(linear)) {
    y0 <- fitted(object$model,
                 newdata = model.frame(object),
                 re_formula = NA,
                 summary = FALSE
    )
  }
  
  out <- vector("list", length(parts))
  
  for (k in seq_along(parts)) {
    w <- (-x) / rowSums(x[, -..k])
    w <- as.data.table(w)
    w[, (parts[k]) := +1]
    
    x2 <- x + (h * w)
    
    expect_equal(rowSums(x2), rowSums(x), tolerance = 1e-3)
    
    rm(w) ## cleanup
    
    d2 <- copy(object$complr$data)
    d2[, (parts) := x2]
    
    rm(x2) ## cleanup
    
    cilr2 <- complr(
      data = d2,
      sbp = object$complr$sbp,
      parts = parts,
      idvar = object$complr$idvar
    )
    
    if  (isFALSE(linear)) {
      switch(level,
             within = {
               bilr2 <- object$complr$between_logratio
               wilr2 <- cilr2$logratio - object$complr$between_logratio
               names(wilr2) <- names(object$complr$within_logratio)
             },
             between = {
               # expect_equal(cilr$within_logratio, cilr2$within_logratio, tolerance = 1e-3)
               bilr2 <- cilr2$between_logratio
               wilr2 <- object$complr$within_logratio
             }
      )
      
      rm(cilr2) ## cleanup
      
      d2 <- cbind(d2, bilr2, wilr2)
      
      rm(bilr2, wilr2)
      
      y2 <- fitted(object$model,
                   newdata = d2,
                   re_formula = NA,
                   summary = FALSE
      )
      
      rm(d2) ## cleanup
      
      out[[k]] <- rowMeans((y2 - y0) )
      
      rm(y2)
    } else if (isTRUE(linear)) {
      switch(level,
             within = {
               wilr2 <- cilr2$logratio - object$complr$between_logratio
               out[[k]] <- fixef(object$model, summary = FALSE)[, colnames(object$complr$within_logratio)] %*% 
                 colMeans(wilr2 - object$complr$within_logratio)
             },
             between = {
               out[[k]] <- fixef(object$model, summary = FALSE)[, colnames(object$complr$between_logratio)] %*% 
                 colMeans(cilr2$between_logratio - object$complr$between_logratio)
             }
      )
    }
  }
  
  if (isTRUE(linear)) rm(x) else rm(x, y0) ## cleanup
  
  finalout <- cbind(
    Part = parts,
    as.data.table(do.call(rbind, lapply(out, posterior_summary)))
  )
  
  return(finalout)
}

#' Estimate pivot balance coordinates 
#' 
#' @param object An object of class \code{brmcoda}.
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
#'   m1 <- brmcoda(complr = cilr,
#'                 formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#'                                    wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#'                 chain = 1, iter = 500,
#'                 backend = "cmdstanr")
#'   
#'   m_pb <- fixef.brmcoda_pivot(m1)
#'   m_pb_raw <- fixef.brmcoda_pivot(m1, summary = FALSE)
#'   brms::posterior_summary(m_pb_raw)
#'   }}
#' @export
#' 
fixef.brmcoda_pivot <- function (object, summary = TRUE, ...) {
  
  out_d <- vector("list")
  
  b_sbp_0 <- fixef(object,
                   # newdata = model.frame(object),
                   summary = FALSE,
                   ...
  )
  
  # what type of model is being estimated
  model_fixef <- rownames(fixef(object))
  model_ranef <- if(dim(object$model$ranef)[1] > 0) (names(ranef(object))) else (NULL)
  
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
    
    pars <- c("Intercept", "bilr", "wilr", "ilr")
    
    b_sbp_target_i <- vector("list", length = length(pars))
    names(b_sbp_target_i) <- c("intercept", "between_logratio", "within_logratio", "logratio")
    
    for (i in seq_along(pars)) {
      if (i == 1) {
        b_sbp_target_i[[i]] <- b_sbp_0[, "Intercept"]
      } else {
        
        ## bilr
        if (i == 2) {
          if (length(grep("bilr", model_fixef, value = T)) > 0) {
            
            b_bilr_sbp_target_j <- vector("list")
            for (j in seq_len(ncol(b_bilr_sbp_0))) {
              
              bilr_sbp_1 <- clr_d$between_logratio
              
              bilr_sbp_prime      <- bilr_sbp_1
              bilr_sbp_prime[, j] <- bilr_sbp_prime[, j] + 1
              
              d_prime_b <- ilrInv(bilr_sbp_prime, V = gsi.buildilrBase(t(clr_d$sbp)))
              d_prime_b <- acomp(d_prime_b, total = object$complr$total)
              
              bilr_prime <- ilr(d_prime_b, V = gsi.buildilrBase(t(object$complr$sbp)))
              colnames(bilr_prime) <- colnames(bilr_sbp_0)
              
              b_bilr_sbp_target <- apply(b_bilr_sbp_0, 1, function(x)
                colMeans(((bilr_prime - bilr_sbp_0) %*% x))
              )
              b_bilr_sbp_target_j[[j]] <- b_bilr_sbp_target
            }
            b_sbp_target_i[[i]] <- do.call(cbind, b_bilr_sbp_target_j)
          } 
        }
        
        ## wilr
        if (i == 3) {
          if (length(grep("wilr", model_fixef, value = T)) > 0 ) {
            b_wilr_sbp_target_j <- vector("list")
            for (j in seq_len(ncol(b_wilr_sbp_0))) {
              
              wilr_sbp_1 <- clr_d$within_logratio
              
              wilr_sbp_prime      <- wilr_sbp_1
              wilr_sbp_prime[, j] <- wilr_sbp_prime[, j] + 1
              
              d_prime_w <- ilrInv(wilr_sbp_prime, V = gsi.buildilrBase(t(clr_d$sbp)))
              d_prime_w <- acomp(d_prime_w, total = 1)
              
              wilr_prime <- ilr(d_prime_w, V = gsi.buildilrBase(t(object$complr$sbp)))
              colnames(wilr_prime) <- colnames(wilr_sbp_0)
              
              # b_sbp_target[[i]] <- sweep(b_wilr_sbp_0, 2, (wilr_prime - wilr_sbp_0)[1,], `%*%`)
              b_wilr_sbp_target <- apply(b_wilr_sbp_0, 1, function(x) 
                colMeans(((wilr_prime - wilr_sbp_0) %*% x))
              )
              
              b_wilr_sbp_target_j[[j]] <- b_wilr_sbp_target
            }
            b_sbp_target_i[[i]] <- do.call(cbind, b_wilr_sbp_target_j)
          }
        }
        
        ## ilr
        if (i == 4) {
          if ((length(grep("ilr", model_fixef, value = T)) > 0) 
              && (length(grep("[b|w]ilr", model_fixef, value = T)) == 0)) {
            
            b_ilr_sbp_target_j <- vector("list")
            for (j in seq_len(ncol(b_ilr_sbp_0))) {
              
              ilr_sbp_1 <- clr_d$logratio
              
              ilr_sbp_prime <- ilr_sbp_1
              ilr_sbp_prime[, j] <- ilr_sbp_prime[, j] + 1
              
              d_prime <- ilrInv(ilr_sbp_prime, V = gsi.buildilrBase(t(clr_d$sbp)))
              d_prime <- acomp(d_prime, total = object$complr$total)
              
              ilr_prime <- ilr(d_prime, V = gsi.buildilrBase(t(object$complr$sbp)))
              colnames(ilr_prime) <- colnames(ilr_sbp_0)
              
              b_ilr_sbp_target <- apply(b_ilr_sbp_0, 1, function(x) 
                colMeans(((ilr_prime - ilr_sbp_0) %*% x))
              )
              b_ilr_sbp_target_j[[j]] <- b_ilr_sbp_target
            }
            b_sbp_target_i[[i]] <- do.call(cbind, b_ilr_sbp_target_j)
          }
        }
      }
    }
    
    # take only non-empty elements (between vs within vs aggregate results)
    b_sbp_target_i <- Filter(Negate(is.null), b_sbp_target_i)
    
    if ("logratio" %in% names(b_sbp_target_i)) {
      level <- "aggregate"
    } else {
      level <- c("between", "within")
    }
    
    # summarise posteriors
    if (isTRUE(summary)) {
      b_sbp_target_summary <- lapply(b_sbp_target_i, posterior_summary)
      b_sbp_target_summary <- do.call(rbind, b_sbp_target_summary)
      
      rownames(b_sbp_target_summary) <- rownames(fixef(object))
      b_sbp_target_summary <- b_sbp_target_summary[rownames(b_sbp_target_summary) %in% c("bilr1", "wilr1", "ilr1"), ]
      
      # assemble output table
      b_sbp_target_summary <- cbind.data.frame(Part = d,
                                               Level = level,
                                               b_sbp_target_summary)
    } else {
      b_sbp_target_summary <- matrix(unlist(b_sbp_target_i), ncol = ndraws(object), byrow = TRUE)
      rownames(b_sbp_target_summary) <- rownames(fixef(object))
      b_sbp_target_summary <- b_sbp_target_summary[rownames(b_sbp_target_summary) %in% c("bilr1", "wilr1", "ilr1"), ]
    }
    b_sbp_summary_d[[d]] <- b_sbp_target_summary
  }
  
  if (isTRUE(summary)) {
    b_sbp_summary_d <- do.call(rbind, b_sbp_summary_d)
    rownames(b_sbp_summary_d) <- NULL
  } else {
    b_sbp_summary_d <- array(unlist(b_sbp_summary_d), 
                             dim = c(length(level),
                                     ndraws(object),
                                     length(object$complr$parts))
    )
    # b_sbp_summary_d <- brms::do_call(abind::abind, c(b_sbp_summary_d, along = 3))
    b_sbp_summary_d <- aperm(b_sbp_summary_d, c(2, 1, 3)) 
    
    dimnames(b_sbp_summary_d)[[1]] <- 1:ndraws(object)
    dimnames(b_sbp_summary_d)[[2]] <- level
    dimnames(b_sbp_summary_d)[[3]] <- object$complr$parts
  }
  b_sbp_summary_d
}

# # #### THIS IS WHAT WE HAVE
# # print(sbp) ## starting SBP
# # 
# clr_sbp <- complr(data = mcompd, sbp = sbp,
#                   parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#                   idvar = "ID")
# 
# m1 <- lm(Stress ~ ilr1 + ilr2 + ilr3 + ilr4,
#          data = cbind(clr_sbp$data, clr_sbp$logratio))
# 
# coef(m1)
# 
# proof_of_concept <- function(d, b_sbp, clr_sbp, sbp_target) {
#   ilr_sbp <- as.matrix(clr_sbp$logratio)
#   sbp  <- clr_sbp$sbp
# 
#   clr_target <- complr(data = d, sbp = sbp_target,
#                        parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#                        idvar = "ID")
# 
#   ilr_sbp_target <- as.matrix(clr_target$logratio)
# 
#   b_sbp_target <- vector("numeric", length(b_sbp))
# 
#   for (i in seq_along(b_sbp)) {
#     if (i == 1) {
#       b_sbp_target[i] <- b_sbp[i]
#     } else {
#       ilr_sbp_target_prime <- ilr_sbp_target
#       ilr_sbp_target_prime[, i - 1] <- ilr_sbp_target_prime[, i - 1] + 1
#       d_prime <- ilrInv(ilr_sbp_target_prime,
#                         V = gsi.buildilrBase(t(sbp_target))
#       )
# 
#       ilr_sbp_prime <- ilr(d_prime, V = gsi.buildilrBase(t(sbp)))
# 
#       b_sbp_target[i] <- as.numeric(((ilr_sbp_prime - ilr_sbp)[1,] %*% b_sbp[-1]))
#     }
#   }
#   b_sbp_target
# }
# 
# 
# ## THIS IS WHAT WE WANT (i.e., the "target")
# sbp_target <- sbp
# sbp_target[1, ] <- c( 1, -1, -1, -1, -1)
# sbp_target[2, ] <- c( 0, 1, -1, -1, -1)
# sbp_target[3, ] <- c( 0, 0,  1, -1, -1)
# sbp_target[4, ] <- c( 0, 0,  0,  1, -1)
# 
# proof_of_concept(d = mcompd, b_sbp = coef(m1),
#                  clr = clr_sbp, sbp_target = sbp_target)
# 
# ## let's compare against refitting the model
# 
# clr_sbp_target <- complr(data = mcompd, sbp = sbp_target,
#                          parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#                          idvar = "ID")
# 
# m_target <- lm(Stress ~ ilr1 + ilr2 + ilr3 + ilr4,
#                data = cbind(clr_sbp_target$data, clr_sbp_target$logratio))
# 
# coef(m_target)
# 
