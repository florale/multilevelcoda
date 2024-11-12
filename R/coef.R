#' #' Calculate "coefficients" based on substitutions for each compositional part
#' #'
#' #' @param object An object of class \code{brmcoda}.
#' #' @param level A character string specifying the level of the coefficients to be calculated.
#' #'   Either \dQuote{between} or \dQuote{within}.
#' #' @param h A numeric value specifying the step size for the substitution.
#' #' @return A data table of results.
#' #' @importFrom data.table as.data.table copy
#' #' @importFrom stats fitted model.frame
#' #' @importFrom testthat expect_equal
#' #' @importFrom brms posterior_summary
#' #' @importFrom nlme fixef
#' #' @export
#' #' @examples
#' #' \donttest{
#' #' if(requireNamespace("cmdstanr")){
#' #' sbp2 <- sbp
#' #' sbp2[1, ] <- c(-1, 1, -1, -1, -1)
#' #' sbp2[2, ] <- c( 1, 0, -1, -1, -1)
#' #' sbp2[3, ] <- c( 0, 0,  1, -1, -1)
#' #' sbp2[4, ] <- c( 0, 0,  0,  1, -1)
#' #' 
#' #' cilr <- complr(data = mcompd, sbp = sbp2, 
#' #'   parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#' #'   idvar = "ID")
#' #' 
#' #' m1 <- brmcoda(complr = cilr,
#' #'               formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
#' #'                                  wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
#' #'               chain = 4, iter = 1000, cores = 4L,
#' #'               backend = "cmdstanr")
#' #' substition_coef(m1, level = "between", h = 10)
#' #' substition_coef(m1, level = "within", h = 10)
#' #' rm(sbp2, cilr, m1) ## cleanup
#' #' }
#' #' }
#' substition_coef <- function(object, level = c("between", "within"), h = 10) {
#'   level <- match.arg(level)
#'   expect_s3_class(object, "brmcoda")
#'   
#'   if (object$model$family$family == "gaussian" && object$model$family$link == "identity") {
#'     linear <- TRUE
#'   } else {
#'     linear <- FALSE
#'   }
#'   
#'   parts <- object$complr$parts
#'   x <- object$complr$data[, ..parts]
#'   
#'   if (isFALSE(linear)) {
#'     y0 <- fitted(object$model,
#'                  newdata = model.frame(object),
#'                  re_formula = NA,
#'                  summary = FALSE
#'     )
#'   }
#'   
#'   out <- vector("list", length(parts))
#'   
#'   for (k in seq_along(parts)) {
#'     w <- (-x) / rowSums(x[, -..k])
#'     w <- as.data.table(w)
#'     w[, (parts[k]) := +1]
#'     
#'     x2 <- x + (h * w)
#'     
#'     expect_equal(rowSums(x2), rowSums(x), tolerance = 1e-3)
#'     
#'     rm(w) ## cleanup
#'     
#'     d2 <- copy(object$complr$data)
#'     d2[, (parts) := x2]
#'     
#'     rm(x2) ## cleanup
#'     
#'     cilr2 <- complr(
#'       data = d2,
#'       sbp = object$complr$sbp,
#'       parts = parts,
#'       idvar = object$complr$idvar
#'     )
#'     
#'     if  (isFALSE(linear)) {
#'       switch(level,
#'              within = {
#'                bilr2 <- object$complr$between_logratio
#'                wilr2 <- cilr2$logratio - object$complr$between_logratio
#'                names(wilr2) <- names(object$complr$within_logratio)
#'              },
#'              between = {
#'                # expect_equal(cilr$within_logratio, cilr2$within_logratio, tolerance = 1e-3)
#'                bilr2 <- cilr2$between_logratio
#'                wilr2 <- object$complr$within_logratio
#'              }
#'       )
#'       
#'       rm(cilr2) ## cleanup
#'       
#'       d2 <- cbind(d2, bilr2, wilr2)
#'       
#'       rm(bilr2, wilr2)
#'       
#'       y2 <- fitted(object$model,
#'                    newdata = d2,
#'                    re_formula = NA,
#'                    summary = FALSE
#'       )
#'       
#'       rm(d2) ## cleanup
#'       
#'       out[[k]] <- rowMeans((y2 - y0) )
#'       
#'       rm(y2)
#'     } else if (isTRUE(linear)) {
#'       switch(level,
#'              within = {
#'                wilr2 <- cilr2$logratio - object$complr$between_logratio
#'                out[[k]] <- fixef(object$model, summary = FALSE)[, colnames(object$complr$within_logratio)] %*% 
#'                  colMeans(wilr2 - object$complr$within_logratio)
#'              },
#'              between = {
#'                out[[k]] <- fixef(object$model, summary = FALSE)[, colnames(object$complr$between_logratio)] %*% 
#'                  colMeans(cilr2$between_logratio - object$complr$between_logratio)
#'              }
#'       )
#'     }
#'   }
#'   
#'   if (isTRUE(linear)) rm(x) else rm(x, y0) ## cleanup
#'   
#'   finalout <- cbind(
#'     Part = parts,
#'     as.data.table(do.call(rbind, lapply(out, posterior_summary)))
#'   )
#'   
#'   return(finalout)
#' }
