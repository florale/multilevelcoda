#' #' @title Margianl Effects of Between-person Substitution
#' #'
#' #' Estimate the difference in outcomes
#' #' when compositional variables are substituted for a specific time period
#' #' at between-person level.
#' #'
#' #' @param data A fitted \code{brms} model data. Required.
#' #' @param substitute A data frame or data table indicating the possible substitution of variables. This dataset can be computed using \code{possub}. Required.
#' #' @param minute A integer or numeric value indicating the minute that compositional variable are substituted to/from. Default is 60L.
#' #' @return A list
#' #' @importFrom data.table as.data.table copy := setnames
#' #' @importFrom compositions acomp ilr clo
#' #' @importFrom extraoperators %snin% %sin%
#' #' @export
#' #' @examples
#' #' ps <- possub(data = mcompd, composition = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' #' 
#' #' bsubmarginstest <- bsubmargins(data = mcmc, substitute = ps, minute = 10)
#' bsubmargins <- function (data, substitute, minute = 60L) {
#' 
#'   if(isFALSE(missing(minute))) {
#'     if (isFALSE(is.integer(minute))) {
#'       if (isFALSE(is.numeric(minute))) {
#'         stop("'minute' must be an integer or a numeric value > 0.")
#'       }
#'     }
#'     
#'   } else if (isTRUE(missing(minute))) {
#'     minute <- 60L
#'   }
#'   
#'   if (isFALSE(identical(ncol(substitute), length(data$CompIlr$composition)))) {
#'     stop(sprintf("The number of columns in 'substitute' (%d) must be the same 
#'   as the compositional variables in 'composition' (%d).",
#'                  ncol(substitute),
#'                  length(data$CompIlr$composition)))
#'   }
#' 
#'   if (isFALSE(identical(colnames(substitute), data$CompIlr$composition))) {
#'     stop(sprintf("The names of compositional variables must be the same
#'   in 'substitute' (%s) and 'composition' (%s).",
#'                  colnames(substitute),
#'                  data$CompIlr$composition))
#'   }
#' 
#'   tmp <- copy(data)
#' 
#'   # Compute between-person composition
#'   b <- tmp$CompIlr$BetweenComp
#'   b <- clo(b, total = 1440)
#'   b <- as.data.table(b)
#' 
#'   ID <- 1 # why re.form = NA but still needs this?
#'   min <- as.integer(minute)
#'   psi <- tmp$CompIlr$psi
#'   
#'   # list to store final output
#'   iout <- vector("list")
#' 
#'   # deal with covariates
#'   ilrn <- c(names(tmp$CompIlr$BetweenILR), names(tmp$CompIlr$WithinILR))
#'   vn <- do.call(rbind, find_predictors(tmp$BrmModel))
#'   
#'   if (isTRUE(identical(length(vn), length(ilrn)))) {
#' 
#'   for(i in colnames(substitute)) { # compostion level
#'     posub <- copy(substitute)
#'     posub <- as.data.table(posub)
#'     posub <- posub[(get(i) != 0)]
#'     posub <- posub[order(-rank(get(i)))]
#'     
#'     # add substitution variable name
#'     subvar <- colnames(posub) %snin% eval(i) # substitute compositional variables
#'     iv <- i # central compositional variable
#'     
#'     # lists to store results - TODO
#'     nd <- vector("list", length = nrow(b)) # a list to store a pairwise substitution for 1 specific period of all participant
#'     kout <- vector("list", length = nrow(substitute))
#'     jout <- vector("list", length = minute)
#'     
#'     for (j in seq_len(min)) { # time level
#'       sub <- posub * j
#'       for (k in seq_len(nrow(sub))) { # substitution level
#'         for (l in seq_len(nrow(b))) { # participant level
#'           newcomp <- b[l, ] + sub[k, ]
#'           names(newcomp) <- colnames(substitute)
#'           misub <- sub[k, get(i)]
#' 
#'           nd[[l]] <- cbind(b[l, ], newcomp, misub)
#'         }
#'         # Here is a dataset 1 possible substitution per 1 min for the entire sample
#'         # nrow = nrow(b), ncol = 10
#'         newd <- as.data.table(do.call(rbind, nd))
#'         
#'         # Add more useful information to the final output
#'         setnames(newd, "misub", "MinSubstituted")
#'         newd[, Substitute := rep(subvar, length.out = nrow(newd))[k]]
#'         newd$Predictor <- iv
#'         newd$MinSubstituted <- as.numeric(newd$MinSubstituted)
#'         
#'         ## remove impossible reallocation that result in negative values - TODO
#'         cols <- colnames(newd) %snin% c("MinSubstituted", "Substitute", "Predictor")
#'         
#'         noneg <- function(x){
#'           res <- ifelse(x < 0, NA, x)
#'           return(res)
#'         }
#'         newd[, (cols) := lapply(.SD, noneg), .SDcols = cols]
#'         newd <- newd[complete.cases(newd), ]
#'         
#'         # Add comp and ilr for predictions        
#'         bcomp <- acomp(newd[, colnames(tmp$CompIlr$BetweenComp), with = FALSE]) 
#'         tcomp <- acomp(newd[, tmp$CompIlr$composition, with = FALSE])
#'         
#'         bilr <- ilr(bcomp, V = psi) 
#'         tilr <- ilr(tcomp, V = psi) 
#'         wilr <- matrix(0, nrow = nrow(newd), ncol = ncol(bilr))
#'         wilr <- as.data.table(wilr)
#'         
#'         colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
#'         colnames(tilr) <- paste0("bilr", seq_len(ncol(tilr)))
#'         colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
#'         
#'         ## Substitution dataset
#'         subd <- cbind(newd, tilr, wilr, ID)
#'         
#'         ## No change dataset
#'         samed <- cbind(newd, bilr, wilr, ID)
#'         
#'         # Generate prediction
#'         ysub <- as.data.table(fitted(tmp$BrmModel, newdata = subd, re.form = NA, summary = FALSE))
#'         ysame <- as.data.table(fitted(tmp$BrmModel, newdata = samed, re.form = NA, summary = FALSE))
#' 
#'         # Difference between substitution and no change
#'         yd <- ysub - ysame
#'         myd <- rowMeans(yd) # ame of 1 substiution at participant level
#'         # nrow(yd) = 1, ncol = posterior draws
#'         # probably ame of covariates follows here
#'         myd <- as.data.table(describe_posterior(myd, centrality = "mean", 
#'                                                  ci = 0.95, ci_method = "eti"))
#'         myd <- myd[, .(Mean, CI_low, CI_high)]
#'         myd$MinSubstituted <- mean(newd[, MinSubstituted])
#' 
#'         kout[[k]] <- myd # ame at substitution level
#'       }
#'       # Save results
#'       jout[[j]] <- do.call(rbind, kout) 
#'       # nrow(myds) = nrow(posub), ncol() = 4
#'     }
#'     jout <- as.data.table(do.call(rbind, jout))
#'     jout[, Substitute := rep(subvar, length.out = nrow(jout))]
#'     jout$Predictor <- iv
#'     names(jout) <- c("Mean", "CI_low", "CI_high", "MinSubstituted", "Substitute", "Predictor")
#'     
#'     # Final results for entire composition
#'     iout[[i]] <- jout
#'   }
#'   return(iout)
#' 
#' } else {
#'   refg <- as.data.table(ref_grid(tmp$BrmModel) @grid)
#'   cv <- colnames(refg) %snin% c(ilrn, ".wgt.")
#'   
#'   refg <- refg[, cv, with = FALSE] # reference grid for covariates
#'   
#'   for(i in colnames(substitute)) { # compostion level
#'     posub <- copy(substitute)
#'     posub <- as.data.table(posub)
#'     posub <- posub[(get(i) != 0)]
#'     posub <- posub[order(-rank(get(i)))]
#'     
#'     # add substitution variable name
#'     subvar <- colnames(posub) %snin% eval(i) # substitute compositional variables
#'     iv <- i # central compositional variable
#'     
#'     # lists to store results - TODO
#'     nd <- vector("list", length = nrow(b)) # a list to store a pairwise substitution for 1 specific period of all participant
#'     kout <- vector("list", length = nrow(substitute))
#'     jout <- vector("list", length = minute)
#'     
#'     for (j in seq_len(min)) { # time level
#'       sub <- posub * j
#'       for (k in seq_len(nrow(sub))) { # substitution level
#'         for (l in seq_len(nrow(b))) { # participant level
#'           newcomp <- b[l, ] + sub[k, ]
#'           names(newcomp) <- colnames(substitute)
#'           misub <- sub[k, get(i)]
#' 
#'           nd[[l]] <- cbind(b[l, ], newcomp, misub)
#'         }
#'         # Here is a dataset 1 possible substitution per 1 min for the entire sample
#'         # nrow = nrow(b), ncol = 10
#'         newd <- as.data.table(do.call(rbind, nd))
#'         
#'         # Add more useful information to the final output
#'         setnames(newd, "misub", "MinSubstituted")
#'         newd[, Substitute := rep(subvar, length.out = nrow(newd))[k]]
#'         newd$Predictor <- iv
#'         newd$MinSubstituted <- as.numeric(newd$MinSubstituted)
#'         
#'         ## remove impossible reallocation that result in negative values - TODO
#'         cols <- colnames(newd) %snin% c("MinSubstituted", "Substitute", "Predictor")
#'         
#'         noneg <- function(x){
#'           res <- ifelse(x < 0, NA, x)
#'           return(res)
#'         }
#'         newd[, (cols) := lapply(.SD, noneg), .SDcols = cols]
#'         newd <- newd[complete.cases(newd), ]
#'         
#'         # Add comp and ilr for predictions        
#'         bcomp <- acomp(newd[, colnames(tmp$CompIlr$BetweenComp), with = FALSE]) 
#'         tcomp <- acomp(newd[, tmp$CompIlr$composition, with = FALSE])
#'         
#'         bilr <- ilr(bcomp, V = psi) 
#'         tilr <- ilr(tcomp, V = psi) 
#'         wilr <- matrix(0, nrow = nrow(newd), ncol = ncol(bilr))
#'         wilr <- as.data.table(wilr)
#'         
#'         colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
#'         colnames(tilr) <- paste0("bilr", seq_len(ncol(tilr)))
#'         colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
#'         
#'         ## Loop through covariates to get AME
#'         if (isFALSE(identical(nrow(refg, 1)))) {
#'           mout <- vector("list", length = nrow(refg))
#'           
#'           for(m in seq_len(nrow(refg))) {
#'           ## Substitution dataset
#'           subd <- cbind(newd, tilr, wilr, ID, refg[m, ])
#'           
#'           ## No change dataset
#'           samed <- cbind(newd, bilr, wilr, ID, refg[m, ])
#'           
#'           # Generate prediction
#'           ysub <- as.data.table(fitted(tmp$BrmModel, newdata = subd, re.form = NA, summary = FALSE))
#'           ysame <- as.data.table(fitted(tmp$BrmModel, newdata = samed, re.form = NA, summary = FALSE))
#' 
#'           # Difference between substitution and no change
#'           yd <- ysub - ysame
#'           myd <- rowMeans(yd) # ame of 1 substiution at participant level
#'           # nrow(yd) = posterior draws, ncol = 1
#'           if (nrow(refg) == 1) {
#'             myd <- myd
#'             } else {
#'               mout[[m]] <- myd
#'               myd <- Reduce(`+`, mout) / length(mout)
#'             }
#'           }
#'           
#'           # AME again across covariates and posteriors
#'           ame.yd <- Reduce(`+`, myd)/length(myd)
#'           # probably ame of covariates follows here
#'           ame.yd <- as.data.table(describe_posterior(ame.yd, centrality = "mean",
#'                                                  ci = 0.95, ci_method = "eti"))
#'           ame.yd <- ame.yd[, .(Mean, CI_low, CI_high)]
#'           ame.yd$MinSubstituted <- mean(newd[, MinSubstituted])
#' 
#'         kout[[k]] <- ame.yd # ame at substitution level
#'       }
#'       # Save results
#'       jout[[j]] <- do.call(rbind, kout)
#'       # nrow(jout) = nrow(posub), ncol() = 4
#'     }
#'     jout <- as.data.table(do.call(rbind, jout))
#'     jout[, Substitute := rep(subvar, length.out = nrow(jout))]
#'     jout$Predictor <- iv
#'     names(jout) <- c("Mean", "CI_low", "CI_high", "MinSubstituted", "Substitute", "Predictor")
#'     
#'     # Final results for entire composition
#'     iout[[i]] <- jout
#'   }
#'   return(iout)
#' }
#' }}