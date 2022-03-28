#' A series of functions for single level compositional models
#'
#' The \code{\link{subcoda}} function wraps a series of 
#' The workflows is as follow:
#' It firstly compute composition and ILR coordinates, 
#' 
#' It fit a Bayesian single-level compositional model to compositional predictors
#' and users' choices of outcomes and other covariates.
#'
#' @param data A composition or dataset of composition. Required.
#' @param composition A string character indicating the names of compositional variables in `data`.
#' @param sbp A signary matrix indicating sequential binary partition. Required.
#' @param ... Further arguments passed to \code{\link{brm}}.
#' @return 
#' @importFrom compositions acomp ilr gsi.buildilrBase
#' @importFrom data.table copy as.data.table :=
#' @importFrom zCompositions lrEM
#' @importFrom brms brm
#' @export
#' @examples
#'
#' data(mcompd)
#' data(sbp)
#' 
#' davg <- mcompd[, .(SB = mean(SB, na.rm = TRUE),
#'                    LPA = mean(LPA, na.rm = TRUE),
#'                    MVPA = mean(MVPA, na.rm = TRUE),
#'                    TST = mean(TST, na.rm = TRUE),
#'                    WAKE = mean(WAKE, na.rm = TRUE),
#'                    TRESS = mean(STRESS, na.rm = TRUE),
#'                    Age = mean (Age, na.rm = TRUE),
#'                    Female = na.omit(Female)[1]),
#'                    by = .(ID)]
#'
#' scoda <- subcoda(data = davg, composition = c("TST", "WAKE", "MVPA", "LPA", "SB"), 
#'                        sbp = sbp, minute = 10, substitute = posubtest,
#'                        formula = STRESS ~ ilr1 + ilr2 + ilr3 + ilr4)
#'  
#'  ## cleanup
#'  rm(mcompd, sbp, davg, scoda)
subcoda <- function(data, composition, sbp, formula, minute, substitute, ...) {

    if (isFALSE(inherits(data, c("data.table", "data.frame", "matrix")))) {
    stop("data must be a data table, data frame or matrix.")
  }
  if (isFALSE(inherits(sbp, "matrix"))) {
    stop(sprintf("sbp is a '%s' but must be a matrix.",
                 paste(class(sbp), collapse = ";")))
  }
  if (isFALSE(identical(length(composition), ncol(sbp)))) {
    stop(sprintf("The number of compositional variables in composition (%d) 
                 must be the same as in sbp (%d).",
                 length(composition),
                 ncol(sbp)))
  }

  # compilr
  psi <- gsi.buildilrBase(t(sbp))
  
  comp <- acomp(data[, composition, with = FALSE])
  ilr <- ilr(comp, V = psi)
  
  colnames(ilr) <- paste0("ilr", seq_len(ncol(ilr)))
  
  tmpd <- cbind(data, ilr)
  
  # brm model
  m <- brm(eval(formula),
           data = tmpd,
           ...)
  
  # Substitution model
  # compute compositional mean
  mcomp <- mean(comp, robust = TRUE)
  mcomp <- clo(mcomp, total = 1440)
  mcomp <- as.data.table(t(mcomp))
  names(mcomp) <- paste0("M", names(mcomp))
  
  # generate input for substitution model
  min <- as.integer(paste0(minute))

  out <- list()
  for(i in colnames(substitute)) {
  posub <- copy(substitute)
  posub <- as.data.table(posub)
  posub <- posub[(get(i) != 0)]
  posub <- posub[order(-rank(get(i)))]

  # Get substitution variable name for substitution model
  subvar <- colnames(posub) %snin% eval(i)
  iv <- i

  # lists to store results - TODO
  nd <- NULL
  newd <- vector("list")

  # substitution dataset
    for (j in seq_len(min)) {
      sub <- posub * j
      for (k in seq_len(nrow(sub))) {
        new <- mcomp + sub[k, ]
        names(new) <- paste0(names(substitute))
        newd[[k]] <- cbind(mcomp, new, sub[k, ][[i]])
        }
        nd[[j]] <- do.call(rbind, newd)
        }
        newd <- as.data.table(do.call(rbind, nd))

  # add names
  colnames(newd)[ncol(newd)] <- "MinSubstituted"
  newd[, Substitute := rep(subvar, length.out = nrow(newd))]
  newd$Predictor <- iv

  # ## remove impossible reallocation that result in negative values - TODO
  cols <- colnames(newd) %snin% c("MinSubstituted", "Substitute", "Predictor")

  noneg <- function(x){
    res <- ifelse(x < 0, NA, x)
    return(res)
  }

  newd[, (cols) := lapply(.SD, noneg), .SDcols = cols]
  newd <- newd[complete.cases(newd), ]

  ## add comp and ilr
  oldvar <- colnames(newd) %sin% names(mcomp)
  newvar <- colnames(newd) %sin% composition

  oldcomp <- acomp(newd[, oldvar, with = FALSE])
  newcomp <- acomp(newd[, newvar, with = FALSE])

  oldilr <- ilr(oldcomp, V = psi)
  newilr <- ilr(newcomp, V = psi)

  names(oldilr) <- c(paste0("ilr", 1:ncol(ilr)))
  names(newilr) <- c(paste0("ilr", 1:ncol(ilr)))

  ## substitution dataset
  subd <- cbind(newd, newilr)

  ## no change dataset
  samed <- cbind(newd, oldilr)

  # prediction
  ## substitution
  ysub <- as.data.table(fitted(m, newdata = subd, re.form = NA, summary = FALSE))

  ## no change
  ysame <- as.data.table(fitted(m, newdata = samed, re.form = NA, summary = FALSE))

  # difference between substitution and no change
  ydiff <- ysub - ysame
  ydiff <- as.data.table(describe_posterior(ydiff, centrality = "mean",
                                              ci = 0.95, ci_method = "eti"))
  ydiff <- ydiff[, .(Mean, CI_low, CI_high)]

  # save results
  result <- do.call(cbind, ydiff)
  result <- cbind(result, newd[, c("MinSubstituted", "Substitute", "Predictor")])
  result <- as.data.table(result)
  names(result) <- c("Mean", "CI_low", "CI_high", "MinSubstituted", "Substitute", "Predictor")

  ## final results for entire composition
  out[[i]] <- result
  
  }
  
  allout <- list(SubstitutionResults = out,
                      BrmsResults = m,
                      ILR = ilr,
                      Comp = comp)
  
  return(allout)
}
