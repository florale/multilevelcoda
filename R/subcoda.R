#' A series of functions for single level models
#'
#' This provides a few sentence description about the example function.
#'
#' @param data A composition or dataset of composition. Required.
#' @param composition A string character indicating the names of compositional variables in `data`.
#' @param sbp A signary matrix indicating sequential binary partition. Required.
#' @param ... Further arguments passed to \code{\link{brm}}.
#' @param 
#' @return 
#' @importFrom 
#' @export
#' @examples
#'
#' # just use total comp for this
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
#' testsubcoda <- subcoda(data = davg, composition = c("TST", "WAKE", "MVPA", "LPA", "SB"), sbp = sbp,
#'  formula = STRESS ~ ilr1 + ilr2 + ilr3 + ilr4, minute = 10, substitute = posubtest)
#'  
#'  ## cleanup
#'  rm(mcompd, sbp, davg, testsubcoda)

subcoda <- function(data, composition, sbp, formula, minute, substitute, ...) {
  
  # compilr
  psi <- gsi.buildilrBase(t(sbp))
  
  compn <- colnames(data) %sin% names(composition)
  
  comp <- acomp(data[, composition, with = FALSE])
  ilr <- ilr(comp, V = psi)
  
  names(ilr) <- c(paste0("ilr", 1:ncol(ilr)))
  
  copyd <- cbind(data, ilr)
  
  # brm model
  m <- brm(eval(formula),
           data = copyd,
           ...)
  
  # Substitution model
  # compute compositional mean
  mcomp <- mean(comp)
  mcomp <- clo(mcomp, total = 1440)
  mcomp <- as.data.table(t(mcomp))
  names(mcomp) <- paste0("M", names(mcomp))
  
  # generate input for substitution model
  vn <- colnames(substitute)
  min <- as.integer(paste0(minute))

  allout <- list()
  for(i in vn) {
  posub <- copy(substitute)
  posub <- as.data.table(posub)
  posub <- posub[(get(i) != 0)]
  posub <- posub[order(-rank(get(i)))]

  # Get substitution variable name for substitution model
  subvar <- colnames(posub) %snin% eval(i)
  iv <- i

  # lists to store results - TODO
  result <- NULL
  new <- vector('list')
  subd <- vector("list")

  # substitution dataset
  for (j in 1:min) {
    sub <- posub * j
    for (k in 1:nrow(posub)) {
      new <- mcomp + sub[k, ]
      names(new) <- paste0(names(substitute))
      subd[[k]] <- cbind(mcomp, new, sub[k, ][[i]])
    }
    result[[j]] <- do.call(rbind, subd)
  }
  subd <- as.data.table(do.call(rbind, result))

  # add names
  colnames(subd)[ncol(subd)] <- "MinSubstituted"
  subd[, Substitute := rep(subvar, length.out = nrow(subd))]
  subd$Predictor <- iv

  # ## remove impossible reallocation that result in negative values - TODO
  cols <- colnames(subd) %snin% c("MinSubstituted", "Substitute", "Predictor")

  noneg <- function(x){
    res <- ifelse(x < 0, NA, x)
    return(res)
  }

  subd[, (cols) := lapply(.SD, noneg), .SDcols = cols]
  subd <- subd[complete.cases(subd), ]

  ## add comp and ilr
  oldvar <- colnames(subd) %sin% names(mcomp)
  newvar <- colnames(subd) %sin% composition

  oldcomp <- acomp(subd[, oldvar, with = FALSE])
  newcomp <- acomp(subd[, newvar, with = FALSE])

  oldilr <- ilr(oldcomp, V = psi)
  newilr <- ilr(newcomp, V = psi)

  names(oldilr) <- c(paste0("ilr", 1:ncol(ilr)))
  names(newilr) <- c(paste0("ilr", 1:ncol(ilr)))

  ## substitution dataset
  subd <- cbind(subd, newilr)

  ## no change dataset
  samed <- cbind(oldilr)

  # prediction
  ## substitution
  predsub <- as.data.table(fitted(m, newdata = subd, re.form = NA, summary = FALSE))

  ## no change
  predsame <- as.data.table(fitted(m, newdata = samed, re.form = NA, summary = FALSE))

  # difference between substitution and no change
  preddif <- predsub - predsame
  preddif <- as.data.table(describe_posterior(preddif, centrality = "mean",
                                              ci = 0.95, ci_method = "eti"))
  preddif <- preddif[, .(Mean, CI_low, CI_high)]

  # save results
  out <- do.call(cbind, preddif)
  out <- cbind(out, subd[, c("MinSubstituted", "Substitute", "Predictor")])
  out <- as.data.table(out)
  names(out) <- c("Mean", "CI_low", "CI_high", "MinSubstituted", "Substitute", "Predictor")

  ## final results for entire composition
  allout[[i]] <- out
  
  }
  
  finalresult <- list(SubstitutionResults = allout,
                      BrmsResults = m,
                      ILR = ilr,
                      Composition = comp)
  
  return(finalresult)
}
