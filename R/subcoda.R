#' @title A series of functions for single level compositional models
#'
#' @description 
#' This function first computes composition and ILR coordinates,
#' then fits a Bayesian single-level compositional model to compositional predictors
#' and users' choices of outcomes and covariates,
#' and finally estimates the average marginal effects of
#' isotemporal compositional substitution model.
#'
#' @param data A \code{data.frame} or \code{data.table}
#' containing data of all variables used in the analysis. 
#' @param sbp A signary matrix indicating sequential binary partition.
#' @param parts A character vector specifying the names of compositional variables.
#' @param substitute A \code{data.frame} or \code{data.table} of the possible substitution of variables.
#' This dataset can be computed using function \code{possub}.
#' @param formula A object of class \code{formula}, \code{brmsformula}:
#' A symbolic description of the model to be fitted. 
#' @param minute A integer or numeric value indicating the maximum minute for which substitution model is desired.
#' Default to \code{60L} (i.e., the model loops through 1:60L minutes).
#' @param total A numeric value of the total amount to which the compositions should be closed. Default to \code{1440}.
#' @param ... Further arguments passed to \code{\link{brm}}.
#' 
#' @return A list with seven elements.
#' \itemize{
#'   \item{\code{BrmsModel}}{ An object of class \code{brmsfit}, which contains the posterior draws 
#'   along with many other useful information about the model.}
#'   \item{\code{SubstitutionResults}}{ A list containing the result of isotemporal substitution model.}
#'   \item{\code{ILR}}{ Isometric log ratio transform of composition.}
#'   \item{\code{Comp}}{A vector of class \code{acomp} representing one composition
#'   or a matrix of class \code{acomp} representing multiple closed compositions each in one row.}
#'   \item{\code{data}}{ The user's dataset or imputed dataset if the input data contains zeros.}
#'   \item{\code{sbp}}{ The user-defined sequential binary partition matrix.}
#'   \item{\code{parts}}{ Names of compositional variables.}
#'   
#' @importFrom compositions acomp ilr gsi.buildilrBase
#' @importFrom data.table copy as.data.table :=
#' @importFrom zCompositions lrEM
#' @importFrom brms brm
#' @importFrom bayestestR describe_posterior
#' @importFrom stats fitted
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
#'                    STRESS = mean(STRESS, na.rm = TRUE),
#'                    Age = mean (Age, na.rm = TRUE),
#'                    Female = na.omit(Female)[1]),
#'                    by = .(ID)]
#' ps <- possub(parts = c("TST", "WAKE", "MVPA", "LPA", "SB")
#' 
#' scoda <- subcoda(data = davg, parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), 
#'                        sbp = sbp, minute = 10, substitute = ps,
#'                        formula = STRESS ~ ilr1 + ilr2 + ilr3 + ilr4)
#'  
#'  rm(mcompd, sbp, davg, scoda)
subcoda <- function(data, sbp, parts, substitute, formula, 
                    minute = 60, total = 1440, ...) {
  
  if (isFALSE(inherits(data, c("data.table", "data.frame", "matrix")))) {
    stop("data must be a data table, data frame or matrix.")
  }
  if (isFALSE(inherits(sbp, "matrix"))) {
    stop(sprintf("sbp is a '%s' but must be a matrix.",
                 paste(class(sbp), collapse = ";")))
  }
  if (isFALSE(identical(length(parts), ncol(sbp)))) {
    stop(sprintf("The number of compositional variables in parts (%d) 
                 must be the same as in sbp (%d).",
                 length(parts),
                 ncol(sbp)))
  }
  tmp <- as.data.table(data)
  
  # Compilr
  # deal with 0s - imputation
  if (isTRUE(any(apply(tmp[, parts, with = FALSE], 2, function(x) x == 0)))) {
    message(paste("This dataset of composition contains zero(s);",
                  "It is now imputed using the Log-ratio EM algorithm.",
                  "For more details, please see ?zCompositions::lrEM",
                  "If you would like to deal with zeros using other methods,",
                  "please do so before running 'compilr'.",
                  sep = "\n"))
    
    dl1 <- rep(eval(total), length(parts))
    impd <- lrEM(tmp[, parts, with = FALSE], label = 0, dl = dl1, ini.cov = "multRepl")
    names(impd) <- parts
    tmp <- as.data.table(cbind(impd, tmp[, !parts, with = FALSE]))
  }
  
  psi <- gsi.buildilrBase(t(sbp))
  comp <- acomp(tmp[, parts, with = FALSE])
  ilr <- ilr(comp, V = psi)
  colnames(ilr) <- paste0("ilr", seq_len(ncol(ilr)))
  
  if(any(colnames(ilr) %in% colnames(tmp))) {
    stop(paste("data should not have any column names starting with 'ilr';",
               "these variables will be used in subsequent models.",
               "Please rename them before running 'subcoda'.",
               sep = "\n"))
  }
    
  # Brm model
  tmpd <- cbind(tmp, ilr)
  m <- brm(eval(formula), data = tmpd,
           ...)
  
  # Substitution model
  # original composition
  b <- comp
  b <- as.data.table(clo(b, total = total))
  
  min <- as.integer(minute)
  
  # Model for no change
  samed <- tmpd
  ysame <- as.data.table(fitted(m, newdata = samed, re.form = NA, summary = FALSE))
  ysame <- rowMeans(ysame)
  
  # Substitution model
  iout <- vector("list")
  for(i in colnames(substitute)) {
    posub <- copy(substitute)
    posub <- as.data.table(posub)
    posub <- posub[(get(i) != 0)]
    posub <- posub[order(-rank(get(i)))]
    
    subvar <- colnames(posub) %snin% eval(i)
    iv <- i
    
    kout <- vector("list", length = nrow(posub))
    jout <- vector("list", length = min)
    
    for (j in seq_len(min)) {
      sub <- posub * j
      for (k in seq_len(nrow(sub))) {
        subk <- sub[k, ]
        subk <- subk[rep(seq_len(nrow(subk)), nrow(b)), ]
        newcomp <- b + subk
        colnames(newcomp) <- paste0("New_", parts)
        MinSubstituted <- subk[, get(i)]
        kout[[k]] <- cbind(newcomp, tmp, MinSubstituted)
      }
      jout[[j]] <- do.call(rbind, kout)
    }
    newd <- as.data.table(do.call(rbind, jout))
    
    # useful information for the final results
    newd[, Substitute := rep(subvar, length.out = nrow(newd))]
    newd$Predictor <- iv
    
    # remove impossible reallocation that result in negative values 
    cols <- colnames(newd) %snin% c("MinSubstituted", "Substitute", "Predictor")
    newd <- newd[rowSums(newd[, ..cols] < 0) == 0]
    
    # comp and ilr
    newcomp <- acomp(newd[, colnames(newcomp), with = FALSE])
    newilr <- ilr(newcomp, V = psi)
    names(newilr) <- c(paste0("ilr", 1:ncol(ilr)))
    
    # prediction
    ## substitution
    subd <- cbind(newd, newilr)
    ysub <- as.data.table(fitted(m, newdata = subd, re.form = NA, summary = FALSE))
    ysub <- cbind(subd[, .(MinSubstituted, Substitute, Predictor)], as.data.table(t(ysub)))
    ysub <- ysub[, lapply(.SD, mean), by = c("MinSubstituted", "Substitute", "Predictor")]
    
    suppl <- ysub[, .(MinSubstituted, Substitute, Predictor)]
    ysub <- ysub[, -c(1:3)]
    
    # difference between substitution and no change
    ydiff <- t(ysub) - ysame
    # posterior means and intervals
    ymean <- apply(ydiff, 2, function(x) {describe_posterior(x, centrality = "mean", ...)})
    ymean <- rbindlist(ymean)
    ymean <- ymean[, .(Mean, CI_low, CI_high)]
    ymean <- cbind(ymean, suppl)
    
    ## final results for entire composition
    iout[[i]] <- ymean
  }
  out <- list(BrmsModel = m,
              SubstitutionResults = iout,
              ILR = ilr,
              Comp = comp,
              data = tmp,
              sbp = sbp,
              parts = parts)
  out
}
