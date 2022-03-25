#' Compute Between-person, Within-person, and Total Composition and Isometric log ratio transform of a (dataset of) composition(s)
#'
#' This function is designed to help calculate sets of compositions and IRLs
#' for Multilevel Compositional Data models
#'
#' @param data A \code{data.frame} or \code{data.table}
#' containing data of all variables used in the analysis. 
#' It must include a composition, a ID variable, and others. Required.
#' @param sbp A signary matrix indicating sequential binary partition. Required.
#' @param composition A character vector specifying the names of compositional variables. Required.
#' @param idvar A character string specifying the name of the variable containing IDs. Default is \code{ID.}
#'
#' @return A list with eleven elements.
#' \itemize{
#'   \item{\code{BetweenComp}}{ A vector of class \code{acomp} representing one closed between-person composition
#'   or a matrix of class \code{acomp} representing multiple closed between-person compositions each in one row.}
#'   \item{\code{WithinComp}}{ A vector of class \code{acomp} representing one closed within-person composition
#'   or a matrix of class \code{acomp} representing multiple closed within-person compositions each in one row.}
#'   \item{\code{TotalComp}}{ A vector of class \code{acomp} representing one closed total composition
#'   or a matrix of class \code{acomp} representing multiple closed total compositions each in one row.}
#'   \item{\code{BetweenILR}}{ Isometric log ratio transform of between-person composition.}
#'   \item{\code{WithinILR}}{ Isometric log ratio transform of within-person composition.}
#'   \item{\code{TotalILR}}{ Isometric log ratio transform of total composition.}
#'   \item{\code{data}}{ The user's dataset or imputed dataset if the input data contains zeros.}
#'   \item{\code{psi}}{ A ILR matrix associated with user-defined partition structure.}
#'   \item{\code{sbp}}{ The user-defined sequential binary partition matrix.}
#'   \item{\code{composition}}{ Names of compositional variables.}
#'   \item{\code{idvar}}{ Name of the variable containing IDs.}
#' }
#' @importFrom compositions ilr acomp gsi.buildilrBase
#' @importFrom data.table copy as.data.table :=
#' @importFrom zCompositions lrEM
#' @export
#' @examples
#' ## Example 1 - Dataset with no 0
#' data(mcompd)
#' data(sbp)
#' cilr1 <- compilr(data = mcompd, sbp = sbp, 
#'                  composition = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#' 
#' ## Example 2 - Dataset with 0s
#' ## Impute a 0 in 'mcompd'
#' mcompd[3, 1] <- 0
#' cilr2 <- compilr(data = mcompd, sbp = sbp, 
#'                  composition = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#'
#' str(cilr1)
#' 
#' ## cleanup
#' rm(cilr1, cilr2, mcompd, sbp)
compilr <- function(data, sbp, composition, idvar = "ID") {
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

  tmp <- copy(data)
  tmp <- as.data.table(tmp)
  
  psi <- gsi.buildilrBase(t(sbp))
  
  ## 0 imputation
  if (isTRUE(any(apply(tmp[, composition, with = FALSE], 2, function(x) x == 0)))) {
    message(paste("This dataset of composition contains zero(s);",
                  "It is now imputed using the Log-ratio EM algorithm.",
                  "For more details, please see ?zCompositions::lrEM",
                  "If you would like to deal with zeros using other methods,",
                  "please do so before using 'compilr'.",
                  sep = "\n"))
    
    dl1 <- rep(1440, length(composition))
    impd <- lrEM(tmp[, composition, with = FALSE], label = 0, dl = dl1, ini.cov = "multRepl")
    names(impd) <- composition
    
    tmp <- as.data.table(cbind(impd, tmp[, !composition, with = FALSE]))
    
  }
  
  ## Total composition
  tcomp <- acomp(tmp[, composition, with = FALSE])
  tilr <- ilr(tcomp, V = psi)
  
  ## Between-person composition
  for (v in composition) {
    tmp[, (v) := mean(get(v), na.rm = TRUE), by = eval(idvar)]
  }
  
  bcomp <- acomp(tmp[, composition, with = FALSE])
  bilr <- ilr(bcomp, V = psi)

  ## Within-person composition
  wcomp <- tcomp - bcomp
  wilr <- ilr(wcomp, V = psi)

  colnames(bcomp) <- paste0("B", composition)
  colnames(wcomp) <- paste0("W", composition)
  colnames(tcomp) <- composition
  
  colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
  colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
  colnames(tilr) <- paste0("ilr", seq_len(ncol(tilr)))
  
  if(any(c(colnames(bilr), colnames(wilr), colnames(tilr))) %in% colnames(data)) {
    stop(sprintf("data should not have any column names starting with 'bilr', 'wilr', or 'ilr'."))
  }
  
  out <- list(
    BetweenComp = bcomp,
    WithinComp = wcomp,
    TotalComp = tcomp,
    BetweenILR = bilr,
    WithinILR = wilr,
    TotalILR = tilr,
    data = tmp,
    psi = psi,
    sbp = sbp,
    composition = composition,
    idvar = idvar)

  return(out)
}
