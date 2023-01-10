#' Compute useful indices from a (dataset of) multilevel composition(s)
#'
#' Computes sets of compositions and IRLs for Multilevel Compositional Data models. 
#'
#' @param data A \code{data.frame} or \code{data.table}
#' containing data of all variables used in the analysis. 
#' It must include a composition and a ID variable. Required.
#' @param sbp A signary matrix indicating sequential binary partition. Required.
#' @param parts A character vector specifying the names of compositional variables to be used.
#' @param idvar A character string specifying the name of the variable containing IDs. 
#' Default to \code{ID}.
#' @param total A numeric value of the total amount to which the compositions should be closed.
#' Default to \code{1440}.
#'
#' @return A list with twelve elements.
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
#'   \item{\code{parts}}{ Names of compositional variables.}
#'   \item{\code{idvar}}{ Name of the variable containing IDs.}
#'   \item{\code{total}}{ Total amount to which the compositions is closed.}
#' }
#' 
#' @importFrom compositions ilr acomp gsi.buildilrBase
#' @importFrom data.table copy as.data.table :=
#' @importFrom zCompositions lrEM
#' @export
#' @examples
#' data(mcompd)
#' data(sbp)
#' cilr1 <- compilr(data = mcompd, sbp = sbp, 
#'                  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")
#' str(cilr1)
#' 
#' ## cleanup
#' rm(cilr1, mcompd, sbp)
compilr <- function(data, sbp, parts, total = 1440, idvar = "ID") {
  if (isFALSE(inherits(data, c("data.table", "data.frame", "matrix")))) {
    stop("data must be a data table, data frame or matrix.")
  }
  if (isFALSE(inherits(sbp, "matrix"))) {
    stop(sprintf("sbp is a '%s' but must be a matrix.",
                 paste(class(sbp), collapse = " ")))
  }
  if (isFALSE(identical(length(parts), ncol(sbp)))) {
    stop(sprintf(
    "The number of compositional variables in parts (%d) 
  must be the same as in sbp (%d).",
                 length(parts),
                 ncol(sbp)))
  }
  tmp <- as.data.table(data)
  psi <- gsi.buildilrBase(t(sbp))
  
  ## Composition and ILRs
  # total
  tcomp <- acomp(tmp[, parts, with = FALSE])
  tilr <- ilr(tcomp, V = psi)
  
  # between-person
  for (v in parts) {
    tmp[, (v) := mean(get(v), na.rm = TRUE), by = eval(idvar)]
    }
  bcomp <- acomp(tmp[, parts, with = FALSE])
  bilr <- ilr(bcomp, V = psi)

  # within-person 
  wcomp <- tcomp - bcomp
  wilr <- ilr(wcomp, V = psi)

  # name them for later use
  colnames(bcomp) <- paste0("B", parts)
  colnames(wcomp) <- paste0("W", parts)
  colnames(tcomp) <- parts
  colnames(bilr) <- paste0("bilr", seq_len(ncol(bilr)))
  colnames(wilr) <- paste0("wilr", seq_len(ncol(wilr)))
  colnames(tilr) <- paste0("ilr", seq_len(ncol(tilr)))
  
  if(any(c(colnames(bilr), colnames(wilr), colnames(tilr)) %in% colnames(tmp))) {
    stop(paste(
      "'data' should not have any column names starting with 'bilr', 'wilr', or 'ilr';",
      "  these variables will be used in subsequent models.",
      "  Please rename them before running 'compilr'.",
               sep = "\n"))
  }
  out <- structure(
    list(
      BetweenComp = bcomp,
      WithinComp = wcomp,
      TotalComp = tcomp,
      BetweenILR = bilr,
      WithinILR = wilr,
      TotalILR = tilr,
      data = data,
      psi = psi,
      sbp = sbp,
      parts = parts,
      idvar = idvar,
      total = total),
    class = "compilr")
  out
}