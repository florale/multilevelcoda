#' Compute Between-person, Within-person, and Total Composition and Isometric log ratio transform of a (dataset of) composition(s)
#'
#' This function is designed to help calculate sets of compositions and IRLs
#' for Multilevel Compositional Data models
#'
#' @param data A composition or dataset of composition. Required.
#' @param sbp A signary matrix indicating sequential binary partition. Required.
#' @param idvar A character string indicating the name of the variable containing IDs.
#'
#' @return A list with six elements.
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
#' }
#' @importFrom compositions ilr acomp gsi.buildilrBase
#' @importFrom extraoperators %snin%
#' @importFrom data.table copy as.data.table :=
#' @export
#' @examples
#'
#' data(mcompd)
#' test <- compilr(data = mcompd[, 1:6], sbp = sbp, idvar = "ID")
#' str(test)
#' ## cleanup
#' rm(test, mcompd)
#' ##TODO - check for 0 in data
compilr <- function(data, sbp, idvar = "ID") {
  if (isFALSE(inherits(data, c("data.table", "data.frame", "matrix")))) {
    stop("data must be a data table, data frame or matrix.")
  }
  if (isFALSE(inherits(sbp, "matrix"))) {
    stop(sprintf("sbp is a '%s' but must be a matrix.",
                 paste(class(sbp), collapse = ";")))
  }
  if (isFALSE(identical(ncol(data) - 1L, ncol(sbp)))) {
    stop(sprintf("The number of columns in data (%d) must be the same as in sbp (%d).",
                 ncol(data),
                 ncol(sbp)))
  }

  b <- copy(data)
  b <- as.data.table(b)
  vn <- colnames(b) %snin% idvar
  for (v in vn) {
    b[, (v) := mean(get(v), na.rm = TRUE), by = eval(idvar)]
  }
  b <- b[, vn, with = FALSE]

  psi <- gsi.buildilrBase(t(sbp))

  ## Between-person composition
  bcomp <- acomp(b)
  bilr <- ilr(bcomp, V=psi)

  ## Total composition
  tcomp <- acomp(data[, vn, with = FALSE])
  tilr <- ilr(tcomp, V = psi)

  ## Within-person composition
  wcomp <- tcomp - bcomp
  wilr <- ilr(wcomp, V=psi)

  out <- list(
    BetweenComp = bcomp,
    WithinComp = wcomp,
    TotalComp = tcomp,
    BetweenILR = bilr,
    WithinILR = wilr,
    TotalILR = tilr)

  return(out)
}
