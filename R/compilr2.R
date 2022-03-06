#' Compute Between-person, Within-person, and Total Composition and Isometric log ratio transform of a (dataset of) composition(s)
#' 
#' This function is designed to help calculate sets of compositions and IRLs
#' for \code{brms} Multilevel Compositional Data models
#'  
#' @param b A compisition or dataset of composition at between-person level. Required.
#' @param t A composition or dataset of composition. Required.
#' @param sbp A signary matrix indicating sequential binary partition. Required.
#' 
#' @return A list with six elements.
#' \itemize{
#'   \item{\code{BetweenComp}}{ A vector of class \code{acomp} representing one closed between-person composition 
#'   or a matrix of class \code{acomp} representing multiple closed between-person compositions each in one row.}
#'   \item{\code{WithinComp}}{ A vector of class \code{acomp}} representing one closed within-person composition 
#'   or a matrix of class \code{acomp} representing multiple closed within-person compositions each in one row.}
#'   \item{\code{TotalComp}}{ A vector of class \code{acomp}} representing one closed total composition 
#'   or a matrix of class \code{acomp} representing multiple closed total compositions each in one row.}
#'   \item{\code{BetweenILR}}{ Isometric log ratio transform of between-person composition.}
#'   \item{\code{WithinILR}}{ Isometric log ratio transform of within-person composition.}
#'   \item{\code{TotalILR}}{ Isometric log ratio transform of total composition.}
#' }
#' 
#' @importFrom compositions ilr
#' @importFrom compositions acomp
#' @importFrom compositions gsi.buildilrBase
#' @importFrom multilevelTools meanDeviations
#' @export
#' @examples
#' 
#' 
compilr <- function(t, sbp) {
  
  bw <- apply(t, 2, meanDeviations) # need to edit
  
  b <- as.data.table(cbind(bw[[1]][[1]], bw[[2]][[1]], bw[[3]][[1]], bw[[4]][[1]], bw[[5]][[1]]))
  b <- b[rep(seq_len(nrow(b)), nrow(t)), ]
  
  psi <- gsi.buildilrBase(t(sbp))
  
  ## Between-person composition
  bcomp <- acomp(b)
  bilr <- ilr(bcomp, V=psi)
  
  ## Total composition
  tcomp <- acomp(t)
  tilr <- ilr(tcomp, V=psi)
  
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