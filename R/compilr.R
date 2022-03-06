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
#'
#' @export
#' @examples
#' 
#' 
compilr <- function(b , t, sbp) {
  
  if(isTRUE(missing(b))){
    stop(paste( "'b is a required argument and cannot be missing,",
                "it should a be data table or data frame representing a between-person composition.",
                "See ? or the website articles (vignettes) for details.",
                sep = "\n"))
  }
  
  if(isTRUE(missing(t))){
    
    stop(paste( "'t' is a required argument and cannot be missing,",
                "it should a be data table or data frame representing a composition.",
                "See ? or the website articles (vignettes) for details.",
                sep = "\n"))
  }
  
  if(isTRUE(missing(sbp))){
    
    stop(paste( "'sbp' is a required argument and cannot be missing,",
                "it should be a matrix indicating sequential binary partition.",
                "See ? or the website articles (vignettes) for details.",
                sep = "\n"))
  }
  
  ## Between-person composition
  bcomp <- acomp(b)
  
  ## make composition ilr
  psi <- gsi.buildilrBase(t(sbp))
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

