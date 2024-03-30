#' Indices from a (dataset of) Multilevel Composition(s)
#'
#' Compute sets of compositions and log ratio transformation for multilevel compositional data
#'
#' @param data A \code{data.frame} or \code{data.table}
#' containing data of all variables used in the analysis. 
#' It must include a composition and a ID variable. Required.
#' @param sbp A signary matrix indicating sequential binary partition. Required.
#' @param parts A character vector specifying the names of compositional variables to be used.
#' @param transform A character value naming a log ratio transformation to be applied on compositional data.
#' Can be either \code{"ilr"} (isometric logratio), \code{"alr"} (additive logratio), or \code{"clr"} (centered logratio).
#' Default is \code{"ilr"}.
#' @param idvar A character string specifying the name of the variable containing IDs. 
#' Default is \code{"ID"}.
#' @param total A numeric value of the total amount to which the compositions should be closed.
#' Default is \code{1}.
#' 
#' @details 
#' The \emph{ilr}-transform maps the D-part compositional data from the simplex into non-overlapping 
#' subgroups in the (D−1)-dimension Euclidean space isometrically by using an orthonormal basis, 
#' thereby preserving the compo- sitional properties and yielding a full-rank covariance matrix.
#' \emph{ilr} transformation should be preferred. 
#' However, the \emph{alr} and \emph{clr} are alternatives.
#' The \emph{alr}-transform maps a D-part composition 
#' in the Aitchison-simplex non-isometrically to a 
#' D − 1-dimension Euclidian vectors, 
#' commonly treating the last part as the common denominator of the others.
#' \emph{alr} transformation does not rely on distance which breaks 
#' the constraint of compositional data.
#' \emph{clr}-transform maps a D-part composition in the Aitchison-simplex 
#' isometrically to a D-dimensional Euclidian vector subspace.
#' \emph{clr} transformation is not injetive, 
#' resulting in singular covariance matrices. 
#'
#' @return A \code{\link{complr}} object with twelve elements.
#'   \item{\code{BetweenComp}}{ A vector of class \code{acomp} representing one closed between-person composition
#'   or a matrix of class \code{acomp} representing multiple closed between-person compositions each in one row.}
#'   \item{\code{WithinComp}}{ A vector of class \code{acomp} representing one closed within-person composition
#'   or a matrix of class \code{acomp} representing multiple closed within-person compositions each in one row.}
#'   \item{\code{Comp}}{ A vector of class \code{acomp} representing one closed composition
#'   or a matrix of class \code{acomp} representing multiple closed  compositions each in one row.}
#'   \item{\code{BetweenILR}}{ Isometric log ratio transform of between-person composition.}
#'   \item{\code{WithinILR}}{ Isometric log ratio transform of within-person composition.}
#'   \item{\code{ILR}}{ Isometric log ratio transform of composition.}
#'   \item{\code{data}}{ The user's dataset or imputed dataset if the input data contains zeros.}
#'   \item{\code{psi}}{ A ILR matrix associated with user-defined partition structure.}
#'   \item{\code{sbp}}{ The user-defined sequential binary partition matrix.}
#'   \item{\code{parts}}{ Names of compositional variables.}
#'   \item{\code{idvar}}{ Name of the variable containing IDs.}
#'   \item{\code{total}}{ Total amount to which the compositions is closed.}
#' 
#' @importFrom compositions ilr alr clr acomp gsi.buildilrBase
#' @importFrom data.table copy as.data.table :=
#' 
#' @examples
#' cilr <- complr(data = mcompd, sbp = sbp,
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                 idvar = "ID", total = 1440)
#' str(cilr)
#' @export
complr <- function(data, parts, 
                   transform = "ilr",
                   sbp = NULL, 
                   total = 1, 
                   idvar = "ID") {
  
  if (isFALSE(inherits(data, c("data.table", "data.frame", "matrix")))) {
    stop("data must be a data table, data frame or matrix.")
  }
  
  tmp <- as.data.table(data)
  
  # check NAs
  if (isTRUE(any(apply(tmp[, parts, with = FALSE], 2, function(x) any(is.na(x)))))) {
    stop(paste(
      "This dataset of composition contains missing data;",
      "  Missind data hinder the application of compositional data analysis",
      "  because the analysis is based on log-ratios",
      "  Please deal with missing data before running 'complr'.",
      sep = "\n"))
  }
  
  # check 0s
  if (isTRUE(any(apply(tmp[, parts, with = FALSE], 2, function(x) x == 0)))) {
    stop(paste(
      "This dataset of composition contains zero(s);",
      "  Zeros hinder the application of compositional data analysis",
      "  because the analysis is based on log-ratios",
      "  Please deal with zeros before running 'complr'.",
      sep = "\n"))
  }
  
  # allow one transform at a time
  if (isFALSE(length(transform) == 1)) {
    stop("only one type of transforms can be done at a time.")
  }
  
  # check transform
  if (isFALSE(transform %in% c("ilr", "alr", "clr"))) {
    stop(" 'transform' should be one of the following: 'ilr', 'alr', 'clr'")
  }
  
  ## Make Composition ---------------
  # combined
  tcomp <- acomp(tmp[, parts, with = FALSE], total = total)
  # between-person
  for (v in parts) {
    tmp[, (v) := mean(get(v), na.rm = TRUE), by = eval(idvar)]
  }
  bcomp <- acomp(tmp[, parts, with = FALSE], total = total)
  
  # within-person 
  wcomp <- tcomp - bcomp
  
  # name them for later use
  colnames(bcomp) <- paste0("B", parts)
  colnames(wcomp) <- paste0("W", parts)
  colnames(tcomp) <- parts
  
  # ILR ---------------
  if (identical(transform, "ilr")) {
    if (isTRUE(missing(sbp))) {
      stop(" 'sbp', i.e., sequential binary partition, is required for ilr transform.")
    }
    if (isFALSE(inherits(sbp, "matrix"))) {
      stop(sprintf("sbp is a '%s' but must be a matrix.",
                   paste(class(sbp), collapse = " ")))
    }
    if (isTRUE(any(apply(sbp, 2, function(x) x %nin% c(-1, 0, 1))))) {
      stop("sbp should only contain 1, -1 and 0 (a partition)")
    }
    if (isFALSE(identical(length(parts), ncol(sbp)))) {
      stop(sprintf(
        "The number of compositional variables in parts (%d) 
  must be the same as in sbp (%d).",
  length(parts),
  ncol(sbp)))
    }
    psi <- gsi.buildilrBase(t(sbp))
    
    # combined
    tilr <- ilr(tcomp, V = psi)
    
    # between-person
    bilr <- ilr(bcomp, V = psi)
    
    # within-person 
    wilr <- ilr(wcomp, V = psi)
    
    # name them for later use
    colnames(bilr)  <- paste0("bilr", seq_len(ncol(bilr)))
    colnames(wilr)  <- paste0("wilr", seq_len(ncol(wilr)))
    colnames(tilr)  <- paste0("ilr", seq_len(ncol(tilr)))
    
    if (any(c(colnames(bilr), colnames(wilr), colnames(tilr)) %in% colnames(tmp))) {
      stop(
        paste(
          "'data' should not have any column names starting with 'bilr', 'wilr', or 'ilr';",
          "  these variables will be used in subsequent models.",
          "  Please rename them before running 'complr'.",
          sep = "\n"))
    }
    
    out <- structure(
      list(
        Comp = tcomp,
        BetweenComp = bcomp,
        WithinComp = wcomp,
        BetweenILR = bilr,
        WithinILR = wilr,
        ILR = tilr,
        data = data,
        transform = "ilr",
        psi = psi,
        sbp = sbp,
        parts = parts,
        idvar = idvar,
        total = total),
      class = "complr"
    )
  }
  
  # ALR ---------------
  if (identical(transform, "alr")) {
    
    # combined
    talr <- alr(tcomp)
    
    # between-person
    balr <- alr(bcomp)
    
    # within-person 
    walr <- alr(wcomp)
    
    # name them for later use
    colnames(balr)  <- paste0("balr", seq_len(ncol(balr)))
    colnames(walr)  <- paste0("walr", seq_len(ncol(walr)))
    colnames(talr)  <- paste0("alr", seq_len(ncol(talr)))
    
    if (any(c(colnames(balr), colnames(walr), colnames(talr)) %in% colnames(tmp))) {
      stop(
        paste(
          "'data' should not have any column names starting with 'balr', 'walr', or 'alr';",
          "  these variables will be used in subsequent models.",
          "  Please rename them before running 'complr'.",
          sep = "\n"))
    }
    
    out <- structure(
      list(
        Comp = tcomp,
        BetweenComp = bcomp,
        WithinComp = wcomp,
        ALR = talr,
        BetweenALR = balr,
        WithinALR = walr,
        data = data,
        transform = "alr",
        parts = parts,
        idvar = idvar,
        total = total),
      class = "complr"
    )
  }
  
  # CLR ---------------
  if (identical(transform, "clr")) {
    
    # combined
    tclr <- clr(tcomp)
    
    # between-person
    bclr <- clr(bcomp)
    
    # within-person 
    wclr <- clr(wcomp)
    
    # name them for later use
    colnames(bclr)  <- paste0("bclr", seq_len(ncol(bclr)))
    colnames(wclr)  <- paste0("wclr", seq_len(ncol(wclr)))
    colnames(tclr)  <- paste0("clr", seq_len(ncol(tclr)))
    
    if (any(c(colnames(bclr), colnames(wclr), colnames(tclr)) %in% colnames(tmp))) {
      stop(
        paste(
          "'data' should not have any column names starting with 'bclr', 'wclr', or 'clr';",
          "  these variables will be used in subsequent models.",
          "  Please rename them before running 'complr'.",
          sep = "\n"))
    }
    
    out <- structure(
      list(
        Comp = tcomp,
        BetweenComp = bcomp,
        WithinComp = wcomp,
        CLR = tclr,
        BetweenCLR = bclr,
        WithinCLR = wclr,
        data = data,
        transform = "clr",
        parts = parts,
        idvar = idvar,
        total = total),
      class = "complr"
    )
  }
  out
}
