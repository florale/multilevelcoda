#' Indices from a (dataset of) Multilevel Composition(s)
#'
#' Compute sets of compositions and log ratio transformation for multilevel compositional data
#'
#' @param data A \code{data.frame} or \code{data.table}
#' containing data of all variables used in the analysis. 
#' It must include a composition and a ID variable. Required.
#' @param transform A character value naming a log ratio transformation to be applied on compositional data.
#' Can be either \code{"ilr"} (isometric logratio), \code{"alr"} (additive logratio), or \code{"clr"} (centered logratio).
#' Default is \code{"ilr"}.
#' @param parts A character vector specifying the names of compositional variables to be used.
#' @param sbp A signary matrix indicating sequential binary partition.
#' @param total A numeric value of the total amount to which the compositions should be closed.
#' Default is \code{1}.
#' @param idvar Only for multilevel data, a character string specifying the name of the variable containing IDs.
#' 
#' @details 
#' The \emph{ilr}-transform maps the D-part compositional data from the simplex into non-overlapping 
#' subgroups in the (D-1)-dimension Euclidean space isometrically by using an orthonormal basis, 
#' thereby preserving the compositional properties and yielding a full-rank covariance matrix.
#' \emph{ilr} transformation should be preferred. 
#' However, the \emph{alr} and \emph{clr} are alternatives.
#' The \emph{alr}-transform maps a D-part composition 
#' in the Aitchison-simplex non-isometrically to a 
#' (D-1)-dimension Euclidian vectors, 
#' commonly treating the last part as the common denominator of the others.
#' \emph{alr} transformation does not rely on distance which breaks 
#' the constraint of compositional data.
#' \emph{clr}-transform maps a D-part composition in the Aitchison-simplex 
#' isometrically to a D-dimensional Euclidian vector subspace.
#' \emph{clr} transformation is not injetive, 
#' resulting in singular covariance matrices. 
#'
#' @return A \code{\link{complr}} object with at least the following elements.
#'   \item{\code{comp}}{ A vector of class \code{acomp} representing one closed composition
#'   or a matrix of class \code{acomp} representing multiple closed  compositions each in one row.}
#'   \item{\code{between_comp}}{ A vector of class \code{acomp} representing one closed between-person composition
#'   or a matrix of class \code{acomp} representing multiple closed between-person compositions each in one row.}
#'   \item{\code{within_comp}}{ A vector of class \code{acomp} representing one closed within-person composition
#'   or a matrix of class \code{acomp} representing multiple closed within-person compositions each in one row.}
#'   \item{\code{logratio}}{ Log ratio transform of composition.}
#'   \item{\code{between_logratio}}{ Log ratio transform of between-person composition.}
#'   \item{\code{within_logratio}}{ Log ratio transform of within-person composition.}
#'   \item{\code{data}}{ The user's dataset or imputed dataset if the input data contains zeros.}
#'   \item{\code{transform}}{ Type of transform applied on compositional data.}
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
#' 
#' calr <- complr(data = mcompd, sbp = sbp,
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                 idvar = "ID",
#'                 transform = "alr")
#' str(calr)
#' 
#' cclr <- complr(data = mcompd, sbp = sbp,
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                 idvar = "ID",
#'                  transform = "clr")
#' str(cclr)
#' 
#' cilr_wide <- complr(data = mcompd[!duplicated(ID)], sbp = sbp,
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' str(cilr_wide)
#' @export
complr <- function(data,
                   parts,
                   sbp = NULL, 
                   total = 1, 
                   idvar = NULL,
                   transform = "ilr"
) {
  
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
  
  # check single level or multilevel
  if (is.null(idvar)) {
    shape <- "wide"
  } else {
    shape <- "long"
  }

  # allow one transform at a time
  if (length(transform) > 1) {
    stop("only one type of transforms can be done at a time.")
  }
  
  # check transform
  if (isFALSE(any(transform %in% c("ilr", "alr", "clr")))) {
    stop(" 'transform' should be one of the following: \"ilr\", \"alr\", \"clr\".")
  }
  
  # specific for ilr
  if (identical(transform, "ilr")) {
    if (missing(sbp)) {
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
  } else {
    psi <- sbp <- NULL
  }
  
  # check var names
  if (isTRUE(any(grep("ilr|alr|clr", colnames(tmp))))) {
    stop(
      paste(
        "'data' should not have any column names with patterns of \"ilr\", \"alr\", or \"clr\";",
        "  these variables will be computed by 'complr' used in subsequent models.",
        "  Please remove or rename them before running 'complr'.",
        sep = "\n"))
  }
  
  # composition and lr
  if (shape == "wide") {
    # make composition
    tcomp <- acomp(tmp[, parts, with = FALSE], total = total)
    bcomp <- wcomp <- NULL
    colnames(tcomp) <- parts
    
    # ILR
    if (identical(transform, "ilr")) {
      tilr <- ilr(tcomp, V = psi)
      bilr <- wilr <- NULL
      colnames(tilr)  <- paste0("ilr", seq_len(ncol(tilr)))
      
    } else if (identical(transform, "alr")) {
      talr <- alr(tcomp)
      balr <- walr <- NULL
      colnames(talr)  <- paste0("alr", seq_len(ncol(talr)))
      
    } else if (identical(transform, "clr")) {
      tclr <- clr(tcomp)
      bclr <- wclr <- NULL
      colnames(tclr)  <- paste0("clr", seq_len(ncol(tclr)))
    }
  }
  
  if (shape == "long") {
    # make composition
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
    colnames(bcomp) <- paste0("b", parts)
    colnames(wcomp) <- paste0("w", parts)
    colnames(tcomp) <- parts
    
    # ILR ---------------
    if (identical(transform, "ilr")) {
      
      tilr <- ilr(tcomp, V = psi)
      bilr <- ilr(bcomp, V = psi)
      wilr <- ilr(wcomp, V = psi)
      
      colnames(bilr)  <- paste0("bilr", seq_len(ncol(bilr)))
      colnames(wilr)  <- paste0("wilr", seq_len(ncol(wilr)))
      colnames(tilr)  <- paste0("ilr", seq_len(ncol(tilr)))
    }
    
    # ALR 
    if (identical(transform, "alr")) {
      
      talr <- alr(tcomp)
      balr <- alr(bcomp)
      walr <- alr(wcomp)
      
      colnames(balr)  <- paste0("balr", seq_len(ncol(balr)))
      colnames(walr)  <- paste0("walr", seq_len(ncol(walr)))
      colnames(talr)  <- paste0("alr", seq_len(ncol(talr)))
    }
    
    # CLR 
    if (identical(transform, "clr")) {
      
      tclr <- clr(tcomp)
      bclr <- clr(bcomp)
      wclr <- clr(wcomp)
      
      colnames(bclr)  <- paste0("bclr", seq_len(ncol(bclr)))
      colnames(wclr)  <- paste0("wclr", seq_len(ncol(wclr)))
      colnames(tclr)  <- paste0("clr", seq_len(ncol(tclr)))
    }
  }
  logratio <-  if (exists("tilr")) (tilr)
  else if (exists("talr")) (talr)
  else if (exists("tclr")) (tclr)
  
  between_logratio <- if (exists("bilr")) (bilr)
  else if (exists("balr")) (balr)
  else if (exists("bclr")) (bclr)
  else (NULL)
  
  within_logratio <- if (exists("wilr")) (wilr)
  else if (exists("walr")) (walr)
  else if (exists("wclr")) (wclr) 
  else (NULL)
  
  out <- structure(
    list(
      comp = tcomp,
      between_comp = bcomp,
      within_comp = wcomp,
      logratio = logratio,
      between_logratio = between_logratio,
      within_logratio = within_logratio,
      data = as.data.table(data),
      transform = transform,
      psi = if(exists("psi")) (psi) else (NULL),
      sbp = if(exists("sbp")) (sbp) else (NULL),
      parts = parts,
      idvar = idvar,
      total = total
    ),
    class = "complr"
  )
  out
}

#' Indices from a (dataset of) Multilevel Composition(s) (deprecated.)
#'
#' @param ... arguments passed to \code{\link{complr}}.
#' @seealso \code{\link{complr}}
#' 
#' @inherit complr return
#'
#' @export
compilr <- function(...) {
  warning("'compilr' is deprecated. Please use 'complr' instead.")
  UseMethod("complr")
}
