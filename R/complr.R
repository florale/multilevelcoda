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
#' For multiple compositions, a list of character vectors.
#' @param sbp A signary matrix indicating sequential binary partition when \code{transform = "ilr"}.
#' If not supplied, a default sequential binary partition (sbp) will be built using function \code{\link{build.sbp}}.
#' For multiple compositions, a list of sbps can be supplied.
#' @param total A numeric value of the total amount to which the compositions should be closed.
#' For multiple compositions, a list of numeric values.
#' Default is \code{1}.
#' @param idvar Only for multilevel data, a character string specifying the name of the variable containing participants IDs.
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
#'   \item{\code{X}}{ A vector of class \code{acomp} representing one closed composition
#'   or a matrix of class \code{acomp} representing multiple closed  compositions each in one row.}
#'   \item{\code{bX}}{ A vector of class \code{acomp} representing one closed between-person composition
#'   or a matrix of class \code{acomp} representing multiple closed between-person compositions each in one row.}
#'   \item{\code{wX}}{ A vector of class \code{acomp} representing one closed within-person composition
#'   or a matrix of class \code{acomp} representing multiple closed within-person compositions each in one row.}
#'   \item{\code{Z}}{ Log ratio transform of composition.}
#'   \item{\code{bZ}}{ Log ratio transform of between-person composition.}
#'   \item{\code{wZ}}{ Log ratio transform of within-person composition.}
#'   \item{\code{data}}{ The user's dataset or imputed dataset if the iiut data contains zeros.}
#'   \item{\code{transform}}{ Type of transform applied on compositional data.}
#'   \item{\code{parts}}{ Names of compositional variables.}
#'   \item{\code{idvar}}{ Name of the variable containing IDs.}
#'   \item{\code{total}}{ Total amount to which the compositions is closed.}
#'
#' @importFrom compositions ilr alr clr acomp gsi.buildilrBase
#' @importFrom data.table copy as.data.table :=
#' @importFrom extraoperators %nin%
#'
#' @examples
#' x1 <- complr(data = mcompd,
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                 idvar = "ID", total = 1440)
#' str(x1)
#'
#' x2 <- complr(data = mcompd,
#'                 parts = list(c("TST", "WAKE"), c("MVPA", "LPA", "SB")),
#'                 total = list(c(480), c(960)),
#'                 idvar = "ID",
#'                 transform = "ilr")
#' str(x2)
#'
#' x3 <- complr(data = mcompd, sbp = sbp,
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
#'                 idvar = "ID",
#'                  transform = "ilr")
#' str(x3)
#'
#' x_wide <- complr(data = mcompd[!duplicated(ID)], sbp = sbp,
#'                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"))
#' str(x_wide)
#' @export
complr <- function(data,
                   parts,
                   sbp = NULL,
                   total = 1,
                   idvar = NULL,
                   transform = "ilr") {
  if (isFALSE(inherits(data, c("data.table", "data.frame", "matrix")))) {
    stop("data must be a data table, data frame or matrix.")
  }
  
  tmp <- as.data.table(data)
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
  
  # check number of composition
  # check if parts is a list
  if (is.list(parts)) {
    if (length(parts) == 0) {
      stop("parts cannot be an empty list.")
    }
    if (isTRUE(any(sapply(parts, function(x)
      ! is.character(x))))) {
      stop("parts should be a character vector or a list of character vectors.")
    }
  } else if (is.character(parts)) {
    parts <- list(parts)
  } else {
    stop("parts should be a character vector or a list of character vectors.")
  }
  
  # loop through list to compute composition and lr
  output <- vector("list", length = length(parts))
  
  for (nx in seq_along(parts)) {
    partsx <- parts[[nx]]
    totalx <- total[[nx]]
    
    if (length(parts) == 1) {
      sbpx <- if (is.list(sbp))
        sbp[[1]]
      else
        sbp
    }
    else {
      if (length(parts) != length(total)) {
        stop("parts and total should have the same length.")
      }
      if (isFALSE(is.null(sbp))) {
        if (length(parts) != length(sbp)) {
          stop("parts and sbp should have the same length.")
        }
      }
      sbpx <- sbp[[nx]]
    }
    
    # check NAs
    if (isTRUE(any(apply(tmp[, partsx, with = FALSE], 2, function(x)
      any(is.na(
        x
      )))))) {
      stop(
        paste(
          "This dataset of composition contains missing data;",
          "  Missind data hinder the application of compositional data analysis",
          "  because the analysis is based on log-ratios",
          "  Please deal with missing data before running 'complr'.",
          sep = "\n"
        )
      )
    }
    
    # check 0s
    if (isTRUE(any(apply(tmp[, partsx, with = FALSE], 2, function(x)
      x == 0)))) {
      stop(
        paste(
          "This dataset of composition contains zero(s);",
          "  Zeros hinder the application of compositional data analysis",
          "  because the analysis is based on log-ratios",
          "  Please deal with zeros before running 'complr'.",
          sep = "\n"
        )
      )
    }
    
    # specific for ilr
    if (identical(transform, "ilr")) {
      if (is.null(sbpx)) {
        # build default sbp
        message(
          " A sequential binary partition (sbp), is required for ilr transform but is not supplied.
 A default sbp, which is a pivot balance, will be applied."
        )
        sbpx <- build.sbp(parts = partsx)
      }
      if (isFALSE(inherits(sbpx, "matrix"))) {
        message(sprintf(
          "sbp is a '%s' but must be a matrix.",
          paste(class(sbpx), collapse = " ")
        ))
        sbpx <- as.matrix(sbpx)
      }
      if (isTRUE(any(apply(sbpx, 2, function(x)
        x %nin% c(-1, 0, 1))))) {
        stop("sbp should only contain 1, -1 and 0 (a partition)")
      }
      if (isFALSE(identical(length(partsx), ncol(sbpx)))) {
        stop(
          sprintf(
            "The number of compositional variables in parts (%d)
             must be the same as in sbp (%d).",
            length(partsx),
            ncol(sbpx)
          )
        )
      }
      psix <- gsi.buildilrBase(t(sbpx))
    } else {
      psix <- sbpx <- NULL
    }
    
    # MAKE COMPOSITION AND LOG RATIO TRANSFORMATIONS ----------------
    if (shape == "wide") {
      # make composition
      tXx <- acomp(tmp[, partsx, with = FALSE], total = totalx)
      bXx <- wXx <- NULL
      colnames(tXx) <- paste0("t", partsx)
      
      # ILR
      if (identical(transform, "ilr")) {
        tilrx <- ilr(tXx, V = psix)
        bilrx <- wilrx <- NULL
        colnames(tilrx)  <- paste0("z", seq_len(ncol(tilrx)), "_", nx)
        
      } else if (identical(transform, "alr")) {
        talrx <- alr(tXx)
        balrx <- walrx <- NULL
        colnames(talrx)  <- paste0("z", seq_len(ncol(talrx)), "_", nx)
        
      } else if (identical(transform, "clr")) {
        tclrx <- clr(tXx)
        bclrx <- wclrx <- NULL
        colnames(tclrx)  <- paste0("z", seq_len(ncol(tclrx)), "_", nx)
      }
    }
    
    if (shape == "long") {
      # make composition
      # combined
      tXx <- acomp(tmp[, partsx, with = FALSE], total = totalx)
      
      # between-person
      for (v in partsx) {
        tmp[, paste0("b", v) := mean(get(v), na.rm = TRUE), by = eval(idvar)]
      }
      bXx <- acomp(tmp[, colnames(tmp) %in% paste0("b", partsx), with = FALSE], total = totalx)
      
      # within-person (notes unclass(x)/unclass(y))
      wXx <- tXx - bXx
      
      # name them for later use
      colnames(bXx) <- paste0("b", partsx)
      colnames(wXx) <- paste0("w", partsx)
      colnames(tXx) <- paste0("t", partsx)
      
      ## ILR ---------------
      if (identical(transform, "ilr")) {
        tilrx <- ilr(tXx, V = psix)
        bilrx <- ilr(bXx, V = psix)
        wilrx <- ilr(wXx, V = psix)
        
        colnames(bilrx)  <- paste0("bz", seq_len(ncol(bilrx)), "_", nx)
        colnames(wilrx)  <- paste0("wz", seq_len(ncol(wilrx)), "_", nx)
        colnames(tilrx)  <- paste0("z", seq_len(ncol(tilrx)), "_", nx)
      }
      
      ## ALR ---------------
      if (identical(transform, "alr")) {
        talrx <- alr(tXx)
        balrx <- alr(bXx)
        walrx <- alr(wXx)
        
        colnames(balrx)  <- paste0("bz", seq_len(ncol(balrx)), "_", nx)
        colnames(walrx)  <- paste0("bz", seq_len(ncol(walrx)), "_", nx)
        colnames(talrx)  <- paste0("z", seq_len(ncol(talrx)), "_", nx)
      }
      
      ## CLR ---------------
      if (identical(transform, "clr")) {
        tclrx <- clr(tXx)
        bclrx <- clr(bXx)
        wclrx <- clr(wXx)
        
        colnames(bclrx)  <- paste0("bz", seq_len(ncol(bclrx)), "_", nx)
        colnames(wclrx)  <- paste0("bz", seq_len(ncol(wclrx)), "_", nx)
        colnames(tclrx)  <- paste0("z", seq_len(ncol(tclrx)), "_", nx)
      }
    }
    
    Zx <-  if (exists("tilrx"))
      (tilrx)
    else if (exists("talrx"))
      (talrx)
    else if (exists("tclrx"))
      (tclrx)
    else
      NULL
    
    bZx <- if (exists("bilrx"))
      (bilrx)
    else if (exists("balrx"))
      (balrx)
    else if (exists("bclrx"))
      (bclrx)
    else
      NULL
    
    wZx <- if (exists("wilrx"))
      (wilrx)
    else if (exists("walrx"))
      (walrx)
    else if (exists("wclrx"))
      (wclrx)
    else
      NULL
    
    # cbind data output
    dataoutx <- cbind(tXx, bXx, wXx, Zx, bZx, wZx)
    
    output[[nx]] <- list(
      X  = if (exists("tXx")) tXx else NULL,
      bX = if (exists("bXx")) bXx else NULL,
      wX = if (exists("wXx")) wXx else NULL,
      Z  = if (exists("Zx")) Zx else NULL,
      bZ = if (exists("bZx")) bZx else NULL,
      wZ = if (exists("wZx")) wZx else NULL,
      
      dataout = dataoutx,
      parts = partsx,
      total = totalx,
      sbp = sbpx,
      psi = psix
    )
  }
  
  # PATCH OUTPUT ----------------
  dataout <- do.call(cbind, lapply(output, function(x)
    x$dataout))
  
  # check any repetitive names between tmp and dataout before cbind
  if (isTRUE(any(colnames(data) %in% colnames(dataout)))) {
    stop(
      sprintf(
        "'data' cannot have any column names the same as logratio variables. 
   Please ensure that the names of the columns in 'data' are not any of the following:
   %s",
        paste(colnames(data)[colnames(data) %in% colnames(dataout)], collapse = ", ")
      )
    )
  }
  structure(
    list(
      output    = output,
      datain    = as.data.table(data),
      dataout   = if (nrow(dataout) == nrow(data))
        cbind(as.data.table(data), dataout)
      else
        cbind(as.data.table(tmp[colnames(data)]), dataout),
      transform = transform,
      idvar     = idvar
    ),
    class = "complr"
  )
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
