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
  
  ## CHECK NUMBER OF COMPOSITION HERE?
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
  
  for (idx in seq_along(parts)) {
    parts_i <- parts[[idx]]
    total_i <- total[[idx]]
    
    if (length(parts) == 1) {
      sbp_i   <- sbp
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
      sbp_i   <- sbp[[idx]]
    }
    
    # check NAs
    if (isTRUE(any(apply(tmp[, parts_i, with = FALSE], 2, function(x)
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
    if (isTRUE(any(apply(tmp[, parts_i, with = FALSE], 2, function(x)
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
      if (is.null(sbp_i)) {
        # build default sbp
        message(
          " A sequential binary partition (sbp), is required for ilr transform but is not supplied.
 A default sbp, which is a pivot balance, will be applied."
        )
        sbp_i <- build.sbp(parts = parts_i)
      }
      if (isFALSE(inherits(sbp_i, "matrix"))) {
        message(sprintf(
          "sbp is a '%s' but must be a matrix.",
          paste(class(sbp_i), collapse = " ")
        ))
        sbp_i <- as.matrix(sbp_i)
      }
      if (isTRUE(any(apply(sbp_i, 2, function(x)
        x %nin% c(-1, 0, 1))))) {
        stop("sbp should only contain 1, -1 and 0 (a partition)")
      }
      if (isFALSE(identical(length(parts_i), ncol(sbp_i)))) {
        stop(
          sprintf(
            "The number of compositional variables in parts (%d)
  must be the same as in sbp (%d).",
            length(parts_i),
            ncol(sbp_i)
          )
        )
      }
      psi_i <- gsi.buildilrBase(t(sbp_i))
    } else {
      psi_i <- sbp_i <- NULL
    }
    
    # MAKE COMPOSITION AND LOG RATIO TRANSFORMATIONS ----------------
    if (shape == "wide") {
      # make composition
      tX_i <- acomp(tmp[, parts_i, with = FALSE], total = total_i)
      bX_i <- wX_i <- NULL
      colnames(tX_i) <- paste0("t", parts_i)
      
      # ILR
      if (identical(transform, "ilr")) {
        tilr_i <- ilr(tX_i, V = psi_i)
        bilr_i <- wilr_i <- NULL
        colnames(tilr_i)  <- paste0("z", seq_len(ncol(tilr_i)), "_", idx)
        
      } else if (identical(transform, "alr")) {
        talr_i <- alr(tX_i)
        balr_i <- walr_i <- NULL
        colnames(talr_i)  <- paste0("z", seq_len(ncol(talr_i)), "_", idx)
        
      } else if (identical(transform, "clr")) {
        tclr_i <- clr(tX_i)
        bclr_i <- wclr_i <- NULL
        colnames(tclr_i)  <- paste0("z", seq_len(ncol(tclr_i)), "_", idx)
      }
    }
    
    if (shape == "long") {
      # make composition
      # combined
      tX_i <- acomp(tmp[, parts_i, with = FALSE], total = total_i)
      # between-person
      for (v in parts_i) {
        tmp[, (v) := mean(get(v), na.rm = TRUE), by = eval(idvar)]
      }
      bX_i <- acomp(tmp[, parts_i, with = FALSE], total = total_i)
      
      # within-person
      wX_i <- tX_i - bX_i
      
      # name them for later use
      colnames(bX_i) <- paste0("b", parts_i)
      colnames(wX_i) <- paste0("w", parts_i)
      colnames(tX_i) <- paste0("t", parts_i)
      
      ## ILR ---------------
      if (identical(transform, "ilr")) {
        tilr_i <- ilr(tX_i, V = psi_i)
        bilr_i <- ilr(bX_i, V = psi_i)
        wilr_i <- ilr(wX_i, V = psi_i)
        
        colnames(bilr_i)  <- paste0("bz", seq_len(ncol(bilr_i)), "_", idx)
        colnames(wilr_i)  <- paste0("wz", seq_len(ncol(wilr_i)), "_", idx)
        colnames(tilr_i)  <- paste0("z", seq_len(ncol(tilr_i)), "_", idx)
      }
      
      ## ALR ---------------
      if (identical(transform, "alr")) {
        talr_i <- alr(tX_i)
        balr_i <- alr(bX_i)
        walr_i <- alr(wX_i)
        
        colnames(balr_i)  <- paste0("bz", seq_len(ncol(balr_i)), "_", idx)
        colnames(walr_i)  <- paste0("bz", seq_len(ncol(walr_i)), "_", idx)
        colnames(talr_i)  <- paste0("z", seq_len(ncol(talr_i)), "_", idx)
      }
      
      ## CLR ---------------
      if (identical(transform, "clr")) {
        tclr_i <- clr(tX_i)
        bclr_i <- clr(bX_i)
        wclr_i <- clr(wX_i)
        
        colnames(bclr_i)  <- paste0("bz", seq_len(ncol(bclr_i)), "_", idx)
        colnames(wclr_i)  <- paste0("bz", seq_len(ncol(wclr_i)), "_", idx)
        colnames(tclr_i)  <- paste0("z", seq_len(ncol(tclr_i)), "_", idx)
      }
    }
    
    Z_i <-  if (exists("tilr_i")) (tilr_i)
    else if (exists("talr_i")) (talr_i)
    else if (exists("tclr_i")) (tclr_i)
    else  (NULL)
    
    bZ_i <- if (exists("bilr_i")) (bilr_i)
    else if (exists("balr_i")) (balr_i)
    else if (exists("bclr_i")) (bclr_i)
    else (NULL)
    
    wZ_i <- if (exists("wilr_i")) (wilr_i)
    else if (exists("walr_i")) (walr_i)
    else if (exists("wclr_i")) (wclr_i)
    else (NULL)
    
    # cbind data output
    dataout_i <- cbind(tX_i, bX_i, wX_i, Z_i, bZ_i, wZ_i)
    
    output[[idx]] <- list(
      X  = if(exists("tX_i")) (tX_i) else (NULL),
      bX = if(exists("bX_i")) (bX_i) else (NULL),
      wX = if(exists("wX_i")) (wX_i) else (NULL),
      
      Z  = if(exists("Z_i")) (Z_i) else (NULL),
      bZ = if(exists("bZ_i")) (bZ_i) else (NULL),
      wZ = if(exists("wZ_i")) (wZ_i) else (NULL),
      
      dataout = dataout_i,
      parts = parts_i,
      total = total_i,
      sbp = sbp_i,
      psi = psi_i
    )
  }
  
  # PATCH OUTPUT ----------------
  dataout <- do.call(cbind, lapply(output, function(x) x$dataout))
  
  # check any repetitive names between tmp and dataout before cbind
  if (isTRUE(any(colnames(tmp) %in% colnames(dataout)))) {
    stop(
      sprintf(
        "'data' cannot have any column names the same with logratio variables",
        "  Please ensure that the names of the columns in 'data' are not any of the following:",
        "  %s",
        paste(colnames(tmp)[colnames(tmp) %in% colnames(dataout)], collapse = ", "),
        sep = "\n"
      )
    )
  }
  structure(
    list(
      output    = output,
      datain    = as.data.table(tmp),
      dataout   = as.data.table(cbind(tmp, dataout)),
      transform = transform,
      idvar     = idvar
    ),
    class = "complr")
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
