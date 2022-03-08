#' Marginal Effects of Between-person Substitution
#'
#' Estimate the difference in outcomes
#' when compositional variables are substituted for a specific time period
#' at between-person level.
#'
#' @param object A fitted \code{brms} model object. Required.
#' @param b A compisition or dataset of composition at between-peson level. Required.
#' @param comp A data frame or data table indicating the possible substitution of variables.
#' @param at Default is 1:60L.
#' @param names A vector representing names of variables in the substitution data.
#' @param substitute A vector representing names of substituted compositional component for each row of substitution dataset.
#'
#' @return
#' @importFrom data.table as.data.table
#' @importFrom compositions acomp ilr
#' @export
#' @examples
#' ## TODO
#' ## subd <- bsubmargin(object = m, b = b, comp = comp)
bsubmargin <- function (object, b, comp, at,
                        names, substitute) {

  if (isTRUE(missing(object))) { ## J to add 
    stop(paste("TODO"
    ))
  }

  if (isTRUE(missing(b))) {
    stop(paste(" 'b' is a required argument and cannot be missing;",
               "it should be a compositional object resulted from a acomp transformation",
               sep = "\n"))
  }

  if (isFALSE(missing(at))) {
    if (isFALSE(is.integer(at) && all(at > 0))) {
      stop("'at' must be integer values > 0.")
    }
  } else if (isTRUE(missing(at))) {
    at <- 1:60L
  }

  if (isFALSE(missing(names))) {
    if (isFALSE(is.character(names))) {
      stop(paste("'names' must be a character string representing the names of the following:",
                 "between-person compostion, total composition, at.",
                 sep = "\n"))
    }
  } else if (isTRUE(missing(names))) {
    names <- c("BTST", "BAWAKE", "BMVPA", "BLPA", "BSB", "TST", "AWAKE", "MVPA", "LPA", "SB", "at")
  }

  if (isTRUE(missing(substitute))) {
    substitute <- names[7:10]
  }

  # possible composition - comp
  comp <- comp[TST != 0] ## fix
  ID <- 1

  out <- NULL
  result <- NULL
  subd <- vector("list")

  for (i in at) {
    sub <- comp * i
    for (j in 1:nrow(comp)) {
      for (k in 1:nrow(b)) {
        newcomp <- b[k, ] + sub[j, ]
        subd[[k]] <- cbind(nodupb[k, ], newcomp, sleepsub[j, 1])
      }
      subd <- do.call(rbind, subd)
      subd <- as.data.table(subd)
      names(subd) <- names

      ## add substitute variable to indicate which variable was substituted for
      subd[, Substitute := rep(substitute, length.out = nrow(subd))][[j]]

      ## remove impossible reallocation that result in negative values
      subd <- subd %>%
        filter(if_all(-c(Substitute, at), ~ . > 0))

      ## add comp and ilr
      bcomp <- acomp(subd[, 1:5])
      tcomp <- acomp(subd[, 6:10])

      bilr <- ilr(bcomp, V = psi)
      tilr <- ilr(tcomp, V = psi)

      subd$bilr1 <- tilr[, 1]
      subd$bilr2 <- tilr[, 2]
      subd$bilr3 <- tilr[, 3]
      subd$bilr4 <- tilr[, 4]

      subd$wilr1 <- 0
      subd$wilr2 <- 0
      subd$wilr3 <- 0
      subd$wilr4 <- 0

      subd$ID <- ID

      ## dataset for no change
      samed <- data.table(bilr1 = bilr[, 1],
                          bilr2 = bilr[, 2],
                          bilr3 = bilr[, 3],
                          bilr4 = bilr[, 4],
                          wilr1 = 0,
                          wilr2 = 0,
                          wilr3 = 0,
                          wilr4 = 0)

      samed$ID <- ID

      # prediction
      ## substitution
      predsub <- as.data.table(fitted(m, newdata = subd, re.form = NA, summary = FALSE))

      ## no change
      predsame <- as.data.table(fitted(m, newdata = samed, re.form = NA, summary = FALSE))

      # calculate difference between substitution and no change
      preddif <- predsub - predsame
      preddif <- rowMeans(preddif)
      preddif <- as.data.table(describe_posterior(preddif, centrality = "mean", ci = 0.95, ci_method = "eti"))
      preddif <- preddif[, .(Mean, CI_low, CI_high)]

      margsub[[j]] <- preddif
    }

    # save results
    margsub.all[[i]] <- do.call(rbind, margsub)
  }
  out <- do.call(rbind, margsub.all)
  out <- cbind(out, subd[, c("at", "Substitute")])
  names(out) <- c("Mean", "CI_low", "CI_high", "at", "Substitute")
}
