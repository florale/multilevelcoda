#' A series of functions for single level models
#'
#' This provides a few sentence description about the example function.
#'
#' @param x A vector.
#' @return The vector coerced to numeric vector in a data table.
#' @importFrom data.table data.table
#' @export
#' @examples
#'
#' # just use total comp for this
#' data(mcompd)
#'  compilrtest <- compilr(data = mcompd[, 1:6], sbp = sbp, idvar = "ID")
#'  bilr <- compilrtest$BetweenILR
#'  wilr <- compilrtest$WithinILR
#'  names(bilr) <- c(paste0("bilr", 1:ncol(bilr)))
#'  names(wilr) <- c(paste0("wilr", 1:ncol(wilr)))
#'  copyd <- cbind(mcompd, bilr, wilr)
#' coda.lmer <- lmer(STRESS ~ bilr1 + bilr2 + bilr3 + bilr4 + wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID), data = copyd)

subcoda <- function(data, sbp, formula, min, substitute, idvar) {
  b <- compilrtest$BetweenComp
  # Compute compositional mean
  mcomp <- rowMeans(b)
  mcomp <- clo(mcomp, total = 1440)
  mcomp <- as.data.table(t(mcomp))
  names(mcomp) <- paste0("B", names(mcomp))

  # generate input for substitution model
  ID <- 1
  vn <- colnames(substitute)
  min <- as.integer(paste0(minute))
  psi <- data$CompIlr$psi
  posub <- copy(substitute)
  posub <- as.data.table(posub)
  posub <- posub[(get(i) != 0)]
  posub <- posub[order(-rank(get(i)))]

  # Get substitution variable name for substitution model
  subvar <- colnames(posub) %snin% eval(i)
  iv <- i

  # lists to store results - TODO
  result <- NULL
  newcomp <- vector('list')
  subd <- vector("list")

  # substitution dataset
  for (j in 1:min) {
    sub <- posub * j
    for (k in 1:nrow(posub)) {
      newcomp <- mcomp + sub[k, ]
      names(newcomp) <- paste0(names(substitute))
      subd[[k]] <- cbind(mcomp, newcomp, sub[k, ][[i]])
    }
    result[[j]] <- do.call(rbind, subd)
  }
  subd <- as.data.table(do.call(rbind, result))

  # add names
  colnames(subd)[ncol(subd)] <- "MinSubstituted"
  subd[, Substitute := rep(subvar, length.out = nrow(subd))]
  subd$Predictor <- iv

  # ## remove impossible reallocation that result in negative values - TODO
  cols <- colnames(subd) %snin% c("MinSubstituted", "Substitute", "Predictor")

  noneg <- function(x){
    res <- ifelse(x < 0, NA, x)
    return(res)
  }

  subd[, (cols) := lapply(.SD, noneg), .SDcols = cols]
  subd <- subd[complete.cases(subd), ]

  ## add comp and ilr

  bn <- colnames(subd) %sin% names(mcomp)
  tn <- colnames(subd) %sin% names(newcomp)

  bcomp <- acomp(subd[, bn, with = FALSE])
  tcomp <- acomp(subd[, tn, with = FALSE])

  bilr <- ilr(bcomp, V = psi)
  tilr <- ilr(tcomp, V = psi)
  wilr <- matrix(0, nrow = nrow(subd), ncol = ncol(bilr))
  wilr <- as.data.table(wilr)

  names(bilr) <- c(paste0("bilr", 1:ncol(bilr)))
  names(tilr) <- c(paste0("bilr", 1:ncol(tilr)))
  names(wilr) <- c(paste0("wilr", 1:ncol(wilr)))

  ## substitution dataset
  subd <- cbind(subd, tilr, wilr)
  subd$ID <- ID

  ## no change dataset
  samed <- cbind(bilr, wilr)
  samed$ID <- ID

  # prediction
  ## substitution
  predsub <- as.data.table(fitted(data$Results, newdata = subd, re.form = NA, summary = FALSE))

  ## no change
  predsame <- as.data.table(fitted(data$Results, newdata = samed, re.form = NA, summary = FALSE))

  # difference between substitution and no change
  preddif <- predsub - predsame
  preddif <- as.data.table(describe_posterior(preddif, centrality = "mean",
                                              ci = 0.95, ci_method = "eti"))
  preddif <- preddif[, .(Mean, CI_low, CI_high)]

  # save results
  out <- do.call(cbind, preddif)
  out <- cbind(out, subd[, c("MinSubstituted", "Substitute", "Predictor")])
  out <- as.data.table(out)
  names(out) <- c("Mean", "CI_low", "CI_high", "MinSubstituted", "Substitute", "Predictor")

  ## final results for entire composition
  allout[[i]] <- out
}