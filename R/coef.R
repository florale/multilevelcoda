if (FALSE} {
sbp2 <- copy(sbp)
sbp2[1, ] <- c(-1, 1, -1, -1, -1)
sbp2[2, ] <- c( 1, 0, -1, -1, -1)
sbp2[3, ] <- c( 0, 0,  1, -1, -1)
sbp2[4, ] <- c( 0, 0,  0,  1, -1)

cilr <- complr(data = mcompd, sbp = sbp2, 
  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
  idvar = "ID")

m1 <- brmcoda(complr = cilr,
              formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
                                 wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
              chain = 4, iter = 1000, cores = 4L,
              backend = "cmdstanr")
# summary(m1)
substition_coef(m1, level = "between", h = .1, scale = 10)
substition_coef(m1, level = "within", h = .1, scale = 10)

substition_coef <- function(object, level = c("between", "within"), h = .1, scale = 1) {
  ## TODO remove after testing
  # object <- m1
  # h <- .1
  # level <- "between"
  # scale <- 60
  level <- match.arg(level)
  expect_s3_class(object, "brmcoda")

  parts <- object$complr$parts
  x <- object$complr$data[, ..parts]

  y0 <- fitted(object$model,
    newdata = model.frame(object),
    re_formula = NA,
    summary = FALSE
  )

  out <- vector("list", length(parts))

  for (k in seq_along(parts)) {
    w <- (-x) / rowSums(x[, -..k])
    w <- as.data.table(w)
    w[, (parts[k]) := +1]

    x2 <- x + (h * w)

    expect_equal(rowSums(x2), rowSums(x), tolerance = 1e-3)

    rm(w) ## cleanup

    d2 <- copy(object$complr$data)
    d2[, (parts) := x2]

    rm(x2) ## cleanup

    cilr2 <- complr(
      data = d2,
      sbp = object$complr$sbp,
      parts = parts,
      idvar = object$complr$idvar
    )

    switch(level,
      within = {
        bilr2 <- object$complr$between_logratio
        wilr2 <- cilr2$logratio - object$complr$between_logratio
        names(wilr2) <- names(object$complr$within_logratio)
      },
      between = {
        expect_equal(cilr$within_logratio, cilr2$within_logratio, tolerance = 1e-3)
        bilr2 <- cilr2$between_logratio
        wilr2 <- object$complr$within_logratio
      }
    )

    rm(cilr2) ## cleanup

    d2 <- cbind(d2, bilr2, wilr2)

    rm(bilr2, wilr2)

    y2 <- fitted(object$model,
      newdata = d2,
      re_formula = NA,
      summary = FALSE
    )

    rm(d2) ## cleanup

    out[[k]] <- rowMeans((y2 - y0) / h)

    rm(y2)
  }

  rm(x, y0) ## cleanup

  finalout <- cbind(
    Component = parts,
    as.data.table(do.call(rbind, lapply(out, posterior_summary)) * scale)
  )

  return(finalout)
}

}