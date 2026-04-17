library(data.table)

diag_distance <- function(n) {
  data.table(
    distance = seq_len(n),
    expected_distance = seq_len(n) / 2,
    deviates = rep(0, n),
    is_extremevalue = rep(c(TRUE, FALSE), length.out = n)
  )
}

diag_list <- function(n = 8, p = 3) {
  x <- matrix(
    seq_len(n * p) + 10,
    nrow = n,
    dimnames = list(NULL, paste0("part", seq_len(p)))
  )

  list(
    x = x,
    distance = diag_distance(n),
    cutoff = 3,
    ev.perc = .1,
    extremevalues = "theoretical",
    levels = "total"
  )
}

small_coda_data <- function() {
  data.table(
    id = rep(seq_len(8), each = 2),
    A = c(20, 21, 30, 28, 25, 26, 40, 39, 35, 36, 45, 44, 50, 49, 55, 56),
    B = c(30, 29, 25, 26, 35, 34, 20, 21, 25, 24, 30, 31, 20, 22, 25, 23),
    C = c(50, 50, 45, 46, 40, 40, 40, 40, 40, 40, 25, 25, 30, 29, 20, 21)
  )
}

small_complr <- function(idvar = "id") {
  dat <- small_coda_data()
  if (is.null(idvar)) {
    dat[, id := NULL]
  }

  complr(
    data = dat,
    parts = c("A", "B", "C"),
    sbp = build.sbp(c("A", "B", "C")),
    idvar = idvar,
    total = 100
  )
}

minimal_complr <- function(coords = matrix(1:12, nrow = 6),
                           raw = matrix(seq_len(18) + 1, nrow = 6),
                           idvar = NULL) {
  dataout <- if (is.null(idvar)) {
    data.table(row_id = seq_len(max(NROW(coords), NROW(raw), 1L)))
  } else {
    data.table(id = seq_len(max(NROW(coords), NROW(raw), 1L)))
  }

  as.complr(list(
    output = list(list(
      Z = coords, X = raw,
      bZ = coords, bX = raw,
      wZ = coords, wX = raw
    )),
    datain = data.table(),
    dataout = dataout,
    transform = "ilr",
    idvar = idvar
  ))
}

test_that("as.diagnostics validates and preserves diagnostics objects", {
  base <- diag_list()
  diag <- as.diagnostics(base)

  expect_true(is.diagnostics(diag))
  expect_false(is.diagnostics(base))
  expect_identical(as.diagnostics(diag), diag)
  expect_named(diag, c("x", "distance", "cutoff", "ev.perc", "extremevalues", "levels"))

  expect_error(as.diagnostics(1), "is.list")

  missing_name <- base
  missing_name$levels <- NULL
  expect_error(as.diagnostics(missing_name), "exactly the following elements")

  extra_name <- base
  extra_name$extra <- TRUE
  expect_error(as.diagnostics(extra_name), "exactly the following elements")

  bad_x <- base
  bad_x$x <- as.data.frame(bad_x$x)
  expect_error(as.diagnostics(bad_x), "x must be a matrix")

  bad_distance <- base
  bad_distance$distance <- as.data.frame(bad_distance$distance)
  expect_error(as.diagnostics(bad_distance), "distance must be a data.table")

  bad_cutoff <- base
  bad_cutoff$cutoff <- c(1, 2)
  expect_error(as.diagnostics(bad_cutoff), "cutoff must be a length 1 numeric")

  bad_ev_perc <- base
  bad_ev_perc$ev.perc <- "0.1"
  expect_error(as.diagnostics(bad_ev_perc), "ev.perc must be a length 1 numeric")

  bad_extremevalues <- base
  bad_extremevalues$extremevalues <- c("theoretical", "empirical")
  expect_error(as.diagnostics(bad_extremevalues), "extremevalues must be a length 1 character")

  bad_levels <- base
  bad_levels$levels <- 1
  expect_error(as.diagnostics(bad_levels), "levels must be a length 1 character")
})

test_that("diagnostics.complr returns expected diagnostics for each level", {
  coda <- small_complr()

  total_diag <- diagnostics(coda, level = "total", ev.perc = .25)
  between_diag <- diagnostics(coda, level = "between", ev.perc = .25)
  within_diag <- diagnostics(coda, level = "within", ev.perc = .25)
  empirical_diag <- diagnostics(
    coda,
    level = "between",
    ev.perc = .25,
    extremevalues = "empirical"
  )

  expect_s3_class(total_diag, "diagnostics")
  expect_equal(total_diag$levels, "total")
  expect_equal(nrow(total_diag$x), nrow(small_coda_data()))
  expect_equal(nrow(total_diag$distance), nrow(small_coda_data()))
  expect_equal(total_diag$cutoff, stats::qchisq(.75, df = 2))
  expect_type(total_diag$distance$is_extremevalue, "logical")

  expect_equal(between_diag$levels, "between")
  expect_equal(nrow(between_diag$x), uniqueN(small_coda_data()$id))
  expect_equal(nrow(between_diag$distance), uniqueN(small_coda_data()$id))

  expect_equal(within_diag$levels, "within")
  expect_equal(nrow(within_diag$x), nrow(small_coda_data()))
  expect_equal(nrow(within_diag$distance), nrow(small_coda_data()))

  expect_equal(empirical_diag$extremevalues, "empirical")
  expect_equal(
    empirical_diag$cutoff,
    stats::quantile(empirical_diag$distance$distance, probs = .75, names = FALSE)
  )
})

test_that("diagnostics.complr validates inputs and available data", {
  coda <- small_complr()
  wide_coda <- small_complr(idvar = NULL)

  expect_error(
    multilevelcoda:::diagnostics.complr(list()),
    "is.complr"
  )
  expect_error(
    diagnostics(coda, ev.perc = NA_real_),
    "between 0 and 1"
  )
  expect_error(
    diagnostics(coda, ev.perc = 1.5),
    "between 0 and 1"
  )
  expect_error(
    diagnostics(coda, level = "other"),
    "'arg' should be one of"
  )
  expect_error(
    diagnostics(wide_coda, level = "between"),
    "require a complr object with an idvar"
  )

  expect_error(
    diagnostics(minimal_complr(coords = NULL), level = "total"),
    "total-level coordinates are not available"
  )
  expect_error(
    diagnostics(minimal_complr(raw = NULL), level = "total"),
    "total-level raw values are not available"
  )
  expect_error(
    diagnostics(
      minimal_complr(
        coords = matrix(seq_len(6), nrow = 2),
        raw = matrix(seq_len(8) + 1, nrow = 2)
      ),
      level = "total"
    ),
    "complete observations are required"
  )
})

test_that("diagnostics.complr handles covariance edge cases", {
  coda <- small_complr()

  near_zero <- with_mocked_bindings(
    diagnostics(coda, level = "total"),
    covMcd = function(x, ...) {
      list(center = colMeans(x), cov = diag(c(1e-12, 2)))
    },
    .package = "multilevelcoda"
  )

  singular <- with_mocked_bindings(
    diagnostics(coda, level = "total"),
    covMcd = function(x, ...) {
      list(center = colMeans(x), cov = matrix(1, ncol(x), ncol(x)))
    },
    .package = "multilevelcoda"
  )

  expect_s3_class(near_zero, "diagnostics")
  expect_s3_class(singular, "diagnostics")
  expect_false(anyNA(near_zero$distance$distance))
  expect_false(anyNA(singular$distance$distance))
})

test_that("plot.diagnostics builds raw and clr scatterplot matrices", {
  plot_diag <- as.diagnostics(list(
    x = matrix(
      c(
        10, 10, 10, 20, 20, 20, 30, 30,
        20, 22, 24, 26, 28, 30, 32, 34,
        70, 68, 66, 54, 52, 50, 38, 36
      ),
      ncol = 3,
      dimnames = list(NULL, c("A", "B", "C"))
    ),
    distance = diag_distance(8),
    cutoff = 7,
    ev.perc = .1,
    extremevalues = "theoretical",
    levels = "total"
  ))

  raw_plot <- plot(plot_diag, transform = "raw")
  clr_plot <- plot(plot_diag, transform = "clr")

  expect_s3_class(raw_plot, "ggplot")
  expect_s3_class(clr_plot, "ggplot")
  expect_length(raw_plot$layers, 3L)
  expect_length(clr_plot$layers, 3L)
})

test_that("plot.diagnostics validates diagnostic plot inputs", {
  base <- diag_list()

  expect_error(
    multilevelcoda:::plot.diagnostics(list()),
    "is.diagnostics"
  )

  row_mismatch <- as.diagnostics(base)
  row_mismatch$distance <- row_mismatch$distance[-1]
  expect_error(
    plot(row_mismatch),
    "must have the same number of rows"
  )

  name_conflict <- base
  colnames(name_conflict$x)[1] <- "obs_id"
  expect_error(
    plot(as.diagnostics(name_conflict)),
    "must not contain a column named"
  )

  expect_error(
    plot(as.diagnostics(base), transform = "other"),
    "'arg' should be one of"
  )
})
