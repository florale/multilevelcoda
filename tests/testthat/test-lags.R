test_that("lag_by_time matches positional shift on a complete equal grid", {
  d <- data.frame(
    id = rep(1:3, each = 4),
    day = rep(1:4, times = 3),
    y = rnorm(12)
  )
  out <- lag_by_time(d, cols = "y", time = "day", group = "id")
  positional <- unlist(lapply(split(d$y, d$id), function(y) c(NA, head(y, -1L))),
                       use.names = FALSE)

  expect_equal(out$y_lag, positional)
  meta <- attr(out, "lag_by_time")
  expect_identical(meta$time_step, 1)
  expect_true(meta$time_step_inferred)
  expect_identical(meta$lag_columns, c(y = "y_lag"))
  expect_identical(meta$n_gaps, 0L)
  expect_identical(meta$n_irregular, 0L)
})

test_that("lag_by_time yields NA after a gap, silently", {
  d <- data.frame(
    id = c(1, 1, 1, 1, 2, 2, 2),
    day = c(1, 2, 4, 5, 1, 2, 3),
    y = c(10, 11, 12, 13, 20, 21, 22)
  )
  expect_silent(out <- lag_by_time(d, cols = "y", time = "day", group = "id"))

  expect_equal(out$y_lag, c(NA, 10, NA, 12, NA, 20, 21))
  meta <- attr(out, "lag_by_time")
  expect_identical(meta$n_gaps, 1L)
  expect_identical(meta$n_irregular, 0L)
})

test_that("lag_by_time handles Date time with skipped days", {
  d <- data.frame(date = as.Date("2026-01-01") + c(0, 1, 3, 4), y = 1:4)
  inferred <- lag_by_time(d, cols = "y", time = "date")
  explicit <- lag_by_time(d, cols = "y", time = "date", time_step = 1)

  expect_equal(inferred$y_lag, c(NA, 1, NA, 3))
  expect_equal(explicit$y_lag, inferred$y_lag)
  expect_identical(attr(inferred, "lag_by_time")$time_step, 1)
  expect_false(attr(explicit, "lag_by_time")$time_step_inferred)
})

test_that("lag_by_time handles POSIXct time with an hourly step and a gap", {
  base <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  d <- data.frame(t = base + 3600 * c(0, 1, 2, 4), y = 1:4)
  out <- lag_by_time(d, cols = "y", time = "t")

  expect_equal(out$y_lag, c(NA, 1, 2, NA))
  meta <- attr(out, "lag_by_time")
  expect_identical(meta$time_step, 3600)
  expect_identical(meta$n_gaps, 1L)
})

test_that("lag_by_time explicit step differing from the grid marks off-grid rows", {
  d <- data.frame(day = 1:6, y = 1:6)
  out <- lag_by_time(d, cols = "y", time = "day", time_step = 2)

  # days 2, 4, 6 are off the step-2 grid anchored at day 1
  expect_equal(out$y_lag, c(NA, NA, 1, NA, 3, NA))
  meta <- attr(out, "lag_by_time")
  expect_identical(meta$n_irregular, 3L)
})

test_that("lag_by_time errors on duplicate times", {
  d <- data.frame(id = c(1, 1), day = c(2, 2), y = 1:2)
  expect_error(
    lag_by_time(d, cols = "y", time = "day", group = "id"),
    "`day` values must be unique within each `id`."
  )
  expect_error(
    lag_by_time(d[, c("day", "y")], cols = "y", time = "day"),
    "`day` values must be unique."
  )
})

test_that("lag_by_time works on an ungrouped series", {
  d <- data.frame(day = c(1, 2, 4), y = c(5, 6, 7))
  out <- lag_by_time(d, cols = "y", time = "day")
  expect_equal(out$y_lag, c(NA, 5, NA))
})

test_that("lag_by_time lags multiple columns of mixed types in one call", {
  d <- data.frame(
    id = c(1, 1, 1),
    day = c(1, 2, 4),
    y = c(1.5, 2.5, 3.5),
    mood = c("low", "high", "low"),
    stringsAsFactors = FALSE
  )
  out <- lag_by_time(d, cols = c("y", "mood"), time = "day", group = "id")

  expect_equal(out$y_lag, c(NA, 1.5, NA))
  expect_equal(out$mood_lag, c(NA, "low", NA))
  expect_identical(
    attr(out, "lag_by_time")$lag_columns,
    c(y = "y_lag", mood = "mood_lag")
  )
})

test_that("lag_by_time supports higher lag orders", {
  d <- data.frame(day = c(1, 2, 3, 5), y = c(1, 2, 3, 4))
  out <- lag_by_time(d, cols = "y", time = "day", n = 2, suffix = "_lag2")

  expect_equal(out$y_lag2, c(NA, NA, 1, 3))
})

test_that("lag_by_time leaves truly irregular times as silent NA", {
  d <- data.frame(t = c(0, 1, 2.5, 4), y = 1:4)
  expect_silent(out <- lag_by_time(d, cols = "y", time = "t"))

  # inferred step is 1; 2.5 is off-grid, and t = 4 has no source at t = 3
  expect_equal(out$y_lag, c(NA, 1, NA, NA))
  meta <- attr(out, "lag_by_time")
  expect_identical(meta$n_irregular, 1L)
  expect_identical(meta$n_gaps, 1L)
})

test_that("lag_by_time errors when the step cannot be inferred", {
  d <- data.frame(id = 1:3, day = c(1, 1, 1), y = 1:3)
  expect_error(
    lag_by_time(d, cols = "y", time = "day", group = "id"),
    "Cannot infer `time_step`"
  )
  expect_silent(
    lag_by_time(d, cols = "y", time = "day", group = "id", time_step = 1)
  )
})

test_that("lag_by_time preserves class and never modifies in place", {
  d <- data.frame(day = 1:3, y = 1:3)
  out <- lag_by_time(d, cols = "y", time = "day")
  expect_s3_class(out, "data.frame")
  expect_false(data.table::is.data.table(out))
  expect_false("y_lag" %in% names(d))

  dt <- data.table::as.data.table(d)
  out_dt <- lag_by_time(dt, cols = "y", time = "day")
  expect_true(data.table::is.data.table(out_dt))
  expect_false("y_lag" %in% names(dt))
})

test_that("lag_by_time validates its inputs", {
  d <- data.frame(id = c(1, 1), day = c(1, 2), y = 1:2)
  expect_error(lag_by_time(1:3, cols = "y", time = "day"), "data frame")
  expect_error(lag_by_time(d, cols = "z", time = "day"), "missing: z")
  expect_error(lag_by_time(d, cols = "y", time = "when"), "`when` is missing")
  expect_error(
    lag_by_time(d, cols = "y", time = "id", group = "missing_group"),
    "`missing_group` is missing"
  )
  expect_error(
    lag_by_time(data.frame(day = c("a", "b"), y = 1:2), cols = "y", time = "day"),
    "numeric, Date, or POSIXct"
  )
  expect_error(
    lag_by_time(data.frame(day = c(1, NA), y = 1:2), cols = "y", time = "day"),
    "must not contain missing values"
  )
  expect_error(
    lag_by_time(d, cols = "y", time = "day", group = "id", time_step = 0),
    "positive number"
  )
  expect_error(
    lag_by_time(d, cols = "y", time = "day", group = "id", n = 1.5),
    "positive whole number"
  )
  expect_error(
    lag_by_time(d, cols = "y", time = "day", group = "id", suffix = ""),
    "non-empty string"
  )
  d$y_lag <- 0
  expect_error(
    lag_by_time(d, cols = "y", time = "day", group = "id"),
    "collide with existing columns"
  )
})

test_that("lag_by_time errors when NA appears in the grouping column", {
  d <- data.frame(id = c(1, NA), day = c(1, 2), y = 1:2)
  expect_error(
    lag_by_time(d, cols = "y", time = "day", group = "id"),
    "must not contain missing values"
  )
})

test_that("lag_by_time is invariant to row order", {
  d <- data.frame(
    id = c(1, 1, 1, 2, 2),
    day = c(1, 2, 4, 1, 2),
    y = c(10, 11, 12, 20, 21)
  )
  shuffled <- d[c(4, 1, 5, 3, 2), ]
  out_sorted <- lag_by_time(d, cols = "y", time = "day", group = "id")
  out_shuffled <- lag_by_time(shuffled, cols = "y", time = "day", group = "id")

  expect_identical(rownames(out_shuffled), rownames(shuffled))
  merged <- out_shuffled[order(out_shuffled$id, out_shuffled$day), ]
  expect_equal(merged$y_lag, out_sorted$y_lag)
})
