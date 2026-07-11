#' Lag Columns by Time
#'
#' Create lagged versions of one or more columns using the actual time values
#' rather than row position. Each row's lag is taken from the row (within the
#' same group, when `group` is supplied) whose time value equals
#' `time - n * time_step`. When no observation exists at that time -- for
#' example, a skipped day in a daily diary -- the lag is `NA` instead of the
#' value from before the gap. This makes the lags robust to missing time
#' points, unlike positional shifts such as [data.table::shift()].
#'
#' @param data A data frame or `data.table` containing the columns to lag.
#' @param cols Character vector of column names to lag. Columns may be of any
#'   atomic type.
#' @param time Character scalar naming the time column. Time values must be
#'   numeric, `Date`, or `POSIXct`, contain no missing values, and be unique
#'   within each group.
#' @param group Optional character scalar naming a grouping column (for
#'   example, a participant identifier). When `NULL`, the default, all rows
#'   are treated as one series.
#' @param time_step Optional positive scalar giving the time difference that
#'   counts as one lag step. When `NULL`, the default, the step is inferred
#'   as the smallest positive within-group difference between consecutive
#'   time values. For `Date` time the unit is days; for `POSIXct` time the
#'   unit is seconds. A `difftime` value is converted with [as.numeric()].
#' @param n Positive whole number giving the lag order; the lag for a row at
#'   time `t` is the observation at `t - n * time_step`. To create several
#'   lag orders, call the function repeatedly with different `n` and
#'   `suffix` values.
#' @param suffix Non-empty character scalar appended to each name in `cols`
#'   to form the new lag column names. The new names must not collide with
#'   existing columns.
#'
#' @return `data` with one new column per entry of `cols`, named
#'   `paste0(cols, suffix)`, appended in the original row order. The input is
#'   never modified in place; a data frame input returns a data frame and a
#'   `data.table` input returns a new `data.table`. The result carries an
#'   attribute `"lag_by_time"`, a list with:
#'   \item{`time_step`}{The (possibly inferred) step used, as a number.}
#'   \item{`time_step_inferred`}{Logical, whether the step was inferred.}
#'   \item{`lag_columns`}{Named character vector mapping `cols` to the new
#'   lag column names.}
#'   \item{`n_gaps`}{Number of rows whose expected source time falls after
#'   the start of the group's series but has no observation (a gap); these
#'   rows receive `NA` lags.}
#'   \item{`n_irregular`}{Number of rows whose time value is not (within
#'   tolerance) a whole number of steps from the start of the group's
#'   series; these rows receive `NA` lags and are never used as lag
#'   sources.}
#'
#' @details
#' Rows are matched on a per-group integer time grid: within each group,
#' time values are converted to step counts from the earliest observation,
#' `(time - min(time)) / time_step`, and a row's lag source is the row `n`
#' grid points earlier. A value is treated as on the grid when it is within
#' `1e-8 * max(1, time_step)` (on the time scale) of a whole number of
#' steps. Off-grid rows -- irregular spacing that is not a multiple of the
#' step -- silently receive `NA` lags and are counted in `n_irregular`; if
#' this affects many rows, supply the intended `time_step` explicitly.
#'
#' Gaps are silent by design: no message or warning is emitted when a lag is
#' `NA` because the previous time point is missing. Inspect the
#' `"lag_by_time"` attribute to see how many rows were affected.
#'
#' When `time_step` is inferred, at least one group must contain two or more
#' time points; otherwise an error asks for an explicit `time_step`.
#' Duplicate time values within a group are an error.
#'
#' @examples
#' # day 3 is missing, so the lag on day 4 is NA rather than the day-2 value
#' d <- data.frame(
#'   id = c(1, 1, 1, 2, 2, 2),
#'   day = c(1, 2, 4, 1, 2, 3),
#'   y = c(10, 11, 12, 20, 21, 22)
#' )
#' lag_by_time(d, cols = "y", time = "day", group = "id")
#'
#' # Date time with an explicit one-day step
#' d2 <- data.frame(date = as.Date("2026-01-01") + c(0, 1, 3), y = 1:3)
#' lag_by_time(d2, cols = "y", time = "date", time_step = 1)
#'
#' @seealso [prep_sim_analysis()], which uses `lag_by_time()` to build its
#'   `ar1()` lag columns.
#'
#' @export
lag_by_time <- function(data, cols, time, group = NULL,
                        time_step = NULL, n = 1L, suffix = "_lag") {
  if (!is.data.frame(data)) {
    .mlsim_stop("`data` must be a data frame or data.table.")
  }
  cols <- .mlsim_check_vars(cols, arg = "cols")
  missing_cols <- cols[!cols %in% names(data)]
  if (length(missing_cols) > 0L) {
    .mlsim_stop(
      "`cols` must name columns of `data`; missing: %s.",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (!is.character(time) || length(time) != 1L || is.na(time) || !nzchar(time)) {
    .mlsim_stop("`time` must be one column name.")
  }
  if (!time %in% names(data)) {
    .mlsim_stop("`time` column `%s` is missing from `data`.", time)
  }
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1L || is.na(group) || !nzchar(group)) {
      .mlsim_stop("`group` must be `NULL` or one column name.")
    }
    if (!group %in% names(data)) {
      .mlsim_stop("`group` column `%s` is missing from `data`.", group)
    }
    if (anyNA(data[[group]])) {
      .mlsim_stop("`group` column `%s` must not contain missing values.", group)
    }
  }
  time_values <- data[[time]]
  if (!(is.numeric(time_values) || inherits(time_values, "Date") ||
        inherits(time_values, "POSIXt"))) {
    .mlsim_stop("`time` column `%s` must contain numeric, Date, or POSIXct values.", time)
  }
  if (anyNA(time_values)) {
    .mlsim_stop("`time` column `%s` must not contain missing values.", time)
  }
  if (!is.null(time_step)) {
    if (inherits(time_step, "difftime")) {
      time_step <- as.numeric(time_step)
    }
    time_step <- .mlsim_check_scalar_number(time_step, "time_step")
    if (time_step <= 0) {
      .mlsim_stop("`time_step` must be a positive number.")
    }
  }
  n <- .mlsim_check_scalar_number(n, "n")
  if (n < 1 || n != round(n)) {
    .mlsim_stop("`n` must be a positive whole number.")
  }
  n <- as.integer(n)
  if (!is.character(suffix) || length(suffix) != 1L || is.na(suffix) || !nzchar(suffix)) {
    .mlsim_stop("`suffix` must be a non-empty string.")
  }
  new_names <- paste0(cols, suffix)
  existing <- new_names[new_names %in% names(data)]
  if (length(existing) > 0L) {
    .mlsim_stop(
      "Lag columns collide with existing columns in `data`: %s.",
      paste(existing, collapse = ", ")
    )
  }

  n_rows <- nrow(data)
  t_num <- as.numeric(time_values)
  row_groups <- if (is.null(group)) {
    list(seq_len(n_rows))
  } else {
    unname(split(seq_len(n_rows), data[[group]], drop = TRUE))
  }

  .mlc_lag_check_unique_time(row_groups, t_num, time, group)
  inferred <- is.null(time_step)
  if (inferred) {
    time_step <- .mlc_lag_infer_step(row_groups, t_num)
  }

  source_rows <- rep(NA_integer_, n_rows)
  n_gaps <- 0L
  n_irregular <- 0L
  for (rows in row_groups) {
    matched <- .mlc_lag_match_positions(t_num[rows], time_step, n)
    source_rows[rows] <- rows[matched$positions]
    n_gaps <- n_gaps + matched$n_gaps
    n_irregular <- n_irregular + matched$n_irregular
  }

  is_dt <- data.table::is.data.table(data)
  out <- if (is_dt) data.table::copy(data) else data
  for (j in seq_along(cols)) {
    lag_values <- out[[cols[[j]]]][source_rows]
    if (is_dt) {
      data.table::set(out, j = new_names[[j]], value = lag_values)
    } else {
      out[[new_names[[j]]]] <- lag_values
    }
  }
  lag_meta <- list(
    time_step = time_step,
    time_step_inferred = inferred,
    lag_columns = stats::setNames(new_names, cols),
    n_gaps = n_gaps,
    n_irregular = n_irregular
  )
  if (is_dt) {
    # setattr() keeps the data.table self-reference valid; attr<- would not
    data.table::setattr(out, "lag_by_time", lag_meta)
  } else {
    attr(out, "lag_by_time") <- lag_meta
  }
  out
}

.mlc_lag_check_unique_time <- function(row_groups, t_num, time, group) {
  for (rows in row_groups) {
    if (anyDuplicated(t_num[rows])) {
      if (is.null(group)) {
        .mlsim_stop("`%s` values must be unique.", time)
      }
      .mlsim_stop("`%s` values must be unique within each `%s`.", time, group)
    }
  }
  invisible(TRUE)
}

.mlc_lag_infer_step <- function(row_groups, t_num) {
  diffs <- unlist(lapply(row_groups, function(rows) diff(sort(t_num[rows]))))
  diffs <- diffs[diffs > 0]
  if (length(diffs) == 0L) {
    .mlsim_stop(
      "Cannot infer `time_step` because no group has at least two time points; supply `time_step`."
    )
  }
  min(diffs)
}

.mlc_lag_match_positions <- function(t_num, time_step, n) {
  if (length(t_num) == 0L) {
    return(list(positions = integer(), n_gaps = 0L, n_irregular = 0L))
  }
  index_real <- (t_num - min(t_num)) / time_step
  index <- round(index_real)
  on_grid <- abs(index_real - index) * time_step <= 1e-8 * max(1, time_step)
  index[!on_grid] <- NA_real_
  if (anyDuplicated(index[on_grid])) {
    # distinct time values collapsing onto one grid point means the step is
    # too coarse to separate them
    .mlsim_stop(
      "Multiple time values map to the same grid point for `time_step = %s`; supply a smaller `time_step`.",
      format(time_step)
    )
  }
  positions <- match(index - n, index, incomparables = NA)
  gaps <- on_grid & is.na(positions) & index - n >= 0
  list(
    positions = positions,
    n_gaps = sum(gaps),
    n_irregular = sum(!on_grid)
  )
}
