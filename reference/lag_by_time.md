# Lag Columns by Time

Create lagged versions of one or more columns using the actual time
values rather than row position. Each row's lag is taken from the row
(within the same group, when `group` is supplied) whose time value
equals `time - n * time_step`. When no observation exists at that time –
for example, a skipped day in a daily diary – the lag is `NA` instead of
the value from before the gap. This makes the lags robust to missing
time points, unlike positional shifts such as
[`data.table::shift()`](https://rdrr.io/pkg/data.table/man/shift.html).

## Usage

``` r
lag_by_time(
  data,
  cols,
  time,
  group = NULL,
  time_step = NULL,
  n = 1L,
  suffix = "_lag"
)
```

## Arguments

- data:

  A data frame or `data.table` containing the columns to lag.

- cols:

  Character vector of column names to lag. Columns may be of any atomic
  type.

- time:

  Character scalar naming the time column. Time values must be numeric,
  `Date`, or `POSIXct`, contain no missing values, and be unique within
  each group.

- group:

  Optional character scalar naming a grouping column (for example, a
  participant identifier). When `NULL`, the default, all rows are
  treated as one series.

- time_step:

  Optional positive scalar giving the time difference that counts as one
  lag step. When `NULL`, the default, the step is inferred as the
  smallest positive within-group difference between consecutive time
  values. For `Date` time the unit is days; for `POSIXct` time the unit
  is seconds. A `difftime` value is converted with
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html).

- n:

  Positive whole number giving the lag order; the lag for a row at time
  `t` is the observation at `t - n * time_step`. To create several lag
  orders, call the function repeatedly with different `n` and `suffix`
  values.

- suffix:

  Non-empty character scalar appended to each name in `cols` to form the
  new lag column names. The new names must not collide with existing
  columns.

## Value

`data` with one new column per entry of `cols`, named
`paste0(cols, suffix)`, appended in the original row order. The input is
never modified in place; a data frame input returns a data frame and a
`data.table` input returns a new `data.table`. The result carries an
attribute `"lag_by_time"`, a list with:

- `time_step`:

  The (possibly inferred) step used, as a number.

- `time_step_inferred`:

  Logical, whether the step was inferred.

- `lag_columns`:

  Named character vector mapping `cols` to the new lag column names.

- `n_gaps`:

  Number of rows whose expected source time falls after the start of the
  group's series but has no observation (a gap); these rows receive `NA`
  lags.

- `n_irregular`:

  Number of rows whose time value is not (within tolerance) a whole
  number of steps from the start of the group's series; these rows
  receive `NA` lags and are never used as lag sources.

## Details

Rows are matched on a per-group integer time grid: within each group,
time values are converted to step counts from the earliest observation,
`(time - min(time)) / time_step`, and a row's lag source is the row `n`
grid points earlier. A value is treated as on the grid when it is within
`1e-8 * max(1, time_step)` (on the time scale) of a whole number of
steps. Off-grid rows – irregular spacing that is not a multiple of the
step – silently receive `NA` lags and are counted in `n_irregular`; if
this affects many rows, supply the intended `time_step` explicitly.

Gaps are silent by design: no message or warning is emitted when a lag
is `NA` because the previous time point is missing. Inspect the
`"lag_by_time"` attribute to see how many rows were affected.

When `time_step` is inferred, at least one group must contain two or
more time points; otherwise an error asks for an explicit `time_step`.
Duplicate time values within a group are an error.

## See also

[`prep_sim_analysis()`](https://florale.github.io/multilevelcoda/reference/prep_sim_analysis.md),
which uses `lag_by_time()` to build its `ar1()` lag columns.

## Examples

``` r
# day 3 is missing, so the lag on day 4 is NA rather than the day-2 value
d <- data.frame(
  id = c(1, 1, 1, 2, 2, 2),
  day = c(1, 2, 4, 1, 2, 3),
  y = c(10, 11, 12, 20, 21, 22)
)
lag_by_time(d, cols = "y", time = "day", group = "id")
#>   id day  y y_lag
#> 1  1   1 10    NA
#> 2  1   2 11    10
#> 3  1   4 12    NA
#> 4  2   1 20    NA
#> 5  2   2 21    20
#> 6  2   3 22    21

# Date time with an explicit one-day step
d2 <- data.frame(date = as.Date("2026-01-01") + c(0, 1, 3), y = 1:3)
lag_by_time(d2, cols = "y", time = "date", time_step = 1)
#>         date y y_lag
#> 1 2026-01-01 1    NA
#> 2 2026-01-02 2     1
#> 3 2026-01-04 3    NA
```
