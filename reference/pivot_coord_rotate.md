# Estimate pivot balance coordinates by rotating sequential binary partition.

This function is an alias of
[`pivot_coord`](https://florale.github.io/multilevelcoda/reference/pivot_coord.md)
to estimates the pivot balance coordinates by `"rotate"` the sequential
binary partition on the same `brmcoda` object.

## Usage

``` r
pivot_coord_rotate(object, summary = TRUE, parts = 1, ...)
```

## Arguments

- object:

  An object of class `brmcoda`.

- summary:

  Should summary statistics be returned instead of the raw values?
  Default is `TRUE`.

- parts:

  A optional character string specifying names of compositional parts
  that should be considered in the substitution analysis. This should
  correspond to a single set of names of compositional parts specified
  in the `complr` object. Default to the first composition in the
  `complr` object.

- ...:

  Further arguments passed to
  [`posterior_summary`](https://paulbuerkner.com/brms/reference/posterior_summary.html).

## Value

Estimated pivot balance coordinates representing the effect of
increasing one compositional part relative to the remaining
compositional parts.

## See also

[`pivot_coord`](https://florale.github.io/multilevelcoda/reference/pivot_coord.md)

## Examples

``` r
# \donttest{
if(requireNamespace("cmdstanr")){
  x <- complr(data = mcompd, sbp = sbp,
                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID",
                 total = 1440)
  
  m <- brmcoda(complr = x,
                formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
                chain = 1, iter = 500,
                backend = "cmdstanr")
  
  m_pivot_coord_rotate <- pivot_coord_rotate(m)
  summary(m_pivot_coord_rotate)
  
  m_pivot_coord_raw <-  pivot_coord_rotate(m, summary = FALSE)
  lapply(m_pivot_coord_raw$output, brms::posterior_summary)
  
  }# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
