# Update [`brmcoda`](https://florale.github.io/multilevelcoda/reference/brmcoda.md) models

This method allows for updating an existing
[`brmcoda`](https://florale.github.io/multilevelcoda/reference/brmcoda.md)
object.

## Usage

``` r
# S3 method for class 'brmcoda'
update(object, formula. = NULL, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted
  [`brmcoda`](https://florale.github.io/multilevelcoda/reference/brmcoda.md)
  object to be updated.

- formula.:

  Changes to the formula; for details see
  [`update.formula`](https://rdrr.io/r/stats/update.formula.html) and
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.html).

- newdata:

  A `data.frame` or `data.table` containing data of all variables used
  in the analysis. It must include a composition and the same ID
  variable as the existing
  [`complr`](https://florale.github.io/multilevelcoda/reference/complr.md)
  object.

- ...:

  Further arguments passed to
  [`brm`](https://paulbuerkner.com/brms/reference/brm.html).

## Value

A
[`brmcoda`](https://florale.github.io/multilevelcoda/reference/brmcoda.md)
with two elements

- `complr`:

  An object of class `complr` used in the `brm` model.

- `model`:

  An object of class `brmsfit`, which contains the posterior draws along
  with many other useful information about the model.

## See also

[`brmcoda`](https://florale.github.io/multilevelcoda/reference/brmcoda.md)

## Examples

``` r
# \donttest{
if(requireNamespace("cmdstanr")){

# model with compositional predictor at between and within-person levels
fit <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
                               parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
                               idvar = "ID"),
                formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + Female + (1 | ID),
                chain = 1, iter = 500,
              backend = "cmdstanr")

# removing the effect of bz1_1
fit1 <- update(fit, formula. = ~ . - bz1_1)

# using only a subset
fit2 <- update(fit, newdata = mcompd[ID != 1])
}# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
