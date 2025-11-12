# Fit Bayesian generalised (non-)linear multilevel compositional model via full Bayesian inference

Fit a `brm` model with multilevel ILR coordinates

## Usage

``` r
brmcoda(complr, formula, ...)
```

## Arguments

- complr:

  A
  [`complr`](https://florale.github.io/multilevelcoda/reference/complr.md)
  object containing data of composition, ILR coordinates, and other
  variables used in the model.

- formula:

  A object of class `formula`, `brmsformula`: A symbolic description of
  the model to be fitted. Details of the model specification can be
  found in
  [`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.html).

- ...:

  Further arguments passed to
  [`brm`](https://paulbuerkner.com/brms/reference/brm.html).

## Value

A `brmcoda` with two elements

- `complr`:

  An object of class `complr` used in the `brm` model.

- `model`:

  An object of class `brmsfit`, which contains the posterior draws along
  with many other useful information about the model.

## Examples

``` r
# \donttest{
if(requireNamespace("cmdstanr")){
  x1 <- complr(data = mcompd, sbp = sbp,
                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"), idvar = "ID")

  # inspect variables before passing to brmcoda
  get_variables(x1)

  ## model with compositional predictor at between and within-person levels
  m1 <- brmcoda(complr = x1,
                formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
                chain = 1, iter = 500,
                backend = "cmdstanr")

  ## model with compositional outcome
  m2 <- brmcoda(complr = x1,
                formula = mvbind(z1_1, z2_1, z3_1, z4_1) ~ Stress + Female + (1 | ID),
                chain = 1, iter = 500,
                backend = "cmdstanr")

  ## model with compositional predictor and outcome
  x2 <- complr(data = mcompd,
                parts = list(c("TST", "WAKE"), c("MVPA", "LPA", "SB")),
                total = list(c(480), c(960)),
                idvar = "ID",
                transform = "ilr")

  m3 <- brmcoda(complr = x2,
                formula = mvbind(z1_2, z2_2) ~ z1_1 + Female + (1 | ID),
                chain = 1, iter = 500,
                backend = "cmdstanr")
  }# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
