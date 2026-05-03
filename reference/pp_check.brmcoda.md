# Posterior Predictive Checks for `brmcoda` Objects

Perform posterior predictive checks with the help of the bayesplot
package.

## Usage

``` r
# S3 method for class 'brmcoda'
pp_check(
  object,
  type = "dens_overlay",
  ndraws = NULL,
  prefix = c("ppc", "ppd"),
  newdata = NULL,
  draw_ids = NULL,
  nsamples = NULL,
  subset = NULL,
  scale = c("linear", "response"),
  parts = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `brmcoda`.

- type:

  Type of the ppc plot as given by a character string. See
  [`PPC`](https://mc-stan.org/bayesplot/reference/PPC-overview.html) for
  an overview of currently supported types. You may also use an invalid
  type (e.g. `type = "xyz"`) to get a list of supported types in the
  resulting error message.

- ndraws:

  Positive integer indicating how many posterior draws should be used.
  If `NULL` all draws are used. If not specified, the number of
  posterior draws is chosen automatically. Ignored if `draw_ids` is not
  `NULL`.

- prefix:

  The prefix of the bayesplot function to be applied. Either \`"ppc"\`
  (posterior predictive check; the default) or \`"ppd"\` (posterior
  predictive distribution), the latter being the same as the former
  except that the observed data is not shown for \`"ppd"\`.

- newdata:

  An optional data.frame for which to evaluate predictions. If `NULL`
  (default), the original data of the model is used. `NA` values within
  factors (excluding grouping variables) are interpreted as if all dummy
  variables of this factor are zero. This allows, for instance, to make
  predictions of the grand mean when using sum coding. `NA` values
  within grouping variables are treated as a new level.

- draw_ids:

  An integer vector specifying the posterior draws to be used. If `NULL`
  (the default), all draws are used.

- nsamples:

  Deprecated alias of `ndraws`.

- subset:

  Deprecated alias of `draw_ids`.

- scale:

  Specifically for models with compositional responses, either
  `"response"` or `"linear"`. If `"linear"`, results are returned on the
  log-ratio scale. If `"response"`, results are returned on the
  compositional scale of the response variable.

- parts:

  Optional names or indices of compositional response parts to include
  when `scale = "response"`. If `NULL`, all parts are used.

- ...:

  Further arguments passed to
  [`predict.brmsfit`](https://paulbuerkner.com/brms/reference/predict.brmsfit.html)
  as well as to the PPC function specified in `type`.

## See also

[`pp_check.brmsfit`](https://paulbuerkner.com/brms/reference/pp_check.brmsfit.html)

## Examples

``` r
# \donttest{
if(requireNamespace("cmdstanr")){
  ## fit a model
  x <- complr(data = mcompd, sbp = sbp,
                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
                 idvar = "ID", total = 1440)

  m1 <- brmcoda(complr = x,
                formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                                   wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
                chain = 1, iter = 500,
                backend = "cmdstanr")

  ## posterior predictive checks
  bayesplot::pp_check(m1, ndraws = 5)

  ## fit a model with compositional outcome
  m2 <- brmcoda(complr = x,
                formula = mvbind(z1_1, z2_1, z3_1, z4_1) ~ 
                          bz1_1 + bz2_1 + bz3_1 + bz4_1 + Female + (1 | ID),
                chain = 1, iter = 500,
                backend = "cmdstanr")

  ## posterior predictive checks for compositional outcome -- linear scale
  bayesplot::pp_check(m2, resp = "z11", ndraws = 10)
  ## posterior predictive checks for compositional outcome -- original response scale
  bayesplot::pp_check(m2, parts = "WAKE", scale = "response", ndraws = 10)
}# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
