# Create a Summary of a Substitution Model represented by a `substitution` object

Create a Summary of a Substitution Model represented by a `substitution`
object

## Usage

``` r
# S3 method for class 'substitution'
summary(object, delta, to, from, ref, level, digits = 2, ...)
```

## Arguments

- object:

  A `substitution` class object.

- delta:

  A integer, numeric value or vector indicating the desired `delta` at
  which substitution results should be summarised. Default to all
  `delta` available in the `substitution` object.

- to:

  A character value or vector specifying the names of the compositional
  parts that were reallocated to in the model.

- from:

  A character value or vector specifying the names of the compositional
  parts that were reallocated from in the model.

- ref:

  Either a character value or vector ((`"grandmean"` and/or
  `"clustermean"` or `"users"`), Default to all `ref` available in the
  `substitution` object.

- level:

  A character string or vector (`"between"` and/or `"within"`). Default
  to all `level` available in the `substitution` object.

- digits:

  A integer value used for number formatting. Default is `2`.

- ...:

  generic argument, not in use.

## Value

A summary of `substitution` object.

- `Mean`:

  Posterior means.

- `CI_low` and `CI_high`:

  95% credible intervals.

- `Delta`:

  Amount substituted across compositional parts.

- `From`:

  Compositional part that is substituted from.

- `To`:

  Compositional parts that is substituted to.

- `Level`:

  Level where changes in composition takes place. Either `between` or
  `within`.

- `Reference`:

  Either `grandmean`, `clustermean`, or `users`.

## Examples

``` r
# \donttest{
if(requireNamespace("cmdstanr")){
  ## fit a model with compositional predictor at between and between-person levels
  m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
                               parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
                               idvar = "ID", total = 1440),
  formula = Stress ~ bz1_1 + bz2_1 + bz3_1 + bz4_1 +
                     wz1_1 + wz2_1 + wz3_1 + wz4_1 + (1 | ID),
  chain = 1, iter = 500,
  backend = "cmdstanr")

  subm <- substitution(object = m, delta = 5)
  summary(subm)
}# }
#> Error: CmdStan path has not been set yet. See ?set_cmdstan_path.
```
