# Generate Normal, Multivariate Normal, and Compositional Variables

Create a generator specification for univariate or multivariate Gaussian
variables, optionally transformed from ILR coordinates to compositional
parts.

## Usage

``` r
gen_mvn(
  vars,
  level = c("single", "level2", "multilevel"),
  fixed_intercept = NULL,
  residual_cov = NULL,
  random_cov = NULL,
  ...,
  scale_fixed_intercept = NULL,
  residual_cor = NULL,
  compositional = FALSE,
  parts = NULL,
  total = 1,
  keep_ilr = TRUE,
  sbp = NULL
)
```

## Arguments

- vars:

  Character vector naming the generated variables. For compositional
  generators these are ILR coordinate names.

- level:

  Simulation level. `"single"` generates one row-level vector per
  observation, `"level2"` generates one group-level vector and expands
  it to each row in the group, and `"multilevel"` uses a
  random-intercept model.

- fixed_intercept:

  Identity-scale location intercept vector.

- residual_cov:

  Residual covariance matrix for Gaussian generators without a
  row-specific scale model. For univariate generators, it may be a
  scalar residual variance.

- random_cov:

  Group-level random-intercept covariance matrix for multilevel
  generators. When `scale_fixed_intercept` is supplied this may be
  either a location-only covariance or a joint location-scale
  covariance.

- ...:

  Removed direct distribution parameters are rejected.

- scale_fixed_intercept:

  Optional log residual standard-deviation intercepts.

- residual_cor:

  Residual correlation matrix used with `scale_fixed_intercept`.
  Defaults to an identity correlation matrix.

- compositional:

  Logical; if `TRUE`, treat `vars` as ILR coordinates and back-transform
  to composition parts.

- parts:

  Character vector naming the composition parts. Must have
  `length(vars) + 1` entries unless supplied through `sbp` column names.

- total:

  Positive scalar total for closed compositions.

- keep_ilr:

  Logical; if `TRUE`, emit both ILR coordinates and parts. If `FALSE`,
  emit only parts.

- sbp:

  Optional sequential binary partition matrix. Columns name parts; rows
  define ILR balances using `-1`, `0`, and `1`.

## Value

An `mlsim_generator_spec` object for use in
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

## Details

Compositional MVN generators simulate ILR coordinates first, then use
the SBP basis to transform the coordinates into positive composition
parts that sum to `total`.

Multilevel generators additionally emit the latent between- and
within-group components of each generated variable as visible columns.
For non-compositional variables these are `<var>_between` (the
group-level mean, fixed plus random intercept) and `<var>_within` (the
row-level residual), with `<var> = <var>_between + <var>_within`
exactly. For compositional generators the same decomposition is emitted
on the ILR scale (`<ilr>_between`, `<ilr>_within`), together with the
true between composition `<part>_between`, the closed back-transform of
the group-level ILR means. All generated variables are labelled with
column roles so that `between()` and
[`within()`](https://rdrr.io/r/base/with.html) terms in
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md)
formulas resolve to the latent components; roles (and hence
`between()`/[`within()`](https://rdrr.io/r/base/with.html)) apply to the
ILR coordinates of a compositional generator, not to its parts. When
`keep_ilr = FALSE` the observed ILR columns are dropped, but the
between/within columns and their roles are still emitted for multilevel
generators.

Multilevel generators also record their group-level random-effect draws
as per-row truth columns: `.mlsim_<var>_random_intercept` for the
location intercept and, when `random_cov` is a joint location-scale
covariance, `.mlsim_<var>_scale_random_intercept` for the scale
intercept. The other multilevel scalar generators
([`gen_negbin()`](https://florale.github.io/multilevelcoda/reference/count-generators.md),
[`gen_gamma()`](https://florale.github.io/multilevelcoda/reference/continuous-generators.md),
[`gen_beta()`](https://florale.github.io/multilevelcoda/reference/continuous-generators.md),
[`gen_poisson()`](https://florale.github.io/multilevelcoda/reference/count-generators.md),
[`gen_binomial()`](https://florale.github.io/multilevelcoda/reference/count-generators.md))
share this machinery and emit the same columns. See the Value section of
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md)
for the full naming contract.

## See also

Other predictor generators:
[`continuous-generators`](https://florale.github.io/multilevelcoda/reference/continuous-generators.md),
[`count-generators`](https://florale.github.io/multilevelcoda/reference/count-generators.md),
[`gen_categorical()`](https://florale.github.io/multilevelcoda/reference/gen_categorical.md),
[`gen_custom()`](https://florale.github.io/multilevelcoda/reference/gen_custom.md),
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md),
[`gen_template()`](https://florale.github.io/multilevelcoda/reference/gen_template.md)

## Examples

``` r
sim <- simulate_data(
  n = 4,
  seed = 2,
  generators = list(
    x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1),
    z = gen_mvn(
      c("z1", "z2"),
      fixed_intercept = c(0, 0),
      residual_cov = diag(2),
      compositional = TRUE,
      parts = c("sleep", "activity", "sedentary")
    )
  )
)
sim$data
#>    obs_id          x         z1          z2      sleep  activity sedentary
#>     <int>      <num>      <num>       <num>      <num>     <num>     <num>
#> 1:      1 -0.8969145 -1.9844739 -0.08025176 0.04207843 0.4518105 0.5061110
#> 2:      2  0.1848492  0.1387870  0.13242028 0.37108945 0.3438136 0.2850969
#> 3:      3  1.5878453 -0.4176508  0.70795473 0.20997858 0.5777381 0.2122833
#> 4:      4 -1.1303757 -0.9817528 -0.23969802 0.12899890 0.3623855 0.5086157
rowSums(as.matrix(sim$data[, c("sleep", "activity", "sedentary"), with = FALSE]))
#> [1] 1 1 1 1
```
