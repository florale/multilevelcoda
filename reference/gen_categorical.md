# Generate Categorical Variables

Create a generator specification for binary or multicategory categorical
variables.

## Usage

``` r
gen_categorical(
  vars,
  level = c("single", "level2", "multilevel"),
  categories = NULL,
  fixed_intercept = NULL,
  random_cov = NULL,
  reference = NULL,
  output = c("factor", "character", "integer"),
  ordered = FALSE
)
```

## Arguments

- vars:

  Character scalar naming the generated variable.

- level:

  Simulation level. `"single"` generates row-level categories,
  `"level2"` generates group-level categories, and `"multilevel"` uses a
  baseline-category logit random-intercept model.

- categories:

  Vector of category values. Defaults to `c(0L, 1L)`.

- fixed_intercept:

  Baseline-category logits. For `k` categories this has length `k - 1`
  and is named or ordered by non-reference category.

- random_cov:

  Multilevel covariance matrix for baseline-category random intercepts.

- reference:

  Reference category for the baseline-category logit model. Defaults to
  the first category.

- output:

  Output type: `"factor"`, `"character"`, or `"integer"`. Integer output
  uses zero-based category codes.

- ordered:

  Logical; when `TRUE`, return an ordered factor. Requires
  `output = "factor"`.

## Value

An `mlsim_generator_spec` object for use in
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

## See also

Other predictor generators:
[`continuous-generators`](https://florale.github.io/multilevelcoda/reference/continuous-generators.md),
[`count-generators`](https://florale.github.io/multilevelcoda/reference/count-generators.md),
[`gen_custom()`](https://florale.github.io/multilevelcoda/reference/gen_custom.md),
[`gen_mvn()`](https://florale.github.io/multilevelcoda/reference/gen_mvn.md),
[`gen_outcome()`](https://florale.github.io/multilevelcoda/reference/gen_outcome.md),
[`gen_template()`](https://florale.github.io/multilevelcoda/reference/gen_template.md)

## Examples

``` r
sim <- simulate_data(
  n = 6,
  seed = 3,
  generators = list(
    arm = gen_categorical(
      "arm",
      categories = c("control", "treatment"),
      fixed_intercept = stats::qlogis(0.5)
    )
  )
)
sim$data
#>    obs_id       arm
#>     <int>    <fctr>
#> 1:      1   control
#> 2:      2 treatment
#> 3:      3   control
#> 4:      4   control
#> 5:      5 treatment
#> 6:      6 treatment
```
