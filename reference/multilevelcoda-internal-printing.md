# Internal Printing and Summary Helpers

Build compact summaries for `mlsim_data` objects and provide small
utility conversions used by print methods.

## Usage

``` r
.mlsim_na_character(x)

.mlsim_na_integer(x)

.mlsim_has_metadata(x)

.mlsim_generated_columns(x)

.mlsim_design_summary(x)

.mlsim_generator_summary_row(name, metadata)

.mlsim_generator_summary(x)

.mlsim_print_table(x, ..., row.names = FALSE)
```

## Arguments

- x:

  An `mlsim_data`, summary, metadata, or table object.

- name:

  Generator name.

- metadata:

  Generator metadata list.

- ...:

  Additional arguments passed to table printing.

- row.names:

  Logical; passed to [`print()`](https://rdrr.io/r/base/print.html).

## Value

Internal helper outputs used by print and summary methods.

## Examples

``` r
multilevelcoda:::.mlsim_na_character(NULL)
#> [1] NA
sim <- simulate_data(
  n = 2,
  generators = list(
    x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1)
  )
)
multilevelcoda:::.mlsim_generated_columns(sim)
#> [1] "x"
```
