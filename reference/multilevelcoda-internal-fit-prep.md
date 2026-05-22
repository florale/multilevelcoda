# Internal Outcome Fit Preparation Helpers

Resolve outcome generator names, transform formulas and helper columns,
and build the internal `mlsim_fit_prep` object used by
[`prepare_outcome_fit()`](https://florale.github.io/multilevelcoda/reference/prepare_outcome_fit.md).

## Usage

``` r
.mlsim_resolve_outcome_name(sim, outcome)

.mlsim_prepare_generic_fit(sim, outcome, metadata, drop_incomplete, ...)

.mlsim_prepare_brmcoda_fit(sim, outcome, metadata, drop_incomplete, ...)

.mlsim_brmcoda_response_map(response)

.mlsim_brmcoda_helper_map(metadata, response_map)

.mlsim_brmcoda_helper_columns(helper_columns, helper_map)

.mlsim_group_lag_values(data, source, idvar)

.mlsim_formula_replace_symbols(formula, replacements)

.mlsim_map_term_string(term, replacements)

.mlsim_complete_formula_rows(data, formula)

.mlsim_fit_prep(...)
```

## Arguments

- sim:

  An `mlsim_data` object.

- outcome:

  Outcome generator name.

- metadata:

  Outcome generator metadata.

- drop_incomplete:

  Logical flag controlling row filtering.

- ...:

  Additional target-specific arguments.

- response:

  Outcome response names.

- response_map:

  Mapping between simulated and target response columns.

- data:

  Data table or data frame.

- source:

  Source column name.

- idvar:

  Optional grouping identifier.

- formula:

  Formula to transform or inspect.

- replacements:

  Named character vector of symbol replacements.

- term:

  Model term string.

## Value

Internal helper outputs used by
[`prepare_outcome_fit()`](https://florale.github.io/multilevelcoda/reference/prepare_outcome_fit.md).

## Examples

``` r
multilevelcoda:::.mlsim_map_term_string("x:y", c(x = "x_fit"))
#> [1] "x_fit:y"
multilevelcoda:::.mlsim_complete_formula_rows(
  data.table::data.table(y = c(1, NA), x = c(2, 3)),
  y ~ x
)
#> [1]  TRUE FALSE
```
