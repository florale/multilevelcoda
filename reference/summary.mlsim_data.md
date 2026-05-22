# Summarize Simulated Data

Build a compact summary of an `mlsim_data` object.

## Usage

``` r
# S3 method for class 'mlsim_data'
summary(object, ...)
```

## Arguments

- object:

  An `mlsim_data` object returned by
  [`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

- ...:

  Ignored.

## Value

A `summary.mlsim_data` object with `design` and `generators` tables.

## Examples

``` r
sim <- simulate_data(n = 3, generators = list(x = gen_normal("x")))
summary(sim)
#> <summary.mlsim_data>
#> 
#> Design:
#>      n n_cols n_generated_cols n_groups group_id group_size_min
#>  <int>  <int>            <int>    <int>   <char>          <int>
#>      3      2                1       NA     <NA>             NA
#>  group_size_median group_size_max obs_id time_id  seed n_generators
#>              <num>          <int> <char>  <char> <int>        <int>
#>                 NA             NA obs_id    <NA>    NA            1
#> 
#> Generators:
#>  generator distribution  level   vars n_vars parameter_level parameter_count
#>     <char>       <char> <char> <char>  <int>          <char>           <int>
#>          x       normal single      x      1             row               3
#>  has_row_parameters has_group_parameters has_fixed_parameters has_random_cov
#>              <lgcl>               <lgcl>               <lgcl>         <lgcl>
#>                TRUE                FALSE                FALSE          FALSE
#>  has_random_effects has_residuals has_scale_model has_formula has_ar_terms
#>              <lgcl>        <lgcl>          <lgcl>      <lgcl>       <lgcl>
#>               FALSE         FALSE           FALSE       FALSE        FALSE
#>  has_composition has_custom_output
#>           <lgcl>            <lgcl>
#>            FALSE             FALSE
```
