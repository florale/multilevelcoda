# Simulate Data from Generator Specifications

Build a single-level, grouped, or longitudinal simulation design and
evaluate generator specifications into one shared data table.

## Usage

``` r
simulate_data(
  generators,
  n = NULL,
  n_groups = NULL,
  n_per_group = NULL,
  group_id = "group_id",
  obs_id = "obs_id",
  time_id = NULL,
  time_values = NULL,
  time_truncate = TRUE,
  seed = NULL
)
```

## Arguments

- generators:

  Non-empty named list of generator specifications created by `gen_*()`
  functions.

- n:

  Total number of observations for a single-level design. For grouped
  designs, `n` may be supplied instead of `n_per_group`; rows are
  distributed as evenly as possible across groups. Supplying both `n`
  and `n_per_group` is an error.

- n_groups:

  Number of groups. When `NULL`, a single-level design is created.

- n_per_group:

  Group sizes. May be a scalar/vector, a function of `n_groups`, or a
  count-distribution list such as
  `list(distribution = "poisson", lambda = 4, minimum = 1)`. Requires
  `n_groups` and cannot be combined with `n`.

- group_id, obs_id:

  Character scalars naming the group and observation index columns.

- time_id:

  Optional character scalar naming a time column.

- time_values:

  Optional time values. May be a vector shared by units or a function
  called per unit/group.

- time_truncate:

  Logical; when `TRUE`, vector `time_values` are truncated to each group
  size. When `FALSE`, every group must have the same length as
  `time_values`.

- seed:

  Optional random seed. When supplied, the caller's existing random seed
  is restored after simulation.

## Value

An object of class `mlsim_data`, a list with:

- `data`:

  A
  [`data.table::data.table()`](https://rdrr.io/pkg/data.table/man/data.table.html)
  containing design columns and generated variables.

- `metadata`:

  Simulation design metadata, seed, index metadata, and generator order.

- `generator_specs`:

  Generator specifications with simulator closures removed.

- `generator_metadata`:

  Per-generator metadata recorded during simulation.

## Details

Generators are evaluated in list order. Each generator can depend on
columns produced by earlier generators through the simulation context.
Generated column names must be unique and must not collide after
[`make.names()`](https://rdrr.io/r/base/make.names.html).

## Examples

``` r
sim <- simulate_data(
  n_groups = 3,
  n_per_group = 2,
  seed = 10,
  generators = list(
    group_x = gen_mvn("group_x", level = "level2", fixed_intercept = 0, residual_cov = 1),
    y = gen_poisson("y", fixed_intercept = log(2))
  )
)
sim$data
#>    group_id obs_id     group_x     y
#>       <int>  <int>       <num> <int>
#> 1:        1      1  0.01874617     1
#> 2:        1      2  0.01874617     1
#> 3:        2      1 -0.18425254     2
#> 4:        2      2 -0.18425254     2
#> 5:        3      1 -1.37133055     2
#> 6:        3      2 -1.37133055     2
summary(sim)
#> <summary.mlsim_data>
#> 
#> Design:
#>      n n_cols n_generated_cols n_groups group_id group_size_min
#>  <int>  <int>            <int>    <int>   <char>          <int>
#>      6      4                2        3 group_id              2
#>  group_size_median group_size_max obs_id time_id  seed n_generators
#>              <int>          <int> <char>  <char> <int>        <int>
#>                  2              2 obs_id    <NA>    10            2
#> 
#> Generators:
#>  generator distribution  level    vars n_vars parameter_level parameter_count
#>     <char>       <char> <char>  <char>  <int>          <char>           <int>
#>    group_x          mvn level2 group_x      1           group               3
#>          y      poisson single       y      1             row               6
#>  has_row_parameters has_group_parameters has_fixed_parameters has_random_cov
#>              <lgcl>               <lgcl>               <lgcl>         <lgcl>
#>               FALSE                 TRUE                 TRUE          FALSE
#>                TRUE                FALSE                 TRUE          FALSE
#>  has_random_effects has_residuals has_scale_model has_composition
#>              <lgcl>        <lgcl>          <lgcl>          <lgcl>
#>               FALSE          TRUE           FALSE           FALSE
#>               FALSE         FALSE           FALSE           FALSE
#>  has_custom_output
#>             <lgcl>
#>              FALSE
#>              FALSE

longitudinal <- simulate_data(
  n_groups = 2,
  n_per_group = 3,
  time_id = "visit",
  generators = list(x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1))
)
longitudinal$metadata$index
#> $longitudinal
#> [1] TRUE
#> 
#> $unit_id
#> [1] "group_id"
#> 
#> $obs_id
#> [1] "obs_id"
#> 
#> $time_id
#> [1] "visit"
#> 
#> $order
#> [1] "group_id" "visit"   
#> 
#> $time_type
#> [1] "integer"
#> 
#> $time_unique_within_unit
#> [1] TRUE
#> 
```
