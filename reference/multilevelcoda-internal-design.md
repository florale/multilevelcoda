# Internal Simulation Design Helpers

Validate generator lists, build base simulation designs, derive
time/index metadata, and remove simulator closures before storing
specifications on `mlsim_data` objects.

## Usage

``` r
.mlsim_drop_simulators(generators)

.mlsim_validate_generators(generators)

.mlsim_build_design(
  n,
  n_groups,
  n_per_group,
  group_id,
  obs_id,
  time_id = NULL,
  time_values = NULL,
  time_truncate = TRUE
)

.mlsim_check_level_possible(spec, design)

.mlsim_time_values_for_unit(time_values, n, group = NULL)

.mlsim_group_time_values(time_values, group_sizes, time_truncate = FALSE)

.mlsim_call_time_values(time_values, n, group)

.mlsim_check_time_values(values, n)

.mlsim_check_unique_time(data, unit_id, time_id)

.mlsim_time_values_metadata(time_values, time_truncate = FALSE)

.mlsim_index_metadata(data, grouped, group_id, obs_id, time_id)
```

## Arguments

- generators:

  Named list of generator specifications.

- n:

  Number of observations or requested generated values.

- n_groups:

  Number of groups.

- n_per_group:

  Group-size specification.

- group_id, obs_id, time_id:

  Design column names.

- time_values:

  Optional vector or function used to create time columns.

- time_truncate:

  Logical flag controlling vector `time_values` recycling.

- spec:

  Generator specification.

- design:

  Internal design object.

- group:

  Optional group identifier passed to time-value functions.

- group_sizes:

  Integer vector of group sizes.

- values:

  Candidate time values.

- data:

  Data table containing design columns.

- unit_id:

  Unit identifier used for uniqueness checks.

- grouped:

  Logical; whether the design is grouped.

## Value

Internal helper outputs used by
[`simulate_data()`](https://florale.github.io/multilevelcoda/reference/simulate_data.md).

## Examples

``` r
design <- multilevelcoda:::.mlsim_build_design(
  n = NULL,
  n_groups = 2,
  n_per_group = 2,
  group_id = "group_id",
  obs_id = "obs_id"
)
design$data
#>    group_id obs_id
#>       <int>  <int>
#> 1:        1      1
#> 2:        1      2
#> 3:        2      1
#> 4:        2      2
multilevelcoda:::.mlsim_time_values_for_unit(NULL, 3)
#> [1] 1 2 3
```
