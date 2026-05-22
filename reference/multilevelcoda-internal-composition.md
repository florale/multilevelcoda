# Internal Composition Helpers

Validate sequential binary partitions, resolve composition part names,
and transform ILR coordinates back to closed composition parts.

## Usage

``` r
.mlsim_ilr_inverse(
  ilr,
  parts = NULL,
  total = 1,
  coordinate_names = NULL,
  sbp = NULL
)

.mlsim_resolve_composition(parts = NULL, sbp = NULL, d)

.mlsim_validate_sbp(sbp, d)

.mlsim_check_composition_output_names(coordinate_names, parts, keep_ilr)

.mlsim_ilr_coordinate_map(sbp, coordinate_names)
```

## Arguments

- ilr:

  Matrix or data frame of ILR coordinates.

- parts:

  Character vector naming composition parts.

- total:

  Positive scalar total used to close compositions.

- coordinate_names:

  Character vector naming ILR coordinates.

- sbp:

  Sequential binary partition matrix with part names in columns.

- d:

  Number of ILR coordinates.

- keep_ilr:

  Logical flag indicating whether ILR coordinates are emitted alongside
  composition parts.

## Value

Internal composition metadata or transformed values.

## Examples

``` r
resolved <- multilevelcoda:::.mlsim_resolve_composition(
  parts = c("a", "b", "c"),
  d = 2
)
resolved$parts
#> [1] "a" "b" "c"

inv <- multilevelcoda:::.mlsim_ilr_inverse(
  matrix(c(0, 0, 0.1, -0.1), nrow = 2, byrow = TRUE),
  parts = c("a", "b", "c")
)
rowSums(inv$values)
#> [1] 1 1
```
