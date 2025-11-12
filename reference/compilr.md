# Indices from a (dataset of) Multilevel Composition(s) (deprecated.)

Indices from a (dataset of) Multilevel Composition(s) (deprecated.)

## Usage

``` r
compilr(...)
```

## Arguments

- ...:

  arguments passed to
  [`complr`](https://florale.github.io/multilevelcoda/reference/complr.md).

## Value

A
[`complr`](https://florale.github.io/multilevelcoda/reference/complr.md)
object with at least the following elements.

- `X`:

  A vector of class `acomp` representing one closed composition or a
  matrix of class `acomp` representing multiple closed compositions each
  in one row.

- `bX`:

  A vector of class `acomp` representing one closed between-person
  composition or a matrix of class `acomp` representing multiple closed
  between-person compositions each in one row.

- `wX`:

  A vector of class `acomp` representing one closed within-person
  composition or a matrix of class `acomp` representing multiple closed
  within-person compositions each in one row.

- `Z`:

  Log ratio transform of composition.

- `bZ`:

  Log ratio transform of between-person composition.

- `wZ`:

  Log ratio transform of within-person composition.

- `data`:

  The user's dataset or imputed dataset if the iiut data contains zeros.

- `transform`:

  Type of transform applied on compositional data.

- `parts`:

  Names of compositional variables.

- `idvar`:

  Name of the variable containing IDs.

- `total`:

  Total amount to which the compositions is closed.

## See also

[`complr`](https://florale.github.io/multilevelcoda/reference/complr.md)
