# Reference Grid for `substitution` model.

Build a dataset for `fitted.brmcoda` used in `substitution` model

## Usage

``` r
build.rg(object, ref, at, parts, level, weight, fill = FALSE)
```

## Arguments

- object:

  A fitted
  [`brmcoda`](https://florale.github.io/multilevelcoda/reference/brmcoda.md)
  object.

- ref:

  Either a character value or vector or a dataset. Can be `"grandmean"`
  and/or `"clustermean"`, or a `data.frame` or `data.table` of user's
  specified reference grid consisting of combinations of covariates over
  which predictions are made. User's specified reference grid is only
  possible for simple substitution. Single level models are default to
  `"grandmean"`.

- at:

  An optional named list of levels for the corresponding variables in
  the reference grid.

- parts:

  A optional character string specifying names of compositional parts
  that should be considered in the substitution analysis. This should
  correspond to a single set of names of compositional parts specified
  in the `complr` object. Default to the first composition in the
  `complr` object.

- level:

  A character string or vector. Should the estimate of multilevel models
  focus on the `"between"` and/or `"within"` or `"aggregate"` variance?
  Single-level models are default to `"aggregate"`.

- weight:

  A character value specifying the weight to use in calculation of the
  reference composition.

- fill:

  Logical value only relevant when `ref` is an user's specified
  reference grid in which information about some, but not all covariates
  is provided (e.g., models including age and sex as covariate but only
  age was provided in the reference grid). If `TRUE`, the unspecified
  covariates are filled with the default reference grid. If `FALSE`,
  users will be asked to provide a full reference grid. Currently only
  support the default to `FALSE`.

## Value

A reference grid consisting of a combination of covariates in `brmcoda`
