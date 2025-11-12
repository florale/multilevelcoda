# Changelog

## multilevelcoda 1.3.3

CRAN release: 2025-11-11

#### New Features

- [`substitution()`](https://florale.github.io/multilevelcoda/reference/substitution.md)
  supports moderated substituion via `at` argument.
- [`complr()`](https://florale.github.io/multilevelcoda/reference/complr.md)
  now supports multiple composition input (see examples).
- [`substitution()`](https://florale.github.io/multilevelcoda/reference/substitution.md)
  supports compositional outcomes.
- [`substitution()`](https://florale.github.io/multilevelcoda/reference/substitution.md)
  output is now weighted across levels of the reference grid.

#### Other Changes

- Syntax for
  [`brmcoda()`](https://florale.github.io/multilevelcoda/reference/brmcoda.md)
  updated to support multiple compositions (see examples).

## multilevelcoda 1.2.2

#### New Features

- Add a new function
  [`multilevelcoda_sim()`](https://florale.github.io/multilevelcoda/reference/multilevelcoda_sim.md)
  to bundle simulation study results in a `shiny` app.
- Add more methods for `brmcoda` (`loo` and `pp_check`).

## multilevelcoda 1.2.1

CRAN release: 2023-08-13

#### Other Changes

- Fix a bug in the `summary` argument of
  [`predict()`](https://florale.github.io/multilevelcoda/reference/predict.brmcoda.md)
  and
  [`fitted()`](https://florale.github.io/multilevelcoda/reference/fitted.brmcoda.md).

## multilevelcoda 1.2.0

CRAN release: 2023-08-12

#### New Features

- Add various methods for `compilr`, `brmcoda`, and `substitution`
  objects.
- Change the default `total` in `compilr` to `1` instead of `1440`.

#### Other Changes

- Fix bugs in vignettes and reformat.
- Add a new vignette to present a case study in `multilevelcoda`
  simulation study.

## multilevelcoda 1.1.0

CRAN release: 2023-06-16

#### New Features

- Add method `update` to support update an existing `compilr` or
  `brmcoda` objects.
- Add method `summary` to `substitution` to summarise output produced by
  function `substitution`.
- Add argument `weight` to `substitution` to allow for different methods
  for calculating reference compostion.
- Add user control over specifying reference composition and reference
  grid in `substitution` model.
- Ensure that `complir` throws errors for missing data and 0s in input
  compositional data.

#### Other Changes

- Deprecate argument `type` of `substitution` and replace with `ref`.
- Deprecate argument `regrid` of `substitution`

## multilevelcoda 1.0.0

CRAN release: 2023-01-13

#### First CRAN release

## multilevelcoda 0.0.9000

#### Initial release
