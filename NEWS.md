# multilevelcoda 1.2.2

### New Features

* Add a new function `multilevelcoda_sim()` to bundle simulation study results in a `shiny` app.
* Add more methods for `brmcoda` (`loo` and `pp_check`).

# multilevelcoda 1.2.1

### Other Changes

* Fix a bug in the `summary` argument of `predict()` and `fitted()`.

# multilevelcoda 1.2.0

### New Features

* Add various methods for `compilr`, `brmcoda`, and `substitution` objects.
* Change the default `total` in `compilr` to `1` instead of `1440`.

### Other Changes

* Fix bugs in vignettes and reformat.
* Add a new vignette to present a case study in `multilevelcoda` simulation study.

# multilevelcoda 1.1.0

### New Features

* Add method `update` to support update an existing `compilr` or `brmcoda` objects.
* Add method `summary` to `substitution` to summarise output produced by function `substitution`.
* Add argument `weight` to `substitution` to allow for different methods for calculating reference compostion.
* Add user control over specifying reference composition and reference grid in `substitution` model.
* Ensure that `complir` throws errors for missing data and 0s in input compositional data.

### Other Changes

* Deprecate argument `type` of `substitution` and replace with `ref`.
* Deprecate argument `regrid` of `substitution`

# multilevelcoda 1.0.0

### First CRAN release

# multilevelcoda 0.0.9000

### Initial release
