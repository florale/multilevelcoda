# multilevelcoda 0.0.9000

### Initial release

# multilevelcoda 1.0.0

### CRAN release

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

