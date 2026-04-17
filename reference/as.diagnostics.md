# Coerce a list to a `diagnostics` object

This is a constructor function for a `diagnostics` object. It checks
that the input list has the necessary structure and components to be
considered a `diagnostics` object, and if so, it assigns the class
"diagnostics" to it. This allows for method dispatch on `diagnostics`
objects.

## Usage

``` r
as.diagnostics(x)
```

## Arguments

- x:

  A list with elements `x`, `distance`, `cutoff`, `ev.perc`,
  `extremevalues`, and `levels`.

## Value

A `diagnostics` object.
