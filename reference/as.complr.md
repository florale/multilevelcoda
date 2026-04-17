# Coerce a list to a `complr` object

This is a constructor function for a `complr` object. It checks that the
input list has the necessary structure and components to be considered a
`complr` object, and if so, it assigns the class "complr" to it. This
allows for method dispatch on `complr` objects.

## Usage

``` r
as.complr(x)
```

## Arguments

- x:

  A list with elements `output`, `datain`, `dataout`, `transform`, and
  `idvar`.

## Value

A `complr` object.
