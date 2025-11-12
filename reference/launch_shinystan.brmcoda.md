# Interface to shinystan

Provide an interface to shinystan for models fitted with brms

## Usage

``` r
# S3 method for class 'brmcoda'
launch_shinystan(object, ...)
```

## Arguments

- object:

  A fitted model object of class `brmcoda`.

- ...:

  Optional arguments to pass to
  [`launch_shinystan`](https://mc-stan.org/shinystan/reference/launch_shinystan.html)
  or [`runApp`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Value

An S4 shinystan object

## See also

[`launch_shinystan`](https://mc-stan.org/shinystan/reference/launch_shinystan.html)
