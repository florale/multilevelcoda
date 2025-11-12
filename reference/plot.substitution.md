# Substitution Plot

Make a plot of
[`substitution`](https://florale.github.io/multilevelcoda/reference/substitution.md)
model results.

## Usage

``` r
# S3 method for class 'substitution'
plot(x, to, ref, level, ...)
```

## Arguments

- x:

  A
  [`substitution`](https://florale.github.io/multilevelcoda/reference/substitution.md)
  class object.

- to:

  An optional character value or vector specifying the names of the
  compositional parts that were reallocated to in the model.

- ref:

  A character value of ((`"grandmean"` or `"clustermean"` or `"users"`),

- level:

  An optional character value of (`"between"`, `"within"`), or
  `"aggregate"`).

- ...:

  Further components to the plot, followed by a plus sign (+).

## Value

A ggplot graph object showing the estimated difference in outcome when
each pair of compositional variables are substituted for a specific
time.
