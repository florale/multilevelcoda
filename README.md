<!-- badges: start -->
[![lifecycle](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

This package has functions to support estimating and presenting
compositional data analysis (CoDA) with multilevel models.
It supports both compositional data as an outcome using 
multivariate models and compositional data as predictors
of other multilevel outcomes.

## Installation

For the latest, development version, you can use the 
`remotes` package to install it from GitHub:

```r
remotes::install_github("florale/multilevelcoda")
```

If you do not have the `remotes` package installed, 
you can install it from CRAN using this code first before
running the above code to install `multilevelcoda`

```r
install.packages("remotes")
```

Release of `multilevelcoda` to CRAN is planned at a later date.

## Learn More

You can see more examples and learn about the package from these vignettes:

- [Introduction to Introduction to Compositional Multilevel Modelling](https://florale.github.io/multilevelcoda/articles/introduction.html)
- [Multilevel Model with Compositional Outcomes](https://florale.github.io/multilevelcoda/articles/comp-outcome.html)
- [Multilevel Model with Compositional Predictors](https://florale.github.io/multilevelcoda/articles/comp-predictor.html)
- [Compositional Multilevel Substitution Model](https://florale.github.io/multilevelcoda/articles/substitution-model.html)