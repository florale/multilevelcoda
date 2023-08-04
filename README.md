

# multilevelcoda
<!-- badges: start -->
[![R-CMD-check](https://github.com/florale/multilevelcoda/workflows/R-CMD-check/badge.svg)](https://github.com/florale/multilevelcoda/actions)
[![Coverage Status](https://codecov.io/gh/florale/multilevelcoda/branch/main/graphs/badge.svg?branch=main)](https://app.codecov.io/gh/florale/multilevelcoda)
[![CRAN Version](https://www.r-pkg.org/badges/version/multilevelcoda)](https://cran.r-project.org/package=multilevelcoda)
[![lifecycle](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

This package provides functions to model compositional data in 
a multilevel framework using full Bayesian inference.
It integrates the principes of Compositional Data Analysis (CoDA) 
and Multilevel Modelling and supports both compositional data as 
an outcome and predictors in a wide range of 
generalized (non-)linear multivariate multilevel models.

## Installation
To install the latest release version from CRAN, run

```r 
install.packages("multilevelcoda")

```

The current developmental version can be downloaded from github via

```r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("florale/multilevelcoda")
```

Because `multilevelcoda` is built on brms, which is based on Stan, a C++ compiler is required. 
The program Rtools (available on https://cran.r-project.org/bin/windows/Rtools/) comes with a C++ compiler for Windows. On Mac, Xcode is required. For further instructions on how to get the compilers running, see the prerequisites section on https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started.

## Resources

You can learn about the package from these vignettes:

- [Introduction to Compositional Multilevel Modelling](https://florale.github.io/multilevelcoda/articles/introduction.html)
- [Multilevel Models with Compositional Outcomes](https://florale.github.io/multilevelcoda/articles/comp-outcome.html)
- [Multilevel Models with Compositional Predictors](https://florale.github.io/multilevelcoda/articles/comp-predictor.html)
- [Compositional Multilevel Substitution Models](https://florale.github.io/multilevelcoda/articles/substitution-model.html)

## Citing `multilevelcoda` and related software 
TBA

