

# multilevelcoda
<!-- badges: start -->
[![R-CMD-check](https://github.com/florale/multilevelcoda/workflows/R-CMD-check/badge.svg)](https://github.com/florale/multilevelcoda/actions)
[![CRAN Version](https://www.r-pkg.org/badges/version/multilevelcoda)](https://cran.r-project.org/package=multilevelcoda)
[![lifecycle](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://cranlogs.r-pkg.org/badges/grand-total/multilevelcoda)](https://cran.r-project.org/package=multilevelcoda)
<!-- [![Coverage Status](https://codecov.io/gh/florale/multilevelcoda/branch/main/graphs/badge.svg?branch=main)](https://app.codecov.io/gh/florale/multilevelcoda)  -->
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

Because multilevelcoda is built on brms, which is based on Stan, a C++ compiler is required. 
The program Rtools (available on https://cran.r-project.org/bin/windows/Rtools/) comes with a C++ compiler for Windows. On Mac, Xcode is required. For further instructions on how to get the compilers running, see the prerequisites section on https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started.

## Resources

You can learn about the package from these vignettes:

- [Introduction to Compositional Multilevel Modelling](https://florale.github.io/multilevelcoda/articles/A-introduction.html)
- [Multilevel Models with Compositional Predictors](https://florale.github.io/multilevelcoda/articles/B-composition-MLM.html)
- [Multilevel Models with Compositional Outcome](https://florale.github.io/multilevelcoda/articles/C-composition-MMLM.html)
- [Compositional Substitution Multilevel Analysis](https://florale.github.io/multilevelcoda/articles/D-substitution.html)

## Citing `multilevelcoda` and related software 
When using multilevelcoda, please cite one or more of the following publications:

-   Le, F., Stanford, T. E., Dumuid, D., & Wiley, J. F. (2025). 
    Bayesian multilevel compositional data analysis: 
    introduction, evaluation, and application. 
    *Psychological Methods*. https://doi.org/10.1037/met0000750
-   Le F., Dumuid D., Stanford T. E., Wiley J. F. (2025). 
    Bayesian multilevel compositional data analysis with the R package multilevelcoda.
    *Multivariate Behavioral Research*. https://doi.org/10.1080/00273171.2025.2565598

As multilevelcoda depends on brms and Stan, please also consider citing:

-   Bürkner P. C. (2017). brms: An R Package for Bayesian Multilevel
    Models using Stan. *Journal of Statistical Software*. 80(1), 1-28.
    doi.org/10.18637/jss.v080.i01
-   Bürkner P. C. (2018). Advanced Bayesian Multilevel Modeling with the
    R Package brms. *The R Journal*. 10(1), 395-411.
    doi.org/10.32614/RJ-2018-017
-   Bürkner P. C. (2021). Bayesian Item Response Modeling in R with brms
    and Stan. *Journal of Statistical Software*, 100(5), 1-54.
    doi.org/10.18637/jss.v100.i05
-   Stan Development Team. YEAR. Stan Modeling Language Users Guide and
    Reference Manual, VERSION. <https://mc-stan.org>
-   Carpenter B., Gelman A., Hoffman M. D., Lee D., Goodrich B.,
    Betancourt M., Brubaker M., Guo J., Li P., and Riddell A. (2017).
    Stan: A probabilistic programming language. *Journal of Statistical
    Software*. 76(1). doi.org/10.18637/jss.v076.i01
