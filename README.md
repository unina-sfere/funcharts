
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/funcharts)](https://CRAN.R-project.org/package=funcharts) -->

[![R-CMD-check](https://github.com/unina-sfere/funcharts/workflows/R-CMD-check/badge.svg)](https://github.com/unina-sfere/funcharts/actions)
[![Build
Status](https://travis-ci.com/unina-sfere/funcharts.svg?branch=main)](https://travis-ci.com/unina-sfere/funcharts)
<!-- badges: end -->

# funcharts

The goal of `funcharts` is to provide control charts for functional data
densely observed on one-dimensional intervals. The objective is to
provide the methodologies proposed in Capezza et al. (2020) and
Centofanti et al. (2020). Moreover, this package provide a new class
`mfd` for multivariate functional data that is a wrapper of the class
`fd` of the package `fda`. See the
[`vignette("mfd", package = "funcharts")`](https://unina-sfere.github.io/funcharts/articles/mfd.html).

In particular:

-   Capezza et al. (2020) proposed control charts for monitoring a
    scalar response variable and functional covariates using
    scalar-on-function regression. See the
    [`vignette("capezza2020", package = "funcharts")`](https://unina-sfere.github.io/funcharts/articles/capezza2020.html).
-   Centofanti et al. (2020) proposed the functional regression control
    chart, i.e. control charts for monitoring a functional response
    variable conditionally on multivariate functional covariates. See
    the
    [`vignette("centofanti2020", package = "funcharts")`](https://unina-sfere.github.io/funcharts/articles/centofanti2020.html).

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("unina-sfere/funcharts")
```

# References

-   Capezza C, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
    Control charts for monitoring ship operating conditions and
    CO<sub>2</sub> emissions based on scalar-on-function regression.
    *Applied Stochastic Models in Business and Industry*, 36(3):477–500.
    <https://doi.org/10.1002/asmb.2507>
-   Centofanti F, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
    Functional Regression Control Chart. *Technometrics*.
    <https://doi.org/10.1080/00401706.2020.1753581>
