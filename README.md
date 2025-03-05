
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/funcharts)](https://CRAN.R-project.org/package=funcharts)
[![R-CMD-check](https://github.com/unina-sfere/funcharts/workflows/R-CMD-check/badge.svg)](https://github.com/unina-sfere/funcharts/actions)
[![Codecov test
coverage](https://codecov.io/gh/unina-sfere/funcharts/branch/main/graph/badge.svg)](https://app.codecov.io/gh/unina-sfere/funcharts?branch=main)
<!-- badges: end -->

# funcharts

The goal of `funcharts` is to provide control charts for the statistical
process monitoring of multivariate functional data densely observed on
one-dimensional intervals. The package is thoroughly illustrated in the
paper of Capezza et al. (2023). The package provides the methodologies
proposed in Colosimo and Pacella (2010), Capezza et al. (2020),
Centofanti et al. (2021), Capezza et al. (2024a), and Capezza et
al. (2024b). Moreover, this package provides a new class `mfd` for
multivariate functional data that is a wrapper of the class `fd` of the
package `fda`. See the
[`vignette("mfd", package = "funcharts")`](https://unina-sfere.github.io/funcharts/articles/mfd.html).

In particular:

- Colosimo and Pacella (2010) propose control charts for monitoring
  functional data based on functional principal component analysis.
  [`vignette("colosimo2010", package = "funcharts")`](https://unina-sfere.github.io/funcharts/articles/colosimo2010.html)
- Capezza et al. (2020) propose control charts for monitoring a scalar
  response variable and functional covariates using scalar-on-function
  regression. See the
  [`vignette("capezza2020", package = "funcharts")`](https://unina-sfere.github.io/funcharts/articles/capezza2020.html).
- Centofanti et al. (2021) propose the functional regression control
  chart (FRCC), i.e. control charts for monitoring a functional response
  variable conditionally on multivariate functional covariates. See the
  [`vignette("centofanti2021", package = "funcharts")`](https://unina-sfere.github.io/funcharts/articles/centofanti2021.html).
- Capezza et al. (2024a) propose the adaptive multivariate functional
  EWMA (AMFEWMA) control chart.
- Capezza et al. (2024b) propose the robust multivariate functional
  control chart (RoMFCC).
- Centofanti et al. (2024) propose the functional real-time monitoring
  (FRTM) control chart.

## Installation

You can install the CRAN version of the R package `funcharts` by doing:

``` r
install.packages("funcharts")
```

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("unina-sfere/funcharts")
```

# References

- Capezza C, Centofanti F, Lepore A, Menafoglio A, Palumbo B,
  Vantini S. (2023) funcharts: control charts for multivariate
  functional data in R. *Journal of Quality Technology*,
  <doi:10.1080/00224065.2023.2219012>.
- Capezza C, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
  Control charts for monitoring ship operating conditions and
  CO<sub>2</sub> emissions based on scalar-on-function regression.
  *Applied Stochastic Models in Business and Industry*, 36(3):477–500,
  <doi:10.1002/asmb.2507>
- Capezza, C., Capizzi, G., Centofanti, F., Lepore, A., Palumbo, B.
  (2024a) An Adaptive Multivariate Functional EWMA Control Chart. To
  appear in *Journal of Quality Technology*,
  <doi:https://doi.org/10.1080/00224065.2024.2383674>.
- Capezza, C., Centofanti, F., Lepore, A., Palumbo, B. (2024b) Robust
  Multivariate Functional Control Charts. *Technometrics*,
  66(4):531–547, <doi:10.1080/00401706.2024.2327346>.
- Centofanti F, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2021)
  Functional Regression Control Chart. *Technometrics*, 63(3):281--294,
  <doi:10.1080/00401706.2020.1753581>.
- Centofanti, F., A. Lepore, M. Kulahci, and M. P. Spooner (2024).
  Real-time monitoring of functional data. Accepted for publication in
  *Journal of Quality Technology*.
- Colosimo BM, Pacella, M. (2010) A comparison study of control charts
  for statistical monitoring of functional data. *International Journal
  of Production Research*, 48(6), 1575-1601,
  <doi:10.1080/00207540802662888>.
