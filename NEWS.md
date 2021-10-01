# funcharts 1.1.0

## Major changes

* `simulate_mfd()` simulates example data for `funcharts`. 
It creates a data set with three functional covariates, a functional response generated as a function of the three functional covariates through a function-on-function linear model, and a scalar response generated as a function of the three functional covariates through a scalar-on-function linear model. This function covers the simulation study in Centofanti et al. (2020) for the function-on-function case and also simulates data in a similar way for the scalar response case.

## Minor changes

* Added a `NEWS.md` file to track changes to the package.
* `inprod_mfd_diag()` calculates the inner product between two multivariate functional data objects observation by observation, avoiding calculating it between all possible couples of observations. Therefore, there are $n$ calculations instead of $n^2$, saving much computational time when calculating the squared prediction error statistic when $n$ is large.
* Code has been improved so that `scale_mfd()` is pre-computed and therefore is not called many times unnecessarily along the different functions.


# funcharts 1.0.0

* Initial release
