# funcharts 1.2.0

* improved backward compatibility, `funcharts` now depends on an older version of R, i.e., >3.6.0 instead of >4.0.0
* `fof_pc()` now is much faster especially when the number of basis functions of the functional coefficient is large since the tensor product has been vectorized.
* the argument `seed` has been deprecated in all functions, so that reproducibility is achieved by setting externally a seed with `set.seed()`, as it is commonly done in R.
* `sim_funcharts()` simulates data sets automatically using the function `simulate_mfd()`. The only input required is the sample size for the Phase I, tuning and Phase II data sets.
* `control_charts_pca()` allows automatic selection of components.
* `get_mfd_list()` and `get_mfd_array()`, with the corresponding real time versions, are now much faster.
* cross-validation in scalar-on-function regression is now much faster, since the for loop is avoided
* inner products are more precise and much faster, because they rely on the pre-computed inner products of the B-spline basis functions, calculated via `inprod.bspline()`.
* argument `seed` is deprecated in all functions. Instead, a seed must be set before calling the functions by using `set.seed()`.

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
