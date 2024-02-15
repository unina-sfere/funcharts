# funcharts (development version)

# funcharts 1.4.0

* `rpca_mfd()` performs multivariate functional principal component analysis as described in Capezza et al. (2022).
* `functional_filter()` performs the functional filtering step of the robust multivariate functional control chart framework of Capezza et al. (2022).
* `RoMFDI()` performs the robust multivariate functional data imputation step as described in Capezza et al. (2022).
* `RoMFCC_PhaseI()` performs Phase I of the robust multivariate functional control chart framework of Capezza et al. (2022).
* `RoMFCC_PhaseII()` performs Phase II of the robust multivariate functional control chart framework of Capezza et al. (2022).

References:

* Capezza, C., Centofanti, F., Lepore, A., Palumbo, B. (2022) Robust Multivariate Functional Control Charts. arXiv:2207.07978v.

# funcharts 1.3.2

* Updated documentation with the newly published paper, which thoroughly illustrates the funcharts package:
Capezza C, Centofanti F, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2023) funcharts: control charts for multivariate functional data in R.
Journal of Quality Technology, doi:10.1002/asmb.2507.

# funcharts 1.3.1

* the default value of the `parametric_limits` argument in `regr_cc_sof()` is now set to `FALSE`.

# funcharts 1.3.0

* all basis function systems that can be used in the `fda` package now can be used also with `funcharts`, which previously it could be used only with B-spline basis.
In particular, Fourier, exponential, monomial, polygonal, power and constant basis function systems are available.
* `get_outliers_mfd()` allows to find outliers among multivariate functional data using the functional boxplot through the `fbplot()` function of the `roahd` package.
* test coverage has been increased
* `control_charts_sof()` and `control_charts_sof_real_time()` have been deprecated.
Instead, use `regr_cc_sof()` and `regr_cc_sof_real_time()`, respectively, with argument `include_covariates = TRUE`. 
This has been done to make more consistent the regression control chart functions for the scalar (`regr_cc_sof()` and `regr_cc_sof_real_time()`) and functional (`regr_cc_fof()` and `regr_cc_fof_real_time()`) response cases.
* `alpha` parameter in all control charting functions, which previously could only be a list with manually specified values of the type-I error probability in each control chart, now can also be a single number between 0 and 1. In this case, Bonferroni correction is automatically applied to take into account the multiplicity problem when more than one control chart is applied.
* `plot_bifd()` now allows to choose to produce also contour or perspective plots of `bifd` objects.
* `simulate_mfd()` is much more general, now it allows to simulate as many covariates as one wants (before the number was fixed to three), it is possible to provide manually the mean and variance function for each variable, it is possible to select the type of correlation function for each variable.
* `plot_mfd()` now relies on patchwork, while the new function `lines_mfd()` allows to add new curve to an existing plot.

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
* `inprod_mfd_diag()` calculates the inner product between two multivariate functional data objects observation by observation, avoiding calculating it between all possible couples of observations. Therefore, there are n calculations instead of squared n, saving much computational time when calculating the squared prediction error statistic when n is large.
* Code has been improved so that `scale_mfd()` is pre-computed and therefore is not called many times unnecessarily along the different functions.


# funcharts 1.0.0

* Initial release
