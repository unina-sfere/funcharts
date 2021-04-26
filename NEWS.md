# funcharts 1.0.0.9000

## Major changes

## Minor changes

* Added a `NEWS.md` file to track changes to the package.
* `inprod_mfd_diag()` calculates the inner product between two multivariate functional data objects observation by observation, avoiding calculating it between all possible couples of observations. Therefore, there are $n$ calculations instead of $n^2$, saving much computational time when calculating the squared prediction error statistic when $n$ is large.


# funcharts 1.0.0

* Initial release
