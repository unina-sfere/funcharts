% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_phaseII.R
\name{control_charts_pca}
\alias{control_charts_pca}
\title{T^2 and SPE control charts for multivariate functional data}
\usage{
control_charts_pca(
  pca,
  components,
  tuning_data = NULL,
  newdata,
  alpha = list(T2 = 0.025, spe = 0.025),
  limits = "standard",
  seed = 0,
  nfold = 5,
  ncores = 1
)
}
\arguments{
\item{pca}{An object of class \code{pca_mfd}
obtained by doing MFPCA on the
training set of multivariate functional data.}

\item{components}{A vector of integers with the components over which
to project the multivariate functional data.}

\item{tuning_data}{An object of class \code{mfd} containing
the tuning set of the multivariate functional data, used to estimate the
T^2 and SPE control chart limits.
If NULL, the training data, i.e. the data used to fit the MFPCA model,
are also used as the tuning data set, i.e. \code{tuning_data=pca$data}.
Default is NULL.}

\item{newdata}{An object of class \code{mfd} containing
the phase II set of the multivariate functional data to be monitored.}

\item{alpha}{A named list with two elements, named \code{T2} and \code{spe},
respectively, each containing
the desired Type I error probability of
the corresponding control chart.
Note that at the moment you have to take into account manually
the family-wise error rate and adjust
the two values accordingly. See Capezza et al. (2020) and
Centofanti et al. (2020)
for additional details. Default value is
\code{list(T2 = 0.025, spe = 0.025)}.}

\item{limits}{A character value.
If "standard", it estimates the control limits on the tuning
data set. If "cv", the function calculates the control limits only on the
training data using cross-validation
using \code{calculate_cv_limits}. Default is "standard".}

\item{seed}{If \code{limits=="cv"},
since the split in the k groups is random,
you can fix a seed to ensure reproducibility.
Otherwise, this argument is ignored.}

\item{nfold}{If \code{limits=="cv"}, this gives the number of groups k
used for k-fold cross-validation.
If it is equal to the number of observations in the training data set,
then we have
leave-one-out cross-validation.
Otherwise, this argument is ignored.}

\item{ncores}{If \code{limits=="cv"}, if you want perform the analysis
in the k groups in parallel,
give the number of cores/threads.
Otherwise, this argument is ignored.}
}
\value{
A \code{data.frame} with as many rows as the number of
multivariate functional observations in the phase II data set and
the following columns:

* one \code{id} column identifying the multivariate functional observation
in the phase II data set,

* one \code{T2} column containing the Hotelling T^2 statistic
calculated for all observations,

* one column per each functional variable,
containing its contribution to the T^2 statistic,

* one \code{spe} column containing the SPE statistic calculated
for all observations,

* one column per each functional variable,
containing its contribution to the SPE statistic,

* \code{T2_lim} gives the upper control limit of
the Hotelling's T^2 control chart,

* one \code{contribution_T2_*_lim} column per each
functional variable giving the
limits of the contribution of that variable
to the Hotelling's T^2 statistic,

* \code{spe_lim} gives the upper control limit of the SPE control chart

* one \code{contribution_spe*_lim} column per each
functional variable giving the
limits of the contribution of that variable to the SPE statistic.
}
\description{
This function builds a data frame needed to plot
the Hotelling's T^2 and squared prediction error (SPE)
control charts
based on multivariate functional principal component analysis
(MFPCA) performed
on multivariate functional data,
proposed in Capezza et al. (2020) together with the scalar control chart
and used also to build the
functional regression control chart proposed in Centofanti et al. (2020)
(this function is used by \code{\link{regr_cc_fof}}).
The training data have already been used to fit the model.
A tuning data set can be provided that is used to estimate
the control chart limits.
A phase II data set contains the observations
to be monitored with the built control charts.
}
\references{
Capezza C, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
Control charts for
monitoring ship operating conditions and CO2 emissions
based on scalar-on-function regression.
\emph{Applied Stochastic Models in Business and Industry},
36(3):477--500.
<doi:10.1002/asmb.2507>

Centofanti F, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
Functional Regression Control Chart.
\emph{Technometrics}. <doi:10.1080/00401706.2020.1753581>

library(funcharts)
data("air")
air <- lapply(air, function(x) x[1:220, , drop = FALSE])
fun_covariates <- c("CO", "temperature")
mfdobj_x <- get_mfd_list(air[fun_covariates],
                         n_basis = 15,
                         lambda = 1e-2)
y <- rowMeans(air$NO2)
y1 <- y[1:100]
y_tuning <- y[101:200]
y2 <- y[201:220]
mfdobj_x1 <- mfdobj_x[1:100]
mfdobj_x_tuning <- mfdobj_x[101:200]
mfdobj_x2 <- mfdobj_x[201:220]
pca <- pca_mfd(mfdobj_x1)
components <- 1:which(cumsum(pca$varprop) >= .90)[1]
cclist <- control_charts_pca(pca = pca,
                             components = components,
                             tuning_data = mfdobj_x_tuning,
                             newdata = mfdobj_x2)
plot_control_charts(cclist)
}
\seealso{
\code{\link{regr_cc_fof}}
}
