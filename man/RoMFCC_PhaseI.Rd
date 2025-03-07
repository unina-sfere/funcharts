% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmfcc.R
\name{RoMFCC_PhaseI}
\alias{RoMFCC_PhaseI}
\title{Robust Multivariate Functional Control Charts - Phase I}
\usage{
RoMFCC_PhaseI(
  mfdobj,
  mfdobj_tuning = NULL,
  functional_filter_par = list(filter = TRUE),
  imputation_par = list(method_imputation = "RoMFDI"),
  pca_par = list(fev = 0.7),
  alpha = 0.05
)
}
\arguments{
\item{mfdobj}{A multivariate functional data object of class mfd.
A functional filter is applied to this data set, then
flagged functional componentwise outliers are imputed in the
robust imputation step.
Finally robust multivariate functional principal component analysis
is applied to the imputed data set for dimension reduction.}

\item{mfdobj_tuning}{An additional functional data object of class mfd.
After applying the filter and imputation steps on this data set,
it is used to robustly estimate the distribution of the Hotelling's T2 and
SPE statistics in order to calculate control limits
to prevent overfitting issues that could reduce the
monitoring performance of the RoMFCC.
Default is NULL, but it is strongly recommended to use a tuning data set.}

\item{functional_filter_par}{A list with an argument \code{filter} that can be TRUE or FALSE depending
on if the functional filter step must be performed or not.
All the other arguments of this list are passed as arguments to the function
\code{functional_filter} in the filtering step.
All the arguments that are not passed take their default values.
See \code{\link{functional_filter}} for all the arguments and their default
values.
Default is \code{list(filter = TRUE)}.}

\item{imputation_par}{A list with an argument \code{method_imputation}
that can be \code{"RoMFDI"} or \code{"mean"} depending
on if the imputation step must be done by means of \code{\link{RoMFDI}} or
by just using the mean of each functional variable.
If \code{method_imputation = "RoMFDI"},
all the other arguments of this list are passed as arguments to the function
\code{RoMFDI} in the imputation step.
All the arguments that are not passed take their default values.
See \code{\link{RoMFDI}} for all the arguments and their default
values.
Default value is \code{list(method_imputation = "RoMFDI")}.}

\item{pca_par}{A list with an argument \code{fev}, indicating a number between 0 and 1
denoting the fraction of variability that must be explained by the
principal components to be selected in the RoMFPCA step.
All the other arguments of this list are passed as arguments to the function
\code{rpca_mfd} in the RoMFPCA step.
All the arguments that are not passed take their default values.
See \code{\link{rpca_mfd}} for all the arguments and their default
values.
Default value is \code{list(fev = 0.7)}.}

\item{alpha}{The overall nominal type-I error probability used to set
control chart limits.
Default value is 0.05.}
}
\value{
A list of the following elements that are needed in Phase II:
\itemize{
\item \code{T2} the Hotelling's T2 statistic values for the Phase I data set,
\item \code{SPE} the SPE statistic values for the Phase I data set,
\item \code{T2_tun} the Hotelling's T2 statistic values for the tuning data set,
\item \code{SPE_tun} the SPE statistic values for the tuning data set,
\item \code{T2_lim} the Phase II control limit of
the Hotelling's T2 control chart,
\item \code{spe_lim} the Phase II control limit of
the SPE control chart,
\item \code{tuning} TRUE if the tuning data set is provided, FALSE otherwise,
\item \code{mod_pca} the final RoMFPCA model fitted on the Phase I data set,
\item \code{K} = K the number of selected principal components,
\item \code{T_T2_inv} if a tuning data set is provided,
it returns the inverse of the covariance matrix
of the first \code{K} scores, needed to calculate the Hotelling's T2
statistic for the Phase II observations.
\item \code{mean_scores_tuning_rob_mean} if a tuning data set is provided,
it returns the robust location estimate of the scores, needed to calculate
the Hotelling's T2 and SPE
statistics for the Phase II observations.
}
}
\description{
It performs Phase I of the Robust Multivariate Functional Control Chart
(RoMFCC) as proposed by Capezza et al. (2024).
}
\examples{
\dontrun{
library(funcharts)
mfdobj <- get_mfd_list(air, n_basis = 5)
nobs <- dim(mfdobj$coefs)[2]
set.seed(0)
ids <- sample(1:nobs)
mfdobj1 <- mfdobj[ids[1:100]]
mfdobj_tuning <- mfdobj[ids[101:300]]
mfdobj2 <- mfdobj[ids[-(1:300)]]
mod_phase1 <- RoMFCC_PhaseI(mfdobj = mfdobj1,
                            mfdobj_tuning = mfdobj_tuning)
phase2 <- RoMFCC_PhaseII(mfdobj_new = mfdobj2,
                         mod_phase1 = mod_phase1)
plot_control_charts(phase2)
}

}
\references{
Capezza, C., Centofanti, F., Lepore, A., Palumbo, B. (2024)
Robust Multivariate Functional Control Chart.
\emph{Technometrics}, 66(4):531--547, \url{doi:10.1080/00401706.2024.2327346}.
}
\author{
C. Capezza, F. Centofanti
}
