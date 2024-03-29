% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_phaseI_training.R
\name{get_sof_pc_outliers}
\alias{get_sof_pc_outliers}
\title{Get possible outliers of a training data set of a
scalar-on-function regression model.}
\usage{
get_sof_pc_outliers(y, mfdobj)
}
\arguments{
\item{y}{A numeric vector containing the observations of the
scalar response variable.}

\item{mfdobj}{A multivariate functional data object of class mfd
denoting the functional covariates.}
}
\value{
A character vector with the ids of functional observations
signaled as possibly anomalous.
}
\description{
Get possible outliers of a training data set of a
scalar-on-function regression model.
It sets the training data set also as tuning data set for the
calculation of control chart limits,
and as phase II data set to compare monitoring statistics
against the limits and identify
possible outliers.
This is only an empirical approach. It is advised to use methods
appropriately designed for phase I monitoring to identify outliers.
}
\examples{
\dontrun{
library(funcharts)
data("air")
air <- lapply(air, function(x) x[1:10, , drop = FALSE])
fun_covariates <- c("CO", "temperature")
mfdobj_x <- get_mfd_list(air[fun_covariates], lambda = 1e-2)
y <- rowMeans(air$NO2)
get_sof_pc_outliers(y, mfdobj_x)
}

}
