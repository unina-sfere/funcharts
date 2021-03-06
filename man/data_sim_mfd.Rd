% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00_mfd.R
\name{data_sim_mfd}
\alias{data_sim_mfd}
\title{Simulate multivariate functional data}
\usage{
data_sim_mfd(nobs = 5, nbasis = 4, nvar = 2, seed = 0)
}
\arguments{
\item{nobs}{Number of functional observations to be simulated.}

\item{nbasis}{Number of basis functions.}

\item{nvar}{Number of functional covariates.}

\item{seed}{Set seed for reproducibility.}
}
\value{
A simulated object of class `mfd`.
}
\description{
Simulate random coefficients and create a multivariate functional data
object of class `mfd`.
}
\examples{
library(funcharts)
data_sim_mfd()
}
