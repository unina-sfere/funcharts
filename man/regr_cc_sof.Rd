% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_phaseII.R
\name{regr_cc_sof}
\alias{regr_cc_sof}
\title{Scalar-on-Function Regression Control Chart}
\usage{
regr_cc_sof(object, y_new, mfdobj_x_new, alpha = 0.05)
}
\arguments{
\item{object}{A list obtained as output from \code{sof_pc},
i.e. a fitted scalar-on-function linear regression model.}

\item{y_new}{A numeric vector containing the observations of
the scalar response variable
in the phase II data set.}

\item{mfdobj_x_new}{An object of class \code{mfd} containing
the phase II data set of the functional covariates observations.}

\item{alpha}{A numeric value indicating the Type I error for
the regression control chart
and such that this function returns the \code{1-alpha} prediction interval
on the response.
Default is 0.05.}
}
\value{
A \code{data.frame} with as many rows as the
number of functional replications in \code{mfdobj_x_new},
with the following columns:

* \code{y_hat}: the predictions of the response variable
corresponding to \code{mfdobj_x_new},

* \code{y}: the same as the argument \code{y_new} given as input
to this function,

* \code{lwr}: lower limit of the \code{1-alpha} prediction interval
on the response,

* \code{pred_err}: prediction error calculated as \code{y-y_hat},

* \code{pred_err_sup}: upper limit of the \code{1-alpha} prediction interval
on the prediction error,

* \code{pred_err_inf}: lower limit of the \code{1-alpha} prediction interval
on the prediction error.
}
\description{
This function builds a data frame needed
to plot the scalar-on-function regression control chart,
based on a fitted function-on-function linear regression model and
proposed in Capezza et al. (2020) together with the Hotelling's T^2 and
squared prediction error control charts.
The training data have already been used to fit the model.
A tuning data set can be provided that is used to estimate
the control chart limits.
A phase II data set contains the observations to be monitored
with the built control charts.
}
\examples{
library(funcharts)
air <- lapply(air, function(x) x[1:100, , drop = FALSE])
fun_covariates <- c("CO", "temperature")
mfdobj_x <- get_mfd_list(air[fun_covariates],
                         n_basis = 15,
                         lambda = 1e-2)
y <- rowMeans(air$NO2)
y1 <- y[1:80]
y2 <- y[81:100]
mfdobj_x1 <- mfdobj_x[1:80]
mfdobj_x2 <- mfdobj_x[81:100]
mod <- sof_pc(y1, mfdobj_x1)
cclist <- regr_cc_sof(object = mod,
                      y_new = y2,
                      mfdobj_x_new = mfdobj_x2)
plot_control_charts(cclist)

}
\references{
Capezza C, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
Control charts for
monitoring ship operating conditions and CO2 emissions
based on scalar-on-function regression.
\emph{Applied Stochastic Models in Business and Industry},
36(3):477--500.
<doi:10.1002/asmb.2507>
}
