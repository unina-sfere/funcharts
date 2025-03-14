% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FRTM_plot_functions.R
\name{plot.mFPCA}
\alias{plot.mFPCA}
\title{Plot the results of the  Mixed Functional Principal Component Analysis (mFPCA)}
\usage{
\method{plot}{mFPCA}(x, ...)
}
\arguments{
\item{x}{The output of   \code{mFPCA}.}

\item{...}{A variable \code{type} could be provided that can assume two values. If \code{type="single"}, the principal components are plotted separately.
If \code{type="all"}, the principal components are plotted together.}
}
\value{
No return value, called for side effects.
}
\description{
This function provides plots of the principal components of the mFPCA.
}
\examples{
library(funcharts)

data <- simulate_data_FRTM(n_obs = 100)
X <- sapply(1:100, function(ii)
  data$x_true[[ii]])
x_fd <-
  fda::smooth.basis(y = X,
                    argvals =  data$grid,
                    fda::create.bspline.basis(c(0, 1), 30))$fd
H <- sapply(1:100, function(ii)
  data$h[[ii]])
h_fd <-
  fda::smooth.basis(y = H,
                    argvals =  data$grid,
                    fda::create.bspline.basis(c(0, 1), 30))$fd
mod_mFPCA <- mFPCA(x_fd, h_fd, ncom = "ptv", par_ncom = 0.95)
plot(mod_mFPCA)

}
