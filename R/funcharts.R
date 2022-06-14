#' \code{funcharts} package
#'
#' Provides functional control charts
#' for statistical process monitoring of functional data,
#' using the methods of Capezza et al. (2020) <doi:10.1002/asmb.2507> and
#' Centofanti et al. (2021) <doi:10.1080/00401706.2020.1753581>.
#'
#'
#' @docType package
#' @name funcharts
#' @import fda
#' @import parallel
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom stats
#' approxfun as.formula hatvalues lm predict quantile setNames var rnorm dnorm approx
#' @importFrom rlang .data
#' @importFrom RSpectra eigs_sym
#' @importFrom matrixStats rowCumsums
#' @references
#' Capezza C, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
#' Control charts for
#' monitoring ship operating conditions and CO2 emissions based
#' on scalar-on-function regression.
#' \emph{Applied Stochastic Models in Business and Industry},
#' 36(3):477--500.
#' <doi:10.1002/asmb.2507>
#'
#' Centofanti F, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2021)
#' Functional Regression Control Chart.
#' \emph{Technometrics}, 63(3), 281--294. <doi:10.1080/00401706.2020.1753581>
NULL

# Quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
