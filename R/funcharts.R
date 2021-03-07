#' \code{funcharts} package
#'
#' Functional control charts.
#'
#'
#' @docType package
#' @name funcharts
#' @import fda
#' @import parallel
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom stats
#' approxfun as.formula hatvalues lm predict quantile setNames var
#' @importFrom rlang .data
#' @importFrom stats rnorm
NULL

# Quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
