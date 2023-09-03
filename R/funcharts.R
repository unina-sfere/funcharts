#' \code{funcharts} package
#'
#' Provides functional control charts
#' for statistical process monitoring of functional data,
#' using the methods of Capezza et al. (2020) <doi:10.1002/asmb.2507> and
#' Centofanti et al. (2021) <doi:10.1080/00401706.2020.1753581>.
#' The package is thoroughly illustrated in the paper of
#' Capezza et al (2023) <doi:10.1080/00224065.2023.2219012>
#'
#'
#' @docType package
#' @name funcharts
#' @importFrom fda
#' bifd create.bspline.basis eval.bifd eval.fd fd fdPar inprod.bspline is.basis
#' @importFrom fda
#' is.fd is.Lfd int2Lfd mean.fd sd.fd times.fd project.basis minus.fd is.fdPar
#' @importFrom fda
#' smooth.basis create.fourier.basis create.constant.basis center.fd
#' @importFrom fda
#' eval.penalty
#' @importFrom parallel
#' mclapply parLapply stopCluster makeCluster clusterExport stopCluster
#' @importFrom dplyr
#' select mutate filter arrange group_by bind_cols bind_rows case_when
#' @importFrom dplyr
#' contains everything inner_join rename "%>%" slice pull n ungroup
#' @importFrom ggplot2
#' aes aes_string after_stat xlab ylab ggplot ggtitle guides xlim ylim
#' @importFrom ggplot2
#' scale_color_continuous scale_colour_continuous scale_color_discrete labs
#' @importFrom ggplot2
#' scale_colour_discrete scale_fill_gradientn scale_color_gradientn
#' @importFrom ggplot2
#' scale_linetype_discrete
#' @importFrom ggplot2
#' scale_colour_gradientn scale_linetype_manual scale_colour_manual
#' @importFrom ggplot2
#' scale_color_manual scale_x_continuous element_rect element_text
#' @importFrom ggplot2
#' element_blank theme theme_bw geom_line geom_tile geom_contour geom_hline
#' @importFrom ggplot2
#' geom_vline geom_point geom_blank geom_text geom_segment geom_col
#' @importFrom ggplot2
#' scale_fill_manual facet_wrap
#' @importFrom stats
#' approxfun as.formula hatvalues lm predict quantile
#' @importFrom stats
#' setNames var rnorm dnorm approx
#' @importFrom stats
#' model.matrix rstandard formula uniroot
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom rlang .data :=
#' @importFrom RSpectra eigs_sym
#' @importFrom matrixStats rowCumsums
#' @importFrom graphics par persp
#' @importFrom stringr str_count
#' @references
#' Capezza C, Centofanti F, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2023)
#' funcharts: control charts for multivariate functional data in R.
#' \emph{Journal of Quality Technology}, \url{doi:10.1002/asmb.2507}
#'
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
