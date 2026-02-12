#' Robust Multivariate Functional Control Charts - Phase I (casewise version)
#'
#' It performs Phase I of the Robust Multivariate Functional Control Chart
#' (RoMFCC), proposed by Capezza et al. (2024), applied to casewise outlier
#' detection.
#'
#' @details
#' Unlike the original RoMFCC implementation, this version assumes that:
#' * functional filter
#' * robust multivariate functional imputation
#' have already been applied to the training and tuning datasets.
#' Therefore, the input data are expected to be multivariate functional
#' data free of cellwise outliers (casewise outliers may still be present).
#'
#' @param mfdobj_imp
#' A multivariate functional data object of class \code{mfd}, already imputed
#' and filtered, so with no cellwise outliers.
#' A robust multivariate principal component functional analysis
#' is applied to the imputed dataset for dimension reduction.
#' @param mfdobj_imp_tuning
#' An additional functional data object of class \code{mfd}, already imputed
#' and filtered, so with no cellwise outliers.
#' It is used to robustly estimate the distribution of the Hotelling's T2 and
#' SPE statistics in order to calculate control limits
#' to prevent overfitting issues that could reduce the
#' monitoring performance of the RoMFCC.
#' @param pca_par
#' A list with an argument \code{fev}, indicating a number between 0 and 1
#' denoting the fraction of variability that must be explained by the
#' principal components to be selected in the RoMFPCA step.
#' All the other arguments of this list are passed as arguments to the function
#' \code{rpca_mfd} in the RoMFPCA step.
#' All the arguments that are not passed take their default values.
#' See \code{\link{rpca_mfd}} for all the arguments and their default
#' values.
#' Default value is \code{list(fev = 0.7)}.
#' @param alpha_casewise
#' The overall nominal type-I error probability used to set
#' control chart limits and to identify functional casewise outliers
#' Default value is 0.0027.
#' @param verbose
#' If TRUE, it prints messages about the steps of the algorithm.
#' Default is FALSE.
#'
#' @return
#' A list of the following elements that are needed in Phase II:
#'
#' * \code{T2} the Hotelling's T2 statistic values for the Phase I data set,
#'
#' * \code{SPE} the SPE statistic values for the Phase I data set,
#'
#' * \code{T2_tun} the Hotelling's T2 statistic values for the tuning data set,
#'
#' * \code{SPE_tun} the SPE statistic values for the tuning data set,
#'
#' * \code{T2_lim} the Phase II control limit of
#' the Hotelling's T2 control chart,
#'
#' * \code{spe_lim} the Phase II control limit of
#' the SPE control chart,
#'
#' * \code{mod_pca} the final RoMFPCA model fitted on the Phase I data set,
#'
#' * \code{K} = K the number of selected principal components,
#'
#' * \code{T_T2_inv} if a tuning data set is provided,
#' it returns the inverse of the covariance matrix
#' of the first \code{K} scores, needed to calculate the Hotelling's T2
#' statistic for the Phase II observations.
#'
#' * \code{mean_scores_tuning_rob_mean} if a tuning data set is provided,
#' it returns the robust location estimate of the scores, needed to calculate
#' the Hotelling's T2 and SPE
#' statistics for the Phase II observations.
#'
#' @references
#' Capezza, C., Centofanti, F., Lepore, A., Palumbo, B. (2024)
#' Robust Multivariate Functional Control Chart.
#' \emph{Technometrics}, 66(4):531--547, <doi:10.1080/00401706.2024.2327346>.
#'
#' @examples
#' \dontrun{
#' library(funcharts)
#' set.seed(0)
#' dat <- simulate_data_RoMFCC(p_cellwise = 0.05,
#'                             p_casewise = 0.05,
#'                             outlier = "outlier_E",
#'                             M_outlier_cell = 0.03,
#'                             M_outlier_case = 0.01,
#'                             max_n_cellwise = 10)
#' mfdobj <- get_mfd_list(dat$X_mat_list, n_basis = 5)
#' mfdobj_training <- mfdobj[1:333, ]
#' mfdobj_tuning <- mfdobj[334:1000, ]
#' ff_training <- functional_filter(mfdobj = mfdobj_training)
#' ff_tuning <- functional_filter(mfdobj = mfdobj_tuning)
#' x_imp_training <- RoMFDI(mfdobj = ff_training$mfdobj)
#' x_imp_tuning <- RoMFDI(mfdobj = ff_tuning$mfdobj)
#' X_imp_training <- x_imp_training[[1]]
#' X_imp_tuning <- x_imp_tuning[[1]]
#' out_phase1_casewise <- RoMFCC_PhaseI_casewise(
#'   mfdobj_imp = X_imp_training,
#'   mfdobj_imp_tuning = X_imp_tuning
#' )
#' mfd_all_imputed <- rbind_mfd(X_imp_training, X_imp_tuning)
#' out_phase2_casewise <- RoMFCC_PhaseII_casewise(
#'   mfdobj_all_imp = mfdobj_all_imputed,
#'   mod_phaseI_casewise = out_phase1_casewise
#' )
#' plot_control_charts(out_phase2_casewise)
#' }
#'
RoMFCC_PhaseI_casewise <- function(mfdobj_imp,
                                   mfdobj_imp_tuning,
                                   pca_par = list(fev = 0.7),
                                   alpha_casewise = 0.0027,
                                   verbose = FALSE) {


  # Basic input controls
  if (anyNA(mfdobj_imp$coefs)) {
    stop("mfdobj_imp contains NA: RoMFCC_PhaseI_casewise assumes fully imputed data.")
  }

  if (anyNA(mfdobj_imp_tuning$coefs)) {
    stop("mfdobj_imp_tuning contains NA: RoMFCC_PhaseI_casewise assumes fully imputed data.")
  }


  # RoMFPCA default arguments
  if (!is.list(pca_par)) {
    stop("pca_par must be a list.")
  }

  if (is.null(pca_par$fev)) {
    pca_par$fev <- 0.7
  }

  nvar <- dim(mfdobj_imp$coefs)[3]
  nobs <- dim(mfdobj_imp$coefs)[2]


  # Robust PCA (training)
  if (verbose) {
    print("Training set: dimension reduction step...")
  }

  # Pack the imputed dataset into a list (for compatibility with the original code)
  X_imp <- list(mfdobj_imp)
  n_dataset <- length(X_imp)

  pca_default <- formals(rpca_mfd)
  pca_args <- pca_par
  pca_args$fev <- NULL
  pca_args <- c(
    pca_args,
    pca_default[!(names(pca_default) %in%
                    names(pca_par))])

  T2_lim_d <- spe_lim_d <- numeric()

  mod_pca_list <- corr_coef_list <- list()

  for (D in 1:n_dataset) {

    X_imp_d <- X_imp[[D]]
    nobs <- dim(X_imp_d$coefs)[2]
    nb <- X_imp_d$basis$nbasis
    nvars <- dim(X_imp_d$coefs)[3]
    nharm <- min(nobs - 1, nb * nvars)

    pca_args$mfdobj <- X_imp_d
    pca_args$nharm <- nharm
    mod_pca <- do.call(rpca_mfd, pca_args)

    cum_var <- cumsum(mod_pca$values / sum(mod_pca$values))
    mod_pca_list[[D]] <- mod_pca

    Gx <- mod_pca$harmonics$coefs
    Gx <- do.call(rbind, lapply(1:nvars, function(jj) Gx[, , jj]))
    corr_coef <- Gx %*% diag(mod_pca$values) %*% t(Gx)
    corr_coef_list[[D]] <- corr_coef

  }

  corr_coef_mean <- apply(simplify2array(corr_coef_list), 1:2, mean)

  # W <- as.matrix(Matrix::bdiag(rep(list(mod_pca$harmonics$basis$B), nvars)))
  W <- kronecker(diag(nvars), mod_pca$harmonics$basis$B)
  W_sqrt <- chol(W)
  W_sqrt_inv <- solve(W_sqrt)
  VCOV <- W_sqrt %*% corr_coef_mean %*% t(W_sqrt)
  e <- eigen(VCOV)
  loadings_coef <- W_sqrt_inv %*% e$vectors

  loadings_list <- list()
  for (jj in 1:nvars) {
    index <- 1:nb + (jj - 1) * nb
    loadings_list[[jj]] <- loadings_coef[index, 1:nharm]
  }

  loadings_coef_array <- simplify2array(loadings_list)
  loadings <- mfd(loadings_coef_array,
                  mod_pca$harmonics$basis,
                  B = mod_pca$harmonics$basis$B)

  values <- e$values[1:nharm]

  meanfd <- mfd(Reduce("+", lapply(mod_pca_list, function(x) {
    x$meanfd$coefs
  })) / n_dataset,
  mod_pca$harmonics$basis)
  scale_fd <- fda::fd(Reduce("+",
                             lapply(mod_pca_list, function(x) {
                               x$scale_fd$coefs
                             })) / n_dataset,
                      mod_pca$harmonics$basis)

  varprop <- values / sum(values)

  mod_pca_final <- list(harmonics = loadings,
                        values = values,
                        meanfd = meanfd,
                        scale_fd = scale_fd,
                        varprop = varprop)


  # T2 and SPE statistics (training)
  if (verbose) {
    print("Training set: calculating monitoring statistics...")
  }

  K <- which(cumsum(mod_pca_final$varprop) >= pca_par$fev)[1]

  alpha_sid <- 1 - (1 - alpha_casewise)^(1/2)

  T2_lim <- stats::qf(1 - alpha_sid, K, nobs - K) *
    (K * (nobs- 1 ) * (nobs + 1)) /
    ((nobs - K) * nobs)

  values_spe <- values[(K+1):length(values)]
  teta_1 <- sum(values_spe)
  teta_2 <- sum(values_spe^2)
  teta_3 <- sum(values_spe^3)
  h0 <- 1 - (2 * teta_1 * teta_3) / (3 * teta_2^2)
  c_alpha <- sign(h0) * abs(stats::qnorm(alpha_sid))

  spe_lim <- teta_1 * (((c_alpha * sqrt(2 * teta_2 * h0^2)) / teta_1) +
                         1 +
                         (teta_2 * h0 * (h0 - 1)) / teta_1^2)^(1 / h0)

  mod_T2_spe_all <- get_T2_spe(
    pca = mod_pca_final,
    components = 1:K,
    newdata_scaled = scale_mfd(mfdobj_imp,
                               center = mod_pca_final$meanfd,
                               scale = mod_pca_final$scale_fd)
  )

  T2_all <- mod_T2_spe_all$T2
  SPE_all <- mod_T2_spe_all$spe


  # T2 and SPE statistics (tuning)
  if (verbose) {
    print("Tuning set: calculating monitoring statistics and control limits...")
  }

  # Pack the imputed dataset into a list (for compatibility with the original code)
  X_imp_tuning <- list(mfdobj_imp_tuning)
  n_dataset_tun <- length(X_imp_tuning)

  mean_scores_tuning_rob <- list()
  S_scores_tuning_rob <- list()
  scores_tuning <- list()
  S_scores_tuning_rob2 <- list()

  for (D in seq_along(X_imp_tuning)) {
    X_imp_tuning_scaled <- scale_mfd(X_imp_tuning[[D]],
                                     center = mod_pca_final$meanfd,
                                     scale = mod_pca_final$scale_fd)
    scores_tuning[[D]] <- get_scores(
      mod_pca_final,
      components = 1:dim(mod_pca_final$harmonics$coefs)[2],
      newdata_scaled = X_imp_tuning_scaled
    )

    mod <- rrcov::CovMcd(scores_tuning[[D]][,1:K], alpha = 0.8)
    mod2 <- rrcov::CovMcd(scores_tuning[[D]][,-(1:K)], alpha = 0.8)

    mean_scores_tuning_rob[[D]] <- c(mod$center, mod2$center)
    S_scores_tuning_rob[[D]] <- mod$cov
    S_scores_tuning_rob2[[D]] <- mod2$cov

  }

  mean_scores_tuning_rob <- do.call(cbind, mean_scores_tuning_rob)
  mean_scores_tuning_rob_mean <- rowMeans(mean_scores_tuning_rob)
  S_scores_tuning_rob <- simplify2array(S_scores_tuning_rob)
  T_mean <- apply(S_scores_tuning_rob, 1:2, mean)

  T_T2 <- T_mean
  T_T2_inv <- solve(T_T2)
  S_scores_tuning_rob2 <- simplify2array(S_scores_tuning_rob2)
  T_mean <- apply(S_scores_tuning_rob2, 1:2, mean)
  T_SPE <- T_mean

  scores_tuning_cen <- t(t(scores_tuning[[D]][, 1:K]) -
                           mean_scores_tuning_rob_mean[1:K])
  T2_tun <- rowSums((scores_tuning_cen %*% T_T2_inv) * scores_tuning_cen)

  alpha_sid <- 1 - (1 - alpha_casewise)^(1/2)
  ntun <- length(T2_tun)

  T2_lim <- stats::qf(1 - alpha_sid, K, ntun - K) *
    (K * (ntun- 1 ) * (ntun + 1)) /
    ((ntun - K) * ntun)

  scores_tun_SPE <- t(t(scores_tuning[[D]][, -(1:K)]) -
                        mean_scores_tuning_rob_mean[-(1:K)])
  SPE_tun <- rowSums(scores_tun_SPE^2)

  e <- eigen(T_SPE)

  values_spe <- e$values
  teta_1 <- sum(values_spe)
  teta_2 <- sum(values_spe^2)
  teta_3 <- sum(values_spe^3)
  h0 <- 1 - (2 * teta_1 * teta_3) / (3 * teta_2^2)
  c_alpha <- sign(h0) * abs(stats::qnorm(alpha_sid))

  spe_lim <- teta_1 * (((c_alpha * sqrt(2 * teta_2 * h0^2)) / teta_1) +
                         1 +
                         (teta_2 * h0 * (h0 - 1)) / teta_1^2)^(1 / h0)

  out <- list(T2 = T2_all,
              SPE = SPE_all,
              T2_tun = T2_tun,
              SPE_tun = SPE_tun,
              T2_lim = T2_lim,
              spe_lim = spe_lim,
              mod_pca = mod_pca_final,
              K = K,
              T_T2_inv = T_T2_inv,
              mean_scores_tuning_rob_mean = mean_scores_tuning_rob_mean)

  return(out)

}





# RoMFCC functions
#' Robust Multivariate Functional Control Charts - Phase II (casewise version)
#'
#' It performs Phase II of the Robust Multivariate Functional Control Chart
#' (RoMFCC) for casewise outlier detection, computing Hotelling's T2 and SPE
#' monitoring statistics according to the methodology proposed by
#' Capezza et al. (2024).

#' @param mfdobj_all_imp
#' A multivariate functional data object of class \code{mfd}, containing the
#' concatenation of the fully imputed training and tuning sets to be
#' monitored for casewise outliers.
#' @param mod_phaseI_casewise
#' Output obtained by applying the function \code{RoMFCC_PhaseI_casewise}
#' to perform Phase I. See \code{\link{RoMFCC_PhaseI_casewise}}
#'
#' @return
#' A \code{data.frame} with as many rows as the number of
#' multivariate functional observations in the phase II data set and
#' the following columns:
#'
#' * one \code{id} column identifying the multivariate functional observation
#' in the phase II data set,
#'
#' * one \code{T2} column containing the Hotelling T2 statistic
#' calculated for all observations,
#'
#' * one column per each functional variable,
#' containing its contribution to the T2 statistic,
#'
#' * one \code{SPE} column containing the SPE statistic calculated
#' for all observations,
#'
#' * \code{T2_lim} gives the upper control limit of
#' the Hotelling's T2 control chart,
#'
#' * \code{SPE_lim} gives the upper control limit of the SPE control chart
#'
#' @references
#' Capezza, C., Centofanti, F., Lepore, A., Palumbo, B. (2024)
#' Robust Multivariate Functional Control Chart.
#' \emph{Technometrics}, 66(4):531--547, <doi:10.1080/00401706.2024.2327346>.
#'
#' @examples
#' \dontrun{
#' library(funcharts)
#' set.seed(0)
#' dat <- simulate_data_RoMFCC(p_cellwise = 0.05,
#'                             p_casewise = 0.05,
#'                             outlier = "outlier_E",
#'                             M_outlier_cell = 0.03,
#'                             M_outlier_case = 0.01,
#'                             max_n_cellwise = 10)
#' mfdobj <- get_mfd_list(dat$X_mat_list, n_basis = 5)
#' mfdobj_training <- mfdobj[1:333, ]
#' mfdobj_tuning <- mfdobj[334:1000, ]
#' ff_training <- functional_filter(mfdobj = mfdobj_training)
#' ff_tuning <- functional_filter(mfdobj = mfdobj_tuning)
#' x_imp_training <- RoMFDI(mfdobj = ff_training$mfdobj)
#' x_imp_tuning <- RoMFDI(mfdobj = ff_tuning$mfdobj)
#' X_imp_training <- x_imp_training[[1]]
#' X_imp_tuning <- x_imp_tuning[[1]]
#' out_phase1_casewise <- RoMFCC_PhaseI_casewise(
#'   mfdobj_imp = X_imp_training,
#'   mfdobj_imp_tuning,X_imp_tuning
#' )
#' mfd_all_imputed <- rbind_mfd(X_imp_training, X_imp_tuning)
#' out_phase2_casewise <- RoMFCC_PhaseII_casewise(
#'   mfdobj_all_imp = mfdobj_all_imputed,
#'   mod_phaseI_casewise = out_phase1_casewise
#' )
#' plot_control_charts(out_phase2_casewise)
#' }
#'
RoMFCC_PhaseII_casewise <- function(mfdobj_all_imp,
                                    mod_phaseI_casewise) {

  mod_pca <- mod_phaseI_casewise$mod_pca
  K <- mod_phaseI_casewise$K
  mfdobj_all_imp_std <- scale_mfd(mfdobj_all_imp,
                                  center = mod_pca$meanfd,
                                  scale = mod_pca$scale_fd)

  mod_T2_spe <- get_T2_spe(mod_pca, 1:K, mfdobj_all_imp_std)
  T2 <- mod_T2_spe$T2
  SPE <- mod_T2_spe$spe

  scores_new <- get_scores(mod_pca,
                           components = 1:dim(mod_pca$harmonics$coefs)[2],
                           newdata_scaled = mfdobj_all_imp_std)

  scores_new_cen <- t(t(scores_new[, 1:K]) -
                        mod_phaseI_casewise$mean_scores_tuning_rob_mean[1:K])
  T2 <- rowSums((scores_new_cen %*% mod_phaseI_casewise$T_T2_inv) * scores_new_cen)
  scores_new_cen_SPE <- t(t(scores_new[, -(1:K)]) -
                            mod_phaseI_casewise$mean_scores_tuning_rob_mean[-(1:K)])
  SPE <- rowSums(scores_new_cen_SPE^2)

  frac_out <- mean(T2 > mod_phaseI_casewise$T2_lim | SPE > mod_phaseI_casewise$spe_lim)
  ARL <- 1 / frac_out

  out <- data.frame(id = mfdobj_all_imp$fdnames[[2]],
                    T2 = T2,
                    SPE = SPE,
                    T2_lim = mod_phaseI_casewise$T2_lim,
                    SPE_lim = mod_phaseI_casewise$spe_lim)

  return(out)

}







# Robust Adaptive Multivariate Functional EWMA Control Chart - Phase I
#'
#' It performs Phase I of a robust version of the Adaptive Multivariate
#' Functional EWMA control chart (RoAMFEWMA).
#' The procedure combines:
#' \enumerate{
#'   \item Functional cellwise outlier detection via
#'   \code{\link{functional_filter}},
#'   \item Robust Multivariate Functional Data Imputation (RoMFDI) via
#'   \code{\link{RoMFDI}},
#'   \item Casewise outliers detection via RoMFCC
#'   (\code{\link{RoMFCC_PhaseI_casewise}} and
#'   \code{\link{RoMFCC_PhaseII_casewise}}),
#'   on the imputed Phase I data,
#'   \item AMFEWMA Phase I (soft-thresholding) calibration
#'   (\code{\link{AMFEWMA_PhaseI_ST}}) on cellwise and casewise clean data.
#' }
#' The resulting object can be directly used as \code{mod_1} argument
#' in \code{\link{RoAMFEWMA_PhaseII}}
#'
#' @details
#' Among the multiple imputed datasets, the first one is used to build
#' the cleaned training and tuning sets for AMFEWMA.
#'
#' @param mfdobj
#' A multivariate functional data object of class \code{mfd}.
#' This dataset is used as the Phase I **training set**.
#' A functional filter is first applied to detect functional cellwise
#' outliers; the flagged components are then imputed through a robust
#' multivariate functional imputation procedure.
#' The imputed data are subsequently processed by the Robust Multivariate
#' Functional Control Chart to detect and remove
#' functional casewise outliers.
#' The resulting clean dataset (free of both cellwise and casewise outliers)
#' is then used as the training set in the Phase I of the Adaptive
#' Multivariate Functional EWMA control chart.
#' @param  mfdobj_tuning
#' A multivariate functional data object of class \code{mfd}
#' representing the Phase I **tuning set**.
#' This dataset undergoes the same functional filtering and robust
#' imputation steps applied to \code{mfdobj}.
#' The filtered and imputed data are used within the Robust Multivariate
#' Functional Control Chart to estimate tuning quantities and control limits.
#' After casewise cleaning, the resulting clean tuning set is employed in the
#' Phase I calibration of the Adaptive Multivariate Functional EWMA control chart.
#' @param lambda
#' Lambda parameter to be used in the scoring function present in AMFEWMA.
#' See equation (7) or (8) of Capezza et al., 2024.
#' If supplied, it must be a number between 0 and 1. This is passed to
#' \code{\link{AMFEWMA_PhaseI_ST}} and used as a fixed input, skipping the internal
#' optimization step.
#' If NULL, it is automatically selected within \code{\link{AMFEWMA_PhaseI_ST}}
#' according to the optimization procedure presented in Section 2.4 of
#' Capezza et al. (2024).
#' The default value is NULL.
#' @param k
#' Parameter k to be used in the scoring function present in AMFEWMA.
#' See equation (7) or (8) of Capezza et al., 2024.
#' If supplied, it must be a number greater than zero. This is passed to
#' \code{\link{AMFEWMA_PhaseI_ST}} and used as a fixed input, skipping the internal
#' optimization step.
#' If NULL, it is automatically selected within \code{\link{AMFEWMA_PhaseI_ST}}
#' according to the optimization procedure in Section 2.4 of
#' Capezza et al. (2024).
#' The default value is \code{NULL}.
#' @param c
#' Non-negative soft-thresholding constant passed to
#' \code{\link{AMFEWMA_PhaseI_ST}}. If NULL (default), it is internally
#' set to 0, corresponding to no soft-thresholding. If a non-negative
#' value is provided, it is used in the AMFEWMA Phase I (soft-thresholding)
#' calibration step.
#' @param functional_filter_par
#' A list with an argument \code{filter} that can be TRUE or FALSE depending
#' on if the functional filter step must be performed or not.
#' All the other arguments of this list are passed as arguments to the function
#' \code{functional_filter} in the filtering step.
#' All the arguments that are not passed take their default values.
#' See \code{\link{functional_filter}} for all the arguments and their default
#' values.
#' Default is \code{list(filter = TRUE)}.
#' @param imputation_par
#' A list with an argument \code{method_imputation}
#' that can be \code{"RoMFDI"} or \code{"mean"} depending
#' on if the imputation step must be done by means of \code{\link{RoMFDI}} or
#' by just using the mean of each functional variable.
#' If \code{method_imputation = "RoMFDI"},
#' all the other arguments of this list are passed as arguments to the function
#' \code{RoMFDI} in the imputation step.
#' All the arguments that are not passed take their default values.
#' See \code{\link{RoMFDI}} for all the arguments and their default
#' values.
#' Default value is \code{list(method_imputation = "RoMFDI")}.
#' @param verbose
#' If TRUE, it prints messages about the steps of the algorithm.
#' Default is FALSE.
#'
#' @return
#' A list of the following elements:
#'
#' * \code{mod_1} object returned by \code{AMFEWMA_PhaseI_ST}, see the value
#' of \code{\link{AMFEWMA_PhaseI_ST}} for a full description of its component;
#'
#' * \code{mfd_clean_training} training data after complete cleaning,
#' containing no outliers at either cellwise or casewise;

#' * \code{mfd_clean_tuning} tuning data after complete cleaning,
#' containing no outliers at either cellwise or casewise;
#'
#' * \code{mfd_all_clean} full Phase I clean data (training + tuning);
#'
#' * \code{idx_casewise_outliers} indices of observations indetified as
#' casewise outliers by RoMFCC Phase II;
#'
#' * \code{ff_training} training set after the functional filter;
#'
#' * \code{ff_tuning} tuning set after the functional filter;
#'
#' * \code{X_imp_training_1} first imputation of the training set
#' after RoMFDI
#'
#' * \code{X_imp_tuning_1} first imputation of the tuning set
#' after RoMFDI
#'
#' * \code{X_all_imputed} training + tuning data after robust multivariate
#' functional imputation;
#'
#' * \code{mod_RoMFCC_phaseI_casewise} object returned by RoMFCC_PhaseI,
#' see the value of \code{\link{RoMFCC_PhaseI_casewise}} for a full description
#' of its component;
#'
#' * \code{mod_RoMFCC_phaseII_casewise} object returned by RoMFCC_PhaseII,
#' see the value of \code{\link{RoMFCC_PhaseII_casewise}} for a full
#' description of its component;
#'
#' @export
#'
#' @references
#' Capezza, C., Centofanti, F., Lepore, A., Palumbo, B. (2024)
#' Robust Multivariate Functional Control Chart.
#' \emph{Technometrics}, 66(4):531--547, <doi:10.1080/00401706.2024.2327346>.
#'
#' Capezza, C., Capizzi, G., Centofanti, F., Lepore, A., Palumbo, B. (2025)
#' An Adaptive Multivariate Functional EWMA Control Chart.
#' \emph{Journal of Quality Technology},  57(1):1--15,
#' doi:https://doi.org/10.1080/00224065.2024.2383674.
#'
#' @examples
#' \dontrun{
#' set.seed(0)
#' dat_phaseI <- simulate_data_RoMFCC(p_cellwise = 0.05,
#'                             p_casewise = 0.05,
#'                             outlier = "outlier_E",
#'                             M_outlier_cell = 0.03,
#'                             M_outlier_case = 0.01,
#'                             max_n_cellwise = 10)
#' dat_phaseII <- simulate_data_RoMFCC(OC = "OC_E",
#'                                     M_OC = 0.01,
#'                                     which_OC = 5)
#' mfdobj_phaseI <- get_mfd_list(dat_phaseI$X_mat_list, n_basis = 5)
#' mfdobj_phaseII <- get_mfd_list(dat_phaseII$X_mat_list, n_basis = 5)
#' mfdobj_training_phaseI <- mfdobj_phaseI[1:333, ]
#' mfdobj_tuning_phaseI <- mfdobj_phaseI[334:1000, ]
#' out_phaseI <- RoAMFEWMA_PhaseI(mfdobj = mfdobj_training_phaseI,
#'                                mfdobj_tuning = mfdobj_tuning_phaseI)
#' out_phaseII <- RoAMFEWMA_PhaseII(mfdobj_2 = mfdobj_phaseII,
#'                                  mod_1 = out_phaseI)
#' plot_control_charts(out_phaseII$cc)
#' }
RoAMFEWMA_PhaseI <- function(mfdobj,
                             mfdobj_tuning,
                             lambda = NULL,
                             k = NULL,
                             c = NULL,
                             functional_filter_par = list(filter = TRUE),
                             imputation_par = list(method_imputation = "RoMFDI"),
                             verbose = FALSE) {

  # -------------------------
  # 0. Basic input controls
  # -------------------------

  # Filter default arguments

  if (!is.list(functional_filter_par)) {
    stop("functional_filter_par must be a list.")
  }

  if (is.null(functional_filter_par$filter)) {
    functional_filter_par$filter <- TRUE
  }

  # Imputation default arguments

  if (!is.list(imputation_par)) {
    stop("imputation_par must be a list.")
  }

  if (is.null(imputation_par$method_imputation)) {
    imputation_par$method_imputation <- "RoMFDI"
  }

  # Lambda e k dafult arguments

  if (!is.null(lambda)) {
    if (!is.numeric(lambda) || length(lambda) != 1 || is.na(lambda) || lambda <= 0 || lambda > 1) {
      stop("lambda must be a single numeric value in (0, 1].")
    }
  }

  if (!is.null(k)) {
    if (!is.numeric(k) || length(k) != 1 || is.na(k) || k <= 0) {
      stop("k must be a single positive numeric value.")
    }
  }

  # c default arguments
  if (!is.null(c)) {
    if (!is.numeric(c) || length(c) != 1 || is.na(c) || c < 0) {
      stop("Parameter 'c' must be a single non-negative numeric value.")
    }
  }

  nvar <- dim(mfdobj$coefs)[3]
  nobs <- dim(mfdobj$coefs)[2]


  # ----------------------------
  # Step 1: CELLWISE OUTLIERS
  # Functional Filter + RoMFDI
  # ----------------------------


  # --- 1A. Functional Filter on training ---

  if (functional_filter_par$filter) {
    if (verbose) {
      print("Training set: functional filtering...")
    }

    functional_filter_default <- formals(functional_filter)
    functional_filter_args <- functional_filter_par
    functional_filter_args$filter <- NULL
    functional_filter_args <- c(
      functional_filter_args,
      functional_filter_default[!(names(functional_filter_default) %in%
                                    names(functional_filter_par))])
    functional_filter_args$mfdobj <- mfdobj

    ff_training_out <- do.call(functional_filter, functional_filter_args)

    idx_out_ff_training <- ff_training_out$flagged_outliers
    mfd_ff_training <- ff_training_out$mfdobj
  }

  else {
    idx_out_ff_training <- NULL
    mfd_ff_training <- mfdobj
    ff_training_out <- NULL
  }


  # --- 1B. Imputation training ---

  if(anyNA(mfd_ff_training$coefs)) {
    if (verbose) {
      print("Training set: imputation step...")
    }

    if (!(imputation_par$method_imputation %in% c("RoMFDI", "mean"))) {
      stop("method_imputation must be either \"RoMFDI\" or \"mean\"")
    }

    if (imputation_par$method_imputation == "RoMFDI") {
      imputation_default <- formals(RoMFDI)
      imputation_args <- imputation_par
      imputation_args$method_imputation <- NULL
      imputation_args <- c(
        imputation_args,
        imputation_default[!(names(imputation_default) %in%
                               names(imputation_par))])
      imputation_args$mfdobj <- mfd_ff_training
      RoMFDI_training <- do.call(RoMFDI, imputation_args)
    }

    if (imputation_par$method_imputation == "mean") {
      RoMFDI_training <- mfd_ff_training
      for(j in 1:nvar) {
        which_na <- sort(unique(which(is.na(RoMFDI_training$coefs[, , j]),
                                      arr.ind = TRUE)[, 2]))
        which_ok <- setdiff(1:dim(RoMFDI_training$coefs)[2], which_na)
        X_fd <- fda::fd(RoMFDI_training$coefs[, which_ok, j], RoMFDI_training$basis)
        X_fdata <- fda.usc::fdata(X_fd)
        mu_fdata <- rofanova::fusem(X_fdata)$mu
        mu_fd <- fda.usc::fdata2fd(mu_fdata, nbasis = RoMFDI_training$basis$nbasis)
        RoMFDI_training$coefs[, which_na, j] <- mu_fd$coefs
      }
      RoMFDI_training <- list(RoMFDI_training)
    }
  }

  else {
    RoMFDI_training <- list(mfd_ff_training)
  }

  n_dataset <- length(RoMFDI_training)


  # --- 1C. Functional Filter on tuning ---

  if (functional_filter_par$filter) {
    if (verbose) {
      print("Tuning set: functional filtering...")
    }

    functional_filter_args$mfdobj <- mfdobj_tuning

    ff_tuning_out <- do.call(functional_filter, functional_filter_args)

    idx_out_ff_tuning <- ff_tuning_out$flagged_outliers
    mfd_ff_tuning <- ff_tuning_out$mfdobj
  }

  else {
    idx_out_ff_tuning <- NULL
    mfd_ff_tuning <- mfdobj_tuning
    ff_tuning_out <- NULL
  }


  # --- 1D. Imputation tuning ---

  nobs_training <- dim(RoMFDI_training[[1]]$coefs)[2]

  if(anyNA(mfd_ff_tuning$coefs)) {
    if (verbose) {
      print("Tuning set: imputation step...")
    }

    if (imputation_par$method_imputation == "RoMFDI") {
      imputation_args$mfdobj <- rbind_mfd(RoMFDI_training[[1]], mfd_ff_tuning)
      RoMFDI_tuning <- do.call(RoMFDI, imputation_args)
      for (D in 1:length(RoMFDI_tuning)) {
        RoMFDI_tuning[[D]] <- RoMFDI_tuning[[D]][-(1:nobs_training)]
      }
    }

    if (imputation_par$method_imputation == "mean") {
      RoMFDI_tuning <- mfd_ff_tuning
      for (j in 1:nvar) {
        which_na <- sort(unique(which(is.na(RoMFDI_tuning$coefs[, , j]),
                                      arr.ind = TRUE)[, 2]))
        which_ok <- setdiff(1:dim(RoMFDI_tuning$coefs)[2], which_na)
        X_fd <- fda::fd(RoMFDI_tuning$coefs[, which_ok, j], RoMFDI_tuning$basis)
        X_fdata <- fda.usc::fdata(X_fd)
        mu_fdata <- rofanova::fusem(X_fdata)$mu
        mu_fd <- fda.usc::fdata2fd(mu_fdata, nbasis = RoMFDI_tuning$basis$nbasis)
        RoMFDI_tuning$coefs[, which_na, j] <- mu_fd$coefs
      }
      RoMFDI_tuning <- list(RoMFDI_tuning)
      n_dataset <- 1
    }
  }

  else {
    RoMFDI_tuning <- list(mfd_ff_tuning)
  }

  n_dataset <- length(RoMFDI_tuning)


  # --- 1E. Extract first imputation ---

  X_mfdimp_training_1 <- RoMFDI_training[[1]]
  X_mfdimp_tuning_1 <- RoMFDI_tuning[[1]]


  # ----------------------------------------
  # STEP 2 - CASEWISE OUTLIERS via RoMFCC
  # ----------------------------------------


  # --- 2A. Phase I of RoMFCC ---

  if (verbose) {
    message("Phase I: casewise detection with RoMFCC...")
  }

  RoMFCC_PhaseI_casewise_default <- formals(RoMFCC_PhaseI_casewise)

  RoMFCC_PhaseI_casewise_args <- list(
    mfdobj_imp = X_mfdimp_training_1,
    mfdobj_imp_tuning = X_mfdimp_tuning_1
  )

  RoMFCC_PhaseI_casewise_args <- c(
    RoMFCC_PhaseI_casewise_args,
    RoMFCC_PhaseI_casewise_default[!(names(RoMFCC_PhaseI_casewise_default) %in%
                                       names(RoMFCC_PhaseI_casewise_args))])

  RoMFCC_PhaseI_casewise_out <- do.call(RoMFCC_PhaseI_casewise,
                                        RoMFCC_PhaseI_casewise_args)


  # --- 2B. Phase II of RoMFCC ---

  if (verbose) {
    message("Phase II: casewise detection with RoMFCC...")
  }

  # Union of the first imputation of training and tuning
  mfd_all_imputed <- rbind_mfd(X_mfdimp_training_1, X_mfdimp_tuning_1)

  RoMFCC_PhaseII_casewise_args <- list(
    mfdobj_all_imp = mfd_all_imputed,
    mod_phaseI_casewise = RoMFCC_PhaseI_casewise_out
  )

  RoMFCC_PhaseII_casewise_out <- do.call(RoMFCC_PhaseII_casewise,
                                         RoMFCC_PhaseII_casewise_args)

  casewise_outliers <- which(
    RoMFCC_PhaseII_casewise_out$T2 > RoMFCC_PhaseII_casewise_out$T2_lim |
      RoMFCC_PhaseII_casewise_out$SPE > RoMFCC_PhaseII_casewise_out$SPE_lim
  )

  if(length(casewise_outliers) > 0) {

    if(verbose) {
      message("Removing", length(casewise_outliers),
              "functional casewise outliers ....")
    }

    mfd_all_clean <- mfd_all_imputed[-casewise_outliers, ]
  }

  else {

    if(verbose) {
      message("No functional casewise outliers detected.")
    }

    mfd_all_clean <- mfd_all_imputed
  }


  # -----------------------------------------
  # STEP 3 - AMFEWMA Phase I on cleaned data
  # -----------------------------------------

  nobs_clean <- dim(mfd_all_clean$coefs)[2]

  # Splitting clean data into training and tuning

  mfd_clean_training <- mfd_all_clean[1:(nobs_clean/2), ]
  mfd_clean_tuning <- mfd_all_clean[-(1:(nobs_clean/2)), ]

  if(verbose) {
    message("Fitting AMFEWMA Phase I on cleaned data")
  }

  AMFEWMA_default <- formals(AMFEWMA_PhaseI_ST)

  AMFEWMA_args <- list(
    mfdobj = mfd_clean_training,
    mfdobj_tuning = mfd_clean_tuning,
    lambda = lambda,
    k = k,
    c = c
  )

  AMFEWMA_args <- c(
    AMFEWMA_args,
    AMFEWMA_default[!(names(AMFEWMA_default) %in% names(AMFEWMA_args))]
  )

  AMFEWMA_PhaseI_ST_out <- do.call(AMFEWMA_PhaseI_ST, AMFEWMA_args)

  # --------
  # OUTPUT
  # --------

  out <- list(
    mod_1 = AMFEWMA_PhaseI_ST_out,
    mfd_clean_training = mfd_clean_training,
    mfd_clean_tuning = mfd_clean_tuning,
    mfd_all_clean = mfd_all_clean,
    idx_casewise_outliers = casewise_outliers,
    ff_training = ff_training_out,
    ff_tuning = ff_tuning_out,
    X_imp_training_1 = X_mfdimp_training_1,
    X_imp_tuning_1 = X_mfdimp_tuning_1,
    X_all_imputed = mfd_all_imputed,
    mod_RoMFCC_phaseI_casewise = RoMFCC_PhaseI_casewise_out,
    mod_RoMFCC_phaseII_casewise = RoMFCC_PhaseII_casewise_out
  )

  return(out)

}





# Robust Adaptive Multivariate Functional EWMA Control Chart - Phase II
#'
#' This function performs Phase II of the robust Adaptive Multivariate
#' Functional EWMA (RoAMFEWMA) control chart.
#'
#' @details
#' This function is conceptually similar to \code{AMFEWMA_PhaseII}, proposed
#' by Capezza et al. (2024), but adapted to the RoAMFEWMA framework.
#' In Phase II, monitoring relies on the RoAMFEWMA model calibrated in Phase I
#' on data cleaned from both cellwise and casewise outliers.
#' The monitoring statistic, control limit, and bootstrap-based ARL estimation
#' remain unchanged, but the input model must be the robust one obtained
#' through \code{RoAMFEWMA_PhaseI}.
#'
#' @param mfdobj_2
#' An object of class \code{mfd} containing the Phase II multivariate
#' functional data set, to be monitored with the RoAMFEWMA control chart.
#' @param mod_1
#' The output of the Phase I achieved through the
#' \code{\link{RoAMFEWMA_PhaseI}} function.
#' @param n_seq_2
#' If it is 1, the Phase II monitoring statistic is calculated on
#' the data sequence.
#' If it is an integer number larger than 1, a number \code{n_seq_2} of
#' bootstrap sequences are sampled with replacement from \code{mfdobj_2}
#' to allow uncertainty quantification on the estimation of the run length.
#' Default value is 1.
#' @param l_seq_2
#' If \code{n_seq_2} is larger than 1, this parameter sets the
#' length of each bootstrap sequence to be generated.
#' Default value is 2000.
#'
#' @return
#' A list with the following elements.
#'
#' * \code{ARL_2}: the average run length estimated over the
#' bootstrap sequences. If \code{n_seq_2} is 1, it is simply the run length
#' observed over the Phase II sequence, i.e., the number of observations
#' up to the first alarm,
#'
#' * \code{RL}: the run length
#' observed over the Phase II sequence, i.e., the number of observations
#' up to the first alarm,
#'
#' * \code{V2}: a list with length \code{n_seq_2}, containing the
#' AMFEWMA monitoring statistic in Equation (8) of Capezza
#' et al. (2024), calculated in each bootstrap sequence, until the first alarm.
#'
#' * \code{cc}: a data frame with the information needed to plot the
#' AMFEWMA control chart in Phase II, with the following columns.
#' \code{id} contains the id of each multivariate functional observation,
#' \code{amfewma_monitoring_statistic} contains the AMFEWMA monitoring
#' statistic values calculated on the Phase II sequence,
#' \code{amfewma_monitoring_statistic_lim} is the upper control limit.
#'
#' @export
#'
#' @references
#' Capezza, C., Capizzi, G., Centofanti, F., Lepore, A., Palumbo, B. (2025)
#' An Adaptive Multivariate Functional EWMA Control Chart.
#' \emph{Journal of Quality Technology},  57(1):1--15,
#' doi:https://doi.org/10.1080/00224065.2024.2383674.
#'
#' @examples
#' \dontrun{
#' set.seed(0)
#' dat_phaseI <- simulate_data_RoMFCC(p_cellwise = 0.05,
#'                             p_casewise = 0.05,
#'                             outlier = "outlier_E",
#'                             M_outlier_cell = 0.03,
#'                             M_outlier_case = 0.01,
#'                             max_n_cellwise = 10)
#' dat_phaseII <- simulate_data_RoMFCC(OC = "OC_E",
#'                                     M_OC = 0.01,
#'                                     which_OC = 5)
#' mfdobj_phaseI <- get_mfd_list(dat_phaseI$X_mat_list, n_basis = 5)
#' mfdobj_phaseII <- get_mfd_list(dat_phaseII$X_mat_list, n_basis = 5)
#' mfdobj_training_phaseI <- mfdobj_phaseI[1:333, ]
#' mfdobj_tuning_phaseI <- mfdobj_phaseI[334:1000, ]
#' out_phaseI <- RoAMFEWMA_PhaseI(mfdobj = mfdobj_training_phaseI,
#'                                mfdobj_tuning = mfdobj_tuning_phaseI)
#' out_phaseII <- RoAMFEWMA_PhaseII(mfdobj_2 = mfdobj_phaseII,
#'                                  mod_1 = out_phaseI)
#' plot_control_charts(out_phaseII$cc)
#' }
RoAMFEWMA_PhaseII <- function(mfdobj_2,
                              mod_1,
                              n_seq_2 = 1,
                              l_seq_2 = 2000) {

  mod_1 <- mod_1$mod_1
  nobs_2 <- dim(mfdobj_2$coefs)[2]
  nvar <- dim(mfdobj_2$coefs)[3]
  grid_points <- mod_1$grid_points
  mean_mfdobj <- mod_1$mean_mfdobj
  vectors <- mod_1$vectors
  values <- mod_1$values
  lambda <- mod_1$lambda
  k <- mod_1$k
  h <- mod_1$h
  huber <- mod_1$huber
  c <- mod_1$c

  RL <- numeric(n_seq_2)

  mfdobj_2_cen <- funcharts::scale_mfd(mfdobj_2,
                                       center = mean_mfdobj,
                                       scale = FALSE)

  mfdobj_2_cen_eval <- fda::eval.fd(grid_points, mfdobj_2_cen)

  X2 <- matrix(aperm(mfdobj_2_cen_eval, c(2, 1, 3)), nrow = nobs_2)

  V2 <- list()

  for (jj in 1:n_seq_2) {
    if (n_seq_2 == 1) {
      idx2 <- 1:nobs_2
    } else {
      idx2 <- sample(1:nobs_2, l_seq_2, replace = TRUE)
    }
    output <- get_RL_cpp(
      X2 = X2,
      X_IC = matrix(),
      idx2 = idx2,
      idx_IC = numeric(),
      lambda = lambda,
      k = k,
      huber = huber,
      c = c,
      h = h,
      Values = values,
      Vectors = vectors
    )
    V2[[jj]] <- output$T2
    RL[jj] <- output$RL
  }

  ARL_2 <- mean(RL, na.rm = TRUE)

  if (mean(is.na(RL)) == 1) {
    ARL_2 <- 1e10
  }

  idx2 <- 1:nobs_2

  output <- get_RL_cpp(
    X2 = X2,
    X_IC = matrix(),
    idx2 = idx2,
    idx_IC = numeric(),
    lambda = lambda,
    k = k,
    huber = huber,
    c = c,
    h = 1e10,
    Values = values,
    Vectors = vectors
  )

  cc <- data.frame(
    id = mfdobj_2$fdnames[[2]],
    amfewma_monitoring_statistic = output$T2[, 1],
    amfewma_monitoring_statistic_lim = mod_1$h
  )

  return(list(
    ARL_2 = ARL_2,
    RL = RL,
    V2 = V2,
    cc = cc
  ))

}





# AMFEWMA function
#' Adaptive Multivariate Functional EWMA control chart - Phase I
#' (soft-thresholding version)
#'
#' This function performs Phase I of a soft-thresholded version of the
#' Adaptive Multivariate Functional EWMA (AMFEWMA) control chart proposed
#' by Capezza et al. (2024). In the RoAMFEWMA framework, it is applied to
#' Phase I data that have already been cleaned from cellwise and casewise
#' outliers.
#' Differently from the original AMFEWMA procedure, after the
#' multivariate functional PCA step (used for representation), the contribution
#' of each retained principal component to the monitoring statistic is filtered
#' via soft-thresholding. In particular, standardized squared PCA scores
#' contribute only through their excess over a threshold c.
#'
#' @param mfdobj
#' An object of class \code{mfd} containing the Phase I multivariate
#' functional training data set. In the robust RoAMFEWMA procedure,
#' this corresponds to the cleaned Phase I training set, obtained after
#' functional filtering, robust imputation, and casewise outlier removal.
#' @param mfdobj_tuning
#' An object of class \code{mfd} containing the Phase I multivariate
#' functional tuning data set, used to estimate the AMFEWMA control chart
#' limit. In the RoAMFEWMA procedure, this corresponds to the cleaned
#' Phase I tuning set.
#' @param lambda
#' lambda parameter to be used in the score function.
#' See Equation (7) or (8) of Capezza et al. (2024).
#' If it is provided, it must be a number between zero and one.
#' If NULL, it is chosen through the selected according to the
#' optimization procedure presented in Section 2.4 of Capezza et al. (2024).
#' In this case, it is chosen among the values of
#' \code{optimization_pars$lambda_grid}.
#' Default value is NULL.
#' @param k
#' k parameter to be used in the score function.
#' See Equation (7) or (8) of Capezza et al. (2024).
#' If it is provided, it must be a number greater than zero.
#' If NULL, it is chosen through the selected according to the
#' optimization procedure presented in Section 2.4 of Capezza et al. (2024).
#' In this case, it is chosen among the values of
#' \code{optimization_pars$k_grid}.
#' Default value is NULL.
#' @param ARL0
#' The nominal in-control average run length.
#' Default value is 200.
#' @param bootstrap_pars
#' Parameters of the bootstrap procedure described in
#' Section 2.4 of Capezza et al. (2024) for the estimation of the
#' control chart limit.
#' It must be a list with two arguments.
#' \code{n_seq} is the number of bootstrap sequences to be generated.
#' \code{l_seq} is the length of each bootstrap sequence, i.e., the number
#' of observations to be sampled with replacement from the tuning set.
#' Default value is \code{list(n_seq = 200, l_seq = 2000)}.
#' @param optimization_pars
#' Parameters to be used in the optimization procedure described in Section
#' 2.4 of Capezza et al. (2024) for the selection of the parameters
#' lambda and k.
#' It must be a list of the following parameters.
#' \code{lambda_grid} contains the possible values of
#' the parameter \code{lambda}.
#' \code{k_grid} contains the possible values of the parameter \code{k}.
#' \code{epsilon} is the parameter used in Equation (10) of
#' Capezza et al. (2024).
#' When performing the parameter optimization,
#' first the parameters lambda and k are selected to minimize the ARL
#' with respect to a large shift,
#' then the same parameters are chosen to minimize the ARL
#' with respect to a small shift, given that the resulting ARL with respect
#' to the previous large shift does not increase, in percentage,
#' more than \code{epsilon}*100.
#' Default value is 0.1.
#' \code{sd_small} is a positive constant that multiplies the
#' standard deviation function to define the small shift delta_1 in Section
#' 2.4 of Capezza et al. (2024). In fact, the small shift is defined as
#' delta_1(t) = mu_0(t) + \code{sd_small} * sigma(t), where
#' mu_0(t) is the estimated in-control mean function and
#' sigma(t) is the estimated standard deviation function.
#' Default value is 0.25.
#' \code{sd_big} is a positive constant that multiplies the
#' standard deviation function to define the large shift delta_2 in Section
#' 2.4 of Capezza et al. (2024). In fact, the large shift is defined as
#' delta_2(t) = mu_0(t) + \code{sd_large} * sigma(t), where
#' mu_0(t) is the estimated in-control mean function and
#' sigma(t) is the estimated standard deviation function.
#' Default value is 2.
#' @param discrete_grid_length
#' The number of equally spaced argument values at which the \code{mfd}
#' objects are discretized.
#' Default value is 25.
#' @param score_function
#' Score function to be used in Equation (7) or (8) of Capezza et al. (2024),
#' to calculate the weighting parameter of the EWMA statistic
#' for each observation of the sequence.
#' Two values are possible.
#' If "huber", it uses the score function (7)
#' inspired by the Huber's function.
#' If "tukey", it uses the score function (8)
#' inspired by the Tukey's bisquare function.
#' @param fev
#' Number between 0 and 1 denoting the fraction
#' of variability that must be explained by the
#' principal components to be selected after
#' applying multivariate functional principal component analysis
#' on \code{mfdobj}. Default is 0.99.
#' @param c
#' Non-negative soft-thresholding constant used to shrink the contribution
#' of PCA scores in the monitoring statistic.
#' After standardization, components whose squared scores are below c
#' are effectively ignored, while components exceeding c contribute
#' only by the amount that surpasses the threshold.
#' Smaller values of c allow more components (and potentially more noise) to
#' enter the statistic; larger values filter out more components but may
#' remove signal.
#' Default value is NULL.
#' @param n_skip
#' The upper control limit of the AMFEWMA control chart is set
#' to achieve a desired in-control ARL, evaluated after the
#' monitoring statisic has reached steady state.
#' A monitoring statistic is in a steady state
#' if the process has been in control long enough
#' for the effect of the starting value to become negligible
#' (Lucas and Saccucci, 1990).
#' In this regard, the first \code{n_skip} observations
#' are excluded from the calculation of the run length.
#' Default value is 100.
#'
#' @return
#' A list with the following elements.
#' \code{lambda} is the selected lambda parameter.
#' \code{k} is the selected k parameter.
#' \code{mod_1} contains the estimated Phase I model. It is a list with
#' the following elements.
#'
#' * \code{mfdobj} the \code{mfdobj} object passed as input to this function,
#'
#' * \code{mfdobj_tuning} the \code{mfdobj_tuning} object
#' passed as input to this function,
#'
#' * \code{inv_sigmaY_reg}: the matrix containing the discretized
#' version of the function K^*(s,t) defined in Equation (9) of
#' Capezza et al. (2024),
#'
#' * \code{mean_mfdobj}: the estimated mean function,
#'
#' * \code{h}: the calculated upper control limit of the AMFEWMA control chart,
#'
#' * \code{ARL0}: the estimated in-control ARL, which should be close to the
#' nominal value passed as input to this function,
#'
#' * \code{lambda}: the lambda parameter selected by the optimization
#' procedure described in Section 2.4 of Capezza et al. (2024).
#'
#' * \code{k}: The function C_j(t)=k sigma_j(t) appearing in the score
#' functions (7) and (8) of Capezza et al. (2024).
#'
#' * \code{grid_points}: the grid containing the points over which
#' the functional data are discretized before computing the AMFEWMA monitoring
#' statistic and estimating all the model parameters.
#'
#' * \code{V2_mat}: the \code{n_seq}X\code{l_seq} matrix containing,
#' in each column, the AMFEWMA monitoring statistic values of each
#' bootstrap sequence.
#' This matrix is used to set the control chart limit \code{h} to
#' ensure that the desired average run length is achieved.
#'
#' * \code{n_skip}: the \code{n_skip} input parameter passed to this function,
#'
#' * \code{huber}: if the input parameter \code{score_function} is
#' \code{"huber"}, this is TRUE, else is FALSE,
#'
#' * \code{vectors}: the discretized eigenfunctions psi_l(t) of
#' the covariance function, appearing in Equation (9) of Capezza et al. (2024).
#'
#' * \code{values}: the eigenvalues rho_l of
#' the covariance function, appearing in Equation (9) of Capezza et al. (2024).
#'
#' * \code{c}: the soft-thresholding constant passed as input to this function.
#'
#' @references
#' Capezza, C., Capizzi, G., Centofanti, F., Lepore, A., Palumbo, B. (2025)
#' An Adaptive Multivariate Functional EWMA Control Chart.
#' \emph{Journal of Quality Technology},  57(1):1--15,
#' doi:https://doi.org/10.1080/00224065.2024.2383674.
#'
#' Lucas, J. M., Saccucci, M. S. (1990)
#' Exponentially weighted moving average control schemes:
#' properties and enhancements. \emph{Technometrics}, 32(1), 1-12.
#'
#' @examples
#' \dontrun{
#' set.seed(0)
#' library(funcharts)
#' dat_I <- simulate_mfd(nobs = 1000,
#'                       correlation_type_x = c("Bessel", "Bessel", "Bessel"),
#'                       sd_x = c(0.3, 0.3, 0.3))
#' dat_tun <- simulate_mfd(nobs = 1000,
#'                         correlation_type_x = c("Bessel", "Bessel", "Bessel"),
#'                         sd_x = c(0.3, 0.3, 0.3))
#' dat_II <- simulate_mfd(nobs = 200,
#'                        correlation_type_x = c("Bessel", "Bessel", "Bessel"),
#'                        shift_type_x = c("C", "C", "C"),
#'                        d_x = c(2, 2, 2),
#'                        sd_x = c(0.3, 0.3, 0.3))
#' mfdobj_I   <- get_mfd_list(dat_I$X_list)
#' mfdobj_tun <- get_mfd_list(dat_tun$X_list)
#' mfdobj_II  <- get_mfd_list(dat_II$X_list)
#' # In the RoAMFEWMA framework, mfdobj_I and mfdobj_tun would correspond
#' # to the cleaned Phase I training and tuning sets (after functional
#' # filtering, robust imputation and RoMFCC-based casewise cleaning).
#' # Here we use simulated in-control data as a proxy for such clean data.
#' p <- plot_mfd(mfdobj_I[1:100])
#' lines_mfd(p, mfdobj_II, col = "red")
#' mod <- AMFEWMA_PhaseI_ST(mfdobj = mfdobj_I,
#'                       mfdobj_tuning = mfdobj_tun)
#' print(mod$lambda)
#' print(mod$k)
#' cc <- AMFEWMA_PhaseII_ST(mfdobj_2 = rbind_mfd(mfdobj_I[1:100], mfdobj_II),
#'                       mod_1 = mod)
#' plot_control_charts(cc$cc, nobsI = 100)
#' }
#'
AMFEWMA_PhaseI_ST <- function(mfdobj,
                              mfdobj_tuning,
                              lambda = NULL,
                              k = NULL,
                              ARL0 = 200,
                              bootstrap_pars = list(n_seq = 200,
                                                    l_seq = 2000),
                              optimization_pars = list(
                                lambda_grid = c(0.1, 0.2, 0.3, 0.5, 1),
                                k_grid = c(1, 2, 3, 4),
                                epsilon = 0.1,
                                sd_small = 0.25,
                                sd_big = 2
                              ),
                              discrete_grid_length = 25,
                              score_function = "huber",
                              fev = 0.99,
                              c = NULL,
                              n_skip = 100) {

  if (!is.null(lambda)) {
    optimization_pars$lambda_grid <- lambda
  }
  if (!is.null(k)) {
    optimization_pars$k_grid <- k
  }

  if (!is.list(optimization_pars)) {
    stop("optimization_pars must be a list.")
  }

  if (is.null(optimization_pars$lambda_grid)) {
    optimization_pars$lambda_grid <- c(0.1, 0.2, 0.3, 0.5, 1)
  }

  if (is.null(optimization_pars$k_grid)) {
    optimization_pars$k_grid <- c(1, 2, 3, 4)
  }

  if (is.null(optimization_pars$epsilon)) {
    optimization_pars$epsilon <- 0.1
  }

  if (is.null(optimization_pars$sd_small)) {
    optimization_pars$sd_small <- 0.25
  }

  if (is.null(optimization_pars$sd_big)) {
    optimization_pars$sd_big <- 2
  }

  if (!is.list(bootstrap_pars)) {
    stop("bootstrap_pars must be a list.")
  }

  if (is.null(bootstrap_pars$n_seq)) {
    bootstrap_pars$n_seq <- 200
  }

  if (is.null(bootstrap_pars$l_seq)) {
    bootstrap_pars$l_seq <- 2000
  }

  lambda_grid <- optimization_pars$lambda_grid
  k_grid <- optimization_pars$k_grid
  epsilon <- optimization_pars$epsilon
  sd_small <- optimization_pars$sd_small
  sd_big <- optimization_pars$sd_big
  n_seq <- bootstrap_pars$n_seq
  l_seq <- bootstrap_pars$l_seq


  # Algorithm to select lambda and k
  sd_mfdobj <- fda::sd.fd(mfdobj)

  p <- dim(mfdobj$coefs)[3]
  n_basis <- dim(mfdobj$coefs)[1]
  nobs_tun <- dim(mfdobj_tuning$coefs)[2]

  # Basis coefficients of out-of-control observations
  # with small shift (0.5*sd_mfdobj)

  coef_small <- array(0, dim = c(n_basis, nobs_tun, p))

  # with large shift (2*sd_mfdobj)

  coef_big <- array(0, dim = c(n_basis, nobs_tun, p))

  sd_mfdobj_coefs <- array(sd_mfdobj$coefs, dim = c(n_basis, 1, p))

  for (hh in 1:nobs_tun) {
    mfdobj_tuning_coefs_hh <- mfdobj_tuning$coefs[, hh, , drop = FALSE]
    coef_small[, hh,] <- sd_mfdobj_coefs * sd_small + mfdobj_tuning_coefs_hh
    coef_big[, hh,] <- sd_mfdobj_coefs * sd_big + mfdobj_tuning_coefs_hh
  }

  fase2_small <- mfd(coef_small, mfdobj_tuning$basis, mfdobj_tuning$fdnames)
  fase2_big <- mfd(coef_big, mfdobj_tuning$basis, mfdobj_tuning$fdnames)

  # Calculate ARL for each combination of lambda and k

  n_lambda <- length(lambda_grid)
  n_k <- length(k_grid)
  matrix_ARL_big <- matrix(0, n_k, n_lambda)
  matrix_ARL_small <- matrix(0, n_k, n_lambda)
  rownames(matrix_ARL_big) <- k_grid
  colnames(matrix_ARL_big) <- lambda_grid
  rownames(matrix_ARL_small) <- k_grid
  colnames(matrix_ARL_small) <- lambda_grid

  mod_1_AMFEWMA_list <- list()
  for (ii in 1:n_lambda) {
    mod_1_AMFEWMA_list[[ii]] <- list()
    for (jj in 1:n_k) {
      lambda <- lambda_grid[ii]
      k <- k_grid[jj]
      mod_1_AMFEWMA_list[[ii]][[jj]] <-
        AMFEWMA_PhaseI_given_pars_ST(mfdobj = mfdobj,
                                     mfdobj_tuning = mfdobj_tuning,
                                     lambda = lambda,
                                     k = k,
                                     ARL0 = ARL0,
                                     bootstrap_pars = bootstrap_pars,
                                     discrete_grid_length = discrete_grid_length,
                                     score_function = score_function,
                                     fev = fev,
                                     c = c,
                                     n_skip = n_skip)
      mod_ii_jj <- list(mod_1 = mod_1_AMFEWMA_list[[ii]][[jj]],
                        lambda = lambda,
                        k = k)
      mod_2_AMFEWMA_big <- AMFEWMA_PhaseII_ST(mfdobj_2 = fase2_big,
                                              mod_1 = mod_ii_jj,
                                              n_seq_2 = n_seq,
                                              l_seq_2 = l_seq)
      mod_2_AMFEWMA_small <- AMFEWMA_PhaseII_ST(mfdobj_2 = fase2_small,
                                                mod_1 = mod_ii_jj,
                                                n_seq_2 = n_seq,
                                                l_seq_2 = l_seq)
      ARL_big <- mod_2_AMFEWMA_big$ARL_2
      ARL_small <- mod_2_AMFEWMA_small$ARL_2

      matrix_ARL_big[jj, ii] <- ARL_big
      matrix_ARL_small[jj, ii] <- ARL_small
    }
  }
  permat_small <- matrix_ARL_small
  mat_big <- matrix_ARL_big

  # Algorithm to select best lambda and k

  ARL_2_min <- min(matrix_ARL_big)
  for (ii in 1:(n_k * n_lambda)) {
    min_1 <- which(matrix_ARL_small == min(matrix_ARL_small), arr.ind = TRUE)
    ARL_2_p <- matrix_ARL_big[min_1[1, 1], min_1[1, 2]]

    if (ARL_2_p <= (1 + epsilon) * ARL_2_min) {
      lambda_opt <- lambda_grid[min_1[1, 2]]
      k_opt <- k_grid[min_1[1, 1]]
      break
    } else {
      matrix_ARL_small[min_1[1, 1], min_1[1, 2]] <- 10000
    }
  }

  # Optimal parameters

  lambda <- lambda_opt
  k <- k_opt

  which_lambda_opt <- which(lambda_grid == lambda_opt)
  which_k_opt <- which(k_grid == k_opt)

  mod_1 <- mod_1_AMFEWMA_list[[which_lambda_opt]][[which_k_opt]]

  ret <- list(mod_1 = mod_1,
              lambda = lambda,
              k = k,
              c = c)

  return(ret)
}


AMFEWMA_PhaseI_given_pars_ST <- function(mfdobj,
                                         mfdobj_tuning,
                                         lambda,
                                         k,
                                         ARL0 = 200,
                                         bootstrap_pars = list(n_seq = 200,
                                                               l_seq = 2000),
                                         discrete_grid_length = 25,
                                         score_function = "huber",
                                         fev = 0.99,
                                         c,
                                         n_skip = 100) {

  nobs <- dim(mfdobj$coefs)[2]
  nvar <- dim(mfdobj$coefs)[3]

  n_seq <- bootstrap_pars$n_seq
  l_seq <- bootstrap_pars$l_seq

  mean_mfdobj <- fda::mean.fd(mfdobj)
  sd_mfdobj <- fda::sd.fd(mfdobj)

  rn <- mfdobj$basis$rangeval
  grid_points <- seq(rn[1], rn[2], length.out = discrete_grid_length)
  sd_fd <- fda::eval.fd(grid_points, sd_mfdobj)
  sd_vec <- as.numeric(sd_fd)
  k_fun <- k * sd_vec
  basis <- mfdobj$basis

  # (TRAINING)

  mfdobj_cen <- funcharts::scale_mfd(mfdobj,
                                     center = mean_mfdobj,
                                     scale = FALSE)
  Xeval_cen <- fda::eval.fd(grid_points, mfdobj_cen)
  X <- matrix(aperm(Xeval_cen, c(2, 1, 3)), nrow = nobs)

  if (!(score_function %in% c("huber", "tukey"))) {
    stop("score_function must be \"huber\" or \"tukey\"")
  }

  if (score_function == "huber") {
    huber <- TRUE
  } else {
    huber <- FALSE
  }

  Y_array_tra <- statisticY_EWMA_cpp(X,
                                     lambda = lambda,
                                     k = k_fun,
                                     huber = huber,
                                     idx = 1:nrow(X))
  sigmaY <- Rfast::cova(Y_array_tra)

  eigen <- eigen(sigmaY, symmetric = TRUE)
  var_spiegata <- cumsum(eigen$values) / sum(eigen$values)
  n_eigenvalues <- which(var_spiegata > fev)[1]
  A <- eigen$values[1:n_eigenvalues]
  B <- eigen$vectors[, 1:n_eigenvalues, drop = FALSE]
  inv_sigmaY_reg <- B %*% diag(1 / A) %*% t(B)

  # (TUNING)

  mfdobj_tuning_cen <- funcharts::scale_mfd(mfdobj_tuning,
                                            center = mean_mfdobj,
                                            scale = FALSE)
  Xeval_tun_cen <- fda::eval.fd(grid_points, mfdobj_tuning_cen)

  nobs_tun <- dim(Xeval_tun_cen)[2]
  nvar <- dim(Xeval_tun_cen)[3]
  Xtun <- matrix(aperm(Xeval_tun_cen, c(2, 1, 3)), nrow = nobs_tun)

  par_fun_tun <- function(qq) {
    idx_tun <- sample(1:nobs_tun, l_seq, replace = TRUE)
    Y_array_tun <- statisticY_EWMA_cpp(
      Xtun,
      lambda = lambda,
      k = k_fun,
      huber = huber,
      idx = idx_tun
    )
    score_tun <- Y_array_tun %*% B
    score_std <- t(t(score_tun[, , drop = FALSE]) / sqrt(A))
    V2 <- rowSums(pmax(score_std^2 - c, 0))
    return(V2)
  }

  V2_mat <- do.call(rbind, lapply(1:n_seq, function(x) par_fun_tun(x)))

  ARL <- -100
  iter <- 0
  hmin <- 0
  hmax <- max(V2_mat)
  while (abs(ARL - ARL0) > 1e-2 & (hmax - hmin) > 1e-2 & iter < 50) {
    iter <- iter + 1
    h <- (hmin + hmax) / 2
    cond <- V2_mat[, -(1:n_skip)] > h
    if (max(V2_mat[, -(1:n_skip)]) < h) {
      ARL <- Inf
      hmax <- h
    } else {
      RL <- apply(cond, 1, function(x) which(x)[1])
      ARL <- mean(RL, na.rm = TRUE)
      if (ARL < ARL0) {
        hmin <- h
      } else {
        hmax <- h
      }
    }
  }
  ARL0 <- ARL

  ret <- list(mfdobj = mfdobj,
              mfdobj_tuning = mfdobj_tuning,
              inv_sigmaY_reg = inv_sigmaY_reg,
              mean_mfdobj = mean_mfdobj,
              h = h,
              ARL0 = ARL0,
              lambda = lambda,
              k = k_fun,
              grid_points = grid_points,
              V2_mat = V2_mat,
              n_skip = n_skip,
              huber = huber,
              vectors = B,
              values = A,
              c = c)
  return(ret)
}
