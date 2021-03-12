#' T^2 and SPE control charts for multivariate functional data
#'
#' This function builds a data frame needed to plot
#' the Hotelling's T^2 and squared prediction error (SPE)
#' control charts
#' based on multivariate functional principal component analysis
#' (MFPCA) performed
#' on multivariate functional data,
#' proposed in Capezza et al. (2020) together with the scalar control chart
#' and used also to build the
#' functional regression control chart proposed in Centofanti et al. (2020)
#' (this function is used by \code{\link{regr_cc_fof}}).
#' The training data have already been used to fit the model.
#' A tuning data set can be provided that is used to estimate
#' the control chart limits.
#' A phase II data set contains the observations
#' to be monitored with the built control charts.
#'
#' @param pca
#' An object of class \code{pca_mfd}
#' obtained by doing MFPCA on the
#' training set of multivariate functional data.
#' @param components
#' A vector of integers with the components over which
#' to project the multivariate functional data.
#' @param tuning_data
#' An object of class \code{mfd} containing
#' the tuning set of the multivariate functional data, used to estimate the
#' T^2 and SPE control chart limits.
#' If NULL, the training data, i.e. the data used to fit the MFPCA model,
#' are also used as the tuning data set, i.e. \code{tuning_data=pca$data}.
#' Default is NULL.
#' @param newdata
#' An object of class \code{mfd} containing
#' the phase II set of the multivariate functional data to be monitored.
#' @param alpha
#' A named list with two elements, named \code{T2} and \code{spe},
#' respectively, each containing
#' the desired Type I error probability of
#' the corresponding control chart.
#' Note that at the moment you have to take into account manually
#' the family-wise error rate and adjust
#' the two values accordingly. See Capezza et al. (2020) and
#' Centofanti et al. (2020)
#' for additional details. Default value is
#' \code{list(T2 = 0.025, spe = 0.025)}.
#' @param limits
#' A character value.
#' If "standard", it estimates the control limits on the tuning
#' data set. If "cv", the function calculates the control limits only on the
#' training data using cross-validation
#' using \code{calculate_cv_limits}. Default is "standard".
#' @param seed
#' If \code{limits=="cv"},
#' since the split in the k groups is random,
#' you can fix a seed to ensure reproducibility.
#' Otherwise, this argument is ignored.
#' @param nfold
#' If \code{limits=="cv"}, this gives the number of groups k
#' used for k-fold cross-validation.
#' If it is equal to the number of observations in the training data set,
#' then we have
#' leave-one-out cross-validation.
#' Otherwise, this argument is ignored.
#' @param ncores
#' If \code{limits=="cv"}, if you want perform the analysis
#' in the k groups in parallel,
#' give the number of cores/threads.
#' Otherwise, this argument is ignored.
#'
#' @return
#' A \code{data.frame} with as many rows as the number of
#' multivariate functional observations in the phase II data set and
#' the following columns:
#'
#' * one \code{id} column identifying the multivariate functional observation
#' in the phase II data set,
#'
#' * one \code{T2} column containing the Hotelling T^2 statistic
#' calculated for all observations,
#'
#' * one column per each functional variable,
#' containing its contribution to the T^2 statistic,
#'
#' * one \code{spe} column containing the SPE statistic calculated
#' for all observations,
#'
#' * one column per each functional variable,
#' containing its contribution to the SPE statistic,
#'
#' * \code{T2_lim} gives the upper control limit of
#' the Hotelling's T^2 control chart,
#'
#' * one \code{contribution_T2_*_lim} column per each
#' functional variable giving the
#' limits of the contribution of that variable
#' to the Hotelling's T^2 statistic,
#'
#' * \code{spe_lim} gives the upper control limit of the SPE control chart
#'
#' * one \code{contribution_spe*_lim} column per each
#' functional variable giving the
#' limits of the contribution of that variable to the SPE statistic.
#'
#' @export
#'
#' @seealso \code{\link{regr_cc_fof}}
#'
#' @references
#' Capezza C, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
#' Control charts for
#' monitoring ship operating conditions and CO2 emissions
#' based on scalar-on-function regression.
#' \emph{Applied Stochastic Models in Business and Industry},
#' 36(3):477--500.
#' <doi:10.1002/asmb.2507>
#'
#' Centofanti F, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
#' Functional Regression Control Chart.
#' \emph{Technometrics}. <doi:10.1080/00401706.2020.1753581>
#'
#' library(funcharts)
#' data("air")
#' air <- lapply(air, function(x) x[1:220, , drop = FALSE])
#' fun_covariates <- c("CO", "temperature")
#' mfdobj_x <- get_mfd_list(air[fun_covariates],
#'                          n_basis = 15,
#'                          lambda = 1e-2)
#' y <- rowMeans(air$NO2)
#' y1 <- y[1:100]
#' y_tuning <- y[101:200]
#' y2 <- y[201:220]
#' mfdobj_x1 <- mfdobj_x[1:100]
#' mfdobj_x_tuning <- mfdobj_x[101:200]
#' mfdobj_x2 <- mfdobj_x[201:220]
#' pca <- pca_mfd(mfdobj_x1)
#' components <- 1:which(cumsum(pca$varprop) >= .90)[1]
#' cclist <- control_charts_pca(pca = pca,
#'                              components = components,
#'                              tuning_data = mfdobj_x_tuning,
#'                              newdata = mfdobj_x2)
#' plot_control_charts(cclist)
#'
control_charts_pca <- function(pca,
                               components,
                               tuning_data = NULL,
                               newdata,
                               alpha = list(T2 = .025, spe = .025),
                               limits = "standard",
                               seed = 0,
                               nfold = 5,
                               ncores = 1) {

  if (!is.list(pca)) {
    stop("pca must be a list produced by pca_mfd.")
  }

  if (!identical(names(pca), c(
    "harmonics",
    "values",
    "scores",
    "varprop",
    "meanfd",
    "pcscores",
    "data",
    "scale"
  ))) {
    stop("pca must be a list produced by pca_mfd.")
  }

  if (is.null(tuning_data)) tuning_data <- pca$data

  if (length(limits) != 1) {
    stop("Only one type of 'limits' allowed.")
  }
  if (!(limits %in% c("standard", "cv"))) {
    stop("'limits' argument can be only 'standard' or 'cv'.")
  }

  if (limits == "standard") lim <- calculate_limits(
    pca = pca,
    tuning_data = tuning_data,
    components = components,
    alpha = alpha)
  if (limits == "cv") lim <- calculate_cv_limits(
    pca = pca,
    components = components,
    alpha = alpha,
    seed = seed,
    nfold = nfold,
    ncores = ncores)

  T2 <- get_T2(pca, components, newdata = newdata)
  spe <- get_spe(pca, components, newdata = newdata)
  id <- data.frame(id = newdata$fdnames[[2]])

  cbind(id, T2, spe, lim)
}


#' Scalar-on-Function Regression Control Chart
#'
#' This function builds a data frame needed
#' to plot the scalar-on-function regression control chart,
#' based on a fitted function-on-function linear regression model and
#' proposed in Capezza et al. (2020) together with the Hotelling's T^2 and
#' squared prediction error control charts.
#' The training data have already been used to fit the model.
#' A tuning data set can be provided that is used to estimate
#' the control chart limits.
#' A phase II data set contains the observations to be monitored
#' with the built control charts.
#'
#' @param object
#' A list obtained as output from \code{sof_pc},
#' i.e. a fitted scalar-on-function linear regression model.
#' @param mfdobj_x_new
#' An object of class \code{mfd} containing
#' the phase II data set of the functional covariates observations.
#' @param y_new
#' A numeric vector containing the observations of
#' the scalar response variable
#' in the phase II data set.
#' @param alpha
#' A numeric value indicating the Type I error for
#' the regression control chart
#' and such that this function returns the \code{1-alpha} prediction interval
#' on the response.
#' Default is 0.05.
#'
#' @return
#' A \code{data.frame} with as many rows as the
#' number of functional replications in \code{mfdobj_x_new},
#' with the following columns:
#'
#' * \code{y_hat}: the predictions of the response variable
#' corresponding to \code{mfdobj_x_new},
#'
#' * \code{y}: the same as the argument \code{y_new} given as input
#' to this function,
#'
#' * \code{lwr}: lower limit of the \code{1-alpha} prediction interval
#' on the response,
#'
#' * \code{pred_err}: prediction error calculated as \code{y-y_hat},
#'
#' * \code{pred_err_sup}: upper limit of the \code{1-alpha} prediction interval
#' on the prediction error,
#'
#' * \code{pred_err_inf}: lower limit of the \code{1-alpha} prediction interval
#' on the prediction error.
#'
#' @export
#'
#' @references
#' Capezza C, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
#' Control charts for
#' monitoring ship operating conditions and CO2 emissions
#' based on scalar-on-function regression.
#' \emph{Applied Stochastic Models in Business and Industry},
#' 36(3):477--500.
#' <doi:10.1002/asmb.2507>
#'
#' @examples
#' library(funcharts)
#' air <- lapply(air, function(x) x[1:100, , drop = FALSE])
#' fun_covariates <- c("CO", "temperature")
#' mfdobj_x <- get_mfd_list(air[fun_covariates],
#'                          n_basis = 15,
#'                          lambda = 1e-2)
#' y <- rowMeans(air$NO2)
#' y1 <- y[1:80]
#' y2 <- y[81:100]
#' mfdobj_x1 <- mfdobj_x[1:80]
#' mfdobj_x2 <- mfdobj_x[81:100]
#' mod <- sof_pc(y1, mfdobj_x1)
#' cclist <- regr_cc_sof(object = mod,
#'                       y_new = y2,
#'                       mfdobj_x_new = mfdobj_x2)
#' plot_control_charts(cclist)
#'
regr_cc_sof <- function(object,
                        y_new,
                        mfdobj_x_new,
                        alpha = .05) {

  if (!is.list(object)) {
    stop("object must be a list produced by sof_pc.")
  }

  if (!identical(names(object), c(
    "mod",
    "pca",
    "beta_fd",
    "components",
    "selection",
    "single_min_variance_explained",
    "tot_variance_explained",
    "gcv",
    "PRESS"
  ))) {
    stop("object must be a list produced by sof_pc.")
  }

  if (!is.numeric(y_new)) {
    stop("y_new must be numeric.")
  }
  if (!is.null(mfdobj_x_new)) {
    if (!is.mfd(mfdobj_x_new)) {
      stop("mfdobj_x_new must be an object from mfd class.")
    }
    if (dim(mfdobj_x_new$coefs)[3] != dim(object$pca$data$coefs)[3]) {
      stop(paste0("mfdobj_x_new must have the same number of variables ",
      "as training data."))
    }
  }
  nobsx <- dim(mfdobj_x_new$coefs)[2]
  nobsy <- length(y_new)
  if (nobsx != nobsy) {
    stop(paste0("y_new and mfdobj_x_new must have ",
                "the same number of observations."))
  }

  y_hat <- predict_sof_pc(object = object,
                          newdata = mfdobj_x_new,
                          alpha = alpha)

  data.frame(
    y_hat,
    y = y_new,
    pred_err = y_new - y_hat$fit,
    pred_err_sup = y_hat$upr - y_hat$fit,
    pred_err_inf = y_hat$lwr - y_hat$fit)

}


#' Control charts for monitoring multivariate functional covariates
#' and a scalar response
#'
#' @description
#' This function builds a data frame needed to
#' plot control charts
#' for monitoring a multivariate functional covariates based on
#' multivariate functional principal component analysis (MFPCA) and
#' a related scalar response variable using the
#' scalar-on-function regression control chart,
#' as proposed in Capezza et al. (2020).
#'
#' In particular, this function provides:
#'
#' * the Hotelling's T^2 control chart,
#'
#' * the squared prediction error (SPE) control chart,
#'
#' * the scalar regression control chart.
#'
#' This function calls \code{control_charts_pca} for the control charts on
#' the multivariate functional covariates and \code{\link{regr_cc_sof}}
#' for the scalar regression control chart.
#'
#' The training data have already been used to fit the model.
#' A tuning data set can be provided that is used to estimate
#' the control chart limits.
#' A phase II data set contains the observations to be monitored
#' with the built control charts.
#'
#'
#' @param mod
#' A list obtained as output from \code{sof_pc},
#' i.e. a fitted scalar-on-function linear regression model.
#' @param y_test
#' A numeric vector containing the observations
#' of the scalar response variable
#' in the phase II data set.
#' @param mfdobj_x_test
#' An object of class \code{mfd} containing
#' the phase II data set of the functional covariates observations.
#' @param mfdobj_x_tuning
#' An object of class \code{mfd} containing
#' the tuning set of the multivariate functional data, used to estimate the
#' T^2 and SPE control chart limits.
#' If NULL, the training data, i.e. the data used to fit the MFPCA model,
#' are also used as the tuning data set, i.e. \code{tuning_data=pca$data}.
#' Default is NULL.
#' @param alpha
#' A named list with three elements, named \code{T2}, \code{spe},
#' and code{y},
#' respectively, each containing
#' the desired Type I error probability of the corresponding control chart
#' (\code{T2} corresponds to the T^2 control chart,
#' \code{spe}  corresponds to the SPE control chart,
#' \code{y} corresponds to the scalar regression control chart).
#' Note that at the moment you have to take into account manually
#' the family-wise error rate and adjust
#' the two values accordingly. See Capezza et al. (2020)
#' for additional details. Default value is
#' \code{list(T2 = 0.0125, spe = 0.0125, y = 0.025)}.
#' @param limits
#' A character value.
#' If "standard", it estimates the control limits on the tuning
#' data set. If "cv", the function calculates the control limits only on the
#' training data using cross-validation
#' using \code{calculate_cv_limits}. Default is "standard".
#' @param seed
#' If \code{limits=="cv"},
#' since the split in the k groups is random,
#' you can fix a seed to ensure reproducibility.
#' Otherwise, this argument is ignored.
#' @param nfold
#' If \code{limits=="cv"}, this gives the number of groups k
#' used for k-fold cross-validation.
#' If it is equal to the number of observations in the training data set,
#' then we have
#' leave-one-out cross-validation.
#' Otherwise, this argument is ignored.
#' @param ncores
#' If \code{limits=="cv"}, if you want perform the analysis
#' in the k groups in parallel,
#' give the number of cores/threads.
#' Otherwise, this argument is ignored.
#'
#' @return
#' A \code{data.frame} with as many rows as the number of
#' multivariate functional observations in the phase II data set and
#' the following columns:
#'
#' * one \code{id} column identifying the multivariate functional observation
#' in the phase II data set,
#'
#' * one \code{T2} column containing the Hotelling T^2 statistic calculated
#' for all observations,
#'
#' * one column per each functional variable, containing its contribution
#' to the T^2 statistic,
#'
#' * one \code{spe} column containing the SPE statistic calculated
#' for all observations,
#'
#' * one column per each functional variable, containing its contribution
#' to the SPE statistic,
#'
#' * \code{T2_lim} gives the upper control limit of the
#' Hotelling's T^2 control chart,
#'
#' * one \code{contribution_T2_*_lim} column per each
#' functional variable giving the
#' limits of the contribution of that variable to the
#' Hotelling's T^2 statistic,
#'
#' * \code{spe_lim} gives the upper control limit of the SPE control chart
#'
#' * one \code{contribution_spe*_lim} column per
#' each functional variable giving the
#' limits of the contribution of that variable to the SPE statistic.
#'
#' * \code{y_hat}: the predictions of the response variable
#' corresponding to \code{mfdobj_x_new},
#'
#' * \code{y}: the same as the argument \code{y_new}
#' given as input to this function,
#'
#' * \code{lwr}: lower limit of the \code{1-alpha}
#' prediction interval on the response,
#'
#' * \code{pred_err}: prediction error calculated as \code{y-y_hat},
#'
#' * \code{pred_err_sup}: upper limit of the \code{1-alpha}
#' prediction interval on the prediction error,
#'
#' * \code{pred_err_inf}: lower limit of the \code{1-alpha}
#'  prediction interval on the prediction error.
#'
#' @export
#'
#' @seealso \code{\link{control_charts_pca}}, \code{\link{regr_cc_sof}}
#'
#' @examples
#' library(funcharts)
#' data("air")
#' air <- lapply(air, function(x) x[201:300, , drop = FALSE])
#' fun_covariates <- c("CO", "temperature")
#' mfdobj_x <- get_mfd_list(air[fun_covariates],
#'                          n_basis = 15,
#'                          lambda = 1e-2)
#' y <- rowMeans(air$NO2)
#' y1 <- y[1:60]
#' y2 <- y[91:100]
#' mfdobj_x1 <- mfdobj_x[1:60]
#' mfdobj_x_tuning <- mfdobj_x[61:90]
#' mfdobj_x2 <- mfdobj_x[91:100]
#' mod <- sof_pc(y1, mfdobj_x1)
#' cclist <- control_charts_sof_pc(mod = mod,
#'                                 y_test = y2,
#'                                 mfdobj_x_test = mfdobj_x2,
#'                                 mfdobj_x_tuning = mfdobj_x_tuning)
#' plot_control_charts(cclist)
#'
control_charts_sof_pc <- function(mod,
                                  y_test,
                                  mfdobj_x_test,
                                  mfdobj_x_tuning = NULL,
                                  alpha = list(
                                    T2  = .0125,
                                    spe = .0125,
                                    y   = .025),
                                  limits = "standard",
                                  seed = 0,
                                  nfold = NULL,
                                  ncores = 1) {

  if (!is.list(mod)) {
    stop("object must be a list produced by sof_pc.")
  }

  if (!identical(names(mod), c(
    "mod",
    "pca",
    "beta_fd",
    "components",
    "selection",
    "single_min_variance_explained",
    "tot_variance_explained",
    "gcv",
    "PRESS"
  ))) {
    stop("object must be a list produced by sof_pc.")
  }

  if (is.null(mfdobj_x_tuning)) mfdobj_x_tuning <- mod$pca$data

  pca <- mod$pca
  components <- which(colnames(pca$pcscores) %in% names(mod$mod$coefficients))
  if (is.null(nfold)) nfold <- length(mod$pca$data$fdnames[[2]])
  out <- control_charts_pca(pca = pca,
                            components = components,
                            newdata = mfdobj_x_test,
                            tuning_data = mfdobj_x_tuning,
                            alpha = alpha,
                            limits = limits,
                            seed = seed,
                            nfold = nfold,
                            ncores = ncores)
  y <- regr_cc_sof(object = mod,
                   y_new = y_test,
                   mfdobj_x_new = mfdobj_x_test,
                   alpha = alpha$y)
  cbind(out, y)

}



#' Functional Regression Control Chart
#'
#' It builds a data frame needed to plot the
#' Functional Regression Control Chart
#' introduced in Centofanti et al. (2020), based on a fitted
#' function-on-function linear regression model.
#' The training data have already been used to fit the model.
#' A tuning data set can be provided that is used to estimate
#' the control chart limits.
#' A phase II data set contains the observations to be monitored
#' with the built control charts.
#'
#' @param object
#' A list obtained as output from \code{fof_pc},
#' i.e. a fitted function-on-function linear regression model.
#' @param mfdobj_y_new
#' An object of class \code{mfd} containing
#' the phase II data set of the functional response
#' observations to be monitored.
#' @param mfdobj_x_new
#' An object of class \code{mfd} containing
#' the phase II data set of the functional covariates
#' observations to be monitored.
#' @param mfdobj_y_tuning
#' An object of class \code{mfd} containing
#' the tuning data set of the functional response observations,
#' used to estimate the control chart limits.
#' If NULL, the training data, i.e. the data used to fit the
#' function-on-function linear regression model,
#' are also used as the tuning data set, i.e.
#' \code{mfdobj_y_tuning=object$pca_y$data}.
#' Default is NULL.
#' @param mfdobj_x_tuning
#' An object of class \code{mfd} containing
#' the tuning data set of the functional covariates observations,
#' used to estimate the control chart limits.
#' If NULL, the training data, i.e. the data used to fit the
#' function-on-function linear regression model,
#' are also used as the tuning data set, i.e.
#' \code{mfdobj_x_tuning=object$pca_x$data}.
#' Default is NULL.
#' @param alpha
#' A named list with two elements,
#' named \code{T2} and \code{spe}, respectively, each containing
#' the desired Type I error probability of the corresponding control chart.
#' Note that at the moment you have to take into account manually
#' the family-wise error rate and adjust
#' the two values accordingly. See Centofanti et al. (2020)
#' for additional details.
#' Default value is \code{list(T2 = 0.025, spe = 0.025)}.
#'
#' @return
#' A \code{data.frame} containing the output of the
#' function \code{control_charts_pca} applied to
#' the prediction errors.
#' @export
#'
#' @seealso \code{\link{control_charts_pca}}
#'
#' @references
#' Centofanti F, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
#' Functional Regression Control Chart.
#' \emph{Technometrics}. <doi:10.1080/00401706.2020.1753581>
#'
#' @examples
#' library(funcharts)
#' data("air")
#' air <- lapply(air, function(x) x[1:100, , drop = FALSE])
#' fun_covariates <- c("CO", "temperature")
#' mfdobj_x <- get_mfd_list(air[fun_covariates],
#'                          n_basis = 15,
#'                          lambda = 1e-2)
#' mfdobj_y <- get_mfd_list(air["NO2"],
#'                          n_basis = 15,
#'                          lambda = 1e-2)
#' mfdobj_y1 <- mfdobj_y[1:60]
#' mfdobj_y_tuning <- mfdobj_y[61:90]
#' mfdobj_y2 <- mfdobj_y[91:100]
#' mfdobj_x1 <- mfdobj_x[1:60]
#' mfdobj_x_tuning <- mfdobj_x[61:90]
#' mfdobj_x2 <- mfdobj_x[91:100]
#' mod_fof <- fof_pc(mfdobj_y1, mfdobj_x1)
#' cclist <- regr_cc_fof(mod_fof,
#'                       mfdobj_y_new = mfdobj_y2,
#'                       mfdobj_x_new = mfdobj_x2,
#'                       mfdobj_y_tuning = NULL,
#'                       mfdobj_x_tuning = NULL)
#' plot_control_charts(cclist)
#'
regr_cc_fof <- function(object,
                        mfdobj_y_new,
                        mfdobj_x_new,
                        mfdobj_y_tuning = NULL,
                        mfdobj_x_tuning = NULL,
                        alpha = list(T2 = 0.025, spe = 0.025)) {

  if (!is.list(object)) {
    stop("object must be a list produced by fof_pc.")
  }

  if (!identical(names(object), c(
    "mod",
    "beta_fd",
    "fitted.values",
    "residuals_original_scale",
    "residuals",
    "type_residuals",
    "pca_x",
    "pca_y",
    "pca_res",
    "components_x",
    "components_y",
    "components_res",
    "y_standardized",
    "tot_variance_explained_x",
    "tot_variance_explained_y",
    "tot_variance_explained_res",
    "get_studentized_residuals"
  ))) {
    stop("object must be a list produced by fof_pc.")
  }

  if (is.null(mfdobj_y_tuning) | is.null(mfdobj_x_tuning)) {
    mfdobj_y_tuning <- object$pca_y$data
    mfdobj_x_tuning <- object$pca_x$data
  }
  tuning <- predict_fof_pc(
    object = object,
    mfdobj_y_new = mfdobj_y_tuning,
    mfdobj_x_new = mfdobj_x_tuning)

  phase_II <- predict_fof_pc(
    object = object,
    mfdobj_y_new = mfdobj_y_new,
    mfdobj_x_new = mfdobj_x_new)

  out <- control_charts_pca(
    pca = object$pca_res,
    components = object$components_res,
    tuning_data = tuning$pred_error,
    newdata = phase_II$pred_error,
    alpha = alpha,
    limits = "standard"
  )
  select(out, -contains("contribution_"))

}



#' Plot control charts
#'
#' This function takes as input a data frame produced
#' with functions such as
#' \code{\link{control_charts_pca}} and \code{\link{control_charts_sof_pc}} and
#' produces a ggplot with the desired control charts, i.e.
#' it plots a point for each
#' observation in the phase II data set against
#' the corresponding control limits.
#'
#' @param cclist
#' A \code{data.frame} produced by
#' \code{\link{control_charts_pca}}, \code{\link{control_charts_sof_pc}}
#' \code{\link{regr_cc_fof}}, or \code{\link{regr_cc_sof}}.
#'
#' @return A ggplot with the functional control charts.
#'
#' @details Out-of-control points are signaled by colouring them in red.
#'
#' @export
#' @examples
#' library(funcharts)
#' data("air")
#' air <- lapply(air, function(x) x[1:100, , drop = FALSE])
#' fun_covariates <- c("CO", "temperature")
#' mfdobj_x <- get_mfd_list(air[fun_covariates],
#'                          n_basis = 15,
#'                          lambda = 1e-2)
#' mfdobj_y <- get_mfd_list(air["NO2"],
#'                          n_basis = 15,
#'                          lambda = 1e-2)
#' mfdobj_y1 <- mfdobj_y[1:60]
#' mfdobj_y_tuning <- mfdobj_y[61:90]
#' mfdobj_y2 <- mfdobj_y[91:100]
#' mfdobj_x1 <- mfdobj_x[1:60]
#' mfdobj_x_tuning <- mfdobj_x[61:90]
#' mfdobj_x2 <- mfdobj_x[91:100]
#' mod_fof <- fof_pc(mfdobj_y1, mfdobj_x1)
#' cclist <- regr_cc_fof(mod_fof,
#'                       mfdobj_y_new = mfdobj_y2,
#'                       mfdobj_x_new = mfdobj_x2,
#'                       mfdobj_y_tuning = NULL,
#'                       mfdobj_x_tuning = NULL)
#' plot_control_charts(cclist)
#'
plot_control_charts <- function(cclist) {

  if (!is.data.frame(cclist)) {
    stop(paste0("cclist must be a data frame ",
                "containing data to produce control charts."))
  }

  df_hot <- NULL
  if (!is.null(cclist$T2)) {
    df_hot <- data.frame(statistic = cclist$T2,
                         lcl = 0,
                         ucl = cclist$T2_lim) %>%
      mutate(id = 1:n(),
             ooc = .data$statistic > .data$ucl,
             type = "HOTELLING~T^2~CONTROL~CHART")

    hot_range <- df_hot %>%
      select(.data$statistic, .data$ucl, .data$lcl) %>%
      range() %>%
      diff()

    df_hot <- df_hot %>%
      mutate(ytext = case_when(
        statistic < lcl ~ statistic - hot_range * 0.2,
        statistic > ucl ~ statistic + hot_range * 0.2,
        TRUE ~ statistic,
      ))
  }

  df_spe <- NULL
  if (!is.null(cclist$spe)) {
    df_spe <- data.frame(statistic = cclist$spe,
                         lcl = 0,
                         ucl = cclist$spe_lim) %>%
      mutate(id = 1:n(),
             ooc = .data$statistic > .data$ucl,
             type = "SPE~CONTROL~CHART")

    spe_range <- df_spe %>%
      select(.data$statistic, .data$ucl, .data$lcl) %>%
      range() %>%
      diff()

    df_spe <- df_spe %>%
      mutate(ytext = case_when(
        statistic < lcl ~ statistic - spe_range * 0.2,
        statistic > ucl ~ statistic + spe_range * 0.2,
        TRUE ~ statistic,
      ))
  }

  df_y <- NULL
  if (!is.null(cclist$pred_err)) {
    df_y <- data.frame(statistic = cclist$pred_err,
                       lcl = cclist$pred_err_inf,
                       ucl = cclist$pred_err_sup)

    y_range <- df_y %>%
      select(c(.data$statistic, .data$ucl, .data$lcl)) %>%
      range() %>%
      diff()

    df_y <- df_y %>%
      mutate(ytext = case_when(
        statistic < lcl ~ statistic - y_range * 0.2,
        statistic > ucl ~ statistic + y_range * 0.2,
        TRUE ~ statistic,
      ))

    df_y <- df_y %>%
      mutate(id = 1:n(),
             ooc = .data$statistic > .data$ucl | .data$statistic < .data$lcl,
             type = "REGRESSION~CONTROL~CHART")

  }

  plot_list <- list()

  if (!is.null(cclist$T2)) {
    plot_list$p_hot <- ggplot(df_hot, aes(x = .data$id, y = .data$statistic)) +
      geom_line() +
      geom_point(aes(colour = .data$ooc)) +
      geom_blank(aes(y = 0)) +
      geom_point(aes(y = .data$ucl),
                 pch = "-", size = 5) +
      scale_x_continuous(
        limits = c(0, max(df_hot$id) + 1),
        breaks = seq(1, max(df_hot$id), by = round(max(df_hot$id) / 50) + 1),
        expand = c(0, 0)) +
      scale_color_manual(values = c("black", "red")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5)) +
      geom_text(aes(y = .data$ytext, label = .data$id),
                data = filter(df_hot, .data$ooc),
                size = 3) +
      xlab("Observation") +
      ylab(expression(T^2~statistic)) +
      ggtitle(expression(HOTELLING~T^2~CONTROL~CHART))
  }

  if (!is.null(cclist$spe)) {
    plot_list$p_spe <- ggplot(df_spe, aes(x = .data$id, y = .data$statistic)) +
      geom_line() +
      geom_point(aes(colour = .data$ooc)) +
      geom_blank(aes(y = 0)) +
      geom_point(aes(y = .data$ucl),
                 pch = "-", size = 5) +
      scale_x_continuous(
        limits = c(0, max(df_spe$id) + 1),
        breaks = seq(1, max(df_spe$id), by = round(max(df_spe$id) / 50) + 1),
        expand = c(0, 0)) +
      scale_color_manual(values = c("black", "red")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5)) +
      geom_text(aes(y = .data$ytext, label = .data$id),
                data = filter(df_spe, .data$ooc),
                size = 3) +
      xlab("Observation") +
      ylab(expression(SPE~statistic)) +
      ggtitle(expression(SPE~CONTROL~CHART))
  }

  if (!is.null(cclist$pred_err)) {
    plot_list$p_y <- ggplot(df_y, aes(x = .data$id, y = .data$statistic)) +
      geom_line() +
      geom_point(aes(colour = .data$ooc)) +
      geom_blank(aes(y = 0)) +
      geom_line(aes(y = .data$lcl), lty = 2) +
      geom_line(aes(y = .data$ucl), lty = 2) +
      scale_x_continuous(
        limits = c(0, max(df_y$id) + 1),
        breaks = seq(1, max(df_y$id), by = round(max(df_y$id) / 50) + 1),
        expand = c(0, 0)) +
      scale_color_manual(values = c("black", "red")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5)) +
      geom_text(aes(y = .data$ytext, label = .data$id),
                data = filter(df_y, .data$ooc),
                size = 3) +
      xlab("Observation") +
      ylab(expression("Prediction error [t]")) +
      ggtitle(expression(REGRESSION~CONTROL~CHART))
  }

  patchwork::wrap_plots(plot_list, ncol = 1)

}
