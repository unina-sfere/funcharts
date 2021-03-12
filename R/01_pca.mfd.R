#' Multivariate functional principal components analysis
#'
#' Multivariate functional principal components analysis (MFPCA)
#' performed on an object of class \code{mfd}.
#' It is a wrapper to \code{fda::\link[fda]{pca.fd}},
#' providing some additional arguments.
#'
#' @param mfdobj
#' A multivariate functional data object of class mfd.
#' @param scale
#' If TRUE, it scales data before doing MFPCA
#' using \code{scale_mfd}. Default is TRUE.
#' @param nharm
#' Number of multivariate functional principal components
#' to be calculated. Default is 20.
#'
#' @return
#' Modified \code{pca.fd} object, with
#' multivariate functional principal component scores summed over variables
#' (\code{fda::\link[fda]{pca.fd}} returns an array of scores
#' when providing a multivariate functional data object).
#' Moreover, the multivariate functional principal components
#' given in \code{harmonics}
#' are converted to the \code{mfd} class.
#' @export
#' @seealso \code{\link{scale_mfd}}
#'
#' @examples
#' library(funcharts)
#' mfdobj <- data_sim_mfd()
#' pca_obj <- pca_mfd(mfdobj)
#' plot_pca_mfd(pca_obj)
#'
pca_mfd <- function(mfdobj, scale = TRUE, nharm = 20) {

  obs_names <- mfdobj$fdnames[[2]]
  variables <- mfdobj$fdnames[[3]]

  data <- scale_mfd(mfdobj, scale = scale)

  data_pca <- if (length(variables) == 1)
    fd(data$coefs[, , 1],
       data$basis,
       data$fdnames) else data

  nobs <- length(obs_names)
  nvar <- length(variables)
  nbasis <- mfdobj$basis$nbasis
  nharm <- min(nobs - 1, nvar * nbasis, nharm)

  pca <- pca.fd(data_pca, nharm = nharm, centerfns = FALSE)
  pca$harmonics$fdnames[c(1, 3)] <- mfdobj$fdnames[c(1, 3)]

  pca$pcscores <- if(length(dim(pca$scores)) == 3) {
    apply(pca$scores, 1:2, sum)
  } else pca$scores

  pc_names <- pca$harmonics$fdnames[[2]]
  rownames(pca$scores) <- rownames(pca$pcscores) <- obs_names
  colnames(
    pca$scores) <- colnames(pca$pcscores) <- pc_names
  if (length(variables) > 1) dimnames(pca$scores)[[3]] <- variables
  pca$data <- mfdobj
  pca$scale <- scale

  if (length(variables) > 1) {
    coefs <- pca$harmonics$coefs
  } else {
    coefs <- array(pca$harmonics$coefs, dim = c(nrow(pca$harmonics$coefs),
                                       ncol(pca$harmonics$coefs),
                                       1))
    dimnames(coefs) <- list(
      dimnames(pca$harmonics$coefs)[[1]],
      dimnames(pca$harmonics$coefs)[[2]],
      variables
    )
  }


  if (length(variables) == 1) {
    pca$scores <- array(pca$scores,
                        dim = c(nrow(pca$scores),ncol(pca$scores), 1))
    dimnames(pca$scores) <-
      list(obs_names, pc_names, variables)

  }

  pca$harmonics <- mfd(coefs,
                       pca$harmonics$basis,
                       pca$harmonics$fdnames,
                       NULL,
                       NULL)

  class(pca) <- c("pca_mfd", "pca.fd")
  pca

}

#' Plot the harmonics of a \code{pca_mfd} object
#'
#' @param pca
#' A fitted multivariate functional principal component analysis
#' (MFPCA) object of class \code{pca_mfd}.
#' @param harm
#' A vector of integers with the harmonics to plot.
#' If 0, all harmonics are plotted. Default is 0.
#' @param scaled
#' If TRUE, eigenfunctions are multiplied by the square root of the
#' corresponding eigenvalues, if FALSE the are not scaled and the
#' all have unit norm.
#' Default is FALSE
#'
#' @return
#' A ggplot of the harmonics/multivariate functional
#' principal components contained in the object \code{pca}.
#' @export
#' @examples
#' library(funcharts)
#' mfdobj <- data_sim_mfd()
#' pca_obj <- pca_mfd(mfdobj)
#' plot_pca_mfd(pca_obj)
#'
plot_pca_mfd <- function(pca, harm = 0, scaled = FALSE) {

  if (harm[1] == 0) harm <- seq_along(pca$harmonics$fdnames[[2]])
  if (scaled) {
    scaled_coefs <- apply(
      pca$harmonics$coefs[, harm, , drop = FALSE],
      3,
      function(x) t(t(x) * sqrt(pca$values[harm])))
    nbasis <- pca$harmonics$basis$nbasis
    nharm <- length(harm)
    nvar <- length(pca$harmonics$fdnames[[3]])
    scaled_coefs <- array(scaled_coefs, dim = c(nbasis, nharm, nvar))
    dimnames(scaled_coefs) <- dimnames(pca$harmonics$coefs[, harm, , drop = FALSE])
    pca$harmonics$coefs <- scaled_coefs
  }

  p_functions <- ggplot() +
    geom_mfd(aes(col = id), mfdobj = pca$harmonics[harm]) +
    ggtitle("Eigenfunctions")

  components <- which(cumsum(pca$varprop) < .99)
  # p_values <- data.frame(eigenvalues = pca$values[components]) %>%
  #   mutate(n_comp = 1:n()) %>%
  #   ggplot() +
  #   geom_col(aes(n_comp, eigenvalues)) +
  #   theme_bw() +
  #   xlab("Number of components")

  p_functions

}

#' Calculate the scores given a fitted \code{pca_mfd} object
#'
#' This function calculates the
#' multivariate functional principal component scores
#' given a multivariate functional principal component analysis (MFPCA) object.
#' Note that scores are already provided with a fitted \code{pca_mfd} object,
#' but they correspond to the
#' multivariate functional observations used to perform MFPCA,
#' With this function scores corresponding to new
#' multivariate functional data can be calculated.
#'
#' @param pca
#' A fitted multivariate functional principal
#' component analysis (MFPCA) object of class \code{pca_mfd}.
#' @param components
#' A vector of integers with the multivariate
#' functional principal components
#' corresponding to the scores to be calculated.
#' @param newdata
#' An object of class \code{mfd} containing new
#' multivariate functional data for which
#' the scores must be calculated.
#' If NULL, it is set to the data used to get \code{pca}, i.e. \code{pca$data}.
#' Default is NULL.
#'
#' @return
#' a NxM matrix, where
#'  N is the number of replications in newdata and
#'  M is the number of components given in the \code{components} argument.
#' @noRd
#' @seealso \code{\link{get_pre_scores}}
#'
get_scores <- function(pca, components, newdata = NULL) {

  inprods <- get_pre_scores(pca, components, newdata)
  apply(inprods, 1:2, sum)

}

#' Calculate the inner products to be summed to get scores given a
#' fitted \code{pca_mfd} object
#'
#' This function calculates inner products needed to calculate
#' the multivariate functional principal component scores
#' given a multivariate functional principal component analysis (MFPCA) object.
#' See details.
#' This function is called by \code{get_scores}.
#'
#' @param pca
#' A fitted MFPCA object of class \code{pca_mfd}.
#' @param components
#' A vector of integers with the components
#' for which to calculate the inner products.
#' @param newdata
#' An object of class \code{mfd} containing
#' new multivariate functional data for which
#' the inner products needed to calculate the scores must be calculated.
#' If NULL, it is set to the data used to get \code{pca}, i.e. \code{pca$data}.
#' Default is NULL.
#'
#' @return
#' A three-dimensional array:
#'
#' * the first dimension is the number of replications in \code{newdata}.
#' * the second dimension is the number of components
#' given in the \code{components} argument.
#' * the third dimension is the number of functional variables.
#'
#' @details
#' The MFPCA scores for a multivariate functional data observation
#' X(t)=(X_1(t),\dots,X_p(t)) in \code{newdata}
#' are calculated as \eqn{<X_1,\psi_1>+\dots+<X_p,\psi_p>}, where
#' <.,.> denotes the inner product in L^2 and
#' \eqn{\psi_m(t)=(\psi_m1(t),\dots,\psi_pm(t))}
#' is the vector of functions of the m-th
#' harmonics/multivariate functional principal component in \code{pca} object.
#' This function calculates the
#' individual inner products \eqn{<X_j,\psi_j>},
#' for all replications in \code{newdata}
#' and all components.
#'
#' @noRd
#' @seealso \code{\link{get_scores}}
#'
get_pre_scores <- function(pca, components, newdata = NULL) {

  fd_std <- scale_mfd(pca$data, scale = pca$scale)

  if (is.null(newdata)) {
    inprods <- pca$scores[, components, , drop = FALSE]
  } else {
    center <- attr(fd_std, "scaled:center")
    scale <- if (pca$scale) attr(fd_std, "scaled:scale") else FALSE
    fd_std <- scale_mfd(newdata, center, scale)
    inprods <- inprod_mfd(fd_std, pca$harmonics[components])
  }

  inprods

}

#' Get the projected data onto a
#' multivariate functional principal component subspace.
#'
#' Get the projection of multivariate functional data
#' onto the multivariate functional principal component subspace
#' defined by a subset of functional principal components.
#' The obtained data can be considered a prediction or
#' a finite-dimensional approximation of the original data,
#' with the dimension being equal to the number of selected components.
#'
#' @param pca
#' A fitted MFPCA object of class \code{pca_mfd}.
#' @param components
#' A vector of integers with the components over which
#' to project the data.
#' @param newdata
#' An object of class \code{mfd} containing
#' new multivariate functional data to be projected onto the desired subspace.
#' If NULL, it is set to the data used to get \code{pca}, i.e. \code{pca$data}.
#' Default is NULL.
#'
#' @return
#' An object of class \code{mfd} with the same variables
#' and observations as in \code{newdata},
#' containing the projection of \code{newdata}
#' onto the the multivariate functional principal component subspace
#' defined by the \code{components} argument.
#'
#' @details
#' In case you want to calculate the projection of the data
#' directly from scores already calculated,
#' see \code{\link{get_fit_pca_given_scores}}.
#'
#' @noRd
#'
#' @seealso \code{\link{get_fit_pca_given_scores}}
#'
get_fit_pca <- function(pca, components, newdata = NULL) {

  scores <- get_scores(pca, components, newdata)
  pca_coefs <- pca$harmonics$coefs[, components, , drop = FALSE]

  get_fit_pca_given_scores(scores, pca$harmonics[components])

}

#' Get the projected data onto the
#' multivariate functional principal component subspace,
#' given scores and harmonics.
#'
#' Get the projection of multivariate functional data
#' onto the multivariate functional principal component subspace
#' defined by a subset of functional principal components.
#' If you want to calculate projected data directly from original data,
#' see \code{\link{get_fit_pca}}.
#' Here you must provide the matrix of scores and the harmonics.
#' This happens for example when the scores are be predicted directly
#' with a regression model.
#' The obtained data can be considered a prediction
#' or a finite-dimensional approximation of the original data,
#' with the dimension being equal to the number of selected components.
#'
#' @param scores
#' A NxM matrix containing the values of the scores, where
#' N is the number of observations for which
#' we want to calculate the projection,
#' M is the number of selected multivariate functional principal components.
#' @param harmonics
#' An object of class \code{mfd}
#' containing the multivariate functional principal components.
#' Its M functional replications must correspond each to the
#' corresponding column of \code{scores}.
#'
#' @return
#' An object of class \code{mfd} with the same functional variables
#' as in \code{harmonics}
#' and the same observations as the rows in \code{scores},
#' containing the projection of \code{newdata} onto the
#' multivariate functional principal component subspace
#' defined by the components given as the argument \code{harmonics}.
#' @noRd
#'
#' @seealso \code{\link{get_fit_pca}}
#'
get_fit_pca_given_scores <- function(scores, harmonics) {

  basis <- harmonics$basis
  nbasis <- basis$nbasis
  basisnames <- basis$names
  variables <- harmonics$fdnames[[3]]
  nvar <- length(variables)
  obs <- rownames(scores)
  nobs <- nrow(scores)

  fit_coefs <- array(NA, dim = c(nbasis, nobs, nvar))
  dimnames(fit_coefs) <- list(basisnames, obs, variables)
  fit_coefs[] <- apply(harmonics$coefs, 3, function(x) x %*% t(scores))

  fdnames <- list(harmonics$fdnames[[1]],
                  obs, variables)

  mfd(fit_coefs, basis, fdnames)
}

#' Calculate the Hotelling's T^2 statistics of multivariate functional data
#'
#' Calculate the Hotelling's T^2 statistics of
#' multivariate functional data projected
#' onto a multivariate functional principal component subspace.
#'
#' @param pca
#' A fitted MFPCA object of class \code{pca_mfd}.
#' @param components
#' A vector of integers with the components
#' over which to project the data.
#' @param newdata
#' An object of class \code{mfd} containing
#' new multivariate functional data to be projected onto the desired subspace.
#' If NULL, it is set to the data used to get \code{pca},
#' i.e. \code{pca$data}.
#' Default is NULL.
#'
#' @return
#' A \code{data.frame} with as many rows as the number of
#' functional replications in \code{newdata}.
#' It has one \code{T2} column containing the Hotelling T^2
#' statistic calculated for all observations, as well as
#' one column per each functional variable, containing its
#' contribution to the T^2 statistic.
#' See Capezza et al. (2020) for definition of contributions.
#' @noRd
#'
#' @references
#' Capezza C, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
#' Control charts for monitoring ship operating conditions and CO2
#' emissions based on scalar-on-function regression.
#' \emph{Applied Stochastic Models in Business and Industry},
#' 36(3):477--500. <doi:10.1002/asmb.2507>
#'
get_T2 <- function(pca, components, newdata = NULL) {

  inprods <- get_pre_scores(pca, components, newdata)
  scores <- apply(inprods, 1:2, sum)
  values <- pca$values[components]
  T2 <- colSums(t(scores^2) / values)
  variables <- dimnames(inprods)[[3]]
  obs <- if (is.null(newdata)) pca$data$fdnames[[2]] else newdata$fdnames[[2]]
  contribution <- sapply(variables, function(variable) {
    rowSums(t(t(inprods[, , variable] * scores) / values))
  })
  contribution <- matrix(contribution,
                         nrow = length(obs),
                         ncol = length(variables))
  rownames(contribution) <- obs
  contribution <- data.frame(contribution)
  names_contribution <- paste0("contribution_T2_", variables)

  out <- data.frame(T2 = T2, contribution)
  names(out) <- c("T2", names_contribution)
  out
}


#' Calculate the Squared Prediction Error statistics of
#' multivariate functional data
#'
#' Calculate the Squared Prediction Error (SPE) statistics
#' of multivariate functional data projected
#' onto a multivariate functional principal component subspace,
#' i.e. the squared norm of the difference between the original data and
#' their projection onto the subspace.
#'
#' @param pca
#' A fitted MFPCA object of class \code{pca_mfd}.
#' @param components
#' A vector of integers with the components over
#' which to project the data.
#' @param newdata
#' An object of class \code{mfd} containing
#' new multivariate functional data to be projected onto the desired subspace.
#' If NULL, it is set to the data used to get \code{pca}, i.e. \code{pca$data}.
#' Default is NULL.
#'
#' @return
#' A \code{data.frame} with as many rows as the number of
#' functional replications in \code{newdata}.
#' It has one \code{spe} column containing the SPE statistic calculated
#' for all observations,
#' as well as one column per each functional variable, containing its
#' contribution to the SPE statistic.
#' See Capezza et al. (2020) for definition of contributions.
#' @noRd
#'
#' @references
#' Capezza C, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
#' Control charts for
#' monitoring ship operating conditions and CO2 emissions based on
#' scalar-on-function regression.
#' \emph{Applied Stochastic Models in Business and Industry},
#' 36(3):477--500. <doi:10.1002/asmb.2507>
#'
get_spe <- function(pca, components, newdata = NULL) {

  fd_std <- scale_mfd(pca$data, scale = pca$scale)

  if (!is.null(newdata)) {
    center <- attr(fd_std, "scaled:center")
    scale <- if (pca$scale) attr(fd_std, "scaled:scale") else FALSE
    fd_std <- scale_mfd(newdata, center, scale)
  }

  obs <- fd_std$fdnames[[2]]
  vars <- fd_std$fdnames[[3]]

  fit <- get_fit_pca(pca, components, newdata)
  res_fd <- minus.fd(fd_std, fit)
  res_fd <- mfd(res_fd$coefs, res_fd$basis, res_fd$fdnames)

  cont_spe <- inprod_mfd(res_fd)
  cont_spe <- apply(cont_spe, 3, diag)
  cont_spe <- matrix(cont_spe,
                     nrow = length(obs),
                     ncol = length(vars))
  rownames(cont_spe) <- obs
  colnames_cont_spe <- paste0("contribution_spe_", vars)
  spe <- rowSums(cont_spe)

  out <- data.frame(spe = spe, cont_spe)
  colnames(out) <- c("spe", colnames_cont_spe)
  out

}

