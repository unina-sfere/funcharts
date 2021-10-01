#' Generate covariance functions
#'
#' This function generates the covariance structure from which data are simulated
#' and precomputes the eigendecomposition of the covariance operators.
#' The output of this function is saved in the file \code{cov_str.Rdata},
#' so that \code{simulate_mfd()} can load the data
#' instead of computing the covariance,
#' reducing the time needed to create the simulated data.
#'
#' @return
#' A list with two objects:
#' \code{e} contains the eigendecomposition of the covariance operator
#' of the multivariate functional covariates.
#' \code{eigy} contains the eigendecomposition of the covariance operator
#' of the functional response.
#'
#' @seealso \code{\link{simulate_mfd}},
#'
generate_cov_str <- function() {

  p <- 3
  P <- 500
  n_comp_x <- 50
  n_comp_y <- 10
  x_seq <- seq(0, 1, l = P)

  get_mat <- function(cov) {
    P <- length(cov)
    covmat <- matrix(0, nrow = P, ncol = P)
    for (ii in 1:P) {
      covmat[ii, ii:P] <- cov[1:(P - ii + 1)]
      covmat[P - ii + 1, 1:(P - ii + 1)] <- rev(cov[1:(P - ii + 1)])
    }
    covmat
  }

  cov1 <- besselJ(x_seq * 4, 0)
  cov2 <- exp(-x_seq^2)
  cov3 <- exp(-sqrt(x_seq))
  covy <- cov1

  covmat1 <- get_mat(cov1)
  covmat2 <- get_mat(cov2)
  covmat3 <- get_mat(cov3)
  covmaty <- covmat1

  covmatList <- list(covmat1, covmat2, covmat3)

  eig1 <- RSpectra::eigs_sym(covmat1, n_comp_x + 10)
  eig2 <- RSpectra::eigs_sym(covmat2, n_comp_x + 10)
  eig3 <- RSpectra::eigs_sym(covmat3, n_comp_x + 10)

  w <- 1 / P
  eig1$vectors <- eig1$vectors / sqrt(w)
  eig2$vectors <- eig2$vectors / sqrt(w)
  eig3$vectors <- eig3$vectors / sqrt(w)
  eig1$values <- eig1$values * w
  eig2$values <- eig2$values * w
  eig3$values <- eig3$values * w
  eigy <- eig1

  eigList <- list(eig1, eig2, eig3)

  eigenvalues <- rowMeans(cbind(eig1$values, eig2$values, eig3$values))

  corr_mat <- vector("list", 3)
  for (ii in 1:3) {
    corr_mat[[ii]] <- vector("list", 3)
  }

  for (ii in 1:3) {
    for (jj in ii:3) {
      if (jj == ii) {
        corr_mat[[ii]][[jj]] <- covmatList[[ii]]
      } else {
        sum <- 0
        for (kk in 1:n_comp_x) {
          sum <- sum +
            eigenvalues[kk] *
            outer(eigList[[ii]]$vectors[, kk],
                  eigList[[jj]]$vectors[, kk]) / 3
        }
        corr_mat[[ii]][[jj]] <- sum
        corr_mat[[jj]][[ii]] <- t(sum)
      }
    }
    corr_mat[[ii]] <- do.call("cbind", corr_mat[[ii]])
  }
  corr_mat <- do.call("rbind", corr_mat)
  e <- RSpectra::eigs_sym(corr_mat, n_comp_x + 10)

  eigy$values <- eigy$values[1:10]
  eigy$vectors <- eigy$vectors[, 1:10]
  list(e = e, eigy = eigy)
}
