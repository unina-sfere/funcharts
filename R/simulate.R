#' Simulate example data for funcharts
#'
#' Function used to simulate a data set to illustrate the use of \code{funcharts}.
#' It creates a data set with three functional covariates,
#' a functional response generated as a function of the three functional covariates
#' through a function-on-function linear model,
#' and a scalar response generated as a function of the three functional covariates
#' through a scalar-on-function linear model.
#' This function covers the simulation study in Centofanti et al. (2020)
#' for the function-on-function case and also simulates data in a similar way
#' for the scalar response case.
#' In the default case, the function generates in-control data.
#' Additional arguments can be used to generate additional data that are out of control,
#' with mean shifts according to the scenarios proposed by Centofanti et al. (2020).
#' Each simulated observation of a functional variable consists of
#' a vector of 150 discrete points, equally spaced between 0 and 1, generated with noise.
#' @param nobs
#' The number of observation to simulate
#' @param R2
#' The desired coefficient of determination in the regression.
#' In both the scalar and functional response cases, only three values are allowed,
#' i.e. 0.97, 0.86, 0.74.
#' Default is 0.97.
#' @param seed Set seed for reproducibility. Default is 0.
#' @param shift_type_y
#' The shift type for the functional response.
#' There are five possibilities: "0" if there is no shift,
#' "A", "B", "C" or "D" for the corresponding shift types
#' shown in Centofanti et al. (2020).
#' Default is "0".
#' @param shift_type_x1
#' #' The shift type for the first functional covariate.
#' There are five possibilities: "0" if there is no shift,
#' "A", "B", "C" or "D" for the corresponding shift types
#' shown in Centofanti et al. (2020).
#' Default is "0".
#' @param shift_type_x2
#' #' The shift type for the second functional covariate.
#' #' There are five possibilities: "0" if there is no shift,
#' "A", "B", "C" or "D" for the corresponding shift types
#' shown in Centofanti et al. (2020).
#' Default is "0".
#' @param shift_type_x3
#' #' The shift type for the third functional covariate.
#' #' There are five possibilities: "0" if there is no shift,
#' "A", "B", "C" or "D" for the corresponding shift types
#' shown in Centofanti et al. (2020).
#' Default is "0".
#' @param d_y
#' A number indicating the severity of the shift type for the functional response.
#' Default is 0.
#' @param d_x1
#' A number indicating the severity of the shift type for the first functional covariate.
#' Default is 0.
#' @param d_x2
#' A number indicating the severity of the shift type for the second functional covariate.
#' Default is 0.
#' @param d_x3
#' A number indicating the severity of the shift type for the third functional covariate.
#' Default is 0.
#' @param d_y_scalar
#' A number indicating the severity of the shift type for the scalar response.
#' Default is 0.
#' @param save_beta
#' If TRUE, the true regression coefficients of both the function-on-function
#' and the scalar-on-function models are saved.
#' Default is FALSE.
#'
#' @return
#' A list with the following elements:
#'
#' * \code{X1} is a \code{nobs}x150 matrix with the simulated observations of the first functional covariate
#'
#' * \code{X2} is a \code{nobs}x150 matrix with the simulated observations of the second functional covariate
#'
#' * \code{X3} is a \code{nobs}x150 matrix with the simulated observations of the third functional covariate
#'
#' * \code{Y} is a \code{nobs}x150 matrix with the simulated observations of the functional response
#'
#' * \code{y_scalar} is a vector of length 150 with the simulated observations of the scalar response
#'
#' * \code{beta_fof}, if \code{save_beta = TRUE}, is a list of three 500x500 matrices
#' with the discretized functional coefficients of the funciton-on-function regression
#'
#' * \code{beta_sof}, if \code{save_beta = TRUE}, is a list of three vectors of length 500
#' with the discretized functional coefficients of the scalar-on-function regression
#'
#'
#' @export
#'
#' @references
#'
#' Centofanti F, Lepore A, Menafoglio A, Palumbo B, Vantini S. (2020)
#' Functional Regression Control Chart.
#' \emph{Technometrics}. <doi:10.1080/00401706.2020.1753581>
simulate_mfd <- function(nobs = 1000,
                         R2 = 0.97,
                         seed = 0,
                         shift_type_y = "0",
                         shift_type_x1 = "0",
                         shift_type_x2 = "0",
                         shift_type_x3 = "0",
                         d_y = 0,
                         d_x1 = 0,
                         d_x2 = 0,
                         d_x3 = 0,
                         d_y_scalar = 0,
                         save_beta = FALSE) {

  if (!(R2 %in% c(0.97, 0.86, 0.74))) {
    stop("only 0.97, 0.86, 0.74 values allowed for R2.")
  }

  if (!(toupper(shift_type_y) %in% c("A", "B", "C", "D", "0"))) {
    stop("only shift types A, B, C, D and 0 are allowed.")
  }
  if ((toupper(shift_type_y) %in% c("A", "B", "C", "D")) & d_y == 0) {
    stop("If there is a shift, the corresponding severity d must be greater than zero.")
  }
  if (toupper(shift_type_y) == "0" & d_y != 0) {
    stop("If there is no shift, the corresponding severity d must be zero.")
  }

  if (!(toupper(shift_type_x1) %in% c("A", "B", "C", "D", "0"))) {
    stop("only shift types A, B, C, D and 0 are allowed.")
  }
  if (toupper(shift_type_x1) %in% c("A", "B", "C", "D") & d_x1 == 0) {
    stop("If there is a shift, the corresponding severity d must be greater than zero.")
  }
  if (toupper(shift_type_x1) == "0" & d_x1 != 0) {
    stop("If there is no shift, the corresponding severity d must be zero.")
  }

  if (!(toupper(shift_type_x2) %in% c("A", "B", "C", "D", "0"))) {
    stop("only shift types A, B, C, D and 0 are allowed.")
  }
  if (toupper(shift_type_x2) %in% c("A", "B", "C", "D") & d_x2 == 0) {
    stop("If there is a shift, the corresponding severity d must be greater than zero.")
  }
  if (toupper(shift_type_x2) == "0" & d_x2 != 0) {
    stop("If there is no shift, the corresponding severity d must be zero.")
  }

  if (!(toupper(shift_type_x3) %in% c("A", "B", "C", "D", "0"))) {
    stop("only shift types A, B, C, D and 0 are allowed.")
  }
  if (toupper(shift_type_x3) %in% c("A", "B", "C", "D") & d_x3 == 0) {
    stop("If there is a shift, the corresponding severity d must be greater than zero.")
  }
  if (toupper(shift_type_x3) == "0" & d_x3 != 0) {
    stop("If there is no shift, the corresponding severity d must be zero.")
  }

  e <- cov_str$e
  eigy <- cov_str$eigy
  P <- nrow(eigy$vectors)
  x_seq <- seq(0, 1, l = P)
  w <- 1 / P
  n_comp_x <- 50
  meig1 <- e$vectors[1:500, 1:n_comp_x] / sqrt(w)
  meig2 <- e$vectors[501:1000, 1:n_comp_x] / sqrt(w)
  meig3 <- e$vectors[1001:1500, 1:n_comp_x] / sqrt(w)
  meigenvalues <- e$values[1:n_comp_x] * w

  set.seed(seed)

  b_perfect <- sqrt(eigy$values / meigenvalues[1:10])
  b_perfect_sof <- rep(1 / sqrt(sum(meigenvalues[1:10])), 10)

  if (R2 == 0.97) a <- 0.01
  if (R2 == 0.86) a <- 0.05
  if (R2 == 0.74) a <- 0.1

  b <- pmax(b_perfect - a, 0)
  var_yexp <- rowSums(t(t(eigy$vectors^2) * (b^2 * meigenvalues[1:10])))
  vary <- rowSums(t(t(eigy$vectors^2) * eigy$values))
  R2_check <- sum(var_yexp / vary) * w
  B <- rbind(diag(b), matrix(0, nrow = 40, ncol = 10))
  eigeps <- diag(diag(eigy$values) - t(B) %*% diag(meigenvalues) %*% B)

  csi_X <- rnorm(n = length(meigenvalues) * nobs,
                 mean = 0,
                 sd = sqrt(rep(meigenvalues, each = nobs)))
  csi_X <- matrix(csi_X, nrow = nobs)
  csi_eps <- rnorm(n = length(eigeps) * nobs, mean = 0,
                   sd = sqrt(rep(eigeps, each = nobs)))
  csi_eps <- matrix(csi_eps, nrow = nobs)
  csi_Y <- csi_X %*% B + csi_eps

  fun <- function(x, par, r, m, s) {
    a <- par[1]
    b <- par[2]
    c <- par[3]
    r <- par[4]
    dnorm_vec <- Vectorize(dnorm, "x")
    function(x) {
      norm_dens <- dnorm_vec(x, mean = m, sd = s)
      sum_term <- if (identical(s, 0)) 0 else colSums(norm_dens)
      a * x^2 + b * x + c + r * sum_term
    }
  }

  par1 <- c(-20, 20, 20, 0.05)
  par2 <- c(0, 4, 0, 0)
  par3 <- c(-10, 14, -8, 0.05)
  pary <- c(0, 30, 0, 0)

  if (shift_type_x1 == "A") par1 <- par1 + d_x1 * c(1, 0, 0, 0)
  if (shift_type_x1 == "B") par1 <- par1 + d_x1 * c(0, 1, 0, 0)
  if (shift_type_x1 == "C") par1 <- par1 + d_x1 * c(0, 0, 1, 0)
  if (shift_type_x1 == "D") par1 <- par1 + d_x1 * c(1, 1, 0, 0)

  if (shift_type_x2 == "A") par2 <- par2 + d_x2 * c(1, 0, 0, 0)
  if (shift_type_x2 == "B") par2 <- par2 + d_x2 * c(0, 1, 0, 0)
  if (shift_type_x2 == "C") par2 <- par2 + d_x2 * c(0, 0, 1, 0)
  if (shift_type_x2 == "D") par2 <- par2 + d_x2 * c(1, 1, 0, 0)

  if (shift_type_x3 == "A") par3 <- par3 + d_x3 * c(1, 0, 0, 0)
  if (shift_type_x3 == "B") par3 <- par3 + d_x3 * c(0, 1, 0, 0)
  if (shift_type_x3 == "C") par3 <- par3 + d_x3 * c(0, 0, 1, 0)
  if (shift_type_x3 == "D") par3 <- par3 + d_x3 * c(1, 1, 0, 0)

  if (shift_type_y == "A") pary <- pary + d_y * c(1, 0, 0, 0)
  if (shift_type_y == "B") pary <- pary + d_y * c(0, 1, 0, 0)
  if (shift_type_y == "C") pary <- pary + d_y * c(0, 0, 1, 0)
  if (shift_type_y == "D") pary <- pary + d_y * c(1, 1, 0, 0)

  f1_m <- fun(par = par1,
              m = c(0.075, 0.100, 0.250, 0.350, 0.500, 0.650, 0.850, 0.900, 0.950),
              s = c(0.050, 0.030, 0.050, 0.050, 0.100, 0.050, 0.100, 0.040, 0.035))
  f2_m <- fun(par = par2, m = 0, s = 0)
  f3_m <- fun(par = par3,
              m = c(0.075, 0.100, 0.150, 0.225, 0.400, 0.525, 0.550,
                    0.600, 0.625, 0.650, 0.850, 0.900, 0.925),
              s = c(0.050, 0.060, 0.050, 0.040, 0.050, 0.035, 0.045,
                    0.045, 0.040, 0.030, 0.015, 0.010, 0.015))
  fy_m <- fun(par = pary, m = 0, s = 0)

  f1_v <- fun(par = c(0, 0, 1, 0.1),
              m = c(0.075, 0.100, 0.125, 0.150, 0.400, 0.650, 0.850, 0.900, 0.925),
              s = c(0.050, 0.060, 0.075, 0.075, 0.075, 0.075, 0.045, 0.045, 0.040))
  f2_v <- fun(par = c(0, 0.02, 1, 0),
              m = 0,
              s = 0)
  f3_v <- fun(par = c(-40, 150, 30, 2),
              m = c(0.075, 0.100, 0.150, 0.225, 0.400, 0.525, 0.550,
                    0.600, 0.625, 0.650, 0.850, 0.900, 0.925),
              s = c(0.050, 0.060, 0.050, 0.040, 0.050, 0.035, 0.045,
                    0.045, 0.040, 0.030, 0.015, 0.010, 0.015))
  fy_v <- fun(par = c(0, 8, 1, 0),
              m = 0,
              s = 0)

  ngrid <- 150
  grid <- seq(0, 1, length.out = ngrid)

  meig1_interpolate <- apply(t(meig1), 1, function(y) approx(x_seq, y, grid)$y)
  X1_scaled <- csi_X %*% t(meig1_interpolate)
  X1 <- t(f1_m(grid) + t(X1_scaled) * sqrt(f1_v(grid)))

  meig2_interpolate <- apply(t(meig2), 1, function(y) approx(x_seq, y, grid)$y)
  X2_scaled <- csi_X %*% t(meig2_interpolate)
  X2 <- t(f2_m(grid) + t(X2_scaled) * sqrt(f2_v(grid)))

  meig3_interpolate <- apply(t(meig3), 1, function(y) approx(x_seq, y, grid)$y)
  X3_scaled <- csi_X %*% t(meig3_interpolate)
  X3 <- t(f3_m(grid) + t(X3_scaled) * sqrt(f3_v(grid)))

  eigy_interpolate <- apply(t(eigy$vectors), 1, function(y) approx(x_seq, y, grid)$y)
  Y_scaled <- csi_Y %*% t(eigy_interpolate)
  Y <- t(fy_m(grid) + t(Y_scaled) * sqrt(fy_v(grid)))

  X1 <- X1 + rnorm(length(X1), sd = 0.3)
  X2 <- X2 + rnorm(length(X1), sd = 0.05)
  X3 <- X3 + rnorm(length(X3), sd = 0.3)
  Y <- Y + rnorm(length(X1), sd = 0.3)

  vary_scalar <- 1
  b_sof <- rep(sqrt(R2 * vary_scalar / sum(meigenvalues[1:10])), 10)
  R2_sof_check <- sum(b_sof^2 * meigenvalues[1:10])

  var_eps_scalar <- vary_scalar - sum(b_sof^2 * meigenvalues[1:10])
  eps_scalar <- rnorm(n = nobs, mean = 0, sd = sqrt(var_eps_scalar))
  y_scalar <- as.numeric(d_y_scalar + csi_X[, 1:10] %*% b_sof + eps_scalar)

  beta_fof <- beta_sof <- NULL
  if (save_beta) {
    beta_fof <- 0
    for (kk in 1:length(diag(B))) {
      beta_fof <- beta_fof + diag(B)[kk] * outer(e$vectors[, 1], eigy$vectors[, 1])
    }
    beta_fof <- list(beta_fof[1:500, ], beta_fof[500 + 1:500, ], beta_fof[1000 + 1:500, ])

    beta_sof <- as.numeric(e$vectors[, 1:10] %*% b_sof)
    beta_sof <- list(beta_sof[1:500], beta_sof[500 + 1:500], beta_sof[1000 + 1:500])
  }

  out <- list(X1 = X1,
       X2 = X2,
       X3 = X3,
       Y = Y,
       y_scalar = y_scalar)
  if (save_beta) {
    out$beta_fof = beta_fof
    out$beta_sof = beta_sof
  }
  out

}

