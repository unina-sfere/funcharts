#' Define a Multivariate Functional Data Object
#'
#' This is the constructor function for objects of the mfd class.
#' It is a wrapper to \code{fda::\link[fda]{fd}},
#' but it forces the coef argument to  be
#' a three-dimensional array of coefficients even if
#' the functional data is univariate.
#' Moreover, it allows to include the original raw data from which
#' you get the smooth functional data.
#'
#' @param coef
#' A three-dimensional array of coefficients:
#'
#' * the first dimension corresponds to basis functions.
#'
#' * the second dimension corresponds to the number of
#' multivariate functional observations.
#'
#' * the third dimension corresponds to variables.
#' @param basisobj
#' A functional basis object defining the basis,
#' as provided to \code{fda::\link[fda]{fd}}, but there is no default.
#' @param fdnames
#' A list of length 3, each member being a string vector
#' containing labels for the levels of the corresponding dimension
#' of the discrete data.
#'
#' The first dimension is for a single character indicating the argument values,
#' i.e. the variable on the functional domain.
#'
#' The second is for replications, i.e. it denotes the functional observations.
#'
#' The third is for functional variables' names.
#' @param raw
#' A data frame containing the original discrete data.
#' Default is NULL, however, if provided, it must contain:
#'
#' a column (indicated by the \code{id_var} argument)
#' denoting the functional observations,
#' which must correspond to values in \code{fdnames[[2]]},
#'
#' a column named as \code{fdnames[[1]]},
#' returning the argument values of each function
#'
#' as many columns as the functional variables,
#' named as in \code{fdnames[[3]]},
#' containing the discrete functional values for each variable.
#' @param id_var
#' A single character value indicating the column
#' in the \code{raw} argument
#' containing the functional observations (as in \code{fdnames[[2]]}),
#' default is NULL.
#'
#' @return
#' A multivariate functional data object
#' (i.e., having class \code{mfd}),
#' which is a list with components named
#' \code{coefs}, \code{basis}, and \code{fdnames},
#' as for class \code{fd},
#' with possibly in addition the components \code{raw} and \code{id_var}.
#'
#' @details
#' To check that an object is of this class, use function is.mfd.
#'
#' @references
#' Ramsay, James O., and Silverman, Bernard W. (2006),
#' \emph{Functional Data Analysis}, 2nd ed., Springer, New York.
#'
#' Ramsay, James O., and Silverman, Bernard W. (2002),
#' \emph{Applied Functional Data Analysis}, Springer, New York.
#'
#' @export
#'
#' @examples
#' library(funcharts)
#' set.seed(0)
#' nobs <- 5
#' nbasis <- 4
#' nvar <- 2
#' coef <- array(rnorm(nobs * nbasis * nvar), dim = c(nbasis, nobs, nvar))
#' bs <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
#' mfdobj <- mfd(coef = coef, basisobj = bs)
#' plot_mfd(mfdobj)
#'
mfd <- function(coef, basisobj, fdnames = NULL, raw = NULL, id_var = NULL) {

  if (is.null(fdnames)) {
    fdnames <- list(
      "arg",
      paste0("rep", 1:dim(coef)[2]),
      paste0("var", 1:dim(coef)[3])
    )
  }
  if (sum(is.null(raw), is.null(id_var)) == 1) {
    stop("If either raw or id_var are not NULL, both must be not NULL")
  }
  if (!is.null(raw)) {
    if (!(fdnames[[1]] %in% names(raw))) {
      stop("fdnames[[1]] must be the name of a column of the raw data.")
    }
    if (!is.data.frame(raw)) {
      stop("raw must be a data frame.")
    }
    if (!(is.character(id_var) & length(id_var) == 1)) {
      stop ("id_var must be a single character")
    }
    if (!(class(raw[[id_var]]) %in% c("character", "factor"))) {
      stop("column id_var of raw data frame must be a character or factor.")
    }
    if (!(setequal(unique(raw[[id_var]]), fdnames[[2]]))) {
      stop(paste0("column id_var of raw data frame must contain",
                  "the same values as fdnames[[2]]."))
    }
  }

  if (!is.array(coef)) {
    stop("'coef' is not array")
  }
  coefd <- dim(coef)
  ndim <- length(coefd)
  if (ndim != 3) {
    stop("'coef' not of dimension 3")
  }

  fdobj <- fd(coef, basisobj, fdnames)
  fdobj$raw <- raw
  fdobj$id_var <- id_var
  class(fdobj) <- c("mfd", "fd")
  fdobj
}


#' Confirm Object has Class \code{mfd}
#'
#' Check that an argument is a multivariate
#' functional data object of class \code{mfd}.
#'
#' @param mfdobj An object to be checked.
#'
#' @return a logical value: TRUE if the class is correct, FALSE otherwise.
#' @export
#'
is.mfd <- function(mfdobj) if (inherits(mfdobj, "mfd")) TRUE else FALSE


#' Simulate multivariate functional data
#'
#' Simulate random coefficients and create a multivariate functional data
#' object of class `mfd`.
#'
#' @param nobs Number of functional observations to be simulated.
#' @param nbasis Number of basis functions.
#' @param nvar Number of functional covariates.
#' @param seed Set seed for reproducibility.
#'
#' @return
#' A simulated object of class `mfd`.
#' @export
#'
#' @examples
#' library(funcharts)
#' data_sim_mfd()
data_sim_mfd <- function(nobs = 5,
                         nbasis = 4,
                         nvar = 2,
                         seed = 0) {
  set.seed(seed)
  coef <- array(stats::rnorm(nobs * nbasis * nvar),
                dim = c(nbasis, nobs, nvar))
  bs <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
  mfd(coef = coef, basisobj = bs)
}



#' Extract observations and/or variables from \code{mfd} objects.
#'
#' @param mfdobj An object of class \code{mfd}.
#' @param i
#' Index specifying functional observations to extract or replace.
#' They can be numeric, character,
#' or logical vectors or empty (missing) or NULL.
#' Numeric values are coerced to integer as by as.integer
#' (and hence truncated towards zero).
#' The can also be negative integers,
#' indicating functional observations to leave out of the selection.
#' Logical vectors indicate TRUE for the observations to select.
#' Character vectors will be matched
#' to the argument \code{fdnames[[2]]} of \code{mfdobj},
#' i.e. to functional observations' names.
#' @param j
#' Index specifying functional variables to extract or replace.
#' They can be numeric, logical,
#' or character vectors or empty (missing) or NULL.
#' Numeric values are coerced to integer as by as.integer
#' (and hence truncated towards zero).
#' The can also be negative integers,
#' indicating functional variables to leave out of the selection.
#' Logical vectors indicate TRUE for the variables to select.
#' Character vectors will be matched
#' to the argument \code{fdnames[[3]]} of \code{mfdobj},
#' i.e. to functional variables' names.
#'
#' @return a \code{mfd} object with selected observations and variables.
#'
#' @details
#' This function adapts the \code{fda::"[.fd"}
#' function to be more robust and suitable
#' for the \code{mfd} class.
#' In fact, whatever the number of observations
#' or variables you want to extract,
#' it always returns a \code{mfd} object with a three-dimensional coef array.
#' In other words, it behaves as you would
#' always use the argument \code{drop=FALSE}.
#' Moreover, you can extract observations
#' and variables both by index numbers and by names,
#' as you would normally do when using
#' \code{`[`} with standard vector/matrices.
#'
#' @export
#' @examples
#' library(funcharts)
#'
#' # In the following, we extract the first one/two observations/variables
#' # to see the difference with `[.fd`.
#' mfdobj <- data_sim_mfd()
#' fdobj <- fd(mfdobj$coefs, mfdobj$basis, mfdobj$fdnames)
#'
#' # The argument `coef` in `fd` objects is converted to a matrix when possible.
#' dim(fdobj[1, 1]$coef)
#' # Not clear what is the second dimension:
#' # the number of replications or the number of variables?
#' dim(fdobj[1, 1:2]$coef)
#' dim(fdobj[1:2, 1]$coef)
#'
#' # The argument `coef` in `mfd` objects is always a three-dimensional array.
#' dim(mfdobj[1, 1]$coef)
#' dim(mfdobj[1, 1:2]$coef)
#' dim(mfdobj[1:2, 1]$coef)
#'
#' # Actually, `[.mfd` works as `[.fd` when passing also `drop = FALSE`
#' dim(fdobj[1, 1, drop = FALSE]$coef)
#' dim(fdobj[1, 1:2, drop = FALSE]$coef)
#' dim(fdobj[1:2, 1, drop = FALSE]$coef)
#'
"[.mfd" <- function(mfdobj, i = TRUE, j = TRUE) {
  if (!(is.mfd(mfdobj))) {
    stop(paste0("First argument must be a ",
                "multivariate functional data object."))
  }

  coefs <- mfdobj$coefs[, i, j, drop = FALSE]
  fdnames <- mfdobj$fdnames
  fdnames[[2]] <- dimnames(coefs)[[2]]
  fdnames[[3]] <- dimnames(coefs)[[3]]
  if (is.null(mfdobj$raw))
    raw_filtered <- id_var <- NULL else {
      raw <- mfdobj$raw
      id_var <- mfdobj$id_var
      raw_filtered <- raw[raw[[id_var]] %in% fdnames[[2]], , drop = FALSE]
    }

  mfd(coef = coefs, basisobj = mfdobj$basis, fdnames = fdnames,
      raw = raw_filtered, id_var = id_var)
}

#' Inner products of functional data contained in \code{mfd} objects.
#'
#' @param mfdobj1
#' A multivariate functional data object of class \code{mfd}.
#' @param mfdobj2
#' A multivariate functional data object of class \code{mfd}.
#' It must have the same functional variables as \code{mfdobj1}.
#' If NULL, it is equal to \code{mfdobj1}.
#'
#' @return
#' a three-dimensional array of \emph{L^2} inner products.
#' The first dimension is the number of functions in argument mfdobj1,
#' the second dimension is the same thing for argument mfdobj2,
#' the third dimension is the number of functional variables.
#' If you sum values over the third dimension,
#' you get a matrix of inner products
#' in the product Hilbert space of multivariate functional data.
#'
#' @details
#' Note that \emph{L^2} inner products are not calculated
#' for couples of functional data
#' from different functional variables.
#' This function is needed to calculate the
#' inner product in the product Hilbert space
#' in the case of multivariate functional data,
#' which for each observation is the sum of the \emph{L^2}
#' inner products obtained for each functional variable.
#'
#' @export
#' @examples
#' library(funcharts)
#' mfdobj1 <- data_sim_mfd(seed = 123)
#' mfdobj2 <- data_sim_mfd(seed = 987)
#' inprod_mfd(mfdobj1)
#' inprod_mfd(mfdobj1, mfdobj2)
inprod_mfd <- function(mfdobj1, mfdobj2 = NULL) {

  if (!(is.mfd(mfdobj1))) {
    stop("First argument is not a multivariate functional data object.")
  }
  if (!is.null(mfdobj2)) {
    if (!(is.mfd(mfdobj2))) {
      stop("Second argument is not a multivariate functional data object.")
    }
    if (length(mfdobj1$fdnames[[3]]) != length(mfdobj2$fdnames[[3]])) {
      stop(paste0("mfdobj1 and mfdobj2 do not have the ",
                  "same number of functional variables."))
    }
  } else mfdobj2 <- mfdobj1

  variables <- mfdobj1$fdnames[[3]]
  n_var <- length(variables)
  ids1 <- mfdobj1$fdnames[[2]]
  ids2 <- mfdobj2$fdnames[[2]]
  n_obs1 <- length(ids1)
  n_obs2 <- length(ids2)
  bs1 <- mfdobj1$basis
  bs2 <- mfdobj2$basis

  inprods <- array(NA, dim = c(n_obs1, n_obs2, n_var),
                   dimnames = list(ids1, ids2, variables))
  inprods[] <- sapply(1:n_var, function(jj) {
    mfdobj1_jj <- fd(mfdobj1$coefs[, , jj], bs1)
    mfdobj2_jj <- fd(mfdobj2$coefs[, , jj], bs2)
    inprod(mfdobj1_jj, mfdobj2_jj)
  })

  inprods

}

#' Norm of Multivariate Functional Data
#'
#' Norm of multivariate functional data contained
#' in a \code{mfd} object.
#'
#' @param mfdobj A multivariate functional data object of class \code{mfd}.
#'
#' @return
#' A vector of length equal to the number of replications
#' in \code{mfdobj},
#' containing the norm of each multivariate functional observation
#' in the product Hilbert space,
#' i.e. the sum of \emph{L^2} norms for each functional variable.
#' @export
#' @examples
#' library(funcharts)
#' mfdobj <- data_sim_mfd()
#' norm.mfd(mfdobj)
#'
norm.mfd <- function(mfdobj) {

  if (!(is.mfd(mfdobj))) {
    stop("Input is not a multivariate functional data object.")
  }

  inprods <- inprod_mfd_diag(mfdobj)
  sqrt(rowSums(inprods))

}



#' Inner product of two multivariate functional data objects, for each observation
#'
#' @param mfdobj1 A multivariate functional data object of class \code{mfd}.
#' @param mfdobj2 A multivariate functional data object of class \code{mfd},
#' with the same number of functional variables and observations as \code{mfdobj1}.
#' If NULL, then \code{mfdobj2=mfdobj1}. Default is NULL.
#'
#' @return It calculates the inner product of two multivariate functional data objects.
#' The main function \code{inprod} of the package \code{fda} calculates inner products among
#' all possible couples of observations.
#' This means that, if \code{mfdobj1} has \code{n1} observations
#' and \code{mfdobj2} has \code{n2} observations,
#' then for each variable \code{n1 X n2} inner products are calculated.
#' However, often one is interested only in calculating the \code{n} inner products
#' between the \code{n} observations of \code{mfdobj1} and the corresponding \code{n}
#' observations of \code{mfdobj2}. This function provides this "diagonal" inner products only,
#' saving a lot of computation with respect to using \code{fda::inprod} and then extracting the
#' diagonal elements.
#' Note that the code of this function calls a modified version of \code{fda::inprod()}.
#' @export
#'
#' @examples
#' mfdobj <- data_sim_mfd()
#' inprod_mfd_diag(mfdobj)
#'
inprod_mfd_diag <- function(mfdobj1, mfdobj2 = NULL) {

  if (!is.mfd(mfdobj1)) {
    stop("Only mfd class allowed for mfdobj1 input")
  }

  if (!(is.fd(mfdobj2) | is.null(mfdobj2))) {
    stop("mfdobj2 input must be of class mfd, or NULL")
  }

  nvar1 <- dim(mfdobj1$coefs)[3]
  nobs1 <- dim(mfdobj1$coefs)[2]

  if (is.fd(mfdobj2)) {
    nvar2 <- dim(mfdobj2$coefs)[3]
    if (nvar1 != nvar2) {
      stop("mfdobj1 and mfdobj2 must have the same number of variables")
    }
    nobs2 <- dim(mfdobj2$coefs)[2]
    if (nobs1 != nobs2) {
      stop("mfdobj1 and mfdobj2 must have the same number of observations")
    }
  }

  inprods <- sapply(1:nvar1, function(jj) {
    fdobj1_jj <- fd(matrix(mfdobj1$coefs[, , jj],
                           nrow = dim(mfdobj1$coefs)[1],
                           ncol = dim(mfdobj1$coefs)[2]),
                    mfdobj1$basis)
    if (is.null(mfdobj2)) {
      out <- inprod_fd_diag_single(fdobj1_jj)
    } else {
      fdobj2_jj <- fd(matrix(mfdobj2$coefs[, , jj],
                             nrow = dim(mfdobj2$coefs)[1],
                             ncol = dim(mfdobj2$coefs)[2]),
                      mfdobj2$basis)
      out <- inprod_fd_diag_single(fdobj1_jj, fdobj2_jj)
    }
    out
  })

  if (nobs1 == 1) inprods <- matrix(inprods, nrow = 1)
  inprods

}

#' @noRd
#'
inprod_fd_diag_single <- function(fdobj1, fdobj2 = NULL) {

  if (!is.fd(fdobj1)) {
    stop("Only fd class allowed for fdobj1 input")
  }

  nrep1 <- dim(fdobj1$coefs)[2]
  coef1 <- fdobj1$coefs
  basisobj1 <- fdobj1$basis
  type1 <- basisobj1$type
  range1 <- basisobj1$rangeval

  if (!(is.fd(fdobj2) | is.null(fdobj2))) {
    stop("fdobj2 input must be of class fd, or NULL")
  }

  if (is.fd(fdobj2)) {
    type_calculated <- "inner_product"
    nrep2 <- dim(fdobj2$coefs)[2]
    if (nrep1 != nrep2) {
      stop("fdobj1 and fdobj2 must have the same number of observations")
    }
    coef2 <- fdobj2$coefs
    basisobj2 <- fdobj2$basis
    type2 <- basisobj2$type
    range2 <- basisobj2$rangeval
    if (!all(range1 == range2)) {
      stop("fdobj1 and fdobj2 must have the same domain")
    }
  }

  if (is.null(fdobj2)) {
    type_calculated <- "norm"
    nrep2 <- nrep1
    coef2 <- coef1
    basisobj2 <- basisobj1
    type2 <- type1
    range2 <- range1
  }

  iter <- 0
  rngvec <- range1
  if ((all(c(coef1) == 0) || all(c(coef2) == 0))) {
    return(numeric(nrep1))
  }

  JMAX <- 25
  JMIN <- 5
  EPS <- 1e-6

  inprodmat <- numeric(nrep1)
  nrng <- length(rngvec)
  for (irng in 2:nrng) {
    rngi <- c(rngvec[irng - 1], rngvec[irng])
    if (irng > 2)
      rngi[1] <- rngi[1] + 1e-10
    if (irng < nrng)
      rngi[2] <- rngi[2] - 1e-10
    iter <- 1
    width <- rngi[2] - rngi[1]
    JMAXP <- JMAX + 1
    h <- rep(1, JMAXP)
    h[2] <- 0.25
    s <- vector(mode = "list", length = JMAXP)
    fx1 <- eval.fd(rngi, fdobj1)
    if (type_calculated == "inner_product") {
      fx2 <- eval.fd(rngi, fdobj2)
      s[[1]] <- width * colSums(fx1 * fx2) / 2
    }
    if (type_calculated == "norm") {
      s[[1]] <- width * colSums(fx1^2) / 2
    }
    tnm <- 0.5
    for (iter in 2:JMAX) {
      tnm <- tnm * 2
      if (iter == 2) {
        x <- mean(rngi)
      } else {
        del <- width/tnm
        x <- seq(rngi[1] + del/2, rngi[2] - del/2, del)
      }
      fx1 <- eval.fd(x, fdobj1)

      if (type_calculated == "inner_product") {
        fx2 <- eval.fd(x, fdobj2)
        chs <- width * colSums(fx1 * fx2) / tnm
      }
      if (type_calculated == "norm") {
        chs <- width * colSums(fx1^2) / tnm
      }

      s[[iter]] <- (s[[iter - 1]] + chs) / 2

      if (iter >= 5) {
        ind <- (iter - 4):iter
        ya <- s[ind]
        xa <- h[ind]
        absxa <- abs(xa)
        absxamin <- min(absxa)
        ns <- min((1:length(absxa))[absxa == absxamin])
        cs <- ya
        ds <- ya
        y <- ya[[ns]]
        ns <- ns - 1
        for (m in 1:4) {
          for (i in 1:(5 - m)) {
            ho <- xa[i]
            hp <- xa[i + m]
            w <- (cs[[i + 1]] - ds[[i]])/(ho - hp)
            ds[[i]] <- hp * w
            cs[[i]] <- ho * w
          }
          if (2 * ns < 5 - m) {
            dy <- cs[[ns + 1]]
          } else {
            dy <- ds[[ns]]
            ns <- ns - 1
          }
          y <- y + dy
        }
        ss <- y
        errval <- max(abs(dy))
        ssqval <- max(abs(ss))
        if (all(ssqval > 0)) {
          crit <- errval/ssqval
        } else {
          crit <- errval
        }
        if (crit < EPS && iter >= JMIN) break
      }
      s[[iter + 1]] <- s[[iter]]
      h[iter + 1] <- 0.25 * h[iter]
      if (iter == JMAX)
        warning("Failure to converge.")
    }
    inprodmat <- inprodmat + ss
  }
  names(inprodmat) <- NULL
  inprodmat
}



#' Get Multivariate Functional Data from a data frame
#'
#' @param dt
#' A \code{data.frame} containing the discrete data.
#' For each functional variable, a single column,
#' whose name is provided in the argument \code{variables},
#' contains discrete values of that variable for all functional observation.
#' The column indicated by the argument \code{id}
#' denotes which is the functional observation in each row.
#' The column indicated by the argument \code{arg}
#' gives the argument value at which
#' the discrete values of the functional variables are observed for each row.
#' @param domain
#' A numeric vector of length 2 defining
#' the interval over which the functional data object
#' can be evaluated.
#' @param arg
#' A character variable, which is the name of
#' the column of the data frame \code{dt}
#' giving the argument values at which the functional variables
#' are evaluated for each row.
#' @param id
#' A character variable indicating
#' which is the functional observation in each row.
#' @param variables
#' A vector of characters of the column names
#' of the data frame \code{dt}
#' indicating the functional variables.
#' @param n_basis
#' An integer variable specifying the number of basis functions;
#' default value is 30.
#' See details on basis functions.
#' @param lambda
#' A non-negative real number.
#' If you want to use a single specified smoothing parameter
#' for all functional data objects in the dataset,
#' this argument is passed to the function \code{fda::fdPar}.
#' Default value is NULL, in this case the smoothing parameter is chosen
#' by minimizing the generalized cross-validation (GCV)
#' criterion over the grid of values given by the argument.
#' See details on how smoothing parameters work.
#' @param lambda_grid
#' A vector of non-negative real numbers.
#' If \code{lambda} is provided as a single number, this argument is ignored.
#' If \code{lambda} is NULL, then this provides the grid of values
#' over which the optimal smoothing parameter is
#' searched. Default value is \code{10^seq(-10,1,l=20)}.
#' @param ncores
#' If you want parallelization, give the number of cores/threads
#' to be used when doing GCV separately on all observations.
#'
#' @details
#' Basis functions are created with
#' \code{fda::create.bspline.basis(domain, n_basis)}, i.e.
#' B-spline basis functions of order 4 with equally spaced knots
#' are used to create \code{mfd} objects.
#'
#' The smoothing penalty lambda is provided as
#' \code{fda::fdPar(bs, 2, lambda)},
#' where bs is the basis object and 2 indicates
#' that the integrated squared second derivative is penalized.
#'
#' Rather than having a data frame with long format,
#' i.e. with all functional observations in a single column
#' for each functional variable,
#' if all functional observations are observed on a common equally spaced grid,
#' discrete data may be available in matrix form for each functional variable.
#' In this case, see \code{get_mfd_list}.
#'
#' @seealso \code{\link{get_mfd_list}}
#'
#' @return An object of class \code{mfd}.
#' See also \code{?mfd} for additional details on the
#' multivariate functional data class.
#' @export
#' @examples
#' library(funcharts)
#'
#' x <- seq(1, 10, length = 25)
#' y11 <- cos(x)
#' y21 <- cos(2 * x)
#' y12 <- sin(x)
#' y22 <- sin(2 * x)
#' df <- data.frame(id = factor(rep(1:2, each = length(x))),
#'                  x = rep(x, times = 2),
#'                  y1 = c(y11, y21),
#'                  y2 = c(y12, y22))
#'
#' mfdobj <- get_mfd_df(dt = df,
#'                      domain = c(1, 10),
#'                      arg = "x",
#'                      id = "id",
#'                      variables = c("y1", "y2"),
#'                      lambda = 1e-5)
#'
get_mfd_df <- function(dt,
                       domain,
                       arg,
                       id,
                       variables,
                       n_basis = 30,
                       lambda = NULL,
                       lambda_grid = 10^seq(-10, 1, length.out = 10),
                       ncores = 1) {

  if (!(is.data.frame(dt))) {
    stop("dt must be a data.frame")
  }
  if (!(arg %in% names(dt))) {
    stop("domain_var must be the name of a column of dt")
  }
  if (!(id %in% names(dt))) {
    stop("id must be the name of a column of dt")
  }
  if (sum(variables %in% names(dt)) < length(variables)) {
    stop("variables must contain only names of columns of dt")
  }
  if (!is.numeric(domain) | length(domain) != 2) {
    stop("domain must be a vector with two numbers.")
  }

  bs <- create.bspline.basis(domain, n_basis)
  ids <- levels(factor(dt[[id]]))
  n_obs <- length(ids)
  n_var <- length(variables)

  lambda_search <- if (!is.null(lambda)) lambda else lambda_grid

  n_lam <- length(lambda_search)
  ncores <- min(ncores, n_obs)

  fun_gcv <- function(ii) {

    rows_i <- dt[[id]] == ids[ii] &
      dt[[arg]] >= domain[1] &
      dt[[arg]] <= domain[2]
    dt_i <- dt[rows_i, , drop = FALSE]

    x <- dt_i[[arg]]
    y <- data.matrix(dt_i[, variables])

    # Generalized cross validation to choose smoothing parameter
    coefList <- list()
    gcv <- matrix(data=NA,n_var,n_lam)
    rownames(gcv) <- variables
    colnames(gcv) <- lambda_search
    for (h in 1:n_lam) {
      fdpenalty <- fdPar(bs, 2, lambda_search[h])
      smoothObj <- smooth.basis(x, y, fdpenalty)
      coefList[[h]] <- smoothObj$fd$coefs
      gcv[, h] <- smoothObj$gcv
    }

    # If only NA in a row,
    # consider as optimal smoothing parameter the first one in the sequence.
    gcvmin <- apply(gcv, 1, function(x) {
      if(sum(is.na(x)) < n_lam) which.min(x) else 1
    })
    opt_lam <- lambda_search[gcvmin]
    names(opt_lam) <- variables

    coefs <- sapply(1:n_var, function(jj) coefList[[gcvmin[jj]]][, jj])
    colnames(coefs) <- variables

    coefs
  }

  # You need to perform gcv separately for each observation
  if (ncores == 1) {
    coefs_list <- lapply(1:n_obs, fun_gcv)
  } else {
    if (.Platform$OS.type == "unix") {
      coefs_list <- mclapply(1:n_obs, fun_gcv, mc.cores = ncores)
    } else {
      cl <- makeCluster(ncores)
      clusterExport(cl,
                    c("dt",
                      "ids",
                      "domain",
                      "variables",
                      "n_var",
                      "n_lam",
                      "lambda_search"),
                    envir = environment())
      coefs_list <- parLapply(cl, 1:n_obs, fun_gcv)
      stopCluster(cl)
    }
  }

  names(coefs_list) <- ids

  coefs <- array(do.call(rbind, coefs_list),
                 dim = c(bs$nbasis, n_obs, n_var),
                 dimnames = list(bs$names, ids, variables))

  fdObj <- mfd(coefs, bs, list(arg, ids, variables),
               dt[, c(id, arg, variables)], id)
  fdObj

}

#' Get Multivariate Functional Data from a list of matrices
#'
#' @param data_list
#' A named list of matrices.
#' Names of the elements in the list denote the functional variable names.
#' Each matrix in the list corresponds to a functional variable.
#' All matrices must have the same dimension, where
#' the number of rows corresponds to replications, while
#' the number of columns corresponds to the argument values at which
#' functions are evaluated.
#' @param grid
#' A numeric vector, containing the argument values at which
#' functions are evaluated.
#' Its length must be equal to the number of columns
#' in each matrix in data_list.
#' Default is NULL, in this case a vector equally spaced numbers
#' between 0 and 1 is created,
#' with as many numbers as the number of columns in each matrix in data_list.
#' @param n_basis
#' An integer variable specifying the number of basis functions;
#' default value is 30.
#' See details on basis functions.
#' @param lambda
#' A non-negative real number.
#' If you want to use a single specified smoothing parameter
#' for all functional data objects in the dataset,
#' this argument is passed to the function \code{fda::fdPar}.
#' Default value is NULL, in this case the smoothing parameter is chosen
#' by minimizing the generalized cross-validation (GCV) criterion
#' over the grid of values given by the argument.
#' See details on how smoothing parameters work.
#' @param lambda_grid
#' A vector of non-negative real numbers.
#' If \code{lambda} is provided as a single number, this argument is ignored.
#' If \code{lambda} is NULL, then this provides
#' the grid of values over which the optimal smoothing parameter is
#' searched. Default value is \code{10^seq(-10,1,l=20)}.
#' @param ncores
#' If you want parallelization, give the number of cores/threads
#' to be used when doing GCV separately on all observations.
#'
#'
#' @details
#' Basis functions are created with
#' \code{fda::create.bspline.basis(domain, n_basis)}, i.e.
#' B-spline basis functions of order 4 with equally spaced knots
#' are used to create \code{mfd} objects.
#'
#' The smoothing penalty lambda is provided as
#' \code{fda::fdPar(bs, 2, lambda)},
#' where bs is the basis object and 2 indicates that
#' the integrated squared second derivative is penalized.
#'
#' Rather than having a list of matrices,
#' you may have a data frame with long format,
#' i.e. with all functional observations in a single column
#' for each functional variable.
#' In this case, see \code{get_mfd_df}.
#'
#' @return
#' An object of class \code{mfd}.
#' See also \code{\link{mfd}} for additional details
#' on the multivariate functional data class.
#'
#' @seealso \code{\link{mfd}},
#' \code{\link{get_mfd_list}},
#' \code{\link{get_mfd_array}}
#'
#' @export
#'
#' @examples
#' library(funcharts)
#' data("air")
#' # Only take first 5 multivariate functional observations
#' # and only two variables from air
#' air_small <- lapply(air[c("NO2", "CO")], function(x) x[1:5, ])
#' mfdobj <- get_mfd_list(data_list = air_small)
#'
get_mfd_list <- function(data_list,
                         grid = NULL,
                         n_basis = 30,
                         lambda = NULL,
                         lambda_grid = 10^seq(-10, 1, length.out = 10),
                         ncores = 1) {

  if (!(is.list(data_list))) {
    stop("data_list must be a list of matrices")
  }
  if (is.null(names(data_list))) {
    stop("data_list must be a named list")
  }
  if (length(unique(lapply(data_list, dim))) > 1) {
    stop("data_list must be a list of matrices all of the same dimensions")
  }
  n_args <- ncol(data_list[[1]])
  if (!is.null(grid) & (length(grid) != n_args)) {
    stop(paste0(
      "grid length, ", length(grid),
      " has not the same length as number of ",
      "observed data per functional observation, ",
      ncol(data_list[[1]])))
  }

  if (is.null(grid)) grid <- seq(0, 1, l = n_args)
  domain <- range(grid)
  bs <- create.bspline.basis(domain, n_basis)

  variables <- names(data_list)
  n_var <- length(variables)

  ids <- rownames(data_list[[1]])
  if (is.null(ids)) ids <- 1:nrow(data_list[[1]])
  n_obs <- length(ids)

  lambda_search <- if (!is.null(lambda)) lambda else lambda_grid
  n_lam <- length(lambda_search)
  ncores <- min(ncores, n_obs)


  df <- cbind(
    data.frame(arg = rep(grid, n_obs)),
    data.frame(id = rep(ids, each = n_args)),
    lapply(seq_along(data_list), function(ii) {
      data_list[[ii]] %>%
        t() %>%
        as.data.frame() %>%
        setNames(ids) %>%
        pivot_longer(everything()) %>%
        mutate(name = factor(.data$name, levels = ids)) %>%
        arrange(.data$name) %>%
        select(.data$value) %>%
        setNames(variables[ii])
    }) %>%
      bind_cols()
  ) %>%
    mutate(id = factor(id, levels = ids))

  get_mfd_df(dt = df, domain = domain, arg = "arg", id = "id",
             variables = variables, n_basis = n_basis,
             lambda = lambda,
             ncores = ncores)

}



#' Get Multivariate Functional Data from a three-dimensional array
#'
#'
#' @param data_array
#' A three-dimensional array.
#' The first dimension corresponds to argument values,
#' the second to replications,
#' and the third to variables within replications.
#' @param grid
#' See \code{\link{get_mfd_list}}.
#' @param n_basis
#' See \code{\link{get_mfd_list}}.
#' @param lambda
#' See \code{\link{get_mfd_list}}.
#' @param lambda_grid
#' See \code{\link{get_mfd_list}}.
#' @param ncores
#' See \code{\link{get_mfd_list}}.
#'
#' @return
#' An object of class \code{mfd}.
#' See also \code{?mfd} for additional details on the
#' multivariate functional data class.
#' @export
#' @seealso
#' \code{\link{get_mfd_list}}, \code{\link{get_mfd_df}}
#'
#' @examples
#' library(funcharts)
#' data("CanadianWeather")
#' mfdobj <- get_mfd_array(CanadianWeather$dailyAv[, 1:10, ],
#'                         lambda = 1e-5)
#' plot_mfd(mfdobj)
#'
get_mfd_array <- function(data_array,
                          grid = NULL,
                          n_basis = 30,
                          lambda = NULL,
                          lambda_grid = 10^seq(- 10, 1, length.out = 10),
                          ncores = 1) {

  n_var <- dim(data_array)[3]
  data_list <- list()
  for (ii in 1:n_var) data_list[[ii]] <- t(data_array[, , ii])
  names(data_list) <- dimnames(data_array)[[3]]
  get_mfd_list(data_list,
               grid = grid,
               n_basis = n_basis,
               lambda = lambda,
               lambda_grid = lambda_grid,
               ncores = ncores)
}


#' Convert a \code{fd} object into a Multivariate Functional Data object.
#'
#'
#' @param fdobj
#' An object of class fd.
#'
#' @return
#' An object of class \code{mfd}.
#' See also \code{?mfd} for additional details on the
#' multivariate functional data class.
#' @export
#' @seealso
#' \code{mfd}
#'
#' @examples
#' library(funcharts)
#' fdobj <- fd()
#' mfdobj <- get_mfd_fd(fdobj)
#' plot_mfd(mfdobj)
get_mfd_fd <- function(fdobj) {

  if (length(fdobj$fdnames[[1]]) > 1) fdobj$fdnames[[1]] <- "time"

  if (!is.fd(fdobj)) {
    stop("fdobj must be an object of class fd.")
  }
  coefs <- fdobj$coefs
  if (length(dim(coefs)) == 2) {
    if (!(identical(colnames(coefs), fdobj$fdnames[[2]]) |
        identical(colnames(coefs), fdobj$fdnames[[3]]))) {
      stop(paste0("colnames(fdobj$coefs) must correspond either to ",
                  "fdobj$fdnames[[2]] (i.e. replication names) or to ",
                  "fdobj$fdnames[[3]] (i.e. variable names)."))
    }
    if (identical(colnames(coefs), fdobj$fdnames[[2]])) {
      coefs <- array(coefs, dim = c(nrow(coefs), ncol(coefs), 1))
      dimnames(coefs) <- list()
      dimnames(coefs)[[1]] <- dimnames(fdobj$coefs)[[1]]
      dimnames(coefs)[[2]] <- fdobj$fdnames[[2]]
      dimnames(coefs)[[3]] <- fdobj$fdnames[[3]]
    }
    if (identical(colnames(coefs), fdobj$fdnames[[3]])) {
      coefs <- array(coefs, dim = c(nrow(coefs), 1, ncol(coefs)))
      dimnames(coefs) <- list()
      dimnames(coefs)[[1]] <- dimnames(fdobj$coefs)[[1]]
      dimnames(coefs)[[2]] <- fdobj$fdnames[[2]]
      dimnames(coefs)[[3]] <- fdobj$fdnames[[3]]
    }
  }
  mfd(coefs,
      fdobj$basis,
      fdobj$fdnames)
}



#' Standardize Multivariate Functional Data.
#'
#' Scale multivariate functional data contained
#' in an object of class \code{mfd}
#' by subtracting the mean function and dividing
#' by the standard deviation function.
#'
#' @param mfdobj
#' A multivariate functional data object of class \code{mfd}.
#' @param center
#' A logical value, or a \code{fd} object.
#' When providing a logical value, if TRUE, \code{mfdobj} is centered,
#' i.e. the functional mean function is calculated and subtracted
#' from all observations in \code{mfdobj},
#' if FALSE, \code{mfdobj} is not centered.
#' If \code{center} is a \code{fd} object, then this function
#' is used as functional mean for centering.
#' @param scale
#' A logical value, or a \code{fd} object.
#' When providing a logical value, if TRUE, \code{mfdobj}
#' is scaled after possible centering,
#' i.e. the functional standard deviation is calculated
#' from all functional observations in \code{mfdobj} and
#' then the observations are divided by this calculated standard deviation,
#' if FALSE, \code{mfdobj} is not scaled.
#' If \code{scale} is a \code{fd} object,
#' then this function is used as standard deviation function for scaling.
#'
#' @return
#' A standardized object of class \code{mfd}, with two attributes,
#' if calculated,
#' \code{center} and \code{scale}, storing the mean and
#' standard deviation functions used for standardization.
#'
#' @details
#' This function has been written to work similarly
#' as the function \code{\link{scale}} for matrices.
#' When calculated, attributes \code{center} and \code{scale}
#' are of class \code{fd}
#' and have the same structure you get
#' when you use \code{fda::\link[fda]{mean.fd}}
#' and \code{fda::\link[fda]{sd.fd}}.
#' @export
#' @examples
#' library(funcharts)
#' mfdobj <- data_sim_mfd()
#' mfdobj_scaled <- scale_mfd(mfdobj)
scale_mfd <- function(mfdobj, center = TRUE, scale = TRUE) {

  if (!is.mfd(mfdobj)) {
    stop("Only mfd class allowed for mfdobj input")
  }

  bs <- mfdobj$basis
  n_obs <- length(mfdobj$fdnames[[2]])

  if (n_obs == 1 & (!is.fd(scale) | !is.fd(center))) {
    stop("There is only one observation in the data set")
  }

  # Center
  if (!(is.logical(center) | is.fd(center))) {
    stop("Only logical or fd classes allowed for center input")
  }

  mean_fd <- NULL
  if (is.logical(center) && center == FALSE) {
    cen_fd <- mfdobj
  } else {
    if (is.logical(center) && center == TRUE) {
      mean_fd <- mean.fd(mfdobj)
      cen_fd <- center.fd(mfdobj)
      cen_fd$fdnames <- mfdobj$fdnames
    }
    if (is.fd(center)) {
      mean_fd <- center
      mean_fd_coefs <- array(
        mean_fd$coefs[, 1, ], dim = c(dim(mean_fd$coefs)[c(1, 3)], n_obs))
      mean_fd_coefs <- aperm(mean_fd_coefs, c(1, 3, 2))
      mean_fd_rep <- fd(mean_fd_coefs, bs, mfdobj$fdnames)
      cen_fd <- minus.fd(mfdobj, mean_fd_rep)
    }
  }

  # Scale
  if (!(is.logical(scale) | is.fd(scale))) {
    stop("Only logical or fd classes allowed for scale input")
  }

  sd_fd <- NULL
  if (is.logical(scale) && scale == FALSE) {
    fd_std <- cen_fd
  } else {
    if (is.logical(scale) && scale == TRUE) {
      sd_fd <- sd.fd(mfdobj)
      sd_fd$coefs <- pmax(sd_fd$coefs, 0)
    }
    if (is.fd(scale)) {
      sd_fd <- scale
    }
    if (sd_fd$basis$nbasis < 15) {
      domain <- mfdobj$basis$rangeval
      x_eval <- seq(domain[1], domain[2], length.out = 1000)
      sd_eval <- eval.fd(evalarg = x_eval, sd_fd)

      bs_sd <- create.bspline.basis(rangeval = domain, nbasis = 20)

      fdpar_more_basis <- fdPar(bs_sd, 2, 0)
      sd_fd_more_basis <- smooth.basis(x_eval, sd_eval, fdpar_more_basis)$fd
      sd_inv <- sd_fd_more_basis^(-1)
      sd_inv_eval <- eval.fd(evalarg = x_eval, sd_inv)
      sd_inv <- smooth.basis(x_eval, sd_inv_eval, fdPar(bs, 2, 0))$fd

    } else {
      sd_inv <- sd_fd^(-1)
    }
    sd_inv_coefs <- array(sd_inv$coefs, dim = c(dim(sd_inv$coefs), n_obs))
    sd_inv_coefs <- aperm(sd_inv_coefs, c(1, 3, 2))
    sd_inv_rep <- fd(sd_inv_coefs, bs, mfdobj$fdnames)
    fd_std <- times.fd(cen_fd, sd_inv_rep, bs)
    fd_std$fdnames <- mfdobj$fdnames
    dimnames(fd_std$coefs) <- dimnames(mfdobj$coefs)
    dimnames(fd_std$fdnames) <- NULL
  }

  attr(fd_std, "scaled:center") <- mean_fd
  attr(fd_std, "scaled:scale") <- sd_fd

  fd_std$raw <- NULL
  fd_std$id_var <- mfdobj$id_var
  class(fd_std) <- c("mfd", "fd")

  fd_std
}


#' De-standardize a standardized Multivariate Functional Data object
#'
#' This function takes a scaled
#' multivariate functional data contained in an object of class \code{mfd},
#' multiplies it by function provided by the argument \code{scale}
#' and adds the function provided by the argument \code{center}.
#'
#' @param scaled_mfd
#' A scaled multivariate functional data object,
#' of class \code{mfd}.
#' @param center
#' A functional data object of class \code{fd},
#' having the same structure you get
#' when you use \code{fda::\link[fda]{mean.fd}}
#' over a \code{mfd} object.
#' @param scale
#' A functional data object of class \code{fd},
#' having the same structure you get
#' when you use \code{fda::\link[fda]{sd.fd}}
#' over a \code{mfd} object.
#'
#' @return
#' A de-standardized object of class \code{mfd},
#' obtained by multiplying \code{scaled_mfd} by \code{scale} and then
#' adding \code{center}.
#' @noRd
#'
descale_mfd <- function(scaled_mfd, center, scale) {

  basis <- scaled_mfd$basis
  nbasis <- basis$nbasis
  nobs <- length(scaled_mfd$fdnames[[2]])
  nvar <- length(scaled_mfd$fdnames[[3]])

  coef_sd_list <- lapply(1:nvar, function(jj) {
    matrix(scale$coefs[, jj], nrow = nbasis, ncol = nobs)
  })
  coef_sd <- simplify2array(coef_sd_list)
  sd_fd <- fd(coef_sd, scale$basis, scaled_mfd$fdnames)

  coef_mean_list <- lapply(1:nvar, function(jj) {
    matrix(center$coefs[, 1, jj], nrow = nbasis, ncol = nobs)
  })
  coef_mean <- simplify2array(coef_mean_list)
  mean_fd <- fd(coef_mean, center$basis, scaled_mfd$fdnames)

  centered <- times.fd(scaled_mfd, sd_fd, basisobj = basis)
  descaled <- centered + mean_fd
  dimnames(descaled$coefs) <- dimnames(scaled_mfd$coefs)
  descaled$fdnames <- scaled_mfd$fdnames

  mfd(descaled$coefs, descaled$basis, descaled$fdnames)
}


#' Tensor product of two Multivariate Functional Data objects
#'
#' This function returns the tensor product of two
#' Multivariate Functional Data objects.
#' Each object must contain only one replication.
#'
#' @param mfdobj1
#' A multivariate functional data object, of class \code{mfd},
#' having only one functional observation.
#' @param mfdobj2
#' A multivariate functional data object, of class \code{mfd},
#' having only one functional observation.
#' If NULL, it is set equal to \code{mfdobj1}. Default is NULL.
#'
#' @return
#' An object of class \code{bifd}.
#' If we denote with x(s)=(x_1(s),\dots,x_p(s))
#' the vector of p functions represented by \code{mfdobj1} and
#' with y(t)=(y_1(t),\dots,y_q(t)) the vector of q functions
#' represented by \code{mfdobj2},
#' the output is the
#' vector of pq bivariate functions
#'
#' f(s,t)=(x_1(s)y_1(t),\dots,x_1(s)y_q(t),
#' \dots,x_p(s)y_1(t),\dots,x_p(s)y_q(t)).
#'
#'
#' @export
#' @examples
#' library(funcharts)
#' mfdobj1 <- data_sim_mfd(nobs = 1, nvar = 3)
#' mfdobj2 <- data_sim_mfd(nobs = 1, nvar = 2)
#' tensor_product_mfd(mfdobj1)
#' tensor_product_mfd(mfdobj1, mfdobj2)
#'
tensor_product_mfd <- function(mfdobj1, mfdobj2 = NULL) {

  if (!is.mfd(mfdobj1)) {
    stop("First argument must be a mfd object.")
  }
  if (is.null(mfdobj2)) mfdobj2 <- mfdobj1
  if (!is.mfd(mfdobj2)) {
    stop("First argument must be a mfd object.")
  }
  obs1 <- mfdobj1$fdnames[[2]]
  obs2 <- mfdobj2$fdnames[[2]]
  nobs1 <- length(obs1)
  nobs2 <- length(obs2)
  if (nobs1 > 1) {
    stop("mfdobj1 must have only 1 observation")
  }
  if (nobs2 > 1) {
    stop("mfdobj2 must have only 1 observation")
  }
  variables1 <- mfdobj1$fdnames[[3]]
  variables2 <- mfdobj2$fdnames[[3]]
  nvar1 <- length(variables1)
  nvar2 <- length(variables2)

  coef1 <- mfdobj1$coefs
  coef2 <- mfdobj2$coefs
  coef1 <- matrix(coef1, nrow = dim(coef1)[1], ncol = dim(coef1)[3])
  coef2 <- matrix(coef2, nrow = dim(coef2)[1], ncol = dim(coef2)[3])

  basis1 <- mfdobj1$basis
  basis2 <- mfdobj2$basis

  nbasis1 <- basis1$nbasis
  nbasis2 <- basis2$nbasis

  coef <- outer(coef1, coef2)
  coef <- aperm(coef, c(1, 3, 4, 2))
  coef <- array(coef, dim = c(nbasis1, nbasis2, 1, nvar1 * nvar2))

  dim_names <- expand.grid(variables1, variables2)
  variables <- paste(dim_names[, 1], dim_names[, 2])

  obs <- paste0("s.", obs1, " t.", obs2)

  fdnames <- list(basis1$names,
                  basis2$names,
                  obs,
                  variables)

  dimnames(coef) <- fdnames

  bifd(coef, basis1, basis2, fdnames)

}


#' Bind variables of two Multivariate Functional Data Objects
#'
#' @param mfdobj1
#' An object of class mfd, with the same number of replications of mfdobj2
#' and different variable names with respect to mfdobj2.
#' @param mfdobj2
#' An object of class mfd, with the same number of replications of mfdobj1,
#' and different variable names with respect to mfdobj1.
#'
#' @return
#' An object of class mfd, whose replications are the same of mfdobj1 and
#' mfdobj2 and whose functional variables are the union of the functional
#' variables in mfdobj1 and mfdobj2.
#' @export
#'
#' @examples
#' library(funcharts)
#' mfdobj1 <- data_sim_mfd(nvar = 3, seed = 1)
#' mfdobj2 <- data_sim_mfd(nvar = 2, seed = 2)
#' dimnames(mfdobj2$coefs)[[3]] <- mfdobj2$fdnames[[3]] <- c("var10", "var11")
#'
#' plot_mfd(mfdobj1)
#' plot_mfd(mfdobj2)
#' mfdobj_cbind <- cbind_mfd(mfdobj1, mfdobj2)
#' plot_mfd(mfdobj_cbind)
#'
cbind_mfd <- function(mfdobj1, mfdobj2) {

  if (!is.null(mfdobj1) & !is.mfd(mfdobj1)) {
    stop("mfdobj1 must be an object of class mfd if not NULL")
  }
  if (!is.null(mfdobj2) & !is.mfd(mfdobj2)) {
    stop("mfdobj2 must be an object of class mfd if not NULL")
  }

  if (is.null(mfdobj2)) {
    return(mfdobj1)
  }
  if (is.null(mfdobj1)) {
    return(mfdobj2)
  }

  if (dim(mfdobj1$coefs)[2] != dim(mfdobj2$coefs)[2]) {
    stop("mfdobj1 and mfdobj2 must have the same number of replications.")
  }
  if (!identical(mfdobj1$basis, mfdobj2$basis)) {
    stop("mfdobj1 and mfdobj2 must have the same basis.")
  }
  mfd(coef = array(c(mfdobj1$coefs,
                     mfdobj2$coefs),
                   dim = c(dim(mfdobj1$coefs)[1],
                           dim(mfdobj1$coefs)[2],
                           dim(mfdobj1$coefs)[3] + dim(mfdobj2$coefs)[3])),
      basisobj = mfdobj1$basis,
      fdnames = list(mfdobj1$fdnames[[1]],
                     mfdobj1$fdnames[[2]],
                     c(mfdobj1$fdnames[[3]],
                       mfdobj2$fdnames[[3]])))
}


#' Bind replications of two Multivariate Functional Data Objects
#'
#' @param mfdobj1
#' An object of class mfd, with the same variables of mfdobj2
#' and different replication names with respect to mfdobj2.
#' @param mfdobj2
#' An object of class mfd, with the same variables of mfdobj1,
#' and different replication names with respect to mfdobj1.
#'
#' @return
#' An object of class mfd, whose variables are the same of mfdobj1 and
#' mfdobj2 and whose replications are the union of the replications
#' in mfdobj1 and mfdobj2.
#' @export
#'
#' @examples
#' library(funcharts)
#' mfdobj1 <- data_sim_mfd(nvar = 3, seed = 1, nobs = 4)
#' mfdobj2 <- data_sim_mfd(nvar = 3, seed = 2, nobs = 5)
#' dimnames(mfdobj2$coefs)[[2]] <-
#'   mfdobj2$fdnames[[2]] <-
#'   c("rep11", "rep12", "rep13", "rep14", "rep15")
#' mfdobj_rbind <- rbind_mfd(mfdobj1, mfdobj2)
#' plot_mfd(mfdobj_rbind)
#'
rbind_mfd <- function(mfdobj1, mfdobj2) {

  if (!is.null(mfdobj1) & !is.mfd(mfdobj1)) {
    stop("mfdobj1 must be an object of class mfd if not NULL")
  }
  if (!is.null(mfdobj2) & !is.mfd(mfdobj2)) {
    stop("mfdobj2 must be an object of class mfd if not NULL")
  }

  if (is.null(mfdobj2)) {
    return(mfdobj1)
  }
  if (is.null(mfdobj1)) {
    return(mfdobj2)
  }

  if (dim(mfdobj1$coefs)[3] != dim(mfdobj2$coefs)[3]) {
    stop("mfdobj1 and mfdobj2 must have the same number of variables")
  }
  if (!identical(mfdobj1$basis, mfdobj2$basis)) {
    stop("mfdobj1 and mfdobj2 must have the same basis.")
  }

  coef1_aperm <- aperm(mfdobj1$coefs, c(1, 3, 2))
  coef2_aperm <- aperm(mfdobj2$coefs, c(1, 3, 2))
  coef_aperm <- array(c(coef1_aperm, coef2_aperm),
                      dim = c(dim(mfdobj1$coefs)[1],
                              dim(mfdobj1$coefs)[3],
                              dim(mfdobj1$coefs)[2] + dim(mfdobj2$coefs)[2]))
  coef <- aperm(coef_aperm, c(1, 3, 2))

  mfd(coef = coef,
      basisobj = mfdobj1$basis,
      fdnames = list(mfdobj1$fdnames[[1]],
                     c(mfdobj1$fdnames[[2]],
                       mfdobj2$fdnames[[2]]),
                     mfdobj2$fdnames[[3]]))
}



# Plots -------------------------------------------------------------------

#' Convert a Multivariate Functional Data Object into a data frame
#'
#' This function extracts the argument \code{raw} of an
#' object of class \code{mfd}.
#'
#' @param mfdobj A multivariate functional data object of class mfd.
#'
#' @return
#' A \code{data.frame} in the long format,
#' which has been used to create the object \code{mfdobj}
#' using the function \code{\link{get_mfd_df}},
#' or that has been created when creating \code{mfdobj}
#' using the function \code{\link{get_mfd_list}}.
#'
#' @seealso \code{\link{get_mfd_df}}, \code{\link{get_mfd_list}}
#' @noRd
mfd_to_df_raw <- function(mfdobj) {

  if (!(is.mfd(mfdobj))) {
    stop("Input must be a multivariate functional data object.")
  }

  dt <- mfdobj$raw
  id_var <- mfdobj$id_var

  arg_var <- mfdobj$fdnames[[1]]
  obs <- mfdobj$fdnames[[2]]
  variables <- mfdobj$fdnames[[3]]
  dt %>%
    select(variables, !!id_var, !!arg_var) %>%
    rename(id = !!id_var) %>%
    filter(id %in% !!obs) %>%
    pivot_longer(variables, names_to = "var") %>%
    arrange(id, var, !!arg_var) %>%
    drop_na()
}


#' Discretize a Multivariate Functional Data Object
#'
#' This function discretizes an object of class \code{mfd}
#' and stores it into a \code{data.frame} object in the long format.
#'
#' @param mfdobj A multivariate functional data object of class mfd.
#'
#' @return
#' A \code{data.frame} in the long format,
#' obtained by discretizing the functional values using
#' \code{fda::\link[fda]{eval.fd}}
#' on a grid of 200 equally spaced argument values
#' covering the functional domain.
#' @noRd
mfd_to_df <- function(mfdobj) {

  if (!(is.mfd(mfdobj))) {
    stop("Input must be a multivariate functional data object.")
  }

  n_obs <- length(mfdobj$fdnames[[2]])
  arg_var <- mfdobj$fdnames[[1]]
  range <- mfdobj$basis$rangeval
  evalarg <- seq(range[1], range[2], l = 200)
  X <- eval.fd(evalarg, mfdobj)
  id <- mfdobj$fdnames[[2]]
  variables <- mfdobj$fdnames[[3]]
  lapply(seq_along(variables), function(jj) {
    variable <- variables[jj]

    .df <- as.data.frame(X[, , jj, drop = FALSE]) %>%
      setNames(id)
    .df[[arg_var]] <- evalarg
    .df$var <- variable
    pivot_longer(.df, -c(!!arg_var, var), names_to = "id")
  }) %>%
    bind_rows() %>%
    mutate(var = factor(var, levels = !!variables),
           id = factor(id, levels = !!id))
}


#' Plot a Multivariate Functional Data Object.
#'
#' Plot an object of class \code{mfd} using \code{ggplot2}.
#'
#' @param mfdobj A multivariate functional data object of class mfd.
#'
#' @return an object of class ggplot, created using
#' \code{ggplot() + geom_mfd(mfdobj = mfdobj)}.
#' @export
#'
#' @seealso \code{\link{geom_mfd}}
#' @examples
#' library(funcharts)
#' mfdobj <- data_sim_mfd()
#' plot_mfd(mfdobj)
#'
plot_mfd <- function(mfdobj) {
  ggplot() + geom_mfd(mfdobj = mfdobj)
}

#' Creates a geom layer to plot a Multivariate Functional Data Object
#' with \code{ggplot}
#'
#' @param mfdobj
#' A multivariate functional data object of class mfd.
#' @param data
#' A \code{data.frame} providing columns
#' to create additional aesthetic mappings.
#' It must contain a column "id" with the replication values
#' as in \code{mfdobj$fdnames[[2]]}.
#' If it contains a column "var", this must contain
#' the functional variables as in \code{mfdobj$fdnames[[3]]}.
#' @param mapping
#' Set of aesthetic mappings additional
#' to \code{x} and \code{y} as passed to the function \code{ggplot2::geom:line}.
#' @param stat
#' See \code{ggplot2::\link[ggplot2]{geom_line}}.
#' @param position
#' See \code{ggplot2::\link[ggplot2]{geom_line}}.
#' @param na.rm
#' See \code{ggplot2::\link[ggplot2]{geom_line}}.
#' @param orientation
#' See \code{ggplot2::\link[ggplot2]{geom_line}}.
#' @param show.legend
#' See \code{ggplot2::\link[ggplot2]{geom_line}}.
#' @param inherit.aes
#' See \code{ggplot2::\link[ggplot2]{geom_line}}.
#' @param ...
#' See \code{ggplot2::\link[ggplot2]{geom_line}}.
#' @param type_mfd
#' A character value equal to "mfd" or "raw".
#' If "mfd", the smoothed functional data are plotted, if "raw",
#' the original discrete data are plotted.
#'
#' @return
#' A geom_line layer to be added to
#' \code{ggplot2::\link[ggplot2]{ggplot}()}
#' in order to plot \code{mfdobj}.
#' @export
#' @examples
#' library(funcharts)
#' mfdobj <- data_sim_mfd()
#' ids <- mfdobj$fdnames[[2]]
#' df <- data.frame(id = ids, first_two_obs = ids %in% c("rep1", "rep2"))
#' ggplot() +
#'   geom_mfd(mapping = aes(colour = first_two_obs),
#'            data = df,
#'            mfdobj = mfdobj)
#'
geom_mfd <- function(mapping = NULL,
                     data = NULL,
                     mfdobj,
                     stat = "identity",
                     position = "identity",
                     na.rm = TRUE,
                     orientation = NA,
                     show.legend = NA,
                     inherit.aes = TRUE,
                     type_mfd = "mfd",
                     ...) {

  if (!(is.mfd(mfdobj))) {
    stop("First argument must be a multivariate functional data object.")
  }

  if (!(type_mfd %in% c("mfd", "raw"))) {
    stop("type_mfd not 'mfd' or 'raw'")
  }
  if (type_mfd == "mfd") df <- mfd_to_df(mfdobj)
  if (type_mfd == "raw") df <- mfd_to_df_raw(mfdobj)

  if (!is.null(data)) {
    join_vars <- "id"
    if ("var" %in% names(data)) join_vars <- c(join_vars, "var")
    df <- inner_join(df, data, by = join_vars)
  }
  variables <- mfdobj$fdnames[[3]]
  df$var <- factor(as.character(df$var), levels = variables)
  arg_var <- mfdobj$fdnames[[1]]
  mapping1 <- aes_string(arg_var, "value", group = "id")
  mapping_tot <- c(mapping1, mapping)
  class(mapping_tot) <- "uneval"

  list(
    geom_line(mapping = mapping_tot,
              data = df,
              stat = stat,
              position = position,
              na.rm = na.rm,
              orientation = orientation,
              show.legend = show.legend,
              inherit.aes = inherit.aes,
              ...),
    ylab(NULL), xlab(arg_var),
    theme_bw(),
    theme(strip.background = element_blank(),
          strip.placement = "outside"),
    facet_wrap(~var, scales = "free_y",
               strip.position = "left",
               # labeller = as_labeller(variable_labels)
    )
  )
}




#' Plot a Bivariate Functional Data Object.
#'
#' Plot an object of class \code{bifd} using
#' \code{ggplot2} and \code{geom_tile}.
#' The object must contain only one single functional replication.
#'
#' @param bifd_obj A bivariate functional data object of class bifd,
#' containing one single replication.
#'
#' @return
#' A ggplot with a geom_tile layer providing a plot of the
#' bivariate functional data object as a heat map.
#' @export
#'
#' @examples
#' library(funcharts)
#' mfdobj <- data_sim_mfd(nobs = 1)
#' tp <- tensor_product_mfd(mfdobj)
#' plot_bifd(tp)
#'
plot_bifd <- function(bifd_obj) {

  if (class(bifd_obj) != "bifd") {
    stop("bifd_obj must be an object of class bifd")
  }
  if (length(dim(bifd_obj$coef)) != 4) {
    stop("length of bifd_obj$coef must be 4")
  }
  if (dim(bifd_obj$coef)[3] != 1) {
    stop("third dimension of bifd_obj$coef must be 1")
  }


  s_eval <- seq(bifd_obj$sbasis$rangeval[1],
                bifd_obj$sbasis$rangeval[2],
                l = 100)
  t_eval <- seq(bifd_obj$tbasis$rangeval[1],
                bifd_obj$tbasis$rangeval[2],
                l = 100)
  X_eval <- eval.bifd(s_eval, t_eval, bifd_obj)

  seq_along(bifd_obj$bifdnames[[4]]) %>%
    lapply(function(ii) {
      X_eval[, , , ii] %>%
        data.frame() %>%
        setNames(t_eval) %>%
        mutate(s = s_eval) %>%
        pivot_longer(-.data$s, names_to = "t", values_to = "value") %>%
        mutate(t = as.numeric(.data$t),
               variable = bifd_obj$bifdnames[[4]][ii])
    }) %>%
    bind_rows() %>%
    mutate(variable = factor(.data$variable, levels = bifd_obj$bifdnames[[4]])) %>%
    ggplot() +
    geom_tile(aes(.data$s, .data$t, fill = .data$value)) +
    facet_wrap(~variable) +
    scale_fill_gradientn(
      colours = c("blue", "white", "red"),
      limits = c(- max(abs(X_eval)), max(abs(X_eval)))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))

}


