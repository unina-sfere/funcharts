#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat score(const arma::mat& e,
                const double lambda,
                const arma::mat& k,
                bool huber=true) {


  arma::mat phi = lambda * e;

  // Huber's Score Function
  if (huber) {

    arma::uvec idx1 = find(e > k);
    arma::uvec idx2 = find(e < -k);

    phi.elem(idx1) = e.elem(idx1) - (1 - lambda) * k.elem(idx1);
    phi.elem(idx2) = e.elem(idx2) + (1 - lambda) * k.elem(idx2);

  }
  // bs score function
  else {

    phi = e;
    arma::uvec idx = find(abs(e) <= k);
    phi.elem(idx) = e.elem(idx) * (1 - (1 - lambda) * pow(1 - square(e.elem(idx) / k.elem(idx)), 2));

  }

  return phi;
}



// [[Rcpp::export]]
arma::vec score2(const arma::vec& e,
                 const double lambda,
                 const arma::vec& k,
                 bool huber = true) {


  arma::vec phi = lambda * e;

  // Huber's Score Function
  if (huber) {

    arma::uvec idx1 = find(e > k);
    arma::uvec idx2 = find(e < -k);

    phi.elem(idx1) = e.elem(idx1) - (1 - lambda) * k.elem(idx1);
    phi.elem(idx2) = e.elem(idx2) + (1 - lambda) * k.elem(idx2);

  }
  // bs score function
  else {

    phi = e;
    arma::uvec idx = find(abs(e) <= k);
    phi.elem(idx) = e.elem(idx) * (1 - (1 - lambda) * pow(1 - square(e.elem(idx) / k.elem(idx)), 2));

  }

  return phi;
}


// [[Rcpp::export]]
arma::vec statisticY_EWMA_vec(const arma::vec& X,
                              const arma::vec& Y_previous,
                              double lambda,
                              const arma::vec& k,
                              bool huber) {

  int nvars = X.n_elem;
  arma::vec e = X - Y_previous;
  arma::vec phi = score2(e, lambda, k, huber);
  arma::vec W = phi / e;
  for (int i = 0; i < nvars; ++i) {
    if (e[i] == 0) {
      W[i] = lambda;
    }
  }
  arma::vec Y = (W % X) + ((1 - W) % Y_previous);
  return Y;
}

// [[Rcpp::export]]
arma::mat statisticY_EWMA_cpp(const arma::mat& X,
                              double lambda,
                              const arma::vec& k,
                              bool huber,
                              const arma::vec& idx) {

  int nvars = X.n_cols;
  arma::mat Y(idx.n_elem, nvars); // Matrix to hold all the vectors
  arma::vec Yvec = arma::zeros<arma::vec>(nvars);

  for (arma::uword kk = 0; kk < idx.n_elem; ++kk) {
    Yvec = statisticY_EWMA_vec(X.row(idx(kk) - 1).t(),
                               Yvec,
                               lambda,
                               k,
                               huber);
    Y.row(kk) = Yvec.t();
  }

  return Y;
}


// [[Rcpp::export]]
double calculate_T2(const arma::vec& Y, const arma::mat& Vectors, const arma::vec& Values) {
  // Check dimensions
  if (Y.n_elem != Vectors.n_rows || Values.n_elem != Vectors.n_cols) {
    throw std::runtime_error("Dimension mismatch");
  }

  // Perform the optimized calculation
  arma::vec temp = Vectors.t() * Y; // This is t(Vectors) %*% Y in R
  temp = temp % temp; // Element-wise squaring, equivalent to square each element of the vector
  temp = temp / Values; // Element-wise division by Values, equivalent to each element of the vector divided by corresponding element in Values

  // Sum up the elements for the final result
  double result = sum(temp);

  // Return the result
  return result;
}


// [[Rcpp::export]]
arma::vec calculate_T2_vec(const arma::mat& Y, const arma::mat& Vectors, const arma::vec& Values) {
  // Check dimensions
  if (Y.n_cols != Vectors.n_rows || Values.n_elem != Vectors.n_cols) {
    throw std::runtime_error("Dimension mismatch");
  }

  arma::mat YVsquared = arma::square(Y * Vectors);
  arma::mat temp = YVsquared.each_row() / Values.t();
  arma::vec result = arma::sum(temp, 1);

  return result;
}


// [[Rcpp::export]]
List get_RL_cpp(const arma::mat& X2,
                const arma::mat& X_IC,
                const arma::vec& idx2,
                const arma::vec& idx_IC,
                double lambda,
                const arma::vec& k,
                bool huber,
                double h,
                const arma::vec& Values,
                const arma::mat& Vectors) {

  arma::uword nvars = X2.n_cols;
  arma::uword nmax = idx2.n_elem;
  arma::uword nwarm = idx_IC.n_elem;
  arma::vec T2_IC(nwarm);
  arma::vec T2(nmax);
  int RL = 0;
  arma::vec T2_out;

  arma::vec Y = arma::zeros<arma::vec>(nvars);
  for (arma::uword kk = 0; kk < nwarm; ++kk) {
    arma::vec Xkk = X_IC.row(idx_IC(kk) - 1).t();
    Y = statisticY_EWMA_vec(Xkk, Y, lambda, k, huber);

    T2_IC(kk) = calculate_T2(Y, Vectors, Values);

  }

  for (arma::uword kk = 0; kk < nmax; ++kk) {
    arma::vec Xkk = X2.row(idx2(kk) - 1).t();
    Y = statisticY_EWMA_vec(Xkk, Y, lambda, k, huber);
    T2(kk) = calculate_T2(Y, Vectors, Values);
    if (T2(kk) > h) {
      T2.resize(kk + 1);
      T2_out = T2;
      RL = kk + 1;
      break;
    }
    if ((kk == (nmax - 1)) & (T2(kk) < h)) {
      RL = NA_INTEGER;
    }
  }

  List resultList;
  resultList["RL"] = RL;
  resultList["T2"] = T2;
  resultList["T2_IC"] = T2_IC;

  return resultList;

}





