// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// score
arma::mat score(const arma::mat& e, const double lambda, const arma::mat& k, bool huber);
RcppExport SEXP _funcharts_score(SEXP eSEXP, SEXP lambdaSEXP, SEXP kSEXP, SEXP huberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type e(eSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type huber(huberSEXP);
    rcpp_result_gen = Rcpp::wrap(score(e, lambda, k, huber));
    return rcpp_result_gen;
END_RCPP
}
// score2
arma::vec score2(const arma::vec& e, const double lambda, const arma::vec& k, bool huber);
RcppExport SEXP _funcharts_score2(SEXP eSEXP, SEXP lambdaSEXP, SEXP kSEXP, SEXP huberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type e(eSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type huber(huberSEXP);
    rcpp_result_gen = Rcpp::wrap(score2(e, lambda, k, huber));
    return rcpp_result_gen;
END_RCPP
}
// statisticY_EWMA_vec
arma::vec statisticY_EWMA_vec(const arma::vec& X, const arma::vec& Y_previous, double lambda, const arma::vec& k, bool huber);
RcppExport SEXP _funcharts_statisticY_EWMA_vec(SEXP XSEXP, SEXP Y_previousSEXP, SEXP lambdaSEXP, SEXP kSEXP, SEXP huberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y_previous(Y_previousSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type huber(huberSEXP);
    rcpp_result_gen = Rcpp::wrap(statisticY_EWMA_vec(X, Y_previous, lambda, k, huber));
    return rcpp_result_gen;
END_RCPP
}
// statisticY_EWMA_cpp
arma::mat statisticY_EWMA_cpp(const arma::mat& X, double lambda, const arma::vec& k, bool huber, const arma::vec& idx);
RcppExport SEXP _funcharts_statisticY_EWMA_cpp(SEXP XSEXP, SEXP lambdaSEXP, SEXP kSEXP, SEXP huberSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type huber(huberSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(statisticY_EWMA_cpp(X, lambda, k, huber, idx));
    return rcpp_result_gen;
END_RCPP
}
// calculate_T2
double calculate_T2(const arma::vec& Y, const arma::mat& Vectors, const arma::vec& Values);
RcppExport SEXP _funcharts_calculate_T2(SEXP YSEXP, SEXP VectorsSEXP, SEXP ValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Vectors(VectorsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Values(ValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_T2(Y, Vectors, Values));
    return rcpp_result_gen;
END_RCPP
}
// calculate_T2_vec
arma::vec calculate_T2_vec(const arma::mat& Y, const arma::mat& Vectors, const arma::vec& Values);
RcppExport SEXP _funcharts_calculate_T2_vec(SEXP YSEXP, SEXP VectorsSEXP, SEXP ValuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Vectors(VectorsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Values(ValuesSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_T2_vec(Y, Vectors, Values));
    return rcpp_result_gen;
END_RCPP
}
// get_RL_cpp
List get_RL_cpp(const arma::mat& X2, const arma::mat& X_IC, const arma::vec& idx2, const arma::vec& idx_IC, double lambda, const arma::vec& k, bool huber, double h, const arma::vec& Values, const arma::mat& Vectors);
RcppExport SEXP _funcharts_get_RL_cpp(SEXP X2SEXP, SEXP X_ICSEXP, SEXP idx2SEXP, SEXP idx_ICSEXP, SEXP lambdaSEXP, SEXP kSEXP, SEXP huberSEXP, SEXP hSEXP, SEXP ValuesSEXP, SEXP VectorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X_IC(X_ICSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type idx2(idx2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type idx_IC(idx_ICSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type huber(huberSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Values(ValuesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Vectors(VectorsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_RL_cpp(X2, X_IC, idx2, idx_IC, lambda, k, huber, h, Values, Vectors));
    return rcpp_result_gen;
END_RCPP
}
// der_c
double der_c(double asn, double smin, double smax, double der_0);
RcppExport SEXP _funcharts_der_c(SEXP asnSEXP, SEXP sminSEXP, SEXP smaxSEXP, SEXP der_0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type asn(asnSEXP);
    Rcpp::traits::input_parameter< double >::type smin(sminSEXP);
    Rcpp::traits::input_parameter< double >::type smax(smaxSEXP);
    Rcpp::traits::input_parameter< double >::type der_0(der_0SEXP);
    rcpp_result_gen = Rcpp::wrap(der_c(asn, smin, smax, der_0));
    return rcpp_result_gen;
END_RCPP
}
// loss_c
double loss_c(double z, double y, double z1, double y1, double alpha, double der_h);
RcppExport SEXP _funcharts_loss_c(SEXP zSEXP, SEXP ySEXP, SEXP z1SEXP, SEXP y1SEXP, SEXP alphaSEXP, SEXP der_hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type z1(z1SEXP);
    Rcpp::traits::input_parameter< double >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type der_h(der_hSEXP);
    rcpp_result_gen = Rcpp::wrap(loss_c(z, y, z1, y1, alpha, der_h));
    return rcpp_result_gen;
END_RCPP
}
// DP3
List DP3(int N, int M, arma::vec l, arma::vec u, arma::mat range_x, arma::mat range_tem, arma::mat grid_t, List x_fd_std, List der_x_fd_std, double delta_x, arma::vec template_eval, arma::vec der_template_eval, double smin, double smax, double alpha, double lambda, double der_0, Rcpp::Function eval_fd_c);
RcppExport SEXP _funcharts_DP3(SEXP NSEXP, SEXP MSEXP, SEXP lSEXP, SEXP uSEXP, SEXP range_xSEXP, SEXP range_temSEXP, SEXP grid_tSEXP, SEXP x_fd_stdSEXP, SEXP der_x_fd_stdSEXP, SEXP delta_xSEXP, SEXP template_evalSEXP, SEXP der_template_evalSEXP, SEXP sminSEXP, SEXP smaxSEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP der_0SEXP, SEXP eval_fd_cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type l(lSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type range_x(range_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type range_tem(range_temSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type grid_t(grid_tSEXP);
    Rcpp::traits::input_parameter< List >::type x_fd_std(x_fd_stdSEXP);
    Rcpp::traits::input_parameter< List >::type der_x_fd_std(der_x_fd_stdSEXP);
    Rcpp::traits::input_parameter< double >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type template_eval(template_evalSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type der_template_eval(der_template_evalSEXP);
    Rcpp::traits::input_parameter< double >::type smin(sminSEXP);
    Rcpp::traits::input_parameter< double >::type smax(smaxSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type der_0(der_0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type eval_fd_c(eval_fd_cSEXP);
    rcpp_result_gen = Rcpp::wrap(DP3(N, M, l, u, range_x, range_tem, grid_t, x_fd_std, der_x_fd_std, delta_x, template_eval, der_template_eval, smin, smax, alpha, lambda, der_0, eval_fd_c));
    return rcpp_result_gen;
END_RCPP
}
// get_path_list1
arma::field<arma::mat> get_path_list1(int N, int M, arma::mat range_x, arma::mat range_tem, arma::vec grid_t, arma::vec ind_end1, List grid_search, arma::mat P);
RcppExport SEXP _funcharts_get_path_list1(SEXP NSEXP, SEXP MSEXP, SEXP range_xSEXP, SEXP range_temSEXP, SEXP grid_tSEXP, SEXP ind_end1SEXP, SEXP grid_searchSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type range_x(range_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type range_tem(range_temSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type grid_t(grid_tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind_end1(ind_end1SEXP);
    Rcpp::traits::input_parameter< List >::type grid_search(grid_searchSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(get_path_list1(N, M, range_x, range_tem, grid_t, ind_end1, grid_search, P));
    return rcpp_result_gen;
END_RCPP
}
// get_path_list2
arma::field<arma::mat> get_path_list2(int N, int M, arma::mat range_x, arma::mat range_tem, arma::vec grid_t, arma::vec ind_end2, List grid_search, arma::mat P);
RcppExport SEXP _funcharts_get_path_list2(SEXP NSEXP, SEXP MSEXP, SEXP range_xSEXP, SEXP range_temSEXP, SEXP grid_tSEXP, SEXP ind_end2SEXP, SEXP grid_searchSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type range_x(range_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type range_tem(range_temSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type grid_t(grid_tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind_end2(ind_end2SEXP);
    Rcpp::traits::input_parameter< List >::type grid_search(grid_searchSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(get_path_list2(N, M, range_x, range_tem, grid_t, ind_end2, grid_search, P));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_funcharts_score", (DL_FUNC) &_funcharts_score, 4},
    {"_funcharts_score2", (DL_FUNC) &_funcharts_score2, 4},
    {"_funcharts_statisticY_EWMA_vec", (DL_FUNC) &_funcharts_statisticY_EWMA_vec, 5},
    {"_funcharts_statisticY_EWMA_cpp", (DL_FUNC) &_funcharts_statisticY_EWMA_cpp, 5},
    {"_funcharts_calculate_T2", (DL_FUNC) &_funcharts_calculate_T2, 3},
    {"_funcharts_calculate_T2_vec", (DL_FUNC) &_funcharts_calculate_T2_vec, 3},
    {"_funcharts_get_RL_cpp", (DL_FUNC) &_funcharts_get_RL_cpp, 10},
    {"_funcharts_der_c", (DL_FUNC) &_funcharts_der_c, 4},
    {"_funcharts_loss_c", (DL_FUNC) &_funcharts_loss_c, 6},
    {"_funcharts_DP3", (DL_FUNC) &_funcharts_DP3, 18},
    {"_funcharts_get_path_list1", (DL_FUNC) &_funcharts_get_path_list1, 8},
    {"_funcharts_get_path_list2", (DL_FUNC) &_funcharts_get_path_list2, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_funcharts(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
