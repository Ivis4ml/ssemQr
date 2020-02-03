#include "MatOp.h"
#include "LR.h"

RcppExport SEXP L2Regression(SEXP X_, SEXP Y_, SEXP S_, SEXP gamma_, SEXP n_, SEXP p_, SEXP k_)
{
BEGIN_RCPP
  MatrixXf X = Rcpp::as<MatrixXf>(X_);
  MatrixXf Y = Rcpp::as<MatrixXf>(Y_);
  std::vector<ArrayXd> S = Rcpp::as<std::vector<ArrayXd> >(S_);
  const double gamma  = Rcpp::as<double>(gamma_);
  const int n = Rcpp::as<int>(n_);
  const int p = Rcpp::as<int>(p_);
  const int k = Rcpp::as<int>(k_);
  MatrixXf B;
  std::vector<MatrixXf> f;
  MatrixXf F;
  MatrixXf mu;
  double error2 = 0, df = 0;
  B.resize(p, p);
  B.setZero();
  f.resize(p);
  error2 = L2lr(X, Y, S, B, f, mu, gamma, n, p, k);
  df = n * p - 1;
  double sigma2 = error2 / df;
  F = get_Fs(f, S, k).transpose();
  return Rcpp::List::create(
    Rcpp::Named("B") = B,
    Rcpp::Named("f") = f,
    Rcpp::Named("F") = F,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("sigma2") = sigma2
  );
END_RCPP
}

RcppExport SEXP ObjL2Regression(SEXP X_, SEXP Y_, SEXP fit_)
{
BEGIN_RCPP
  MatrixXf X = Rcpp::as<MatrixXf>(X_);
  MatrixXf Y = Rcpp::as<MatrixXf>(Y_);
  Rcpp::List fit(fit_);
  MatrixXf B = Rcpp::as<MatrixXf>(fit["B"]);
  MatrixXf F = Rcpp::as<MatrixXf>(fit["F"]);
  MatrixXf mu = Rcpp::as<MatrixXf>(fit["mu"]);
  double error2 = 0;
  error2 += LR_Objerr(X, Y, B, F, mu);
  double err = std::sqrt(error2);
  return Rcpp::wrap(err);
  //return Rcpp::wrap(error2);
END_RCPP
}

RcppExport SEXP L2lamax(SEXP X_, SEXP Y_, SEXP S_, SEXP n_, SEXP p_, SEXP k_)
{
BEGIN_RCPP
  MatrixXf X = Rcpp::as<MatrixXf>(X_);
  MatrixXf Y = Rcpp::as<MatrixXf>(Y_);
  std::vector<ArrayXd>  S  = Rcpp::as<std::vector<ArrayXd> >(S_);
  const int n = Rcpp::as<int>(n_);
  const int p = Rcpp::as<int>(p_);
  const int k = Rcpp::as<int>(k_);
  double lambda = 0;
  lambda = L2lamax(X, Y, S, n, p, k);
  return Rcpp::wrap(lambda);
END_RCPP
}

