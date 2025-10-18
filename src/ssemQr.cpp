// src/ssemQr.cpp
#include <Rcpp.h>
#include <RcppEigen.h>
#include "MatOp.h"
#include "LR.h"
#include <vector>
#include <cmath>

using namespace Rcpp;
using Eigen::ArrayXi;
using Eigen::Map;
using Eigen::VectorXi;
using Eigen::MatrixXd;

// [[Rcpp::export]]
SEXP L2Regression(SEXP X_, SEXP Y_, SEXP S_, SEXP gamma_, SEXP lambda_, SEXP n_, SEXP p_, SEXP k_)
{
  BEGIN_RCPP
  // Convert inputs
  MatrixXd X = Rcpp::as<MatrixXd>(X_);
  MatrixXd Y = Rcpp::as<MatrixXd>(Y_);
  Rcpp::List Slist(S_);
  std::vector<ArrayXi> S;
  S.reserve(Slist.size());
  for (int i = 0; i < Slist.size(); ++i) {
    IntegerVector vi = Slist[i];
    // map IntegerVector to Eigen::VectorXi and then to ArrayXi
    Map<VectorXi> m(vi.begin(), vi.size());
    ArrayXi ai = m.array();
    S.push_back(ai);
  }

  const double gamma  = Rcpp::as<double>(gamma_);
  const double lambda = Rcpp::as<double>(lambda_);
  const int n = Rcpp::as<int>(n_);
  const int p = Rcpp::as<int>(p_);
  const int k = Rcpp::as<int>(k_);

  // outputs / temporaries
  MatrixXd B;
  std::vector<MatrixXd> f;
  MatrixXd F;
  MatrixXd mu;
  double error2 = 0.0, df = 0.0;

  B.resize(p, p);
  B.setZero();
  f.resize(p);

  // Call templated routine (explicit template args to help deduction)
  error2 = L2lr<MatrixXd, std::vector<MatrixXd>, std::vector<ArrayXi>>(
    X, Y, S, B, f, mu, gamma, lambda, n, p, k
  );

  df = static_cast<double>(n) * static_cast<double>(p) - 1.0;
  double sigma2 = error2 / df;

  // get_Fs expects the fi vector and S; returns a p x k matrix (MatType)
  F = get_Fs<MatrixXd, std::vector<ArrayXi>>(f, S, k).transpose();

  return Rcpp::List::create(
    Rcpp::Named("B") = B,
    Rcpp::Named("f") = f,
    Rcpp::Named("F") = F,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("sigma2") = sigma2
  );
  END_RCPP
}

// [[Rcpp::export]]
SEXP ObjL2Regression(SEXP X_, SEXP Y_, SEXP fit_)
{
  BEGIN_RCPP
  MatrixXd X = Rcpp::as<MatrixXd>(X_);
  MatrixXd Y = Rcpp::as<MatrixXd>(Y_);
  Rcpp::List fit(fit_);
  MatrixXd B = Rcpp::as<MatrixXd>(fit["B"]);
  MatrixXd F = Rcpp::as<MatrixXd>(fit["F"]);
  MatrixXd mu = Rcpp::as<MatrixXd>(fit["mu"]);
  double error2 = 0.0;
  error2 += LR_Objerr<MatrixXd>(X, Y, B, F, mu);
  double err = std::sqrt(error2);
  return Rcpp::wrap(err);
  END_RCPP
}

// [[Rcpp::export]]
SEXP L2lamax(SEXP X_, SEXP Y_, SEXP S_, SEXP n_, SEXP p_, SEXP k_)
{
  BEGIN_RCPP
  MatrixXd X = Rcpp::as<MatrixXd>(X_);
  MatrixXd Y = Rcpp::as<MatrixXd>(Y_);

  Rcpp::List Slist(S_);
  std::vector<ArrayXi> S;
  S.reserve(Slist.size());
  for (int i = 0; i < Slist.size(); ++i) {
    IntegerVector vi = Slist[i];
    Map<VectorXi> m(vi.begin(), vi.size());
    ArrayXi ai = m.array();
    S.push_back(ai);
  }

  const int n = Rcpp::as<int>(n_);
  const int p = Rcpp::as<int>(p_);
  const int k = Rcpp::as<int>(k_);
  double lambda = 0.0;
  lambda = L2lamax<MatrixXd, std::vector<ArrayXi>>(X, Y, S, n, p, k);
  return Rcpp::wrap(lambda);
  END_RCPP
}
