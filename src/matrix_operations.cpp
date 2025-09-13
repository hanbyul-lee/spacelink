// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Matrix-vector-matrix multiplication: x^T * A * B * y (returns scalar)
// [[Rcpp::export]]
double vMMv_arma(const arma::vec& x, const arma::mat& A, const arma::mat& B, const arma::vec& y) {
  arma::mat result = x.t() * A * B * y;
  return result(0,0);  // Extract scalar from 1x1 matrix
}

// Matrix-vector multiplication: x^T * A * y (returns scalar)
// [[Rcpp::export]]  
double vMv_arma(const arma::vec& x, const arma::mat& A, const arma::vec& y) {
  arma::mat result = x.t() * A * y;
  return result(0,0);  // Extract scalar from 1x1 matrix
}

// Matrix inversion
// [[Rcpp::export]]
arma::mat invM_arma(const arma::mat& mat) {
  return arma::inv(mat);
}