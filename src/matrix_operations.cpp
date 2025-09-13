// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Matrix-vector-matrix multiplication: x^T * A * B * y
// [[Rcpp::export]]
arma::vec vMMv_arma(arma::vec& x, arma::mat& A, arma::mat& B, arma::vec& y) {
  return x.t() * A * B * y;
}

// Matrix-vector multiplication: x^T * A * y  
// [[Rcpp::export]]
arma::vec vMv_arma(arma::vec& x, arma::mat& A, arma::vec& y) {
  return x.t() * A * y;
}

// Matrix inversion
// [[Rcpp::export]]
arma::mat invM_arma(const arma::mat& mat) {
  return arma::inv(mat);
}