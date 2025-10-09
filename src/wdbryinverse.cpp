#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

 //' @title Compute the Woodbury inverse of B, sigma squared and Omega
 //' @param sigma2  Vector of idiosyncratic variances
 //' @param Omega  Matrix of the cluster covariance
 //' @param B  Matrix of factor loadings
 //' @noRd
 // [[Rcpp::export]]

arma::mat wdbryinverse(arma::vec sigma2, arma::mat B, arma::mat Omega) {
  
  vec vinv = 1./sigma2;
  mat vinvB = B;
  vinvB.each_col() %= vinv;
  mat Aux = (-vinvB) * inv_sympd(inv_sympd(Omega) + B.t() * vinvB) * vinvB.t();
  Aux.diag() += vinv;
  return Aux;
  
}


