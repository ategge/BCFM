#include "misc.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' @title Simulate sigma squared
//' @param attributes Attributes of the model and the data
//' @param hyp_parm Hyperparameters
//' @param parm Current estimated parameters
//' @param data A 3-dimensional array dataset
//' @noRd
// [[Rcpp::export]]
arma::vec simulated_sigma(List attributes, List hyp_parm, List parm, arma::cube data){

  int S = attributes["S"];
  int times = attributes["times"];

  double n_sigma = hyp_parm["n.sigma"];
  double n_s2_sigma = hyp_parm["n.s2.sigma"];
  double n_aux = n_sigma + S*times;

  mat B = parm["B"];
  cube X = parm["X"];
  // vec next_sigma(R);

  mat X_matrix = X.slice(0);
  mat data_matrix = data.slice(0);
  for(int tt = 1; tt < times; tt++){
    X_matrix = join_cols(X_matrix, X.slice(tt));
    data_matrix = join_cols(data_matrix, data.slice(tt));
  }

  mat difference = data_matrix - X_matrix * trans(B);
  vec n_s2_aux = n_s2_sigma + csum(square(difference));
  vec next_sigma = 1 / mult_rgamma(n_aux/2, 2/n_s2_aux);

  return next_sigma;

}


