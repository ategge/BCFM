#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @title Simulate taus
//' @param attributes Attributes of the model and the data
//' @param hypparm Hyperparameters
//' @param parm Current estimated parameters
//' @noRd
// [[Rcpp::export]]
NumericVector simulated_tau(List attributes, List hypparm, List parm) {
  
  int L = attributes["L"];
  int R = attributes["R"];
  double n_tau = hypparm["n.tau"];
  double n_s2_tau = hypparm["n.s2.tau"];
  double n_aux = n_tau;
  double n_s2_aux = n_s2_tau;
  NumericMatrix B = parm["B"];
  NumericVector next_tau(L);
  
  for(int i = 0; i < L; i++){
    n_aux = n_tau + R - i - 1;
    n_s2_aux = n_s2_tau;
    for(int j = i+1; j < R; j++){
      double b_temp = B(j,i);
      n_s2_aux += pow(b_temp, 2.0);
    }
    NumericVector tau_temp = rgamma(1, n_aux/2, 2/n_s2_aux);
    next_tau(i) = 1/tau_temp(0);
  }
  
  return next_tau;
}
