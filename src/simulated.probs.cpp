#include "misc.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
//' @title Simulate probabilities
//' @param attributes Attributes of the model and the data
//' @param parm Current estimated parameters
//' @param hyp_parm Hyperparameters
//' @noRd
// [[Rcpp::export]]
arma::vec simulated_probs(List attributes, List parm, List hyp_parm){

  int G = attributes["G"];
  NumericVector Z = parm["Z"];
  vec Zvec = parm["Z"];
  IntegerVector Ztable = table(Z);
  vec p_exponent_G = hyp_parm["p.exponent"];
  vec prob_aux(G);
  mat pnewmat;
  arma::vec pnew;
  arma::vec prob_aux_table = as<arma::vec>(wrap(Ztable));
  prob_aux(span(0,prob_aux_table.n_elem - 1)) = prob_aux_table;
  int prob_aux_n = prob_aux_table.n_elem;
  if(prob_aux_n < G){
    prob_aux = table0(Zvec, linspace(1, G, G));
  }

  pnewmat = rdirichlet_cpp(1, prob_aux + p_exponent_G);
  pnew = arma::conv_to<arma::vec>::from(pnewmat);

  return pnew;

}
