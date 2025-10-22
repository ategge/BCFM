#include "misc.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
//' @title Simulate Zs
//' @param attributes Attributes of the model and the data
//' @param parm Current estimated parameters
//' @noRd
// [[Rcpp::export]]
arma::mat simulated_Z(List attributes, List parm){

  int G = attributes["G"];
  int times = attributes["times"];
  int S = attributes["S"];
  int L = attributes["L"];
  arma::cube X = parm["X"];
  arma::cube Omega = parm["omega"];
  arma::mat mu = parm["mu"];
  arma::vec probs = parm["probs"];
  arma::vec logdx(G);
  NumericVector onetoG = as<NumericVector>(wrap(arma::linspace(1, G, G)));
  arma::mat newZ(S, times);

  for(int tt = 0; tt < times; tt++){
    for(int i = 0; i < S; i++){
      arma::mat Xt = X.slice(tt);
      arma::vec Xnow = arma::conv_to<arma::vec>::from(Xt.row(i));
      for(int l = 0; l < G; l++){
        if(L > 1){
          arma::vec mu_l = arma::conv_to<arma::vec>::from(mu.row(l));
          arma::mat Omega_l = Omega(arma::span(l,l), arma::span(0,L-1), arma::span(0,L-1));
          logdx(l) = logdmvnorm(Xnow, mu_l, Omega_l);
        }
        if(L == 1){
          double Xnow_d = arma::as_scalar(Xt.row(i));
          double mu_l_d = arma::as_scalar(mu.row(l));
          double Omega_l_d = arma::as_scalar(Omega(arma::span(l,l), arma::span(0,L-1), arma::span(0,L-1)));
          logdx(l) = R::dnorm(Xnow_d, mu_l_d, Omega_l_d, TRUE);
        }
      }
      arma::vec logZit_aux = logdx + log(probs);
      double maxlogZit_aux = max(logZit_aux);
      arma::vec Zit_aux = exp(logZit_aux-maxlogZit_aux);
      Zit_aux = Zit_aux / sum(Zit_aux);
      NumericVector Zsamp = Rcpp::sample(onetoG, 1, false, as<NumericVector>(wrap(Zit_aux)));
      newZ(i, tt) = Zsamp(0);
    }
  }

  return newZ;

}
