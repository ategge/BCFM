#include "misc.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//' @title Simulate mus
//' @param attributes Attributes of the model and the data
//' @param hyp_parm Hyperparameters
//' @param parm Current estimated parameters
//' @noRd
// [[Rcpp::export]]
arma::mat simulated_mu(List attributes, List hyp_parm, List parm) {


  int G = attributes["G"];
  int L = attributes["L"];
  int times = attributes["times"];
  arma::cube X = parm["X"];
  NumericVector Z = parm["Z"];
  arma::mat Zmat = parm["Z"];
  vec Zmatvec = parm["Z"];
  IntegerVector Ztable = table(Z);
  arma::cube omega = parm["omega"];
  arma::cube mu_C = hyp_parm["mu.C"];
  arma::mat mu_m = hyp_parm["mu.m"];
  arma::mat munew(G,L);

  vec Ztable_vec = as<arma::vec>(wrap(Ztable));
  vec Zvec(G);
  Zvec(span(0, Ztable_vec.n_elem - 1)) = Ztable_vec;
  int Ztable_vec_n_elem = Ztable_vec.n_elem;
  if(Ztable_vec_n_elem < G){
    Zvec = table0(Zmatvec, linspace(1, G, G));
  }


  for(int l = 0; l < G; l++){
    int n_assigned = Zvec(l);
    if(n_assigned > 0){
      arma::mat X1 = X.slice(0);
      arma::vec Z1 = Zmat.col(0);
      arma::mat x_selected = X1.rows(arma::find(Z1 == (l+1)));

      if(times > 1){
        for(int tt = 1; tt < times; tt++){
          arma::mat Xnow = X.slice(tt);
          arma::vec Znow = Zmat.col(tt);
          x_selected = join_cols(x_selected, Xnow.rows(arma::find(Znow == (l+1))));
        }
      }

      arma::mat mu_Cnow = mu_C(arma::span(l,l), arma::span(0,L-1), arma::span(0,L-1));
      arma::mat omega_now = omega(arma::span(l,l), arma::span(0,L-1), arma::span(0,L-1));
      if(!mu_Cnow.is_symmetric()){mu_Cnow = symmatu(mu_Cnow);} // Force symmetry if not symmetric
      mat mu_C_inv = inv_sympd(mu_Cnow);
      if(!omega_now.is_symmetric()){omega_now = symmatu(omega_now);} // Force symmetry if not symmetric
      mat omega_inv = inv_sympd(omega_now);
      mat var_aux0 = mu_C_inv + n_assigned * omega_inv;
      if(!var_aux0.is_symmetric()){var_aux0 = symmatu(var_aux0);}
      arma::mat mu_C_star = inv_sympd(var_aux0);

      arma::colvec mu_m_star = arma::conv_to<arma::colvec>::from(mu_m.row(l));
      mu_m_star = mu_C_star * (mu_C_inv * mu_m_star + n_assigned * omega_inv * cmean(x_selected));

      munew.row(l) = arma::conv_to<arma::rowvec>::from(mu_m_star + chol(mu_C_star) * arma::randn(L));

    }else{
      arma::mat mu_C_star = mu_C(arma::span(l,l), arma::span(0,L-1), arma::span(0,L-1));
      arma::colvec mu_m_star = arma::conv_to<arma::colvec>::from(mu_m.row(l));
      //// arma::vec mu_m_star = mu_m.row(l);

      munew.row(l) = arma::conv_to<arma::rowvec>::from(mu_m_star + trans(chol(mu_C_star)) * arma::randn(L));
    }

  }

  return(munew);

}
