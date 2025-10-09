
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//' @title Simulate Xs
//' @param attributes Attributes of the model and the data
//' @param parm Current estimated parameters
//' @param data A 3-dimensional array dataset
//' @noRd
// [[Rcpp::export]]
arma::cube simulated_X(List attributes, List parm, arma::cube data){
  
  int G = attributes["G"];
  int L = attributes["L"];
  int times = attributes["times"];
  int S = attributes["S"];
  int R = attributes["R"];
  arma::mat Z = parm["Z"];
  arma::mat B = parm["B"];
  arma::mat mu = parm["mu"];
  arma::cube Omega = parm["omega"];
  arma::vec sigma2 = parm["sigma2"];
  arma::mat SigmaV_inv = diagmat(1./sigma2);
  arma::cube Xnew(S, L, times);
  arma::cube temp;
  
  cube Omegainv(L, L, G);
  cube var_aux(L, L, G);
  cube chol_var_aux(L, L, G);
  
  for(int j = 0; j < G; j++){
    mat OmegaG = Omega(arma::span(j,j), arma::span(0,L-1), arma::span(0, L-1));
    if(!OmegaG.is_symmetric()){OmegaG = symmatu(OmegaG);} // Force symmetric if not symmetric
    Omegainv.slice(j) = inv_sympd(OmegaG);
    mat var_aux0 = Omegainv.slice(j) + trans(B) * SigmaV_inv * B;
    if(!var_aux0.is_symmetric()){var_aux0 = symmatu(var_aux0);}
    var_aux.slice(j) = inv_sympd(var_aux0);
    chol_var_aux.slice(j) = chol(var_aux.slice(j));
  }
  
  for(int tt = 0; tt < times; tt++){
    for(int i = 0; i < S; i++){
      mat data_now = data(arma::span(i,i), arma::span(0,R-1), arma::span(tt,tt));
      colvec data_nowcol = arma::conv_to<arma::colvec>::from(data_now);
      int Znow = Z(i,tt) - 1;
      
      colvec muselect = conv_to<colvec>::from(mu.row(Znow));
      mat mean_aux = var_aux.slice(Znow) * (trans(B) * SigmaV_inv * data_nowcol + Omegainv.slice(Znow) * muselect);
      Xnew(arma::span(i,i), arma::span(0,L-1), arma::span(tt,tt)) = mean_aux + trans(chol_var_aux.slice(Znow)) * arma::randn(L);
    }
  }
  
  return Xnew;
  
}
