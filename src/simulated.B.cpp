#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' @title Simulate Bs
//' @param attributes Attributes of the model and the data
//' @param parm Current estimated parameters
//' @param data Dataset of a 3-dimensional array
//' @noRd
// [[Rcpp::export]]

arma::mat simulated_B(List attributes, List parm, arma::cube data) {
  
  // int S = attributes["S"];
  int R = attributes["R"];
  int L = attributes["L"];
  int times = attributes["times"];
  // NumericVector X = parm["X"];
  arma::cube X = parm["X"];
  arma::mat nextB(R,L);
  nextB.eye(R,L);
  arma::vec sigma2 = parm["sigma2"];
  arma::vec tau = parm["tau"];
  arma::mat XtX(L,L);
  arma::mat XtY(L,R);
  arma::mat identityL;
  identityL.eye(L,L);
  arma::mat temp;
  
  for(int tt = 0; tt < times; tt++){
    // arma::mat Xtt(S,L);
    arma::mat datat = data.slice(tt);
    arma::mat Xtt = X.slice(tt);
    XtX += Xtt.t() * Xtt;
    XtY += Xtt.t() * datat;
  }
  
  for(int j = L; j < R; j++){
    arma::mat var_aux0 = XtX / sigma2(j) + diagmat(1/tau);
    if(!var_aux0.is_symmetric()){var_aux0 = symmatu(var_aux0);} // Force symmetry if not symmetric
    arma::mat var_aux = inv_sympd(var_aux0);
    arma::vec mean_aux = var_aux * XtY.col(j)/ sigma2(j);
    arma::mat chol_var_aux = trans(chol(var_aux));
    arma::vec randnorm = rnorm(L);
    nextB.row(j) = trans(mean_aux + chol_var_aux * randnorm);
  }
  
  // Simulate values from 2nd to Lth row
  if(L > 1){
    for(int j = 1; j < L; j++){
      
      arma::mat XtXlim(j,j);
      arma::mat XtYlim(j,1);
      arma::vec nseq = arma::linspace(0, j-1);
      for(int tt = 0; tt < times; tt++){
        arma::mat Xcurrent = X.slice(tt);
        arma::mat datat = data.slice(tt);
        XtXlim += trans(Xcurrent.cols(0, j-1)) * Xcurrent.cols(0, j-1);
        XtYlim += trans(Xcurrent.cols(0, j-1)) * (datat.col(j) - Xcurrent.col(j));
      }
      
      arma::mat var_aux0 = XtXlim/sigma2(j) + arma::diagmat(1/tau.subvec(0,j-1));
      if(!var_aux0.is_symmetric()){var_aux0 = symmatu(var_aux0);} // Force symmetry if not symmetric
      arma::mat var_aux = arma::inv_sympd(var_aux0);
      arma::mat mean_aux = var_aux * XtYlim / sigma2(j);
      arma::mat chol_var_aux = chol(var_aux);
      arma::vec randnorm = rnorm(j);
      nextB(j, arma::span(0,j-1)) = conv_to<rowvec>::from(mean_aux) + conv_to<rowvec>::from(chol_var_aux.t() * randnorm);
    }
  }
  
  // Simulate values from L+1th row to Rth row
  
  for(int j = L; j < R; j++){
    
    arma::mat var_aux0 = XtX / sigma2(j) + diagmat(1/tau);
    if(!var_aux0.is_symmetric()){var_aux0 = symmatu(var_aux0);} // Force symmetry if not symmetric
    mat var_aux = inv_sympd(var_aux0);
    vec mean_aux = var_aux * XtY.col(j) / sigma2(j);
    mat chol_var_aux = chol(var_aux);
    vec randnorm = rnorm(L);
    nextB.row(j) = conv_to<rowvec>::from(mean_aux) + conv_to<rowvec>::from(chol_var_aux.t() * randnorm);
  }
  
  
  return nextB;
  
}
