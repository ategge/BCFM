#include <RcppArmadillo.h>
#include "misc.h"

using namespace Rcpp;
using namespace arma;

arma::vec table0(arma::vec v, arma::vec ctgry){
    int n_ctgry = ctgry.n_elem;
    arma::vec counts(n_ctgry);
    for(int i = 0; i < n_ctgry; i++){
      counts(i) = sum(v == ctgry(i));
    }
    return(counts);
  }

arma::vec cmean(arma::mat X){
    int nCols = X.n_cols;
    arma::vec out(nCols);
    for(int i = 0; i < nCols; i++){
      out(i) = mean(X.col(i));
    }
    return out;
}

vec mattovec(mat A){

  int nrow = A.n_rows;
  int ncol = A.n_cols;
  vec a(nrow * ncol);

  for(int i = 0; i < ncol; i++){
    a(span(nrow * i, nrow *i + nrow - 1)) = A.col(i);
  }

  return a;
}

arma::vec csum(arma::mat X){
  int nCols = X.n_cols;
  arma::vec out(nCols);
  for(int i = 0; i < nCols; i++){
    out(i) = accu(X.col(i));
  }
  return out;
}

arma::vec rsum(arma::mat X){
  int nRows = X.n_rows;
  arma::vec out(nRows);
  for(int i = 0; i < nRows; i++){
    out(i) = accu(X.row(i));
  }
  return out;
}

vec mult_rgamma(double a, vec b){

  int R = b.n_elem;
  vec outvec(R);
  for(int i = 0; i < R; i++){
    double bnow = b(i);
    Rcpp::NumericVector nextval = Rcpp::rgamma(1, a, bnow);
    outvec(i) = nextval(0);
  }

  return outvec;
}

arma::mat rdirichlet_cpp(int num_samples,
                         arma::vec alpha_m) {
  int distribution_size = alpha_m.n_elem;
  // each row will be a draw from a Dirichlet
  arma::mat distribution = arma::zeros(num_samples, distribution_size);

  for (int i = 0; i < num_samples; ++i) {
    double sum_term = 0;
    // loop through the distribution and draw Gamma variables
    for (int j = 0; j < distribution_size; ++j) {
      double cur = R::rgamma(alpha_m[j],1.0);
      distribution(i,j) = cur;
      sum_term += cur;
    }
    // now normalize
    for (int j = 0; j < distribution_size; ++j) {
      distribution(i,j) = distribution(i,j)/sum_term;
    }
  }
  return(distribution); //outputs a row-wise matrix
} // Got the function from: https://www.mjdenny.com/blog.html

double dmvnorm(arma::vec X, arma::vec mu, arma::mat Sigma){

  arma::vec Yvec;
  arma::colvec Xmat = arma::conv_to<arma::colvec>::from(X);
  arma::colvec mumat = arma::conv_to<arma::colvec>::from(mu);
  double n = Sigma.n_rows;

  Yvec = pow(2 * M_PI, -n/2) * pow(arma::det(Sigma), -0.5) * exp(-0.5 * arma::trans(Xmat - mumat) * arma::inv_sympd(Sigma) * (Xmat - mumat));
  double Y = arma::conv_to<double>::from(Yvec);
  return Y;

}

double logdmvnorm(arma::vec X, arma::vec mu, arma::mat Sigma){

  arma::vec Yvec;
  arma::colvec Xmat = arma::conv_to<arma::colvec>::from(X);
  arma::colvec mumat = arma::conv_to<arma::colvec>::from(mu);
  double n = Sigma.n_rows;
  if(!Sigma.is_symmetric()){Sigma = symmatu(Sigma);} // Force symmetry if not symmetric
  double logdens = arma::conv_to<double>::from((-n/2)*log(2 * M_PI) - 0.5 * (arma::log_det(Sigma)) - 0.5 * (arma::trans(Xmat - mumat) * arma::inv_sympd(Sigma) * (Xmat - mumat)));
  return logdens;

}

// Function to extract row i from Q
arma::vec Qextract1(const List& attributes, const NumericVector& Q, int i) {
  int G = attributes["G"];
  arma::vec q(G);
  for (int j = 0; j < G; ++j) {
    q(j) = Q[i * G + j];
  }
  return q;
}

// Function to extract column j from Q
arma::vec Qextract2(const List& attributes, const NumericVector& Q, int j) {
  int G = attributes["G"];
  arma::vec q(G);
  for (int i = 0; i < G; ++i) {
    q(i) = Q[i * G + j];
  }
  return q;
}
