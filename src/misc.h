#ifndef MISC_H
#define MISC_H

#include <RcppArmadillo.h>


arma::vec table0(arma::vec v, arma::vec ctgry);
arma::vec cmean(arma::mat X);
arma::vec mattovec(arma::mat A);
arma::vec csum(arma::mat X);
arma::vec rsum(arma::mat X);
arma::vec mult_rgamma(double a, arma::vec b);
arma::mat rdirichlet_cpp(int num_samples, arma::vec alpha_m);
double dmvnorm(arma::vec X, arma::vec mu, arma::mat Sigma);
double logdmvnorm(arma::vec X, arma::vec mu, arma::mat Sigma);
arma::vec Qextract1(const Rcpp::List& attributes, const Rcpp::NumericVector& Q, int i);
arma::vec Qextract2(const Rcpp::List& attributes, const Rcpp::NumericVector& Q, int j);

#endif
