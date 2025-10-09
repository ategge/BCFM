
#include <RcppArmadillo.h>
#include "misc.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' @title Simulate Omegas
//' @param attributes Attributes of the model and the data
//' @param hyp_parm Hyperparameters
//' @param parm Current estimated parameters
//' @noRd
// [[Rcpp::export]]
arma::cube simulated_omega(List attributes, List hyp_parm, List parm) {

  // Extract model attributes
  int G = as<int>(attributes["G"]);
  int L = as<int>(attributes["L"]);
  int times = attributes["times"];
  arma::cube X = parm["X"];
  arma::mat mu = parm["mu"];
  arma::mat Z = parm["Z"];
  NumericVector Zmat = parm["Z"];
  IntegerVector Ztable = rep(0,G);

  // Extract hyperparameters
  arma::cube omega_psi = hyp_parm["omega.psi"];
  double omega_nu = hyp_parm["omega.nu"];
  double omega_diag_nu = hyp_parm["omega.diag.nu"];

  // Initialize Omega sum
  arma::cube Omegasum(L, L, G);
  // Store updated Omega
  arma::cube Omeganew(G, L, L);

  for(int l = 0; l < G; l++){
    Ztable(l) = sum(Zmat == (l+1));
  }


  mat X_matrix = X.slice(0);
  for(int tt = 1; tt < times; tt++){
    X_matrix = join_cols(X_matrix, X.slice(tt));
  }
  vec Zvec = mattovec(Z);

  // Simulate Omega for the first group (diagonal matrix)
  int g = 0;  // C++ index starts from 0
  mat X_matrix1 = X_matrix.rows(find(Zvec == (g + 1)));// Find indices for group g+1

  mat t_difference1 = X_matrix1.each_row() - conv_to<rowvec>::from(mu.row(0));

  vec sum_squares = csum(square(t_difference1));

  mat omega_psi_g = omega_psi(arma::span(g,g), arma::span(0,L-1), arma::span(0,L-1));

  if (L == 1) {

    Omeganew(arma::span(g,g), arma::span(0,L-1), arma::span(0,L-1)) =  1 / mult_rgamma(omega_diag_nu + Ztable(g) / 2.0,
             1/(omega_diag_nu * omega_psi_g + sum_squares / 2.0));
  } else {
    Omeganew(arma::span(g,g), arma::span(0,L-1), arma::span(0,L-1)) = diagmat(1 / mult_rgamma(omega_diag_nu + Ztable(g) / 2.0,
                        1/(omega_diag_nu * omega_psi_g.diag() + sum_squares / 2.0)));
  }

  // Simulate Omega for remaining groups
  for (g = 1; g < G; g++) {

    mat X_matrixg = X_matrix.rows(find(Zvec == (g + 1)));

    mat t_difference = X_matrixg.each_row() - conv_to<rowvec>::from(mu.row(g));

    Omegasum.slice(g) = trans(t_difference) * t_difference;

    mat omega_psi_g = omega_psi(arma::span(g,g), arma::span(0,L-1), arma::span(0,L-1));

    if ( G > 1 && L > 1) {
      Omeganew(arma::span(g,g), arma::span(0,L-1), arma::span(0,L-1)) = iwishrnd(omega_psi_g + Omegasum.slice(g), Ztable(g) + omega_nu);
    }
    if ( G > 1 && L == 1) {
      Omeganew(arma::span(g,g), arma::span(0,L-1), arma::span(0,L-1)) = 1 / mult_rgamma(omega_nu + Ztable(g) / 2.0,
               1/(omega_nu * omega_psi_g + Omegasum.slice(g) / 2.0));
    }
  }

  // Return updated parameters
  return Omeganew;
}
