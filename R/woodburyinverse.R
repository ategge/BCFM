#' @title  Inverse by Woodbury lemma
#' @description Inverse by Woodbury lemma. It calculates faster than original method.
#'
#' @param v  Idiosyncratic variances
#' @param B  Factor loadings
#' @param Omega  Group covariance
#'
#' @return  Inverse by Woodbury lemma
#' @noRd

woodburyinverse <- function(v,B,Omega){
  O <- diag(1/v) - ((1/v)*B)%*%solve(solve(Omega) + t(B)%*%((1/v)*B))%*%(t(B)%*%diag(1/v))
  return(O)
}