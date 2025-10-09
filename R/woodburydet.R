#' @title Determinant by Woodbury lemma
#' @description Determinant by Woodbury lemma. It calculates faster than original method.
#'
#' @param v  Idiosyncratic variance
#' @param B  Factor loadings
#' @param Omega  Group covariance
#'
#' @return  Determinant by Woodbury lemma
#' @noRd

woodburydet <- function(v,B,Omega){
  detresult <- prod(v)*det(Omega)*det(solve(Omega) + t(B)%*%((1/v)*B))
  return(detresult)
}