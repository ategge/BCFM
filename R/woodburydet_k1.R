#' @title Determinant by Woodbury lemma, but with one factor
#' @description Determinant by Woodbury lemma, but with one factor. It calculates faster than original method.
#'
#' @param v  Idiosyncratic variance
#' @param B  Factor loadings
#' @param Omega  Group covariance
#'
#' @return  Determinant by Woodbury lemma
#' @noRd

woodburydet_k1 <- function(v,B,Omega){
  detresult <- prod(v)*Omega*(1/Omega + t(B)%*%((1/v)*B))
  return(detresult)
}
