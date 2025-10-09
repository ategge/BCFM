#' @title  Simulate tau squared, variance of factor loadings
#' @description  Simulate tau squared, variance of factor loadings. We assume factors are independent, therefore has a inverse gamma posterior.
#'
#' @param model.attributes  Model attributes from \code{initialize.model.attributes}
#' @param hyp.parm  Hyperparameters from \code{initialize.hyp.parm}
#' @param parm  Parameters of current iteration
#'
#' @importFrom stats rgamma
#'
#' @return  A updated \code{parm} with new tau squared
#' @noRd

simulated.tau <- function(model.attributes, hyp.parm, parm){
  n.aux = hyp.parm$n.tau + model.attributes$R - c(1:model.attributes$L)
  # the minus 1 is to remove the main diagonal elements that are fixed at 1
  n.s2.aux = hyp.parm$n.s2.tau + colSums(parm$B^2) - 1
  parm$tau2 = 1 / rgamma(model.attributes$L,shape=n.aux/2,rate = n.s2.aux/2)
  return(parm)
}
