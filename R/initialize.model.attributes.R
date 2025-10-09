#' @title  Build model attributes from the dataset
#' @description  Basic setting and information of the dataset: number of subjects, timepoints, variables, factors and groups.
#'
#' @param S  Number of subjects
#' @param times  Number of timepoints
#' @param R  Number of covariates
#' @param L  Number of factors
#' @param G  Number of groups
#'
#' @return  A list of model attributes
#' @export

initialize.model.attributes <- function(S=216,times=5,R=9,L=3,G=4){

  model.attributes = NULL
  
  # number of participants
  model.attributes$S = S
  # number of timepoints
  model.attributes$times = times
  # number of variables
  model.attributes$R = R
  # Number of factors
  model.attributes$L = L
  # Number of groups
  model.attributes$G = G
  
  return(model.attributes)
}