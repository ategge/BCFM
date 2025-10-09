#' @title  Get the mode of a vector
#' @description  The function returns the mode of a vector. 
#'
#' @param v  The vector to find the mode.
#'
#' @return  The mode of \code{v}.
#' @export

getmode <- function(v){
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v,uniqv)))]
}