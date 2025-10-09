#' @title  Permute the dataset by the largest absolute value in each eigenvector, and scale
#' @description It finds the vector of permutation to permute data by its largest absolute value in each eigenvector. It sets the order by specified number of factors, and the rest is ordered as they were. The data is permuted, and if needed, scaled.
#'
#' @param data  The dataset
#' @param permutation  Vector with the permutation of the data
#' @param covariance  Logic variable indicating whether the analysis uses covariance or correlation matrix
#' @param return.array  Return the data as 3-dimensional array
#' @param num.layers  Number of timepoints
#'
#' @return  The dataset that is permuted, either in matrix or array 
#' @export

permutation.scale <- function(data, permutation = NA, covariance = FALSE, return.array = TRUE, num.layers = 1){
  if(length(dim(data)) == 3){data.return <- apply(data, 2, c) ; num.layers = dim(data)[3]}
  if(length(dim(data)) == 2){data.return <- as.matrix(data)}
  if(class(data.return)[1]!="matrix"){data.return <- as.matrix(data.return)}
  r <- dim(data)[2]
  if(is.na(permutation[1])){
    permutation <- 1:r
  }  
  
  data.return <- data.return[,permutation]
  if(covariance == TRUE){data.return <- apply(data.return, 2, function(x) x - mean(x))}else{data.return <- apply(data.return, 2, scale)}
    
  data.2d <- data.return
  if(return.array){
    num.rows <- dim(data)[1]
    data.return <- array(NA, c(num.rows, r, num.layers))
    for(i in 1:num.layers){
      data.return[,,i] <- data.2d[seq(num.rows*(i-1)+1, num.rows*i),]
    }
  }
  return(data.return)
}