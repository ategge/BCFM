#' @title  Order of permutation by the largest absolute value in each eigenvector
#' @description It finds the vector of permutation to permute data by its largest absolute value in each eigenvector. It sets the order by specified number of factors, and the rest is ordered as they were.
#'
#' @param data  The dataset
#' @param covariance  Logic variable indicating whether the analysis uses covariance or correlation matrix
#' @param fa  Use factor analysis to sort the variables
#' 
#' @return  The vector of permutation
#' @export

permutation.order <- function(data, covariance = FALSE, L = 1, fa = TRUE){
  if(length(dim(data)) == 3){data.return <- apply(data, 2, c)}
  if(length(dim(data)) == 2){data.return <- data}
  if(class(data.return)[1]!="matrix"){data.return <- as.matrix(data.return)}
  
  if(covariance){ # covariance
    A <- eigen(cov(data.return))$vectors
  } else{ # correlation
    A <- eigen(cor(data.return))$vectors
  }
  r <- dim(data)[2]
  permutation = 1:r
  
  if(!fa){
    for(i in 1:(r-1))
    {
      indice = which.max(abs(A[,i]))
      if(i == 1)
      {
        permutation = c(permutation[indice],permutation[-indice])
        A = A[permutation,]
      }else
      {
        permutation = c(permutation[1:(i-1)], permutation[indice], permutation[-c(1:(i-1),indice)])
        remaining = i:r
        rem = which(remaining == indice)
        A[i:r,] = A[c(indice,remaining[-rem]),]
      }
      A[,i] = A[,i] / A[i,i]
      for(j in (i+1):r) A[,j] = A[,j] - A[i,j] * A[,i]
    }
  }
  
  if(fa){
    if(!covariance){
      data.return = scale(data.return)
    }
    A <- factanal(data.return, factors = L, nstart = 4, lower = 0.05)$loadings
#    A <- fa(data.return, nfactors = L, fm="minres", rotate = "none", covar=FALSE, warnings=FALSE)$loadings
    
    if(L == 1){
      indice = which.max(abs(A[,1]))
      permutation = c(permutation[indice],permutation[-indice])
      A = A[permutation,]
    }
    if(L > 1){
      for(i in 1:(L)){
        indice = which.max(abs(A[,i]))
        if(i == 1){
          permutation = c(permutation[indice],permutation[-indice])
          A = A[permutation,]
        }else{
          permutation = c(permutation[1:(i-1)], permutation[indice], permutation[-c(1:(i-1),indice)])
          remaining = i:r
          rem = which(remaining == indice)
          A[i:r,] = A[c(indice,remaining[-rem]),]
        }
        A[,i] = A[,i] / A[i,i]
        if (i<L){
          for(j in (i+1):L) A[,j] = A[,j] - A[i,j] * A[,i]
        }
      }
    }
  }
  
  
  return(permutation)
  
}

