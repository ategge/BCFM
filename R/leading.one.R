
#' @title  A transformation matrix that turns a matrix satisfying structural hierarchical constraint
#' @description A transformation matrix that turns a matrix satisfying structural hierarchical constraint. Needs to be computed when deriving the cluster hyperparmeters. 
#' @param X  A matrix to be transformed in hierarchical constraint. It comes from a eigenvectors of a dataset after the variables are sorted.
#' @noRd
leading.one <- function(X){
  
  Xmat <- X
  k <- ncol(X)
  l <- 1
  # diag matrix of 1s. dim ncol(X) by ncol(X)
  Imat <- diag(ncol(X))
  Imat.list <- vector("list", length = sum(1:ncol(X)))
  for(i in 1:k){
    for(j in 1:i){
      if(j == i){Imat.list[[l]] <- Imat
      Imat.list[[l]][j,j] <- 1/Xmat[j,j]
      Xmat <- Xmat %*% Imat.list[[l]]
      l <- l + 1
      }
      if(j < i){
        Imat.list[[l]] <- Imat
        Imat.list[[l]][j,i] <- -Xmat[j,i]
        Xmat <- Xmat %*% Imat.list[[l]]
        l <- l + 1
      }
    }
  }
  
  result <- Imat
  for(i in 1:length(Imat.list)){result <- result %*% Imat.list[[i]]}
  return(result)
}

