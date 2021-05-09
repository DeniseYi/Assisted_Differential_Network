networkER <- function(p.){
  # Input  : p. - total dim
  # Output : a matrix with Erdos-Renyi structure
  A <- matrix(0,p.,p.)
  for (i in 1:p.){
    for (j in 1:p.) A[i,j] <- ifelse(rbinom(1, size=1, 0.3), ifelse(rbinom(1, size=1, 0.5), 
                                                                    runif(1, min = -0.4, max = -0.1), 
                                                                    runif(1, min = 0.1, max = 0.4)),0)
  }
  A <- (A+t(A))
  for (i in 1:p.) A[i,i] <- sum(abs(A[i,]))+0.1
  return(A)
}

networkBD <- function(p., rho.){ #this function needs to be updated according to the dim change
  # Input  : p.   - total dim
  #          rho. - the parameter to generate band correlation
  # Output : a matrix with banded structure
  A <- diag(0.5, p.)
  for (i in 2:p.){
    for (j in 1:(i-1)) A[i,j] <- rho.^(abs(i-j))
  }
  A <- A + t(A)
  return(A)
}

sigmaGenerator1 <- function(p.=50, block.=10){ 
  # Output: a covariance matrix with #block. blocks and dim #p.
  # submatrix is in ER structure
  mat = matrix(0,nrow=p., ncol=p.)
  sub.p = p./block.
  
  for(i in 1:block.){
    index1 = 1 + (i-1)*sub.p
    index2 = sub.p + (i-1)*sub.p
    mat[index1:index2 , index1:index2] <- networkER(sub.p)
  }
  print(is.positive.semi.definite(mat))
  return(mat)
}

sigmaGenerator2 <- function(p.=100, bs.=rep(c(5,10,10),4)){ 
  # Output: a covariance matrix with #block. blocks and dim #p.
  # submatrix is in ER structure
  mat = matrix(0,nrow=p., ncol=p.)
  
  index1 = index2 <- 0
  for(i in bs.){
    index1 = 1 + index2
    index2 = index2 + i
    mat[index1:index2 , index1:index2] <- networkER(i)
  }
  print(is.positive.semi.definite(mat))
  return(mat)
}

sigmaChange1 <- function(sigma., change., sub.p.){
  # Input : sigma.  - a matrix needs to be changed
  #         change. - a vector of the order of blocks that need to be changed
  #         sub.p.  - the size of the block
  # Output: a covariance modified matrix based on sigma.
  # submatrix is in ER structure
  for (i in change.){
    index1 = 1 + (i-1)*sub.p.
    index2 = sub.p. + (i-1)*sub.p.
    
    sigma.[index1:index2 , index1:index2] <- networkER(sub.p.)
  }
  print(is.positive.semi.definite(sigma.))
  return(sigma.)
}

sigmaChange2 <- function(sigma., change.){
  # Input : sigma.  - a matrix needs to be changed
  #         change. - a vector of the order of blocks that need to be changed
  #         sub.p.  - the size of the block
  # Output: a covariance modified matrix based on sigma.
  # submatrix is in ER structure
  sigma.[change., change.] <- networkER(length(change.))
  print(is.positive.semi.definite(sigma.))
  return(sigma.)
}

betaGenerator <- function(p.=50, q.=100, block., blockX., blockY., type.="A1C1"){
  # Output : a coeff beta matrix in Y = X %*% Beta + W. It is block diag, with value 1
  mat = matrix(0, nrow=q., ncol=p.)
  sub.p = p./block.
  sub.q = q./block.
  
  if (type. == "A1C1"){
    for (i in 1:block.){
      index1 = 1 + (i-1)*sub.p
      index2 = sub.p + (i-1)*sub.p
      index3 = 1 + (i-1)*sub.q
      index4 = sub.q + (i-1)*sub.q   
      mat[index3:index4 , index1:index2] <- runif(sub.q*sub.p, min = 0.9, max = 1)
    }
  }
  
  if (type. == "A2C1" | type. == "A3C1"){
    index1 = index2 = index3 = index4 <- 0
    for(i in 1:block.){
      index1 = 1 + index2
      index2 = index2 + blockX.[i]
      index3 = 1 + index4
      index4 = index4 + blockY.[i]
      mat[index1:index2 , index3:index4] <- runif(blockX.[i]*blockY.[i], min = 0.9, max = 1)
    }
  }
  
  if (type. ==""){
    
  }
  return(mat)
}

betaChange <- function(Beta., r.){
  p <- dim(Beta.)[1]
  q <- dim(Beta.)[2]
  mat <- matrix(nrow = p, ncol = q)
  for (i in 1:p){
    for (j in 1:q) mat[i,j] <- ifelse(rbinom(1, size=1, r.), runif(1, min = 0.9, max = 1),0)
  }
  res <- Beta.*(Beta.!=0) + mat*(Beta.==0)
  return(res)
}

noiseGenerator <- function(p.=50, bs. = rep(5,10), rho.=0.3){
  mat = matrix(0, nrow=p., ncol=p.)
  index1 = index2 <- 0
  for (i in bs.){
    index1 = 1 + index2
    index2 = index2 + i
    if (i>1){
      mat[index1:index2 , index1:index2] <- networkBD(i, rho.)
    }
    else{
      mat[index1:index2 , index1:index2] <- runif(1,0.9, 1.1)
    }
  }
  print(is.positive.semi.definite(mat))
  return(mat)
} # Revised noiseGenerator() and made it more general

