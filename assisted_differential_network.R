penSpMcp <- function(x,lambda,gamma){
  x        <- abs(x)
  sign.1   <- 1*(x<(lambda*gamma))
  val.1    <- lambda*x -x^2/(2*gamma)
  sign.2   <- 1-sign.1
  val.2    <- (lambda*gamma)^2 /(2*gamma)
  rlt      <- val.1 * sign.1 + val.2* sign.2
  return(rlt)
}

penSpMcp.Deriv <- function(x,lambda,gamma){
  x         <- abs(x)
  sign.1    <- 1*(x<(lambda*gamma))
  val.1     <- lambda -x/gamma
  rlt       <- val.1 * sign.1 
  return(rlt)
}


NormDiff <- function(xnew,xold,measure="Fnorm2",type="R"){
  
  if (measure=="Fnorm2"){
    num1 <- sum((xnew-xold)^2)
    num2 <- max(c(1, sum(xnew^2), sum(xold^2) ) )
  }
  if (measure=="Fnorm"){
    num1 <-  sqrt(sum((xnew-xold)^2))
    num2 <- max(c(1, sqrt(sum(xnew^2)), sqrt(sum(xold^2))))
  }
  if (measure=="abs"){
    num1 <-  sum(abs(xnew-xold))
    num2 <- sum(abs(xold))
  }
  val <- ifelse(type=="R",num1/num2,num1)
  return(val)
}

ssvd_fxn <- function(X., lambda., maxIter.=1000, tol.=10^-6, type.){
  # This is the function for Alt.3 - the bi-clustering sparse SVD approach
  # Step1: Initialization. Apply the standard SVD
  res1 <- svd(X.)
  u <- res1$u[,1] # sqrt(res1$d[1]) * 
  v <- res1$v[,1] #res1$v[,1] # sqrt(res1$d[1]) * 
  
  iter <- 1
  diff <- c(10^6)
  
  # Step2: Update
  while(diff[iter] > tol. & iter < maxIter.){
    
    temp1 <- abs(t(X.) %*% u)-lambda./2
    v.new <- sign(t(X.) %*% u) * temp1 * (temp1>0)
    s.new <- sqrt(sum(v.new^2)) 
    if(sum(s.new^2)!=0){
      v.new <- v.new/s.new
    } 
    
    temp2 <- abs(X. %*% v.new)-lambda./2
    u.new <- sign(X. %*% v.new) * temp2 * (temp2>0)
    s.new <- sqrt(sum(u.new^2)) 
    if(sum(s.new^2)!=0){
      u.new <- u.new/s.new
    } 
    
    diff <- c(diff, max(NormDiff(v,v.new), NormDiff(u,u.new)))
    
    u <- u.new
    v <- v.new
    iter <- iter + 1
  }
  if (iter == maxIter.) {
    print(c(type.,lambda.))
    return(NA)
  }
  else return(list(u=u, v=v, diff=diff))
}

sssvd_fxn <- function(G., R., lbd1., lbd2., lbd3., cor., type.="pre",
                      tol.=10^-6, maxIter.=1000, a.=3, tau.=0.2){
  # This is the main function for the proposed assisted differential network analysis
  # sssvd_fxn() is a function to estimate the left singular vectors of G. and R. using 
  # the proposed penalization method
  # G.      : Diff Network of GE 
  # R.      : Diff Network of the regulator
  # lbd1.   : tunining parameter lambda_1
  # lbd2.   : tunining parameter lambda_2
  # lbd3.   : tunining parameter lambda_3
  # a., tau.: other parameters
  # tol.    : convergent criterion
  # maxIter.: the max number of iterations
  
  # Initialization 
  p = dim(G.)[1]
  q = dim(R.)[1]
  s1 <- svd(G.)
  v  <- s1$u[,1]
  w1 <- s1$v[,1]
  
  s2 <- svd(R.)
  u  <- s2$u[,1]
  w2 <- s2$v[,1]
  
  diff <- 10^5
  obj  <- sum((G.- v %*% t(w1))^2) + sum((R.-u %*% t(w2))^2) + sum(penSpMcp(v, lbd1., a.)) +
    sum(penSpMcp(w1, lbd1., a.)) + sum(penSpMcp(u, lbd2., a.)) + 
    sum(penSpMcp(w2, lbd2., a.)) - drop(lbd3.* t(v!=0) %*% cor. %*% (u!=0))
  iter = 1
  
  # Iteration
  # v.result = u.result = matrix(nrow = maxIter, ncol = p)
  while(diff[iter] > tol. & iter < maxIter.){
    v.pre = v
    w1.pre = w1
    u.pre = u
    w2.pre = w2
    
    # update v, w1
    D2 = matrix(0, nrow = p, ncol = p)
    diag(D2) <- exp(-(v.pre)^2/tau.)*2*v.pre/tau.
    temp1 = 1-exp(-u^2 /tau.)
    temp2 = abs( 2* G.%*%w1.pre + lbd3.* D2 %*% cor. %*% temp1)- penSpMcp.Deriv(v.pre, lbd1., a.)
    v = 0.5 * sign(2*G. %*% w1.pre + lbd3. * D2 %*% cor. %*% temp1) * (temp2*(temp2>0))
    if (sum(v^2) != 0 & sum(w1^2) != 0) {
      v = v/sqrt(sqrt(sum(v^2))*sqrt(sum(w1.pre^2)))
    } 
    
    
    temp5 = abs(2*G. %*% v) - penSpMcp.Deriv(w1.pre, lbd1., a.)
    w1 = 0.5 * sign(G. %*% v)* temp5 * (temp5>0)
    if (sum(w1^2) != 0 & sum(v^2) != 0) {
      w1 = w1/sqrt(sqrt(sum(v^2))*sqrt(sum(w1^2)))
    } 
    
    # update u, w2
    D4 = matrix(0, nrow = q, ncol = q)
    diag(D4) <- exp(-(u.pre)^2/tau.)*2*u.pre/tau.
    temp3 = 1-exp(-v^2 /tau.)
    temp4 = abs(2* R. %*% w2.pre + lbd3.* t(D4) %*% t(cor.) %*% temp3) - penSpMcp.Deriv(u.pre, lbd2., a.)
    #temp4 = abs(2* R. %*% w2.pre + lbd3.* D4 %*% cor. %*% temp3) - penSpMcp.Deriv(u.pre, lbd2., a.)
    u = 0.5 * sign(2* R. %*% w2.pre + lbd3. * t(D4) %*% t(cor.) %*% temp3) * (temp4*(temp4>0))
    if (sum(u^2) != 0 & sum(w2.pre^2) !=0) {
      u = u/sqrt(sqrt(sum(u^2))*sqrt(sum(w2.pre^2)))
    } 
    
    
    temp6 = abs(2*R. %*% u) - penSpMcp.Deriv(w2.pre, lbd1., a.)
    w2 = 0.5 * sign(R. %*% u)* temp6 * (temp6>0)
    if (sum(w2^2) != 0 & sum(u^2) != 0) {
      w2 = w2/sqrt(sqrt(sum(u^2))*sqrt(sum(w2^2)))
    } 
    
    # update objective func
    obj.new <- sum((G.- v %*% t(w1))^2) + sum((R.-u %*% t(w2))^2) + sum(penSpMcp(v, lbd1., a.)) +
      sum(penSpMcp(w1, lbd1., a.)) + sum(penSpMcp(u, lbd2., a.)) + 
      sum(penSpMcp(w2, lbd2., a.)) - drop(lbd3.* t(v!=0) %*% cor. %*% (u!=0))
    obj <- c(obj, obj.new)
    
    # save result
    #v.result[iter,] = v
    #u.result[iter,] = u
    
    # Check convergence
    diff = c(diff,max(NormDiff(v,v.pre), NormDiff(u,u.pre)))
    iter = iter+1
  }
  if (iter == maxIter.) {
    return(NA)
  }
  else {
    return(list(v=v, u=u, obj=obj, diff=diff))
  }
}

solve_svd_fxn <- function(X., thres.=10^-2){
  p <- dim(X.)[1]
  res <- svd(X.)
  Dr <- matrix(0, nrow = p, ncol = p)
  d <- 1/res$d[abs(res$d) > thres.]
  diag(Dr)[1:length(d)] <- d
  out <- solve(t(res$v)) %*% Dr %*% solve(res$u)
  return(out)
}
