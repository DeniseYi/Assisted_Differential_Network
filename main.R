#---------------------------------------------------------------------------------#
# Simulation main file 
#---------------------------------------------------------------------------------#
rm(list=ls())
#dev.off()
setwd("/Users/hy283/Desktop/Research/Differential Network/Publication-code/")

# Step 0: load functions
source("assisted_differential_network.R")
source("networks.R")
library(MASS)
library(matrixcalc)
library(irlba)

set.seed(1)
p = 50
q = 50
n = 200

#set.seed(1)
# Step 1: Generate regulators
sigmaX1 <- sigmaGenerator1(q, block. = 5) # Generate block diag cov matrix
sigmaX2 <- sigmaChange1(sigmaX1, change.=c(3), sub.p.=10) # Make change to sigmaX1
# test = (sigmaX1-sigmaX2)*((sigmaX1-sigmaX2)!=0)

# Step 2: Generate GE
# 2.a: Generate beta (fixed & sparse): block diag or muilky way
Beta1 <- betaGenerator(p.=p, q. = q, block.=5) # block diag

# 2.b: Generate cov matrix of error
sigmaW <- matrix(0, p, p)
diag(sigmaW) <- 1 + rnorm(p, 0, 0.1)
is.positive.semi.definite(sigmaW)

#rm(networkER, betaGenerator, networkBD, noiseGenerator, sigmaChange1, sigmaGenerator1)
#save.image("/Users/hy283/Desktop/Research/Differential Network/DNSVD_Code/submit/AA1B1C1.RData")


# Replicate begins here
# Set tuning parameters for different methods including 3 alternatives and the proposed method
l = 10 ##
lambda1 = seq(0.05, 1.5, length.out = l)
lambda2 = seq(0.5, 24, length.out = l)
lambda3 = seq(0.05, 2.5, length.out = l)
lambda4 = seq(0.05, 1, length.out = l)

lambda5 = seq(1, 24, length.out = l) #1
lambda6 = seq(0.01, 0.5, length.out = l)  #2
lambda7 = seq(0, 1, length.out = l) #3
lambda8 = seq(0.1, 1.5, length.out = l) 
lambda9 = seq(0.05, 2.5, length.out = l) #8
lambda10 = seq(0, 1, length.out = l)

seq1 = seq(1,q,2)
seq2 = seq(1,p,2)
res_covX_svd = res_covY_svd = res_preX_svd = res_preY_svd <- list() 
res_covX_ssvd = res_covY_ssvd = res_preX_ssvd = res_preY_ssvd <- list() 
res_covX_myssvd = res_covY_myssvd = res_preX_myssvd = res_preY_myssvd <- list()
res_cov_sssvd = res_pre_sssvd <- list()

# simulation
id = as.numeric(Sys.time())
nid = (trunc(id)%%10000 + id%%1) *100000
nid = as.integer(nid)
set.seed(nid)

rep <- 1
replicate = 100

while(rep <= replicate){
  print(paste("Replicate =", rep))
  X1 <- mvrnorm(n,rep(0,q),sigmaX1)
  X2 <- mvrnorm(n,rep(0,q),sigmaX2)
  Y1 <- X1 %*% Beta1 + mvrnorm(n,rep(0,p),sigmaW)
  Y2 <- X2 %*% Beta1 + mvrnorm(n,rep(0,p),sigmaW)
  
  covXdiff <- cov(X2)-cov(X1)
  covYdiff <- cov(Y2)-cov(Y1)
  
  Rdiff <- solve_svd_fxn(cov(X2)) - solve_svd_fxn(cov(X1))
  Gdiff <- solve_svd_fxn(cov(Y2)) - solve_svd_fxn(cov(Y1))
  
  corYX <- cor(Y1, X1) + cor(Y2, X2)
  
  # Alternative 1: SVD
  res_covX_svd[[rep]] <- svd(covXdiff, nu=1, nv=1)
  res_covY_svd[[rep]] <- svd(covYdiff, nu=1, nv=1)
  res_preX_svd[[rep]] <- svd(Rdiff, nu=1, nv=1)
  res_preY_svd[[rep]] <- svd(Gdiff, nu=1, nv=1)
  
  # Alternative 2: irlba package
  covX_ssvd = covY_ssvd = preX_ssvd = preY_ssvd <- list() 
  for (i in seq1){ #X
    covX_ssvd[[length(covX_ssvd)+1]] <- ssvd(covXdiff, k = 1, n = i)
    preX_ssvd[[length(preX_ssvd)+1]] <- ssvd(Rdiff, k = 1, n = i)
  }
  for (i in seq2){
    covY_ssvd[[length(covY_ssvd)+1]] <- ssvd(covYdiff, k = 1, n = i)
    preY_ssvd[[length(preY_ssvd)+1]] <- ssvd(Gdiff, k = 1, n = i)
  }
  res_covX_ssvd[[rep]] <- covX_ssvd
  res_covY_ssvd[[rep]] <- covY_ssvd
  res_preX_ssvd[[rep]] <- preX_ssvd
  res_preY_ssvd[[rep]] <- preY_ssvd
  
  # Alternative 3: SSVD_fun
  covX_myssvd = covY_myssvd = preX_myssvd = preY_myssvd = list()
  for (i in 1:length(lambda1)){
    covX_myssvd[[i]] <- ssvd_fxn(X.=covXdiff, lambda. = lambda1[i], type. = "covX")
    covY_myssvd[[i]] <- ssvd_fxn(X.=covYdiff, lambda. = lambda2[i], type. = "covY") 
    preX_myssvd[[i]] <- ssvd_fxn(X.=Rdiff, lambda. = lambda3[i], type. = "preX")
    preY_myssvd[[i]] <- ssvd_fxn(X.=Gdiff, lambda. = lambda4[i], type. = "preY")
  }
  res_covX_myssvd[[rep]] <- covX_myssvd
  res_covY_myssvd[[rep]] <- covY_myssvd
  res_preX_myssvd[[rep]] <- preX_myssvd
  res_preY_myssvd[[rep]] <- preY_myssvd
  print("Alternative 3 finished.")
  
  # Proposed method
  
  
  cov_sssvd = pre_sssvd = list()
  for (i in 1:length(lambda1)){
    for (j in 1:length(lambda2)){
      for (k in 1:length(lambda3)){
        cov_sssvd[[length(cov_sssvd)+1]] <- sssvd_fxn(G.=covYdiff, R.=covXdiff, lbd1.=lambda5[i], lbd2.=lambda6[j], 
                                                      lbd3.=lambda7[k], cor.=corYX, type. = "cov") ## remark: lbd1.=lambda2[i]
        pre_sssvd[[length(pre_sssvd)+1]] <- sssvd_fxn(G.=Gdiff, R.=Rdiff, lbd1.=lambda8[i], lbd2.=lambda9[j], 
                                                      lbd3.=lambda10[k], cor.=corYX)
      }
    }
  }
  res_cov_sssvd[[rep]] <- cov_sssvd
  res_pre_sssvd[[rep]] <- pre_sssvd
  
  print(paste("Simulation finished", rep, "replicates."))
  rep <- rep + 1
}

save.image(file = paste("A1B1C1-", nid, ".RData", sep = ""))

print("Completed")