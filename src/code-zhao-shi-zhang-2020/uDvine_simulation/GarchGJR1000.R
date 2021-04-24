rm(list=ls())
set.seed(7)
options(digits=10)
source("./GarchUtil.R")
library(copula)
library(fpp)
library(forecast)
library(CDVine)
library(vars)
library(foreach)
library(doParallel)
# registerDoParallel(20)

#############################
# Parameter setting
totaltrain <- 1000
totaltest <- 100
total <- totaltrain + totaltest
B <- 500
beta <- c(0.1, 0.05)

# Simulation starts
finalVaR <- matrix(0, nrow=B, ncol=length(beta))
finalVaR <- foreach(i = 1:B, .combine=rbind, .errorhandling='remove', .packages=c('copula','fpp','forecast','CDVine','vars')) %dopar% {
  set.seed(7+i)

  # Simulate GJR-GARCH(1,1) data
  data0 <- simuGjrGARCH(n=total)
  traindata0 <- data0[1:totaltrain]
  testdata0 <- data0[-c(1:(totaltrain))]
  
  # Estimate and Select uDvine 
  pseudo_dvine <- rank(traindata0)/(totaltrain+1)
  pseudo_tree1 <- cbind(pseudo_dvine[1:(totaltrain-1)], pseudo_dvine[2:totaltrain])
  #Select Tree1 from various bivariate copulas
  tree1 <- BiCopSelect(pseudo_tree1[,1], pseudo_tree1[,2], selectioncrit='BIC', indeptest=F)
  
  #Estimate and Select Tree2
  #Construct pseudo observations for tree2
  pseudo_tree2 <- cbind(BiCopHfunc(pseudo_dvine[1:(totaltrain-2)], pseudo_dvine[2:(totaltrain-1)], family=tree1$family,
                                   par=tree1$par, par2=tree1$par2)$hfunc2,
                        BiCopHfunc(pseudo_dvine[3:(totaltrain)], pseudo_dvine[2:(totaltrain-1)], family=tree1$family,
                                   par=tree1$par, par2=tree1$par2)$hfunc2)
  #Select Tree2 from various bivariate copulas
  tree2 <- BiCopSelect(pseudo_tree2[,1], pseudo_tree2[,2], selectioncrit='BIC', indeptest=T)

  # Predictive distribution based on estimated uDvine
  vaR <- matrix(0, ncol=length(beta), nrow=(totaltest-2))
  for(index in 3:totaltest){
    pred <- bootstrapDvine(traindata0, ob1=testdata0[index-2], 
                           ob2=testdata0[index-1], tree1, tree2, eps=1e-7, B=1000)
    vaR[index-2,] <- testdata0[index] <= quantile(pred, beta)
  }
  print(apply(vaR,2,mean))
  apply(vaR,2,mean)
}

meanVaR <- apply(finalVaR, 2, mean)
sdVaR <- apply(finalVaR, 2, sd)
p_value <- 2*(1-pnorm(abs(sqrt(B)*(meanVaR-beta)/sdVaR)))
meanVaR
p_value
