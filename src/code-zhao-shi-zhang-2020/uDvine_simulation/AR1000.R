rm(list=ls())
set.seed(7)
options(digits=10)
source("./GarchUtil.R")
source("./MDvineUtil.R")
library(copula)
library(foreach)
library(doParallel)
library(VineCopula)
# registerDoParallel(100)

###########################
# Parameter setting
totaltrain0 <- 1000
totaltest <- 50
total <- totaltrain0 + totaltest
B <- 500 # repeat time

# beta <- c(0.6,0.5,-0.3,-0.5,0.4,-0.3) # AR(6) parameters
beta <- c(0.6,0.5,-0.3) # AR(3) parameters
ARp <- length(beta)
abs(polyroot(c(-rev(beta),1))) # test stationarity of AR(p)
familyset <- c(1,2,3,4,5,6) # candidate set of bivariate copulas 
independent_test_significance_level <- 0.1 # type I error level of the independent test

###########################
# Simulation starts
finalResult <- foreach(i=1:B, .errorhandling='remove', .packages=c('copula', 'VineCopula')) %dopar% {
  set.seed(7+i)

  # Simulate AR(p) data
  data0 <- arima.sim(model=list(ar=beta),total)
  traindata0 <- data0[1:totaltrain0]
  testdata0 <- data0[-c(1:(totaltrain0))]
  
  # Estimate and Select uDvine 
  pseudo_dvine <- rank(traindata0)/(totaltrain0+1) # transform by ECDF
  independent <- FALSE
  treeLevel <- 1 # level of tree
  individualResult <- list() # store result for trees
  totalfam <- totalpara1 <- totalpara2 <- c()
  while(!independent){
    if(treeLevel==1){
      totaltrain <- length(pseudo_dvine)
      pseudo_tree <- cbind(pseudo_dvine[1:(totaltrain-treeLevel)], pseudo_dvine[(1+treeLevel):totaltrain])
    } else {
      totaltrain <- dim(pseudo_tree)[1]
      pseudo_tree <- cbind(BiCopHfunc(pseudo_tree[1:(totaltrain-1),1], pseudo_tree[1:(totaltrain-1),2], family=tree$family,
                                      par=tree$par, par2=tree$par2)$hfunc2,
                           BiCopHfunc(pseudo_tree[2:totaltrain,2], pseudo_tree[2:totaltrain,1], family=tree$family,
                                      par=tree$par, par2=tree$par2)$hfunc2)
    }
    # Selecting a bivariate copula from all copulas in familyset 
    tree <- BiCopSelect(pseudo_tree[,1], pseudo_tree[,2], familyset=familyset, selectioncrit='BIC', indeptest=F)
    if(tree$p.value.indeptest>independent_test_significance_level){
      independent <- TRUE
      break
    }
    totalfam <- c(totalfam, tree$family)
    totalpara1 <- c(totalpara1, tree$par)
    totalpara2 <- c(totalpara2, tree$par2)
    individualResult[[treeLevel]] <- tree
    treeLevel <- treeLevel + 1
  }
  Dvine_p <- length(individualResult)
  totalpara <- c(totalpara1, totalpara2)
  print(Dvine_p) 
  
  # Estimate and Select AR(p) by AIC
  ar_est <- ar(traindata0, aic=T)
  ar_est_p <- ar_est$order
  
  ###########################
  # Predictive distribution based on estimated uDvine
  predErr_Dvine <- predErr_oracle <- predErr_arest <- predY_Dvine <- predY_oracle <- predY_arest <- c()
  p <- max(Dvine_p, ARp, ar_est_p)
  for(index in (p+1):totaltest){
    pastOb <- sapply(testdata0[(index-Dvine_p):(index-1)], empCDF, data=traindata0)
    tmp_pred <- bootstrapLongDvine(totalfam, totalpara, loc=1:Dvine_p, pastOb, B=1000, eps=1e-7) # This can be time consuming if B is large
    tmp_predY_Dvine <- mean(quantile(traindata0, tmp_pred))
    tmp_predY_oracle <- sum(testdata0[(index-ARp):(index-1)]*rev(beta))
    tmp_predY_arest <- sum(testdata0[(index-ar_est_p):(index-1)]*rev(ar_est$ar)) + ar_est$x.mean
    
    predY_Dvine <- c(predY_Dvine, tmp_predY_Dvine)
    predY_oracle <- c(predY_oracle, tmp_predY_oracle)
    predY_arest <- c(predY_arest, tmp_predY_arest)
    predErr_Dvine <- c(predErr_Dvine, tmp_predY_Dvine-testdata0[index])
    predErr_oracle <- c(predErr_oracle, tmp_predY_oracle-testdata0[index])
    predErr_arest <- c(predErr_arest, tmp_predY_arest-testdata0[index])
  }
  
  simulation_result <- list(MSE_Dvine=mean(predErr_Dvine^2), MSE_oracle=mean(predErr_oracle^2), MSE_arest=mean(predErr_arest^2),
                            predErr_Dvine=predErr_Dvine, predErr_oracle=predErr_oracle, predErr_arest=predErr_arest,
                            predY_Dvine=predY_Dvine, predY_oracle=predY_oracle, predY_arest=predY_arest,
                            totalfam=totalfam, totalpara=totalpara, Dvine_p=Dvine_p, ar_est_para=ar_est$ar, ar_est_p=ar_est_p)
  simulation_result
}

# save.image(paste0("AR", length(beta), "_", totaltrain0, ".RData"))