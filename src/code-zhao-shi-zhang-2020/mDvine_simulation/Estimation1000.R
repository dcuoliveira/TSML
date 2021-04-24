rm(list=ls())
set.seed(7)
source("./MDvineUtil.R")
source("./JointLogLikelihood.R")
library(foreach)
library(doParallel)
library(copula)
library(VineCopula)
library(mvtnorm)
# registerDoParallel(20)

B <- 1000 #length of the time series
reptime <- 500 #number of experiments

###########################
# Simulate the Multivariate D-vine (First Level tau=0.5, second Level tau=0.2)
# uDvine1 - Gaussian + Gumbel
# uDvine2 - t + Clayton
# uDvine3 - Gaussian + Gaussian
totalfam <- list(c(1, 4), c(2, 3), c(1, 1))  
totalpara <- list(c(0.7, 1.25), c(0.7, 0.5, 3), c(0.7, 0.3)) #c(par1, par2)
loc <- list(NULL, c(1), NULL) #location of  par2 (which level of tree has par2)
lb <- list(c(0.1, 1.1), c(0.1, 0.1, 2.1), c(0.1, 0.1))
ub <- list(c(0.9, 4), c(0.9, 4, 6), c(0.9, 0.9))

jointpara <- c(0.2, 0.5, 0.8) #CJ copula Gaussian
eps <- 10^(-7)
totalD <- 3
totalTree <- 2

finale <- list()
finale <- foreach(i=1:reptime, .packages=c('copula','VineCopula','mvtnorm')) %dopar% {
  # Simulate continuous data with Normal(0,1) marginals
  totaldata <- qnorm(simulateMDvine(jointpara, totalfam, totalpara, loc, B, eps))
  
  # Estimation of each uDvine
  estpara <- list()
  for(d in 1:totalD) {
    unitf <- runif(1,0,1)
    init <- unitf*lb[[d]] + (1-unitf)*ub[[d]]
    temp <- nlminb(init, loglikUDvine, totaldata=totaldata[d,], fam=totalfam[[d]], loc=loc[[d]],
                   lower=lb[[d]], upper=ub[[d]])
    estpara[[d]] <- temp$par
  }

  # Estimate cross-sectional copula CJ
  # Transform each uDvine into F(y_t|y_(t-1),...,y_(t-p))
  totalConData <- matrix(0, nrow=totalD, ncol=(B-totalTree)) #Store the conditional CDF
  for(d in 1:totalD) {
    totalOb <- rank(totaldata[d,])/(B+1)  #scaled-ECDF
    for(t in (totalTree+1):B){
      totalConData[d, (t-totalTree)] <-
        predConCDF(totalOb[t], totalOb[(t-totalTree):(t-1)], fam=totalfam[[d]], 
                   para=estpara[[d]], loc=loc[[d]])
    }
  }
  
  # Optimize on total random initial points
  lower <- c(rep(0.1, 2), 0.1)
  upper <- c(rep(0.9, 2), 0.9)
  total <- 2
  zop <- list()
  obj <- rep(0, total)
  for(index in 1:total){
    init <- runif(3,0.1,0.9)
    zop[[index]] <-
      nlminb(init, loglikCJ, totalConData=totalConData, lower=lower, upper=upper)
    obj[index] <- zop[[index]]$objective
  }
  estpara[[(totalD+1)]] <- recoverRho(zop[[which.min(obj)]]$par)
  print(estpara)
  estpara
}


###########################
# Aggregate result
result <- lapply(finale, unlist)
result <- do.call(rbind, result)
estpara <- apply(result, 2, mean)
sd <- apply(result, 2, sd)

truepara <- c(unlist(totalpara), jointpara)
round(estpara, 3)
truepara-round(estpara,3)
round(sd, 3)