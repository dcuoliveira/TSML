library(CDVine)
library(expm)
##################################################
### Conditional distribution of Two-tree Dvine ###
##################################################
condCdf_T1_T2 <- function(u1, u2, u3, tree1, tree2) {
  BiCopHfunc(BiCopHfunc(u3, u2, family=tree1$family, tree1$par, tree1$par2)$hfunc2, 
             BiCopHfunc(u1, u2, family=tree1$family, tree1$par, tree1$par2)$hfunc2, 
             family=tree2$family, tree2$par, tree2$par2)$hfunc2
}

invCondCdf_T1_T2 <- function(u, u1, u2, tree1, tree2, eps) {
  if (condCdf_T1_T2(u1, u2, 1-eps, tree1, tree2)<u) {return(1-eps)}
  if (condCdf_T1_T2(u1, u2, eps, tree1, tree2)>u) {return(eps)}
  return(uniroot(function(x)(condCdf_T1_T2(u1, u2, x, tree1, tree2)-u), 
                 lower=eps, upper=1-eps)[1]$root)
}

##############################
### Rescaled empirical CDF ###
##############################
empCDF <- function(observation, data) {
  n <- length(data)
  sum(data<=observation)/(n+1)
}

###############################
### Bootstrap MCMC Marginal ###
###############################
#Bootstrap AR(1) + GARCH(1,1) Semiparametrically/Parametrically
bootstrapGarch <- function(etaErr, ob_old, sigma_old, B, est) {
  mean <- est[1] + est[2]*ob_old
  sigma <- sqrt(est[3]+est[4]*ob_old^2+est[5]*sigma_old^2)
  c(as.numeric(mean+sigma*quantile(etaErr, runif(B,0,1))), sigma)
}

#Bootstrap AR(p) Semiparametrically/Parametrically
bootstrapAR <- function(testdata, pos, arEst, arErr, B) {
  order <- length(arEst)-1
  mean <- arEst[1] + sum(arEst[2:(order+1)]*testdata[(pos-1):(pos-order)])
  as.numeric(mean+quantile(arErr, runif(B,0,1)))
  #mean
}

############################
### Bootstrap MCMC Joint ###
############################
#Bootstrap two-tree Dvine Semiparametrically
bootstrapDvine <- function(data, ob1, ob2, tree1, tree2, eps, B) {
  as.numeric(quantile(data, sapply(runif(B,0,1), invCondCdf_T1_T2, empCDF(ob1, data), 
                                   empCDF(ob2, data), tree1, tree2, eps)))
}

bootstrapDvineJoint <- function(traindata, ob1, ob2, margDvineResult, jointCopula, eps, B) {
  ob1 <- as.numeric(ob1); ob2 <- as.numeric(ob2)
  totalD <- length(ob1)
  rand <- BiCopSim(B, family=jointCopula$family, par=jointCopula$par, par2=jointCopula$par2)
  
  result <- vector()
  for(index in 1:totalD) {
    tmp <- as.numeric(quantile(traindata[,index], sapply(rand[,index], invCondCdf_T1_T2, u1=empCDF(ob1[index], traindata[,index]), 
                                                     u2=empCDF(ob2[index], traindata[,index]), tree1=margDvineResult[[index]]$tree1, 
                                                     tree2=margDvineResult[[index]]$tree2, eps=eps)))
    
    result <- cbind(result, tmp)
  }
  result
}

# Trivariate Copula Simulation
bootstrapDvineJointTrivariate <- function(traindata, ob1, ob2, margDvineResult, jointCopula, eps, B, dayAhead=1) {
  if(dayAhead==1){
    ob1 <- as.numeric(ob1); ob2 <- as.numeric(ob2); 
    totalD <- length(ob1)
    len <- length(jointCopula@estimate)
    cross_copula <- tCopula(param=jointCopula@estimate[1:(len-1)], df=jointCopula@estimate[len], dim=totalD, dispstr='un')
    rand <- rCopula(B, cross_copula)
    
    result <- vector()
    for(index in 1:totalD) {
      tmp <- as.numeric(quantile(traindata[,index], sapply(rand[,index], invCondCdf_T1_T2, u1=empCDF(ob1[index], traindata[,index]), 
                                                           u2=empCDF(ob2[index], traindata[,index]), tree1=margDvineResult[[index]]$tree1, 
                                                           tree2=margDvineResult[[index]]$tree2, eps=eps)))
      
      result <- cbind(result, tmp)
    }
  }
  return(result)
}

# Trivariate Copula Simulation with Gaussian DCC
bootstrapDvineJointGauDCC <- function(traindata, ob1, ob2, margDvineResult, DynCor, eps, B, dayAhead=1) {
  if(dayAhead==1){
    ob1 <- as.numeric(ob1); ob2 <- as.numeric(ob2); 
    totalD <- length(ob1)
    cross_copula <- normalCopula(param=DynCor, dim=totalD, dispstr='un')
    rand <- rCopula(B, cross_copula)
    
    result <- vector()
    for(index in 1:totalD) {
      tmp <- as.numeric(quantile(traindata[,index], sapply(rand[,index], invCondCdf_T1_T2, u1=empCDF(ob1[index], traindata[,index]), 
                                                           u2=empCDF(ob2[index], traindata[,index]), tree1=margDvineResult[[index]]$tree1, 
                                                           tree2=margDvineResult[[index]]$tree2, eps=eps)))
      
      result <- cbind(result, tmp)
    }
  }
  return(result)
}

# Trivariate Copula Simulation with Gaussian DCC
bootstrapDvineJointtDCC <- function(traindata, ob1, ob2, margDvineResult, DynCor, nu, eps, B, dayAhead=1) {
  if(dayAhead==1){
    ob1 <- as.numeric(ob1); ob2 <- as.numeric(ob2); 
    totalD <- length(ob1)
    cross_copula <- tCopula(param=DynCor, df=nu, dim=totalD, dispstr='un')
    rand <- rCopula(B, cross_copula)
    result <- vector()
    for(index in 1:totalD) {
      tmp <- as.numeric(quantile(traindata[,index], sapply(rand[,index], invCondCdf_T1_T2, u1=empCDF(ob1[index], traindata[,index]), 
                                                           u2=empCDF(ob2[index], traindata[,index]), tree1=margDvineResult[[index]]$tree1, 
                                                           tree2=margDvineResult[[index]]$tree2, eps=eps)))
      
      result <- cbind(result, tmp)
    }
  }
  return(result)
}

#Bootstrap VAR(p) Semiparametrically/Parametrically
bootstrapVAR <- function(testdata, pos, varErr, varResult, B, type="Parametric") {
  if(type=="Parametric") {
    totalD <- dim(varErr)[2]
    order <- varResult$p
    covar <- (cov(varErr))
    
    result <- vector()
    for(index in 1:totalD) {
      arEst <- varResult$varresult[[index]]$coefficients
      mean <- arEst[(totalD*order+1)] + sum(arEst[1:(totalD*order)]*as.vector(t(testdata[(pos-1):(pos-order),])))
      tmp <- as.numeric(mean + sqrt(covar[index,index])*rnorm(B,0,1))
      result <- cbind(result, tmp)
    }
  } else {
    totalD <- dim(varErr)[2]
    order <- varResult$p
    
    result <- vector()
    for(index in 1:totalD) {
      arEst <- varResult$varresult[[index]]$coefficients
      mean <- arEst[(totalD*order+1)] + sum(arEst[1:(totalD*order)]*as.vector(t(testdata[(pos-1):(pos-order),])))
      tmp <- as.numeric(mean + quantile(varErr[,index], runif(B,0,1)))
      result <- cbind(result, tmp)
    }
  }
  result
}

#Bootstrap Multivariate AR Semiparametrically/Parametrically
bootstrapSemiAR <- function(testdata, pos, margARResult, jointCopula, B) {
  totalD <- length(margARResult)
  rand <- BiCopSim(B, family=jointCopula$family, par=jointCopula$par, par2=jointCopula$par2)
  
  result <- vector()
  for(index in 1:totalD) {
    arEst <- margARResult[[index]]$arEst
    arErr <- margARResult[[index]]$arErr
    order <- length(arEst)-1
    mean <- arEst[1] + sum(arEst[2:(order+1)]*testdata[(pos-1):(pos-order), index])
    tmp <- as.numeric(mean+quantile(arErr, rand[,index]))    
    result <- cbind(result, tmp)
  }
  result
}

#Simulate Trivariate AR
bootstrapSemiARTrivariate <- function(testdata, pos, margARResult, jointCopula, B) {
  totalD <- length(margARResult)
  len <- length(jointCopula@estimate)
  cross_copula <- normalCopula(param=jointCopula@estimate[1:len], dim=totalD, dispstr='un')
  rand <- rCopula(B, cross_copula)
  
  result <- vector()
  for(index in 1:totalD) {
    arEst <- margARResult[[index]]$arEst
    arErr <- margARResult[[index]]$arErr
    order <- length(arEst)-1
    mean <- arEst[1] + sum(arEst[2:(order+1)]*testdata[(pos-1):(pos-order), index])
    tmp <- as.numeric(mean+quantile(arErr, rand[,index]))    
    result <- cbind(result, tmp)
  }
  result
}

#crps calculation
crpsComp <- function(X, X1, ob) {
  ob <- as.numeric(ob)
  B <- dim(X)[1] 
  totalD <- dim(X)[2]
  crps <- rep(0, totalD)
  for(index in 1:totalD) {
    crps[index] <- sum(abs(X[,index]-ob[index]))/B - sum(abs(X[,index]-X1[,index]))/2/B
  }
  crps
}

crpsCompDim1 <- function(X, X1, ob){
  ob <- as.numeric(ob)
  B <- length(X)
  sum(abs(X-ob))/B - sum(abs(X-X1))/2/B
}

#qrps calculation
qrpsComp <- function(quantile, pred, ob){
  ob <- as.numeric(ob)
  totalD <- length(ob)
  qrps <- rep(0, totalD)
  for(index in 1:totalD) {
    qrps[index] <- (pred[index]-ob[index])*((ob[index]<=pred[index])-quantile)
  }
  qrps
}

#Extract List
extractList <- function(list, position){
  list[[position]]
}

# Simulate VAR(1) based on estimated VAR(1) model
bootVAR <- function(varResult, varErr, bootlength=10000){
  totalD <- length(varResult$varresult)
  # VAR parameters
  VAR_A <- c()
  for(index in 1:totalD){
    VAR_A <- rbind(VAR_A, varResult$varresult[[index]]$coefficients)
  }
  VAR_const <- VAR_A[,totalD+1]
  VAR_A <- VAR_A[,-c(totalD+1)]
  
  # Covariance
  covar <- cov(varErr)
  totalY <- c()
  lastY <- rep(0, totalD)
  for(index in 1:bootlength){
    # Simulate error terms
    error <- sqrtm(covar)%*%rnorm(totalD)
    lastY <- VAR_const + VAR_A%*%lastY + error
    totalY <- rbind(totalY, t(lastY))
  }
  return(totalY)
}

# Simulate VAR-DCC based on estimated model
bootVARDCC <- function(varResult, dcc_stage1, dcc_stage2, bootlength=10000){
  totalD <- length(varResult$varresult)
  # VAR parameters
  VAR_A <- c()
  for(index in 1:totalD){
    VAR_A <- rbind(VAR_A, varResult$varresult[[index]]$coefficients)
  }
  VAR_const <- VAR_A[,totalD+1]
  VAR_A <- VAR_A[,-c(totalD+1)]
  
  # DCC marginal GARCH parameters
  max_variance <- 10*apply(dcc_stage1$marVol,2,max)^2 # to make sure simulated VAR-DCC is stable
  last_train <- dim(dcc_stage1$marVol)[1]
  last_margvar <- diag(dcc_stage1$marVol[last_train,])^2 # variance scale not std scale
  garch_omega <- diag(dcc_stage1$est[,1])
  garch_alpha <- diag(dcc_stage1$est[,2])
  garch_beta <- diag(dcc_stage1$est[,3])
  
  # DCC correlation matrix parameter
  last_dccR <- last_dccQ <- matrix(dcc_stage2$rho.t[last_train,], totalD, totalD)
  dcc_S <- cor(dcc_stage1$sresi)
  dcc_alpha <- dcc_stage2$estimates[2]
  dcc_beta <- dcc_stage2$estimates[1]
  
  # recover the sequence
  totalY <- c()
  lastY <- rep(0, totalD)
  for(index in 1:bootlength){
    # Simulate error terms
    std_error <- sqrtm(last_dccR)%*%rnorm(totalD)
    error <- sqrtm(last_margvar)%*%std_error
    lastY <- VAR_const + VAR_A%*%lastY + error
    totalY <- rbind(totalY, t(lastY))
    
    # Update variance in GARCH
    last_margvar <- garch_omega + garch_beta%*%last_margvar + garch_alpha%*%diag(diag(error%*%t(error)))
    last_margvar <- diag(pmin(diag(last_margvar), max_variance))
    # Update correlation in DCC
    last_dccQ <- dcc_S*(1-dcc_alpha-dcc_beta) + dcc_alpha*(std_error%*%t(std_error)) + dcc_beta*last_dccQ
    last_dccR <- diag(diag(last_dccQ)^(-0.5))%*%last_dccQ%*%diag(diag(last_dccQ)^(-0.5))
  }
  return(totalY)
}

bootDvine <- function(margDvineResult, jointCopulaDvinetDCC, trans_traindata, traindata, eps, bootlength) {
  # DCC parameters
  DCC_para <- jointCopulaDvinetDCC$par
  nu <- DCC_para[1]
  alpha <- DCC_para[2]
  beta <- DCC_para[3]
  
  # Initial correlation
  trans_data_tscale <- qt(trans_traindata, df=nu)
  S <- cov(trans_data_tscale)
  Q <- S
  DynCor_last <- P2p(cov2cor(Q))
  
  # Initial Y
  totaltrain <- dim(traindata)[1]
  Ylag1 <- traindata[totaltrain,]
  Ylag2 <- traindata[totaltrain-1,]
  
  # Simulate
  totalY <- vector()
  for(t in 1:bootlength){
    cross_copula <- tCopula(param=DynCor_last, df=nu, dim=totalD, dispstr='un')
    rand <- rCopula(1, cross_copula)
    
    newY <- c()
    for(index in 1:totalD) {
      tmp <- as.numeric(quantile(traindata[,index], sapply(rand[,index], invCondCdf_T1_T2, u1=empCDF(Ylag2[index], traindata[,index]), 
                                                           u2=empCDF(Ylag1[index], traindata[,index]), tree1=margDvineResult[[index]]$tree1, 
                                                           tree2=margDvineResult[[index]]$tree2, eps=eps)))
      newY <- c(newY, tmp)
    }
    totalY <- rbind(totalY, newY)
    Ylag2 <- Ylag1
    Ylag1 <- newY
    
    Q <- (1-alpha-beta)*S + alpha*(t(qt(rand, df=nu))%*%qt(rand, df=nu)) + beta*Q
    DynCor_last <- P2p(cov2cor(Q))
  }
  return(totalY)
}