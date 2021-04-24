DCCtCopulaLoglik <- function(para, trans_traindata){
  # print(para)
  nu <- para[1]
  alpha <- para[2]
  beta <- para[3]
  if(alpha+beta>0.99999){
    return(1e10)
  }
  trans_traindata_tscale <- qt(trans_traindata, df=nu)
  S <- cov(trans_traindata_tscale)
  totalT <- dim(trans_traindata_tscale)[1]
  totalD <- dim(trans_traindata_tscale)[2]
  Q <- S
  currentCopula <- tCopula(param=P2p(cov2cor(Q)), dim=totalD, df=nu, dispstr='un')
  loglik <- log(dCopula(trans_traindata[1,], copula=currentCopula))
  for(t in 2:totalT){
    Q <- (1-alpha-beta)*S + alpha*(trans_traindata_tscale[(t-1),]%*%t(trans_traindata_tscale[(t-1),])) + beta*(Q)
    currentCopula <- tCopula(param=P2p(cov2cor(Q)), dim=totalD, df=nu, dispstr='un')
    loglik <- loglik + log(dCopula(trans_traindata[t,], copula=currentCopula))
  }
  return(-loglik)
}

DCCtCopulaLoglikWithNu <- function(para, trans_traindata, nu){
  # print(para)
  alpha <- para[1]
  beta <- para[2]
  if(alpha+beta>0.99999){
    return(1e10)
  }
  trans_traindata_tscale <- qt(trans_traindata, df=nu)
  S <- cov(trans_traindata_tscale)
  totalT <- dim(trans_traindata_tscale)[1]
  totalD <- dim(trans_traindata_tscale)[2]
  Q <- S
  currentCopula <- tCopula(param=P2p(cov2cor(Q)), dim=totalD, df=nu, dispstr='un')
  loglik <- log(dCopula(trans_traindata[1,], copula=currentCopula))
  for(t in 2:totalT){
    Q <- (1-alpha-beta)*S + alpha*(trans_traindata_tscale[(t-1),]%*%t(trans_traindata_tscale[(t-1),])) + beta*(Q)
    currentCopula <- tCopula(param=P2p(cov2cor(Q)), dim=totalD, df=nu, dispstr='un')
    loglik <- loglik + log(dCopula(trans_traindata[t,], copula=currentCopula))
  }
  return(-loglik)
}

DCCGaussianCopulaLoglik <- function(para, trans_traindata){
  # print(para)
  alpha <- para[1]
  beta <- para[2]
  if(alpha+beta>0.99999){
    return(1e10)
  }
  trans_traindata_gauscale <- qnorm(trans_traindata)
  S <- cov(trans_traindata_gauscale)
  totalT <- dim(trans_traindata_gauscale)[1]
  totalD <- dim(trans_traindata_gauscale)[2]
  Q <- S
  currentCopula <- normalCopula(param=P2p(cov2cor(Q)), dim=totalD, dispstr='un')
  loglik <- log(dCopula(trans_traindata[1,], copula=currentCopula))
  for(t in 2:totalT){
    Q <- (1-alpha-beta)*S + alpha*(trans_traindata_gauscale[(t-1),]%*%t(trans_traindata_gauscale[(t-1),])) + beta*(Q)
    currentCopula <- normalCopula(param=P2p(cov2cor(Q)), dim=totalD, dispstr='un')
    loglik <- loglik + log(dCopula(trans_traindata[t,], copula=currentCopula))
  }
  return(-loglik)
}

recoverDCCGau <- function(para, trans_data, S=NULL){
  alpha <- para[1]
  beta <- para[2]
  trans_data_gauscale <- qnorm(trans_data)
  if(is.null(S)){
    S <- cov(trans_data_gauscale)
  }
  totalT <- dim(trans_data_gauscale)[1]
  totalD <- dim(trans_data_gauscale)[2]
  result <- c()
  Q <- S
  result <- cbind(result, P2p(cov2cor(Q)))
  for(t in 2:totalT){
    Q <- (1-alpha-beta)*S + alpha*(trans_data_gauscale[(t-1),]%*%t(trans_data_gauscale[(t-1),])) + beta*(Q)
    result <- cbind(result, P2p(cov2cor(Q)))
  }
  return(result)
}

recoverDCCt <- function(para, trans_data, S=NULL, position=NULL){
  nu <- para[1]
  alpha <- para[2]
  beta <- para[3]
  if(alpha+beta>0.99999){
    return(1e10)
  }
  trans_data_tscale <- qt(trans_data, df=nu)
  if(is.null(S)){
    if(is.null(position)){
      position <- dim(trans_data)[1]
    }
    S <- cov(trans_data_tscale[1:position,])
  }
  totalT <- dim(trans_data_tscale)[1]
  totalD <- dim(trans_data_tscale)[2]
  result <- c()
  Q <- S
  result <- cbind(result, P2p(cov2cor(Q)))
  for(t in 2:totalT){
    Q <- (1-alpha-beta)*S + alpha*(trans_data_tscale[(t-1),]%*%t(trans_data_tscale[(t-1),])) + beta*(Q)
    result <- cbind(result, P2p(cov2cor(Q)))
  }
  return(result)
}

recoverVARDCC <- function(testdata, traindata, varResult, dcc_stage1, dcc_stage2){
  # function to recover VAR(1) + DCC
  total_test <- dim(testdata)[1]
  total_train <- dim(traindata)[1]
  totalD <- dim(testdata)[2]
  
  # VAR(1) parameters
  VAR_A <- c()
  for(index in 1:totalD){
    VAR_A <- rbind(VAR_A, varResult$varresult[[index]]$coefficients)
  }
  VAR_const <- VAR_A[,totalD+1]
  VAR_A <- VAR_A[,-c(totalD+1)]
  
  # DCC marginal GARCH parameters
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
  totaldata <- rbind(traindata, testdata)
  margvar_seq <- dccR_seq <- list()
  for(index in (total_train+1):(total_train+total_test)){
    error <- as.numeric(totaldata[(index-1),]-VAR_A%*%as.numeric(totaldata[(index-2),]))-VAR_const
    std_error <- error/sqrt(diag(last_margvar))
    last_margvar <- garch_omega + garch_beta%*%last_margvar + garch_alpha%*%diag(diag(error%*%t(error)))
    margvar_seq[[index-total_train]] <- last_margvar
    
    last_dccQ <- dcc_S*(1-dcc_alpha-dcc_beta) + dcc_alpha*(std_error%*%t(std_error)) + dcc_beta*last_dccQ
    last_dccR <- diag(diag(last_dccQ)^(-0.5))%*%last_dccQ%*%diag(diag(last_dccQ)^(-0.5))
    dccR_seq[[index-total_train]] <- last_dccR
  }
  result <- list(margvar_seq=margvar_seq, dccR_seq=dccR_seq)  
  return(result)
}

bootstrapVARDCC <- function(testdata, pos, varErr, varResult, B, VARDCC_cov, type="Parametric_Marg") {
  if(type=="Parametric_Marg") {
    totalD <- dim(varErr)[2]
    order <- varResult$p
    covar <- VARDCC_cov$margvar_seq[[pos]]
    
    result <- vector()
    for(index in 1:totalD){
      arEst <- varResult$varresult[[index]]$coefficients
      mean <- arEst[(totalD*order+1)] + sum(arEst[1:(totalD*order)]*as.vector(t(testdata[(pos-1):(pos-order),])))
      tmp <- as.numeric(mean + sqrt(covar[index,index])*rnorm(B,0,1))
      result <- cbind(result, tmp)
    }
  } else {
    totalD <- dim(varErr)[2]
    order <- varResult$p
    covar <- VARDCC_cov$margvar_seq[[pos]]
    covar <-  diag(diag(covar)^(0.5))%*%VARDCC_cov$dccR_seq[[pos]]%*%diag(diag(covar)^(0.5))
    error_term <- sqrtm(covar)%*%matrix(rnorm(B*totalD), totalD, B)
    
    result <- vector()
    for(index in 1:totalD){
      arEst <- varResult$varresult[[index]]$coefficients
      mean <- arEst[(totalD*order+1)] + sum(arEst[1:(totalD*order)]*as.vector(t(testdata[(pos-1):(pos-order),])))
      tmp <- as.numeric(mean + error_term[index,])
      result <- cbind(result, tmp)
    }
  }
  result
}