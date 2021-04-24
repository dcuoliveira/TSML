#Conditional distribution of Two-tree Dvine
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

#Bootstrap two-tree Dvine Semiparametrically
bootstrapDvine <- function(data, ob1, ob2, tree1, tree2, eps, B) {
  as.numeric(quantile(data, sapply(runif(B,0,1), invCondCdf_T1_T2, empCDF(ob1, data), 
                                   empCDF(ob2, data), tree1, tree2, eps)))
}

#Rescaled empirical CDF
empCDF <- function(observation, data) {
  n <- length(data)
  sum(data<=observation)/(n+1)
}

#Simulate AR + GARCH
simuARGARCH<-function(n=1000, phi0=0.01, phi1=0.05, omega0=0.05, omega1=0.85, omega2=0.1, burnin=1000){
  n <- n+burnin
  eta=rnorm(n)
  #AR + Garch Simulation
  y=rep(0,n)
  sigma=rep(0,n)
  sigma[1]=omega0  
  y[1]=phi0+sigma[1]*(eta[1])  
  for(i in 2:n){
    sigma[i]=sqrt(omega0+omega1*sigma[(i-1)]^2+omega2*sigma[(i-1)]^2*eta[(i-1)]^2)
    y[i]=phi0+phi1*y[(i-1)]+sigma[i]*eta[i]
  }
  return(y[-(1:burnin)])
}

#Simulate pure GARCH
simuGARCH<-function(n=1000, omega0=0.05, omega1=0.85, omega2=0.1, burnin=1000){
  n <- n+burnin
  eta=rnorm(n)
  #GARCH Simulation
  y=rep(0,n)
  sigma=rep(0,n)
  sigma[1]=omega0  
  y[1]=sigma[1]*(eta[1])  
  for(i in 2:n){
    sigma[i]=sqrt(omega0+omega1*sigma[(i-1)]^2+omega2*sigma[(i-1)]^2*eta[(i-1)]^2)
    y[i]=sigma[i]*eta[i]
  }
  return(y[-(1:burnin)])
}

#Simulate GJR-GARCH
simuGjrGARCH<-function(n=1000, omega0=0.05, omega1=0.85, omega2=0.1, gamma=0.05, burnin=1000){
  n <- n+burnin
  eta=rnorm(n)
  #GARCH Simulation
  y=rep(0,n)
  sigma=rep(0,n)
  sigma[1]=omega0  
  y[1]=sigma[1]*(eta[1])  
  for(i in 2:n){
    sigma[i]=sqrt(omega0+omega1*sigma[(i-1)]^2+(omega2+gamma*(y[i-1]<0))*sigma[(i-1)]^2*eta[(i-1)]^2)
    y[i]=sigma[i]*eta[i]
  }
  return(y[-(1:burnin)])
}

simuBlockFactorGarch <- function(para, groupsize=groupsize, n=1000, burnin=100){
  # Parameter setting
  omega <- para$omega
  beta <- para$beta
  gamma <- para$gamma
  nu <- para$nu # multivariate t with df = 6
  
  # Number of Group
  p <- length(groupsize)
  
  # Recover the correlation matrix from beta and gamma
  subRho <- matrix(0, nrow=p, ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j){
        subRho[i,j]=(beta[i]^2+gamma[i]^2)/(beta[i]^2+gamma[i]^2+1)
      }
      else{
        subRho[i,j]=(beta[i]*beta[j])/(sqrt(beta[i]^2+gamma[i]^2+1)*sqrt(beta[j]^2+gamma[j]^2+1))
      }
    }
  }
  
  wholeRho=matrix(0, nrow=sum(groupsize), ncol=sum(groupsize))
  tracker=c(1,1)
  for(i in 1:p){
    for(j in 1:p){
      wholeRho[tracker[1]:(tracker[1]-1+groupsize[i]),tracker[2]:(tracker[2]-1+groupsize[j])]=subRho[i,j]
      tracker[2]=sum(groupsize[1:j])+1
    }
    tracker[1]=sum(groupsize[1:i])+1
    tracker[2]=1
  }
  diag(wholeRho)=1
  totalD <- sum(groupsize)
  
  # Multivariate Garch simulation
  n <- n+burnin
  sigma <- diag(omega[1], totalD)
  covar <- sqrtm(sigma)%*%wholeRho%*%sqrtm(sigma)*(nu-2)/nu
  lastY <- rmvt(1, sigma=covar, df=nu)

  totalY <- c()
  covar_sets <- c()
  for(i in 1:n){
    # sigma <- diag(omega[1], totalD) + diag(omega[2], totalD)%*%sigma + diag(omega[3], totalD)%*%diag(c(lastY)^2) # give same result
    sigma <- diag(omega[1], totalD) + diag(omega[2], totalD)%*%sigma + diag(omega[3], totalD)%*%diag(diag(t(lastY)%*%lastY))
    covar <- sqrtm(sigma)%*%wholeRho%*%sqrtm(sigma)*(nu-2)/nu
    lastY <- rmvt(1, sigma=covar, df=nu)
    totalY <- rbind(totalY, lastY)
    covar_sets <- rbind(covar_sets, c(covar))
  }
  
  return(list(totalY=totalY[-c(1:burnin),], covar_sets=covar_sets[-c(1:burnin),]))
}

copulaBlockFactorLoglik <- function(para, groupsize, trans_traindata, nu=6){
  p <- length(groupsize) # Number of Group
  beta <- para[1:p]
  gamma <- para[(p+1):(2*p)]
  if(length(para)>2*p){
    nu <- para[2*p+1] # estimate nu
  }
  
  # Recover the correlation matrix from beta and gamma
  subRho <- matrix(0, nrow=p, ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j){
        subRho[i,j]=(beta[i]^2+gamma[i]^2)/(beta[i]^2+gamma[i]^2+1)
      }
      else{
        subRho[i,j]=(beta[i]*beta[j])/(sqrt(beta[i]^2+gamma[i]^2+1)*sqrt(beta[j]^2+gamma[j]^2+1))
      }
    }
  }
  
  wholeRho=matrix(0, nrow=sum(groupsize), ncol=sum(groupsize))
  tracker=c(1,1)
  for(i in 1:p){
    for(j in 1:p){
      wholeRho[tracker[1]:(tracker[1]-1+groupsize[i]),tracker[2]:(tracker[2]-1+groupsize[j])]=subRho[i,j]
      tracker[2]=sum(groupsize[1:j])+1
    }
    tracker[1]=sum(groupsize[1:i])+1
    tracker[2]=1
  }
  diag(wholeRho)=1
  totalD <- sum(groupsize)
  
  tcopula <- tCopula(param=P2p(wholeRho), dim=totalD, df=nu, dispstr='un')
  loglik <- sum(log(dCopula(trans_traindata, copula=tcopula)))
  
  return(-loglik)
}

bootstrapCuDvine <- function(data, ob1, ob2, margDvineResult, para_est, groupsize, eps, B, nu=6) {
  p <- length(groupsize)
  beta <- para_est[1:p]
  gamma <- para_est[(p+1):(2*p)]
  if(length(para_est)>(2*p)){
    nu <- para_est[2*p+1]
  }
  
  # Recover the correlation matrix from beta and gamma
  subRho <- matrix(0, nrow=p, ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j){
        subRho[i,j]=(beta[i]^2+gamma[i]^2)/(beta[i]^2+gamma[i]^2+1)
      }
      else{
        subRho[i,j]=(beta[i]*beta[j])/(sqrt(beta[i]^2+gamma[i]^2+1)*sqrt(beta[j]^2+gamma[j]^2+1))
      }
    }
  }
  
  totalD <- sum(groupsize)
  wholeRho=matrix(0, nrow=totalD, ncol=totalD)
  tracker=c(1,1)
  for(i in 1:p){
    for(j in 1:p){
      wholeRho[tracker[1]:(tracker[1]-1+groupsize[i]),tracker[2]:(tracker[2]-1+groupsize[j])]=subRho[i,j]
      tracker[2]=sum(groupsize[1:j])+1
    }
    tracker[1]=sum(groupsize[1:i])+1
    tracker[2]=1
  }
  diag(wholeRho)=1
  
  cross_copula <- tCopula(param=P2p(wholeRho), df=nu, dim=totalD, dispstr='un')
  rand <- rCopula(B, cross_copula)
  
  boot_result <- c()
  for(index in 1:totalD){
    tmp <- as.numeric(quantile(data[,index], sapply(rand[,index], invCondCdf_T1_T2, empCDF(ob1[index], data[,index]), 
                                                    empCDF(ob2[index], data[,index]), margDvineResult[[index]]$tree1,
                                                    margDvineResult[[index]]$tree2, eps)))
    boot_result <- cbind(boot_result, tmp)
  }
  return(boot_result)
}

bootstrapLongCuDvine <- function(data, pastOb, margDvineResult, para_est, groupsize, eps, B, nu=6) {
  p <- length(groupsize)
  beta <- para_est[1:p]
  gamma <- para_est[(p+1):(2*p)]
  if(length(para_est)>(2*p)){
    nu <- para_est[2*p+1]
  }
  
  # Recover the correlation matrix from beta and gamma
  subRho <- matrix(0, nrow=p, ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j){
        subRho[i,j]=(beta[i]^2+gamma[i]^2)/(beta[i]^2+gamma[i]^2+1)
      }
      else{
        subRho[i,j]=(beta[i]*beta[j])/(sqrt(beta[i]^2+gamma[i]^2+1)*sqrt(beta[j]^2+gamma[j]^2+1))
      }
    }
  }
  
  totalD <- sum(groupsize)
  wholeRho=matrix(0, nrow=totalD, ncol=totalD)
  tracker=c(1,1)
  for(i in 1:p){
    for(j in 1:p){
      wholeRho[tracker[1]:(tracker[1]-1+groupsize[i]),tracker[2]:(tracker[2]-1+groupsize[j])]=subRho[i,j]
      tracker[2]=sum(groupsize[1:j])+1
    }
    tracker[1]=sum(groupsize[1:i])+1
    tracker[2]=1
  }
  diag(wholeRho)=1
  
  cross_copula <- tCopula(param=P2p(wholeRho), df=nu, dim=totalD, dispstr='un')
  rand <- rCopula(B, cross_copula)
  
  boot_result <- c()
  for(index in 1:totalD){
    tmp_fam <- margDvineResult[[index]]$totalfam # t copula
    tmp_loc <- margDvineResult[[index]]$totalloc # t copula
    tmp_para <- margDvineResult[[index]]$totalpara
    
    tmp_pastOb <- c()
    for(tree_index in 1:length(tmp_fam)){
      tmp_pastOb <- c(tmp_pastOb, empCDF(pastOb[tree_index,index], data[,index]))
    }
    tmp <- c()
    for(i in 1:B){
      tmp <- c(tmp, quantile(data[,index], invPredConCDF(rand[i,index], pastOb=tmp_pastOb,
                                                         fam=tmp_fam, para=tmp_para, loc=tmp_loc, eps=eps)))
    }
    boot_result <- cbind(boot_result, tmp)
  }
  return(boot_result)
}