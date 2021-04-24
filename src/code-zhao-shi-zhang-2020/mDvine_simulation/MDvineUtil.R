options(digits=10)
library(GenSA)
library(numDeriv)

# conditional CDF
condCDF1 <- function(u1, u2, family, par, par2=4) {
  BiCopHfunc2(u1, u2, family, par, par2)
}

# conditional CDF
condCDF2 <- function(u1, u2, family, par, par2=4) {
  BiCopHfunc1(u1, u2, family, par, par2)
}

#evaluate F(x|pastOb) [Can be used for prediction and all kinds of F(|) calculation]
#need length(pastOb)>=2
predConCDF <- function(x, pastOb, fam, para, loc) {
    #The main property used is that D-vine is closed among time shifting
    #All the observations in currentPol are in F(Y) scale (i.e. pastOb is in F() scale)
    currentPol <- c(pastOb, x) # history is (Y(T-p),...,Y(T)) #p is the no of trees
    
    para1 <- para[1:length(fam)]
    para2=rep(0, length(fam))
    if(!is.null(loc)){
      para2[loc] <- para[-c(1:length(fam))]
    }
    
    sizeOb <- length(currentPol) # usually should be of size p+1
    if(sizeOb>=3){
      temp1 <- matrix(0, ncol=(sizeOb-2), nrow=(sizeOb-2)) #store s|s+1,...,s+t-1
      temp2 <- matrix(0, ncol=(sizeOb-2), nrow=(sizeOb-2)) #store s+t|s+1,...,s+t-1
      
      basis <- currentPol #store [F(y(T-p)),..., F(y(T)), x] #For t=2
      for(s in 1:(sizeOb-2)){
        temp1[1,s] <- condCDF1(u1=basis[s], u2=basis[(s+1)], family=fam[1], par=para1[1], par2=para2[1]) #store s|s+1  
        temp2[1,s] <- condCDF2(u1=basis[(s+1)], u2=basis[(s+2)], family=fam[1], par=para1[1], par2=para2[1]) #store s+2|s+1
      }    
      if(sizeOb>3){
        for(t in 3:(sizeOb-1)){
          for(s in 1:(sizeOb-t)) {
            temp1[(t-1),s] <- condCDF1(u1=temp1[(t-2),s], u2=temp2[(t-2),s], family=fam[t-1], par=para1[t-1], par2=para2[t-1])  
            temp2[(t-1),s] <- condCDF2(u1=temp1[(t-2),(s+1)], u2=temp2[(t-2),(s+1)], family=fam[t-1], par=para1[t-1], par2=para2[t-1])
          }
        }
      }
      #store the estimated tempConCDF
      tempConCDF <- condCDF2(u1=temp1[(sizeOb-2),1], u2=temp2[(sizeOb-2),1], family=fam[sizeOb-1], 
                             par=para1[sizeOb-1], par2=para2[sizeOb-1])
    } else {
      basis <- currentPol #store [F(y(T-p)),..., F(y(T)), x] #For t=2
      tempConCDF <- condCDF2(u1=basis[1], u2=basis[2], family=fam[sizeOb-1], par=para1[sizeOb-1],
                             par2=para2[sizeOb-1])
    }
    tempConCDF
}

#Inverse prediction CDF
invPredConCDF <- function(u, pastOb, fam, para, loc, eps) {
  if (predConCDF(1-eps, pastOb, fam, para, loc)<u) {return(1-eps)}
  if (predConCDF(eps, pastOb, fam, para, loc)>u) {return(eps)}
  return(uniroot(function(x)(predConCDF(x, pastOb, fam, para, loc)-u), 
                 lower=eps, upper=1-eps)[1]$root)
}

#evaluate f(x|pastOb) #need length(pastOb)>=2
predConPDF <- function(x, pastOb, fam, para, loc) {
  #The main property used is that D-vine is closed among time shifting
  #All the observations in currentPol are in F(Y) scale (i.e. pastOb is in F() scale)
  currentPol <- c(pastOb, x) # history is (Y(T-p),...,Y(T)) #p is the no of trees
  
  para1 <- para[1:length(fam)]
  para2=rep(0, length(fam))
  if(!is.null(loc)){
    para2[loc] <- para[-c(1:length(fam))]
  }
  
  sizeOb <- length(currentPol) # usually should be of size p+1
  if(sizeOb>=3){
    temp1 <- matrix(0, ncol=(sizeOb-2), nrow=(sizeOb-2)) #store s|s+1,...,s+t-1
    temp2 <- matrix(0, ncol=(sizeOb-2), nrow=(sizeOb-2)) #store s+t|s+1,...,s+t-1
    
    basis <- currentPol #store [F(y(T-p)),..., F(y(T)), x] #For t=2
    for(s in 1:(sizeOb-2)){
      temp1[1,s] <- condCDF1(u1=basis[s], u2=basis[(s+1)], family=fam[1], par=para1[1], par2=para2[1]) #store s|s+1  
      temp2[1,s] <- condCDF2(u1=basis[(s+1)], u2=basis[(s+2)], family=fam[1], par=para1[1], par2=para2[1]) #store s+2|s+1
    }    
    if(sizeOb>3){
      for(t in 3:(sizeOb-1)){
        for(s in 1:(sizeOb-t)) {
          temp1[(t-1),s] <- condCDF1(u1=temp1[(t-2),s], u2=temp2[(t-2),s], family=fam[t-1], par=para1[t-1], par2=para2[t-1])  
          temp2[(t-1),s] <- condCDF2(u1=temp1[(t-2),(s+1)], u2=temp2[(t-2),(s+1)], family=fam[t-1], par=para1[t-1], par2=para2[t-1])
        }
      }
    }
  } else {
    basis <- currentPol #store [F(y(T-p)),..., F(y(T)), x] #For t=2
  }
  #Calculate conditional pdf of f(y_T|y_(T-1),...,y(T-p) based on the conditional CDF's
  f <- BiCopPDF(currentPol[sizeOb-1], currentPol[sizeOb], family=fam[1], par=para1[1], par2=para2[1])
  if(sizeOb>=3){
    for(index in 2:(sizeOb-1)){
      f <- f*BiCopPDF(u1=temp1[(index-1),(sizeOb-index)], u2=temp2[(index-1),(sizeOb-index)], 
                      family=fam[index], par=para1[index], par2=para2[index])
    }
  }
  f 
}

#Simulate the joint Dvine with uniform distribution as marginal (continuous marginal)
simulateMDvine <- function(jointpara, totalfam, totalpara, loc, B, eps) {
  totalTree <- length(totalfam[[1]]) #p
  totalD <- length(totalfam)
  gauC <- normalCopula(param=jointpara, dim=totalD, dispstr='un')
  simuTs <- matrix(0, ncol=B, nrow=totalD) #store the simulated time series
  
  simuTs[,1] <- rCopula(1, gauC) #for t=1
  
  for(t in 2:B){
    temp <- rCopula(1, gauC)
    for(d in 1:totalD){
      simuTs[d,t] <- invPredConCDF(temp[d], pastOb=simuTs[d,(max((1),(t-totalTree)):(t-1))],
                                   fam=totalfam[[d]], para=totalpara[[d]], loc=loc[[d]], eps=eps)
    }
  }
  simuTs
}

#Recover rho from partial correlation JointParaC
recoverRho <- function(jointCpara) {
  rho12 <- jointCpara[1]; 
  rho23 <- jointCpara[2]
  rho13 <- rho12*rho23 + jointCpara[3]*sqrt(1-rho12^2)*sqrt(1-rho23^2)
  rho <- c(rho12, rho13, rho23)
  rho
}

#Return partial correlation based on rho
reverseRho <- function(rho) {
  partial1 <- rho[1] #rho12
  partial2 <- rho[3] #rho23
  partial3 <- (rho[2]-rho[1]*rho[3])/sqrt(1-rho[1]^2)/sqrt(1-rho[3]^2) #partial rho13|2
  return(c(partial1, partial2, partial3))
}