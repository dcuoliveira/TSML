#Function to calculate -Loglikelihood for uDvine
loglikUDvine <- function(para, totaldata, fam, loc) {
  totalTree <- length(fam)
  
  totalOb <- length(totaldata)
  totaldata <- rank(totaldata)/(totalOb+1)  #scaled-ECDF
  loglik <- 0
  for(index in (totalTree+1):totalOb){
    loglik <- loglik + log(predConPDF(totaldata[index], totaldata[(index-totalTree):(index-1)], 
                                      fam=fam, para=para, loc=loc))
  }
  -loglik
}

#Function to calculate -Loglikelihood for joint loglikCJ given F(|...)'s and jointCpara
loglikCJ <- function(jointCpara, totalConData) {
  rho <- recoverRho(jointCpara)
  gauC <- normalCopula(param=rho, dim=3, dispstr='un')

  loglik <- 0
  for(t in 1:dim(totalConData)[2]){
    #conditional CDF transformed observations
    ob1 <- totalConData[1,t] 
    ob2 <- totalConData[2,t]
    ob3 <- totalConData[3,t] 
    
    ob1<-ifelse (ob1 < 1e-9, 1e-9, ob1)
    ob1<-ifelse (ob1 > 1-1e-9, 1-1e-9, ob1)
    ob2<-ifelse (ob2 < 1e-9, 1e-9, ob2)
    ob2<-ifelse (ob2 > 1-1e-9, 1-1e-9, ob2)
    ob3<-ifelse (ob3 < 1e-9, 1e-9, ob3)
    ob3<-ifelse (ob3 > 1-1e-9, 1-1e-9, ob3)    
    
    loglik <- loglik + log(dCopula(c(ob1, ob2, ob3), gauC))
  }
  -loglik
}
