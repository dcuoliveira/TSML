rm(list=ls())
set.seed(7)
options(digits=10)
source("./ElectricityUtil.R")
source("./DCCUtil.R")
library(copula)
library(fpp)
library(forecast)
library(CDVine)
library(vars)
library(foreach)
library(doParallel)
library(mvtnorm)
# registerDoParallel(24)
par(mfrow=c(2,2))

##################################
# Train/Test Separation (LOG(data) has already been day-effect adjusted by dummy variable regression
# and seasonality-trend adjusted by STL)
data0 <- read.csv("./electricity(5dim).csv")  #NSW	VIC	QLD	TAS	SA

start <- 1
trainyear <- 4
testyear <- 2
freq <- 365
traindata0 <- data0[start:(start+trainyear*freq-1),-1]
testdata0 <- data0[(trainyear*freq+start):(start-1+(trainyear+testyear)*freq),-1]

totaltrain <- dim(traindata0)[1]
totaltest <- dim(testdata0)[1]
totalD <- dim(traindata0)[2]
ts.plot(rbind(traindata0, testdata0), col=1:totalD)
ts.plot(apply(rbind(traindata0, testdata0),1,sum))
abline(v=365*(1:7), col='red')

##################################
# Marginal time series (Fit two-tree uDvine)
margDvineResult <- list()
for(index in 1:totalD){
  #subtraindata <- traindata[,index]
  subtraindata <- traindata0[,index]
  #Estimate and Select Tree1
  pseudo_dvine <- rank(subtraindata)/(totaltrain+1)
  pseudo_tree1 <- cbind(pseudo_dvine[1:(totaltrain-1)], pseudo_dvine[2:totaltrain])
  plot(qnorm(pseudo_tree1), main="Dvine Tree1")
  #Select Tree1 from various bivariate copulas
  tree1 <- BiCopSelect(pseudo_tree1[,1], pseudo_tree1[,2], selectioncrit='BIC', indeptest=F)
  
  #Estimate and Select Tree2
  #Construct pseudo observations for tree2
  pseudo_tree2 <- cbind(BiCopHfunc(pseudo_dvine[1:(totaltrain-2)], pseudo_dvine[2:(totaltrain-1)], family=tree1$family,
                                   par=tree1$par, par2=tree1$par2)$hfunc2,
                        BiCopHfunc(pseudo_dvine[3:(totaltrain)], pseudo_dvine[2:(totaltrain-1)], family=tree1$family,
                                   par=tree1$par, par2=tree1$par2)$hfunc2)
  plot(qnorm(pseudo_tree2), main="Dvine Tree2")
  #Select Tree2 from various bivariate copulas
  tree2 <- BiCopSelect(pseudo_tree2[,1], pseudo_tree2[,2], selectioncrit='BIC', indeptest=T)
  
  individualResult <- list(tree1=tree1, tree2=tree2)
  margDvineResult[[index]] <- individualResult
}

##################################
# Cross-sectional conditional F(|) relationship
# Fit a bivariate copula to each pair
#pairs <- list(c(1,2), c(1,3), c(2,3)) #NSW-VIC| NSW-QLD| VIC-QLD
#Joint copula estimation
notree <- 2 #Each marginal is a 2-tree Dvine
trans_traindata <- matrix(0, ncol=totalD, nrow=(totaltrain-notree)) #Store F(|...)
for(index in 1:totalD){
  subtraindata <- traindata0[,index]
  trans_subtraindata <- sapply(subtraindata, empCDF, data=subtraindata)
  trans_traindata[,index]<-
    condCdf_T1_T2(u1=trans_subtraindata[1:(totaltrain-2)], u2=trans_subtraindata[2:(totaltrain-1)],
                  u3=trans_subtraindata[3:totaltrain], tree1=margDvineResult[[index]]$tree1, 
                  tree2=margDvineResult[[index]]$tree2)
}

pairs <- combn(5,2)
for(index in 1:dim(pairs)[2]){
  plot(qnorm(trans_traindata[,pairs[,index]]), main="Cross-sectional copula Dvine")
  biCopulaPair<- BiCopSelect(trans_traindata[,pairs[1,index]], trans_traindata[,pairs[2,index]], 
                             selectioncrit='AIC', indeptest=T)
  print(biCopulaPair)
}

# Use a constant t-copula to fit the result 
jointCopulaDvinetConst <- fitCopula(tCopula(dim=totalD, dispstr='un'), trans_traindata)
        
# Fit a time-varying Gaussian-copula via DCC
ub <- c(0.99, 0.99)
lb <- c(0.01, 0.01)
start <- c(0.1, 0.8)
jointCopulaDvineGaussianDCC <- nlminb(start, DCCGaussianCopulaLoglik, trans_traindata=trans_traindata,
                                      upper=ub, lower=lb)                      
DCCGau <- recoverDCCGau(jointCopulaDvineGaussianDCC$par, trans_traindata)
ts.plot(apply(DCCGau, 2, mean))

# Fit a time-varying t-copula via DCC
ub <- c(20, 0.99, 0.99)
lb <- c(4, 0.01, 0.01)
start <- c(jointCopulaDvinetConst@estimate[11], jointCopulaDvineGaussianDCC$par)
jointCopulaDvinetDCC <- nlminb(start, DCCtCopulaLoglik, trans_traindata=trans_traindata,
                               upper=ub, lower=lb)
DCCt <- recoverDCCt(jointCopulaDvinetDCC$par, trans_traindata)
ts.plot(apply(DCCt, 2, mean))

##################################
#VaR estimation, CRPS, RMSE based on Bootstrap(time consuming)
B <- 2000; eps <- 10^(-7); 
# Recover the DCC correlation matrix for DCCGaussian in the test set.
notree <- 2
trans_data <- matrix(0, ncol=totalD, nrow=(totaltest+totaltrain-notree)) #Store F(|...)
totaldata0 <- rbind(traindata0, testdata0)
for(index in 1:totalD){
  subdata <- totaldata0[,index]
  trans_subdata <- sapply(subdata, empCDF, data=subdata)
  trans_data[,index]<-
    condCdf_T1_T2(u1=trans_subdata[1:(totaltest+totaltrain-2)], u2=trans_subdata[2:(totaltest+totaltrain-1)],
                  u3=trans_subdata[3:(totaltest+totaltrain)], tree1=margDvineResult[[index]]$tree1, 
                  tree2=margDvineResult[[index]]$tree2)
}

orderVAR <- 1 # Maximum order of VAR model
orderAR <- 18 # Maximum order AR + Copula model
teststart <- max(orderVAR, notree, orderAR)+1
DCCt_total <- recoverDCCt(jointCopulaDvinetDCC$par, trans_data)
DCCt_test <- DCCt_total[,(totaltrain+teststart-notree):(totaltrain+totaltest-notree)]

# Plot time-varying correlation
plot(apply(DCCt_total[,c(1:(5*365))],2,mean), type='l', xlab='', ylab='',
     main='The time-varying average correlation among 5 regions', xaxt='n')
abline(v=365*4, col='red')
xtick<-seq(0, 2000, by=365)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3]-0.005,
     labels=c(2009:2014), pos=1, xpd=TRUE)


##################################
#1-day ahead prediction for each method in the test set
testResult <- foreach(i=teststart:totaltest, .packages='copula') %dopar% {
  set.seed(i+7)
  estDvine <- bootstrapDvineJointtDCC(traindata0, testdata0[(i-2),], testdata0[(i-1),],
                                      margDvineResult, DCCt_test[,(i-teststart+1)], 
                                      jointCopulaDvinetDCC$par[1], eps, 2*B)

  #### Univariate Result
  # Prediction Error and CRPS for each univariate series X and Y
  predErrDvine <- as.numeric(apply(estDvine, 2, mean)-testdata0[i,])
  crpsDvine <- crpsComp(estDvine[1:B,], estDvine[(B+1):(2*B),], ob=testdata0[i,])
  # Value at Risk for each univariate series X and Y
  varDvine95 <- as.numeric(apply(estDvine, 2, quantile, 0.95))>testdata0[i,]
  varDvine99 <- as.numeric(apply(estDvine, 2, quantile, 0.99))>testdata0[i,]

  ##### Bivariate Result
  combin <- combn(5, 2)
  pairs <- list()
  for(index in 1:dim(combin)[2]){
    pairs[[index]] <- combin[,index]
  }
  crpsDiffDvine <- crpsDiffAR <- crpsDiffVAR <- vector()
  crpsAddDvine <- crpsAddAR <- crpsAddVAR <- vector()
  varDiffDvine95 <- varDiffDvine99 <- varDiffAR95 <- varDiffAR99 <- varDiffVAR95 <- varDiffVAR99 <- vector()
  varAddDvine95 <- varAddDvine99 <- varAddAR95 <- varAddAR99 <- varAddVAR95 <- varAddVAR99 <- vector()
  qrpsDiffDvine95 <- qrpsDiffDvine99 <- qrpsDiffAR95 <- qrpsDiffAR99 <- qrpsDiffVAR95 <- qrpsDiffVAR99 <- vector()
  qrpsAddDvine95 <- qrpsAddDvine99 <- qrpsAddAR95 <- qrpsAddAR99 <- qrpsAddVAR95 <- qrpsAddVAR99 <- vector()
  
  for(pair in pairs){
    #CRPS for X-Y
    crpsDiffDvine <- c(crpsDiffDvine,
                       crpsCompDim1(estDvine[1:B,pair[2]]-estDvine[1:B,pair[1]], estDvine[(B+1):(2*B),pair[2]]-estDvine[(B+1):(2*B),pair[1]], testdata0[i,pair[2]]-testdata0[i,pair[1]]))
    #CRPS for X+Y
    crpsAddDvine <- c(crpsAddDvine,
                      crpsCompDim1(estDvine[1:B,pair[2]]+estDvine[1:B,pair[1]], estDvine[(B+1):(2*B),pair[2]]+estDvine[(B+1):(2*B),pair[1]], testdata0[i,pair[2]]+testdata0[i,pair[1]]))

    #Value at Risk for X-Y
    varDiffDvine95 <- c(varDiffDvine95,
                        quantile(estDvine[,pair[2]]-estDvine[,pair[1]], 0.95)>testdata0[i,pair[2]]-testdata0[i,pair[1]])
    varDiffDvine99 <- c(varDiffDvine99,
                        quantile(estDvine[,pair[2]]-estDvine[,pair[1]], 0.99)>testdata0[i,pair[2]]-testdata0[i,pair[1]])

    #QRPS for X-Y
    qrpsDiffDvine95 <- c(qrpsDiffDvine95,
                         qrpsComp(0.95, quantile(estDvine[,pair[2]]-estDvine[,pair[1]], 0.95), testdata0[i,pair[2]]-testdata0[i,pair[1]]))
    qrpsDiffDvine99 <- c(qrpsDiffDvine99,
                         qrpsComp(0.99, quantile(estDvine[,pair[2]]-estDvine[,pair[1]], 0.99), testdata0[i,pair[2]]-testdata0[i,pair[1]]))

    #Value at Risk for X+Y
    varAddDvine95 <- c(varAddDvine95,
                       quantile(estDvine[,pair[2]]+estDvine[,pair[1]], 0.95)>testdata0[i,pair[2]]+testdata0[i,pair[1]])
    varAddDvine99 <- c(varAddDvine99,
                       quantile(estDvine[,pair[2]]+estDvine[,pair[1]], 0.99)>testdata0[i,pair[2]]+testdata0[i,pair[1]])

    #QRPS for X+Y
    qrpsAddDvine95 <- c(qrpsAddDvine95,
                        qrpsComp(0.95, quantile(estDvine[,pair[2]]+estDvine[,pair[1]], 0.95), testdata0[i,pair[2]]+testdata0[i,pair[1]]))
    qrpsAddDvine99 <- c(qrpsAddDvine99,
                        qrpsComp(0.99, quantile(estDvine[,pair[2]]+estDvine[,pair[1]], 0.99), testdata0[i,pair[2]]+testdata0[i,pair[1]]))
  }
  
  ##### 5-dimensional Result (CRPS and CI coverage)
  # Price Index for Electricity Market, weighted sum where weights by Demand
  demandWeight <- c(8235.034, 5476.738, 5913.029, 1120.563, 1441.425)
  demandWeight <- demandWeight/sum(demandWeight)
  weightSum <- function(v) {sum(v*demandWeight)}
  
  crpsSumDvine <- 
    crpsCompDim1(apply(estDvine[1:B,], 1, weightSum), apply(estDvine[(B+1):(2*B),], 1, weightSum), weightSum(testdata0[i,]))

  #VAR for weightSum
  varSumDvine95 <- 
    (quantile(apply(estDvine[1:B,], 1, weightSum), 0.95)>weightSum(testdata0[i,]))
 
  varSumDvine99 <- 
    (quantile(apply(estDvine[1:B,], 1, weightSum), 0.99)>weightSum(testdata0[i,]))
 
  #CI for weightSum
  ciDvine95 <- 
    (quantile(apply(estDvine[1:B,], 1, weightSum), 0.975)>weightSum(testdata0[i,]))*(quantile(apply(estDvine[1:B,], 1, weightSum), 0.025)<weightSum(testdata0[i,]))
 
  ciDvine99 <- 
    (quantile(apply(estDvine[1:B,], 1, weightSum), 0.995)>weightSum(testdata0[i,]))*(quantile(apply(estDvine[1:B,], 1, weightSum), 0.005)<weightSum(testdata0[i,]))
  
  #QRPS for weightSum
  qrpsSumDvine95 <- 
    qrpsComp(0.95, quantile(apply(estDvine[1:B,], 1, weightSum), 0.95), weightSum(testdata0[i,]))
 
  qrpsSumDvine99 <- 
    qrpsComp(0.99, quantile(apply(estDvine[1:B,], 1, weightSum), 0.99), weightSum(testdata0[i,]))

  ##### Combine all results
  #Univariate 
  var95 <- list(varDvine95) 
  var99 <- list(varDvine99) 
  crps <- list(crpsDvine)
  predErr <- list(predErrDvine)
  #Bivariate Diff
  varDiff95 <- list(varDiffDvine95)
  varDiff99 <- list(varDiffDvine99)
  qrpsDiff95 <- list(qrpsDiffDvine95)
  qrpsDiff99 <- list(qrpsDiffDvine99)
  crpsDiff <- list(crpsDiffDvine)
  #Bivariate Add
  varAdd95 <- list(varAddDvine95)
  varAdd99 <- list(varAddDvine99)
  qrpsAdd95 <- list(qrpsAddDvine95)
  qrpsAdd99 <- list(qrpsAddDvine99)
  crpsAdd <- list(crpsAddDvine)
  #5-dimensional Sum
  crpsSum <- list(crpsSumDvine)
  varSum95 <- list(varSumDvine95)
  varSum99 <- list(varSumDvine99)
  qrpsSum95 <- list(qrpsSumDvine95)
  qrpsSum99 <- list(qrpsSumDvine99)
  ciSum95 <- list(ciDvine95)
  ciSum99 <- list(ciDvine99)
  
  #Combine final result
  result <- list(var95, var99, crps, predErr, #1
                 varDiff95, varDiff99, crpsDiff, #5
                 varAdd95, varAdd99, crpsAdd, #8
                 crpsSum, varSum95, varSum99, ciSum95, #11
                 qrpsDiff95, qrpsDiff99, qrpsAdd95, qrpsAdd99, #15
                 qrpsSum95, qrpsSum99, #19
                 ciSum99)
  print(i)
  result
}
