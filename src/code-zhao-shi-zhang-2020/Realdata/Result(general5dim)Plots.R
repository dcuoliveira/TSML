rm(list=ls())
options(digits=10)
source("./ElectricityUtil.R")
source("./DCCUtil.R")
library(copula)
library(fpp)
library(forecast)
library(CDVine)
library(mvtnorm)
library(ks)
library(vars)
library(kdecopula)

################################
# Train/Test Separation (LOG(data) has already been day-effect adjusted by dummy variable regression
# and seasonality-trend adjusted by STL)
data0 <- read.csv("electricity(5dim).csv")  #NSW	VIC	QLD	TAS	SA

start <- 1
trainyear <- 4
testyear <- 2
freq <- 365
traindata0 <- data0[start:(start+trainyear*freq-1),-1]
testdata0 <- data0[(trainyear*freq+start):(start-1+(trainyear+testyear)*freq),-1]

totaltrain <- dim(traindata0)[1]
totaltest <- dim(testdata0)[1]
totalD <- dim(traindata0)[2]
bootlength <- 1e4
ts.plot(rbind(traindata0, testdata0), col=1:totalD)

################################
# Competitor: Parametric VAR Model 
# VAR with AIC p selected = 1
library(vars)
varResult <- VAR(traindata0, p=1)
varErr <- cbind(varResult$varresult[[1]]$residuals, varResult$varresult[[2]]$residuals,
                varResult$varresult[[3]]$residuals, varResult$varresult[[4]]$residuals,
                varResult$varresult[[5]]$residuals)
orderVAR <- varResult$p
simu_VAR1 <- bootVAR(varResult, varErr, bootlength=bootlength)

# VAR with DCC
library(MTS)
dcc_stage1 <- dccPre(traindata0, include.mean=T, p=1, cond.dist='norm')
dcc_stage2 <- dccFit(dcc_stage1$sresi, cond.dist='norm', type='Engle')
simu_VARDCC <- bootVARDCC(varResult, dcc_stage1, dcc_stage2, bootlength=bootlength)

################################
# CuDvine model
#Marginal time series (Fit two-tree Dvine)
margDvineResult <- list()
for(index in 1:totalD){
  subtraindata <- traindata0[,index]
  #Estimate and Select Tree1
  pseudo_dvine <- rank(subtraindata)/(totaltrain+1)
  pseudo_tree1 <- cbind(pseudo_dvine[1:(totaltrain-1)], pseudo_dvine[2:totaltrain])
  # plot(qnorm(pseudo_tree1), main="Dvine Tree1")
  #Select Tree1 from various bivariate copulas
  tree1 <- BiCopSelect(pseudo_tree1[,1], pseudo_tree1[,2], selectioncrit='BIC', indeptest=F)
  
  #Estimate and Select Tree2
  #Construct pseudo observations for tree2
  pseudo_tree2 <- cbind(BiCopHfunc(pseudo_dvine[1:(totaltrain-2)], pseudo_dvine[2:(totaltrain-1)], family=tree1$family,
                                   par=tree1$par, par2=tree1$par2)$hfunc2,
                        BiCopHfunc(pseudo_dvine[3:(totaltrain)], pseudo_dvine[2:(totaltrain-1)], family=tree1$family,
                                   par=tree1$par, par2=tree1$par2)$hfunc2)
  # plot(qnorm(pseudo_tree2), main="Dvine Tree2")
  #Select Tree2 from various bivariate copulas
  tree2 <- BiCopSelect(pseudo_tree2[,1], pseudo_tree2[,2], selectioncrit='BIC', indeptest=T)
  individualResult <- list(tree1=tree1, tree2=tree2)
  margDvineResult[[index]] <- individualResult
}

################################
# Cross-sectional conditional F(|) relationship
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
  # plot(qnorm(trans_traindata[,pairs[,index]]), main="Cross-sectional copula Dvine")
  biCopulaPair<- BiCopSelect(trans_traindata[,pairs[1,index]], trans_traindata[,pairs[2,index]], 
                             selectioncrit='AIC', indeptest=T)
  # print(biCopulaPair)
}

# Use a t-copula to fit the result 
jointCopulaDvinetConst <- fitCopula(tCopula(dim=totalD, dispstr='un'), trans_traindata)
ub <- c(0.99, 0.99)
lb <- c(0.01, 0.01)
start <- c(0.1, 0.8)
jointCopulaDvineGaussianDCC <- nlminb(start, DCCGaussianCopulaLoglik, trans_traindata=trans_traindata,
                                      upper=ub, lower=lb)                      

# Fit a time-varying t-copula via DCC
ub <- c(20, 0.99, 0.99)
lb <- c(4, 0.01, 0.01)
start <- c(jointCopulaDvinetConst@estimate[11], jointCopulaDvineGaussianDCC$par)
jointCopulaDvinetDCC <- nlminb(start, DCCtCopulaLoglik, trans_traindata=trans_traindata, upper=ub, lower=lb)
simu_Dvine <- bootDvine(margDvineResult, jointCopulaDvinetDCC, trans_traindata, traindata0, eps=10^(-7), bootlength=bootlength)

################################
# Plot of copula
burnin <- 1000
simu_VAR1 <- simu_VAR1[-c(1:burnin),]
simu_Dvine <- simu_Dvine[-c(1:burnin),]
simu_VARDCC <- simu_VARDCC[-c(1:burnin),]

lags <- 1 # time lags
region_names <- colnames(traindata0)
method_names <- c('Empirical', 'VAR(1)', 'VAR(1)-DCC', 'CuDvine')

# Plot self-lag
for(region_index1 in 1:2){ # NSW and VIC
  par(mfrow=c(1,4), mai=c(0.5, 0.6, 0.4, 0))
  for(method in 1:4){
    if(method==1){
      plot_data <- traindata0 # empirical
    }
    if(method==2){
      plot_data <- simu_VAR1 # VAR1
    }
    if(method==3){
      plot_data <- simu_VARDCC # VAR-DCC
    }
    if(method==4){
      plot_data <- simu_Dvine # CuDvine
    }
    data_length <- dim(plot_data)[1]
    
    region_train1 <- rank(plot_data[,region_index1])/(data_length+1)
    tmp_data <- cbind(region_train1[-c((data_length-lags+1):data_length)], region_train1[-c(1:lags)])
    tmp_copula <- kdecop(tmp_data)
    print(paste(region_names[region_index1], ':', method_names[method]))
    print(dep_measures(tmp_copula))
    contour(tmp_copula, main=paste(method_names[method], '\ntau:', round(dep_measures(tmp_copula)[1],2), ' rho:', round(dep_measures(tmp_copula)[2], 2)),
            ylab=paste0(region_names[region_index1], '(lagged)'), xlab=paste0(region_names[region_index1]), cex.main=1.6, cex.lab=1.6)
  }
}

# Plot cross lag
for(region_index1 in 1){ # NSW v.s. VIC
  for(region_index2 in 2){
    par(mfrow=c(1,4), mai=c(0.5, 0.6, 0.4, 0))
    for(method in 1:4){
      if(method==1){
        plot_data <- traindata0 # empirical
      }
      if(method==2){
        plot_data <- simu_VAR1 # VAR1
      }
      if(method==3){
        plot_data <- simu_VARDCC # VAR-DCC
      }
      if(method==4){
        plot_data <- simu_Dvine # CuDvine
      }
      data_length <- dim(plot_data)[1]
      
      region_train1 <- rank(plot_data[,region_index1])/(data_length+1)
      region_train2 <- rank(plot_data[,region_index2])/(data_length+1)
      
      tmp_data <- cbind(region_train1[-c((data_length-lags+1):data_length)], region_train2[-c(1:lags)])
      tmp_copula <- kdecop(tmp_data)
      print(paste(region_names[region_index1], '(lagged) v.s', region_names[region_index2], ':', method_names[method]))
      print(dep_measures(tmp_copula))

      contour(tmp_copula, main=paste(method_names[method], '\ntau:', round(dep_measures(tmp_copula)[1],2), ' rho:', round(dep_measures(tmp_copula)[2], 2)),
              ylab=paste0(region_names[region_index2], '(lagged)'), xlab=paste0(region_names[region_index1]), cex.main=1.6, cex.lab=1.6)
    }
  }
}
