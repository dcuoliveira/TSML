rm(list=ls())
set.seed(7)
source("./MDvineUtil.R")
source("./JointLogLikelihood.R")
library(foreach)
library(doParallel)
library(copula)
library(VineCopula)
library(mvtnorm)
# registerDoParallel(30)

B <- 1000 #length of the time series
reptime <- 500 #number of experiments

###########################
# Simulate the Multivariate D-vine (First Level Kendall tau=0.5, second Level Kendall tau=0.2)
# uDvine1 - Gaussian + Gumbel
# uDvine2 - t + Clayton
# uDvine3 - Gaussian + Gaussian
totalfam <- list(c(1, 4), c(2, 3), c(1, 1))  
totalpara <- list(c(0.7, 1.25), c(0.7, 0.5, 3), c(0.7, 0.3)) #c(par1, par2)
loc <- list(NULL, c(1), NULL) #location of  par2 (which level of tree has par2)
lb <- list(c(0.1, 1.1), c(0.1, 0.1, 2.1), c(0.1, 0.1))
ub <- list(c(0.9, 4), c(0.9, 4, 6), c(0.9, 0.9))

familyset <- c(1,2,3,4,5,6) # candidate set of bivariate copulas 
#(Gaussian, t, Clayton, Gumbel, Frank, Joe)

jointpara <- c(0.2, 0.5, 0.8)
eps <- 10^(-7)
totalD <- 3
totalTree <- 2

finale <- list()
finale <- foreach(i=1:reptime, .packages=c('copula','VineCopula','mvtnorm')) %dopar% {
  # Simulate continuous data with Normal(0,1) marginals
  totaldata <- qnorm(simulateMDvine(jointpara, totalfam, totalpara, loc, B, eps))
  
  # Tree by tree BIC selection with Independent Test
  estpara <- list()
  for(d in 1:totalD) {
    subtraindata <- totaldata[d,]
    totaltrain <- length(subtraindata)
    
    pseudo_dvine <- rank(subtraindata)/(totaltrain+1) # transform by ECDF
    independent <- FALSE
    treeLevel <- 1 # level of tree
    individualResult <- list() # store result for trees
    while(!independent){
      if(treeLevel==1){
        pseudo_tree <- cbind(pseudo_dvine[1:(totaltrain-treeLevel)], pseudo_dvine[(1+treeLevel):totaltrain])
      } else {
        totaltrain <- dim(pseudo_tree)[1]
        pseudo_tree <- cbind(BiCopHfunc(pseudo_tree[1:(totaltrain-1),1], pseudo_tree[1:(totaltrain-1),2], family=tree$family,
                                         par=tree$par, par2=tree$par2)$hfunc2,
                             BiCopHfunc(pseudo_tree[2:totaltrain,2], pseudo_tree[2:totaltrain,1], family=tree$family,
                                         par=tree$par, par2=tree$par2)$hfunc2)
      }
      #Selecting from all the list in Vinecopula R package
      tree <- BiCopSelect(pseudo_tree[,1], pseudo_tree[,2], familyset=familyset, 
                          selectioncrit='BIC', indeptest=T, rotations=F)
      if(tree$p.value.indeptest>0.01){
        independent <- TRUE
        break
      }
      individualResult[[treeLevel]] <- tree
      treeLevel <- treeLevel + 1
    }
    estpara[[d]] <- individualResult
  }
  print(estpara)
  estpara
}


###########################
# Aggregate result
uDvine1 <- uDvine2 <- uDvine3 <- list()
for(index in 1:reptime){
  uDvine1[[index]] <- finale[[index]][[1]]
  uDvine2[[index]] <- finale[[index]][[2]]
  uDvine3[[index]] <- finale[[index]][[3]]
}
uDvine <- list(uDvine1, uDvine2, uDvine3)

# Right order of uDvine
orderP <- rep(0,3)
for(index in 1:3) {
  orderP[index] <- sum(unlist(lapply(uDvine[[index]], length))==2)/reptime
}

# Family recovery
famRate <- cbind(rep(0,3), rep(0,3))
fam <- list()
for(d in 1:3){
  temp <- matrix(0, ncol=2, nrow=reptime)
  temp1 <- matrix(0, ncol=2, nrow=reptime)
  for(index in 1:reptime){
    if(uDvine[[d]][[index]][[1]]$family==totalfam[[d]][1]){
      temp[index, 1] <- 1
    }
    if(uDvine[[d]][[index]][[2]]$family==totalfam[[d]][2]){
      temp[index, 2] <- 1
    }
    temp1[index, 1] <- uDvine[[d]][[index]][[1]]$family
    temp1[index, 2] <- uDvine[[d]][[index]][[2]]$family
  }
  fam[[d]] <- temp1
  famRate[d,] <- apply(temp, 2, mean)
}

orderP
famRate
summary(as.factor((fam[[1]])[,1]))
summary(as.factor((fam[[1]])[,2]))
summary(as.factor((fam[[2]])[,1]))
summary(as.factor((fam[[2]])[,2]))
summary(as.factor((fam[[3]])[,1]))
summary(as.factor((fam[[3]])[,2]))
