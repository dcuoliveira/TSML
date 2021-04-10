rm(list=ls())
library('SparseTSCGM')
library('dplyr')

## Function signature (Multivariate time series simulation with chain graphical models)

 #sim_tscgm = sim.data(model=c("ar1"), # VAR order
 #                     time=1000, # numer of time points
 #                     n.obs=1, # number of observations or replicates
 #                     n.var=5, # number of variables
 #                     seed=221994, # 
 #                     prob0=0.1, # initial sparsity level
 #                     network=c("random"), # type of the network
 #                     prec=NULL, # precision matrix (inverse of the covariance)
 #                     gamma1=NULL, # autoregressive coefficient matrix at time lag 1
 #                     gamma2=NULL) #  autoregressive coefficient matrix at time lag 2

## Functions output

# data1 - Multivariate time series data in longitudinal format
# theta - Sparse precision matrix
# gamma - Sparse autoregressive coefficients
# sigma - Covariance matrix

datas = sim.data(model="ar1", 
                  time=2,
                  n.obs=100,
                  n.var=5,
                  seed=221994,
                  prob0=0.35,
                  network="random")
simu_data = datas$data1
precision_matrix = datas$theta
coef_matrix = datas$gamma

ts.plot(datas$data1, col=c("blue","green","red","yellow", "black"))

ts.plot(cbind(simu_data[,3], simu_data[,2]), col=c("blue", "red"))
