rm(list=ls())
library('SparseTSCGM')

?sim.data

 sim_tscgm = sim.data(model=c("ar1","ar2"), # VAR order
                      time=time, # numer of time points
                      n.obs=n.obs, # number of observations or replicates
                      n.var=n.var, # number of variables
                      seed=NULL, # 
                      prob0=NULL, # initial sparsity level
                      network=c("random","scale-free","hub","user_defined"), # type of the network
                      prec=NULL, # precision matrix (inverse of the covariance)
                      gamma1=NULL, # autoregressive coefficient matrix at time lag 1
                      gamma2=NULL) #  autoregressive coefficient matrix at time lag 2