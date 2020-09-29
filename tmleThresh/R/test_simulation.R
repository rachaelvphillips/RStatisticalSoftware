# library(simcausal)
# library(dplyr)
#
#
#
# # create DAG object
# setD=generate_covariates()
#
# # Y follows a logit model
# dat_logit=simu_logit_cov(0.1,n=1000,p0=0.2,p1=1,seed=1,setD=setD)
#
# # Y follows a probit model
# dat_probit=simu_probit_cov(0.1,n=1000,p0=0.2,p1=1,seed=1,setD=setD)
#
# # Y follows a step model
# dat_step=simu_step(0.1,n=1000,p0=0.2,p1=1,seed=1,setD)
#
# # Y follows a scaled logit model
# dat_scale=simu_logit_cov(0.1,lambda=0.5,n=1000,p0=0.2,p1=1,seed=1,setD=setD)
#
# head(dat_scale)
