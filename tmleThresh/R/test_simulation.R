# library(simcausal)
# library(dplyr)
#
#
# # variable: covariates W1,W2,W3. Immuresponse S. Disease Y. Sampling G.
# ## 3 cases for covariates: cov in {'non_covariate','non_confound','confound'}
# ## 4 models for P(Y|S,W): model in {'logit','probit','step','scale'}
#
# # Examples:
# ## 1. non-covariate (model:'logit')
# dat_logit<-simu(Risk=0.1,p0=1,p1=1,cov='non_covariate',model='logit')
# data1<-sim_thre(dat_logit,n = 1000,rndseed=1,p=c(0.01,0.03,0.05,0.07,0.09))
# mean(data1$simmed_data$Y)
# plotDAG(dat_logit)
#
# ## 2. non_confound (model: 'logit')
# dat_logit2<-simu(Risk=0.1,p0=1,p1=1,cov='non_confound',model='logit')
# data2<-sim_thre(dat_logit2,n = 1000,rndseed=1,p=c(0.01,0.03,0.05,0.07,0.09))
# mean(data2$simmed_data$Y)
# plotDAG(dat_logit2)
#
# ## 3. confound (model: 'probit') confounder W3
# dat_probit<-simu(Risk=0.1,p0=1,p1=1,cov='confound',model='probit')
# data3<-sim_thre(dat_probit,n = 1000,rndseed=1,p=c(0.01,0.03,0.05,0.07,0.09))
# mean(data3$simmed_data$Y)
# plotDAG(dat_probit)
#
# dat_probit
#







