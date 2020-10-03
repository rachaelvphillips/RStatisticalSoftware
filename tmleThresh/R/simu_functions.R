
library(simcausal)
#'@export
rmults <- function(n, size, prob) {
  size = size[1]
  prob = prob[1]
  apply(rmultinom(n = n, 1:size, rep(prob, size))==1, 2, which)

}
#'@export
generate_covariates<-function(){
  # create a DAG object
  D <- DAG.empty()
  D <- D +
    node("W1", distr = "rmults", size = 10 , prob = 1/10) +
    node("W2", distr = "rmults", size = 10 , prob = 1/10) +
    node("W3", distr = "rmults", size = 10 , prob = 1/10) +
    node("W4", distr = "rconst", const = W1 -3) +
    node("W5", distr = "rconst", const = W2 -3) +
    node("W6", distr = "rconst", const = W3 -3) +
    node("S", distr = "rgamma",  shape=4,rate=1)
  set.DAG(D)
}

expit <- function(x){1/(1+exp(-x))}

expit_combi<-function(lambda,beta0,beta1,w4,w5,w6,s){lambda*expit(beta0+beta1*s+w4/2+w5/2+max(w4*w5,w6)/4-max(abs(w4*w6),abs(w5*w6))/4)}

probit_combi<-function(beta0,beta1,w4,w5,w6,s){pnorm(beta0+beta1*s+w4/2+w5/2+max(w4*w5,w6)/4-max(abs(w4*w6),abs(w5*w6))/4)}

##############--------------------------------------------------------------------------------
#'@export
simu_logit_cov<-function(Risk,lambda=1,beta1=-5,n,p0,p1,seed,setD){
  # Risk: disease risk
  # beta1: slope parameter
  # lambda: scale parameter for scaled logit model
  # n: sample size
  # p0: immune response sampling rate for control group case-control study design
  # p1: immune response sampling rate for case group case-control study design
  # seed: seed for simualtion
  # setD: A DAG subject

  data_ws <- sim(setD, n = n,rndseed=seed)
  beta0=seq(-10,10,0.1)
  risk=rep(0,length(beta0))

  for (i in 1:length(beta0)){
    pp=apply(data_ws,1,function(x) expit_combi(lambda=lambda,beta0[i],beta1,x[5],x[6],x[7],x[8]))
    y=sapply(pp,function (x) rbern(1,prob=x))
    risk[i]=mean(y)
  }
  mydiff <- abs(risk-Risk)
  mymin <- min(mydiff)
  select <- head(beta0[mydiff==mymin],1)
  PY=apply(data_ws,1,function(x) expit_combi(lambda=lambda,select,beta1,x[5],x[6],x[7],x[8]))
  Y=sapply(PY,function (x) rbern(1,prob=x))
  S <- data_ws$S
  ## Creating matricies for observed data under case control, incorporating resulting weights ##
  vacdata <- cbind(data_ws,S, Y,PY)
  vacy1 <- vacdata[Y==1, ]
  vacy0 <- vacdata[Y==0, ]

  v0obs <- sample(dim(vacy0)[1], trunc(p0*(dim(vacy0)[1])), replace=TRUE)
  v1obs <- sample(dim(vacy1)[1], trunc(p1*(dim(vacy1)[1])), replace=TRUE)

  # Creating subset for subjects observed and not observed respectively under case-cohort subsampling
  vacy0obs <- data.frame(vacy0[v0obs,]) %>% mutate(S_data=1)
  vacy1obs <- data.frame(vacy1[v1obs,]) %>% mutate(S_data=1)

  vacy0not <- data.frame(vacy0[-v0obs,]) %>% mutate(S_data=0)
  vacy1not <- data.frame(vacy1[-v1obs,]) %>% mutate(S_data=0)

  vacobs <- rbind(vacy0obs, vacy0not, vacy1obs, vacy1not)

  data_obs <- vacobs %>% mutate(S_obs=ifelse(S_data==1, S, NA)) %>% select(W1,W2,W3,W4,W5,W6,S_obs, Y,S)
}



##############---------------------------------------------------------------------------
#'@export
simu_probit_cov<-function(Risk,beta1=-5,n,p0,p1,seed,setD){

  # Risk: disease risk
  # beta1: slope parameter
  # n: sample size
  # p0: immune response sampling rate for control group case-control study design
  # p1: immune response sampling rate for case group case-control study design
  # seed: seed for simualtion
  # setD: A DAG subject

  data_ws <- sim(setD, n = n,rndseed=seed)
  beta0=seq(-10,10,0.1)
  risk=rep(0,length(beta0))

  for (i in 1:length(beta0)){
    pp=apply(data_ws,1,function(x) probit_combi(beta0[i],beta1,x[5],x[6],x[7],x[8]))
    y=sapply(pp,function (x) rbern(1,prob=x))
    risk[i]=mean(y)
  }
  mydiff <- abs(risk-Risk)
  mymin <- min(mydiff)
  select <- head(beta0[mydiff==mymin],1)
  PY=apply(data_ws,1,function(x) probit_combi(select,beta1,x[5],x[6],x[7],x[8]))
  Y=sapply(PY,function (x) rbern(1,prob=x))
  S <- data_ws$S
  ## Creating matricies for observed data under case control, incorporating resulting weights ##
  vacdata <- cbind(data_ws,S, Y,PY)
  vacy1 <- vacdata[Y==1, ]
  vacy0 <- vacdata[Y==0, ]

  v0obs <- sample(dim(vacy0)[1], trunc(p0*(dim(vacy0)[1])), replace=FALSE)
  v1obs <- sample(dim(vacy1)[1], trunc(p1*(dim(vacy1)[1])), replace=FALSE)

  # Creating subset for subjects observed and not observed respectively under case-cohort subsampling
  vacy0obs <- data.frame(vacy0[v0obs,]) %>% mutate(S_data=1)
  vacy1obs <- data.frame(vacy1[v1obs,]) %>% mutate(S_data=1)

  vacy0not <- data.frame(vacy0[-v0obs,]) %>% mutate(S_data=0)
  vacy1not <- data.frame(vacy1[-v1obs,]) %>% mutate(S_data=0)

  vacobs <- rbind(vacy0obs, vacy0not, vacy1obs, vacy1not)

  data_obs <- vacobs %>% mutate(S_obs=ifelse(S_data==1, S, NA)) %>% select(W1,W2,W3,W4,W5,W6,S_obs, Y,S)
}


##############---------------------------------------------------------------------------
#'@export
simu_step<-function(Risk,n,p0,p1,seed,setD){

  # Risk: disease risk
  # n: sample size
  # p0: immune response sampling rate for control group case-control study design
  # p1: immune response sampling rate for case group case-control study design
  # seed: seed for simualtion
  # setD: A DAG subject

  data_ws <- sim(setD, n = n,rndseed=seed)
  p=seq(0,1,0.01)
  risk=rep(0,length(p))

  for (i in 1:length(p)){
    pp=apply(data_ws,1,function(x) p[i])
    y=sapply(pp,function (x) rbern(1,prob=x))
    risk[i]=mean(y)
  }
  mydiff <- abs(risk-Risk)
  mymin <- min(mydiff)
  select <- head(p[mydiff==mymin],1)
  PY=apply(data_ws,1,function(x) select)
  Y=sapply(PY,function (x) rbern(1,prob=x))
  S <- data_ws$S
  ## Creating matricies for observed data under case control, incorporating resulting weights ##
  vacdata <- cbind(data_ws,S, Y,PY)
  vacy1 <- vacdata[Y==1, ]
  vacy0 <- vacdata[Y==0, ]

  v0obs <- sample(dim(vacy0)[1], trunc(p0*(dim(vacy0)[1])), replace=FALSE)
  v1obs <- sample(dim(vacy1)[1], trunc(p1*(dim(vacy1)[1])), replace=FALSE)

  # Creating subset for subjects observed and not observed respectively under case-cohort subsampling
  vacy0obs <- data.frame(vacy0[v0obs,]) %>% mutate(S_data=1)
  vacy1obs <- data.frame(vacy1[v1obs,]) %>% mutate(S_data=1)

  vacy0not <- data.frame(vacy0[-v0obs,]) %>% mutate(S_data=0)
  vacy1not <- data.frame(vacy1[-v1obs,]) %>% mutate(S_data=0)

  vacobs <- rbind(vacy0obs, vacy0not, vacy1obs, vacy1not)

  data_obs <- vacobs %>% mutate(S_obs=ifelse(S_data==1, S, NA)) %>% select(W1,W2,W3,W4,W5,W6,S_obs, Y,S)
}




