

simu<-function(Risk,p0,p1,cov,model){
  # Arguments:
  # Risk: disease risk
  # n: sample size
  # p0: immune response sampling rate for control group case-control study design
  # p1: immune response sampling rate for case group case-control study design
  # cov: whether inclduing covariates and how to use covariates
  # model: for P(Y|S=v,W=w)
  #------------------------------------------------------------------------
  # output: A DAG subject
  #----------------------------------
  D<-generate_covariates()
  D<-generate_S(D,cov)
  if(cov=='non_covariate'){
    fun_y<-eval(parse(text='generate_Y_noncov'))
    D<-fun_y(Risk,D,model,p0,p1)
  }
  if(cov=='non_confound'){
    fun_y<-eval(parse(text='generate_Y_nonconfound'))
    D<-fun_y(Risk,D,model,p0,p1)
  }
  if(cov=='confound'){
    fun_y<-eval(parse(text='generate_Y_confound'))
    D<-fun_y(Risk,D,model,p0,p1)
  }
  D
}



rmults <- function(n, size, prob) {
  size = size[1]
  prob = prob[1]
  apply(rmultinom(n = n, 1:size, rep(prob, size))==1, 2, which)

}

logit <- function(x){1/(1+exp(-x))}
scale <- function(x){0.7*1/(1+exp(-x))}
probit<-function(x){pnorm(x)}


generate_covariates<-function(){
  # create a DAG object
  D <- DAG.empty()
  D <- D +
    node("A", distr = "rmults", size = 10 , prob = 1/10) +
    node("B", distr = "rmults", size = 10 , prob = 1/10) +
    node("C", distr = "rmults", size = 10 , prob = 1/10) +
    node("W1", distr = "rconst", const = A -5) +
    node("W2", distr = "rconst", const = B -5) +
    node("W3", distr = "rconst", const = C -5)
  D
}


generate_S<-function(D,cov){
  if(cov=='non_covariate') {D<-D+node("S", distr = "rgamma",  shape=4,rate=1)}
  else if(cov=='non_confound') {D<-D+node("S", distr = "rgamma",shape=abs(W3)+1,rate=1)}
  else if(cov=='confound') {D<-D+node("S", distr = "rgamma",shape=abs(W3)+1,rate=1)} ##
  D
}

generate_Y_nonconfound<-function(Risk,D,model,p0,p1){
  fun=eval(parse(text=model))
  beta0<-seq(-5,5,0.1)
  if(model=='step') beta0<-seq(0.01,0.99,0.1)
  risk<-rep(0,length(beta0))
  for(i in 1:length(beta0)){
    new_D<-D
    if(model!='step'){
      new_D<-D+node("PY", distr = "rconst", const = fun(eval(beta0[i])+(S-3)/6))+
        node("Y",distr='rbern',prob=PY)
    }
    else{
      new_D<-D+node("PY",distr = "rconst", const = eval(beta0[i]))+
        node("Y",distr='rbern',prob=PY)
    }
    setD<-set.DAG(new_D)
    data<- sim(setD, n = 5000,rndseed=1)
    risk[i]<-mean(data$Y)
  }
  mydiff <- abs(risk-Risk)
  mymin <- min(mydiff)
  select <- head(beta0[mydiff==mymin],1)
  if(model!='step'){
    D<-D+node("PY", distr = "rconst", const = fun(eval(beta0[i])+(S-3)/6))+
      node("Y",distr='rbern',prob=PY)
  }
  else{
    D<-D+node("PY",distr = "rconst", const = eval(select))+
      node("Y",distr='rbern',prob=PY)
  }
  D<-D+node('G',distr='rbern',prob=(Y==1)*eval(p1)+(Y==0)*eval(p0))
  set.DAG(D)
}

generate_Y_noncov<-function(Risk,D,model,p0,p1){
  fun=eval(parse(text=model))
  beta0<-seq(-10,10,1)
  if(model=='step') beta0<-seq(0.01,0.99,0.1)
  risk<-rep(0,length(beta0))
  for(i in 1:length(beta0)){
    new_D<-D
    if(model!='step'){
      new_D<-D+node("PY", distr = "rconst", const = fun(eval(beta0[i])-5*S))+
        node("Y",distr='rbern',prob=PY)
    }
    else{
      new_D<-D+node("PY",distr = "rconst", const = eval(beta0[i]))+
        node("Y",distr='rbern',prob=PY)
    }
    setD<-set.DAG(new_D)
    data<- sim(setD, n = 50000,rndseed=1)
    risk[i]<-mean(data$Y)
  }
  mydiff <- abs(risk-Risk)
  mymin <- min(mydiff)
  select <- head(beta0[mydiff==mymin],1)
  if(model!='step'){
    D<-D+node("PY", distr = "rconst", const = fun(eval(select)-5*S))+
      node("Y",distr='rbern',prob=PY)
  }
  else{D<-D+node("PY",distr = "rconst", const = eval(select))+
    node("Y",distr='rbern',prob=PY)
  }

  D<-D+node('G',distr='rbern',prob=(Y==1)*eval(p1)+(Y==0)*eval(p0))
  setD<-set.DAG(D)
}

generate_Y_confound<-function(Risk,D,model,p0,p1){
  fun=eval(parse(text=model))
  beta0<-seq(-10,10,1)
  if(model=='step') beta0<-seq(0.01,0.99,0.1)
  risk<-rep(0,length(beta0))
  for(i in 1:length(beta0)){
    new_D<-D
    if(model!='step'){
      new_D<-D+node("PY", distr = "rconst", const = fun(eval(beta0[i])-5*S+W1/2+W2/2+W3/2-S^2*(W3>0)-S^2*(W1>0)-S^2*(W2>0)))+
        node("Y",distr='rbern',prob=PY)
    }
    else{
      new_D<-D+node("PY",distr = "rconst", const = eval(beta0[i]))+
        node("Y",distr='rbern',prob=PY)
    }
    setD<-set.DAG(new_D)
    data<- sim(setD, n = 50000,rndseed=1)
    risk[i]<-mean(data$Y)
  }
  mydiff <- abs(risk-Risk)
  mymin <- min(mydiff)
  select <- head(beta0[mydiff==mymin],1)
  if(model!='step'){
    D<-D+node("PY", distr = "rconst", const = fun(eval(select)-5*S+W1/2+W2/2+W3/2-S^2*(W3>0)-S^2*(W1>0)-S^2*(W2>0)))+
      node("Y",distr='rbern',prob=PY)
  }
  else{
    D<-D+node("PY",distr = "rconst", const = eval(select))+
      node("Y",distr='rbern',prob=PY)
  }
  D<-D+node('G',distr='rbern',prob=(Y==1)*eval(p1)+(Y==0)*eval(p0))
  set.DAG(D)
}


sim_thre<-function(object,n,rndseed=1,p){
  data<-sim(object,n=n ,rndseed=rndseed)
  truethres <- rep(NA, length=length(p))
  v <- seq(0,8,by=.01)
  nv <- length(v)
  risks <- vector(, length=nv)
  for(i in 1:nv){
    risks[i]<-mean(data[data$S>=v[i],]$PY)
  }

  for(j in 1:length(p)){
    mydiff <- abs(risks-p[j])
    mymin <- min(mydiff)
    index <- which(mydiff==mymin)
    truethres[j] <- min(v[index])
  }
  output<-list(data, truethres)
  names(output)[1] <- "simmed_data"
  names(output)[2] <- "true_thresholds"
  names(output[[2]]) <- as.character(p)
  return(output)
}
