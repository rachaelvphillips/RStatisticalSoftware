
rmults <- function(n, size, prob) {
  size = size[1]
  prob = prob[1]
  apply(rmultinom(n = n, 1:size, rep(prob, size))==1, 2, which)

}

adjusted_sim <- function(n = 1000){
  library(simcausal)
  D <- DAG.empty()
  D <- D +
    node("W1", distr = "rmults", size = 15 , prob = 1/15) +
    node("W2", distr = "rmults", size = 15 , prob = 1/15) +
    node("W3", distr = "rmults", size = 15 , prob = 1/15) +
    node("Wa", distr = "rconst", const = W1 -7.5) +
    node("Wb", distr = "rconst", const = W2 -7.5) +
    node("Wc", distr = "rconst", const = W3 -7.5) +
    node("S", distr = "runif",  min = -5, max = 5) +
    node("Y", distr = "rbinom", size =1 , prob =plogis(S + Wa/4 + Wb/4 - max(min(Wc*Wb/10, abs(Wc)), -abs(Wc)))) +
    node("pY", distr = "rconst",const = plogis(S + Wa/4 + Wb/4 - max(min(Wc*Wb/10, abs(Wc)), -abs(Wc))))
  setD <- set.DAG(D)
  dat <- sim(setD, n = n)
  dat

  return(dat)
}

basic_sim <- function(n = 1000){
  library(simcausal)
  D <- DAG.empty()
  D <- D +
    node("W1", distr = "rmults", size = 15 , prob = 1/15) +
    node("W", distr = "rconst", const = W1 + 1) +
    node("S", distr = "runif",  min = -5, max = 5) +
    node("Y", distr = "rbinom", size =1 , prob = plogis(S)) +
    node("pY", distr = "rconst",const = plogis((S)))


  setD <- set.DAG(D)
  dat <- sim(setD, n = n)
  dat

  return(dat)
}

compute_estimates <- function(n,sim_func, cutoffs) {
  get_val <- function(cutoff) {
    dat <- sim_func(n)
    return(mean(dat$Y[(dat$S >= cutoff)]))
  }
  ests <- lapply(cutoffs, get_val)
  names(ests) <- cutoffs
  return(ests)
}
