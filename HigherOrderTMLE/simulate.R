
truth <- function() {
  library(simcausal)


  D <- DAG.empty()
  D <- D +
    node("W", distr = "runif", min = -0.8, max = 0.8) +
    node("W1", distr = "runif", min = -1, max = 1) +
    node("A", distr = "rbinom", size = 1,  prob = plogis(W1)) +
    node("g1", distr = "rconst", const = plogis(W1)) +
    node("Y", distr = "rbinom", size =1 , prob = plogis(( -1.5 + 1 + A + W - A*W1/2 ))) +
    node ("EY1", distr = "rconst", const = plogis(( -1.5 + 1 + 1 + W - 1*W1/2 )))
  setD <- set.DAG(D)
  data <- sim(setD, n = 10000)
  data <- as.data.table(data)
  truth <- mean(data$EY1)
}
simulate <- function(n = 1000){

verbose = T


D <- DAG.empty()
D <- D +
  node("W", distr = "runif", min = -0.8, max = 0.8) +
  node("W1", distr = "runif", min = -1, max = 1) +
  node("A", distr = "rbinom", size = 1,  prob = plogis(W1)) +
  node("g1", distr = "rconst", const = plogis(W1)) +
  node("Y", distr = "rbinom", size =1 , prob = plogis(( -1.5 + 1 + A + W - A*W1/2 ))) +
  node ("EY1", distr = "rconst", const = plogis(( -1.5 + 1 + 1 + W - 1*W1/2 )))
setD <- set.DAG(D)
data <- sim(setD, n = n)
data <- as.data.table(data)


D <- DAG.empty()
D <- D +
  node("W", distr = "runif", min = -0.8, max = 0.8) +
  node("W1", distr = "runif", min = -1, max = 1) +
  node("A", distr = "rbinom", size = 1,  prob = plogis(W1)) +
  node("g1", distr = "rconst", const = plogis(W1)) +
  node("Z", distr = "rnorm", mean = 1) +
  node("Y", distr = "rbinom", size =1 , prob = plogis(( -1.5 + 1 + A + W - A*W1/2 ) + Z/(n)^(0.25))) +
  node ("EY1", distr = "rconst", const = plogis(( -1.5 + 1 + 1 + W - 1*W1/2 )))
setD <- set.DAG(D)
data_bias <- sim(setD, n = n)
data_bias <- as.data.table(data_bias)

npsem <- list(
  define_node("W", c("W", "W1"), c()),
  define_node("A", "A",  c("W")),
  define_node("Y",  "Y", c("A", "W"))
)

task <- tmle3_Task$new(data, npsem)
task_bias <- tmle3_Task$new(data_bias, npsem)

lrnr_glm <- Lrnr_glm$new()
lrnr_mean <- Lrnr_xgboost$new(max_depth = 4)
stopping_criterion <- function(n, residual) {return(min(log(n)/n, 1/ sqrt(n) / log(n)))}

lrnr_undersmooth <- Lrnr_undersmooth$new(stopping_criterion = stopping_criterion)

lrnr_undersmoothA <- lrnr_undersmooth
lrnr_undersmoothY <- lrnr_undersmooth

factor_list <- list(
  LF_emp$new("W"),
  LF_fit$new("A",lrnr_glm),
  LF_fit$new("Y",lrnr_mean, type = "mean")
)

factor_list_tilde <- list(
  LF_emp$new("W"),
  LF_fit$new("A",lrnr_undersmoothA ),
  LF_fit$new("Y",lrnr_undersmoothY, type = "mean" )
)

task$get_regression_task("A")$data
task$get_regression_task("Y")$data

likelihood <- Likelihood$new(factor_list)
likelihood_tilde <- Likelihood$new(factor_list_tilde)

if(verbose) {
  print("Training likelihoods")
}
likelihood <- likelihood$train(task)

likelihood_tilde <- likelihood_tilde$train(task)

if(verbose) {
  print("Done")
}

tlik <- Targeted_Likelihood$new(likelihood, updater = list(maxit = 300, convergence_type = "sample_size", cvtmle = F, constrain_step = T, one_dimensional = T, delta_epsilon = c( 0.01)))
param <- Param_TSM_higher_order$new(tlik, likelihood_tilde)
print(length(tlik$cache$tasks))

if(verbose) {
  print("Second order updating")
}

est <- param$estimates(task)
if(verbose) {
  print(  mean(est$IC))
  print(   mean(est$psi))
  print(  mean(data$EY1))
}

tlik$updater$update(tlik, task)
for(i in 1:5) {
  tlik$updater$update_step(tlik, task)
}
if(verbose) {
  print("Done")
}

print(tlik$updater$step_number)
print(tlik$updater$check_convergence(task))

#lapply(param$clever_covariates(task, for_fitting = T), data.table)
est <- param$estimates(task)
if(verbose) {
  print(  mean(est$IC))
  print(   mean(est$psi))
  print(  mean(data$EY1))
}


tlik2 <- Targeted_Likelihood$new(tlik, updater = list(cvtmle = F))
param <- Param_TSM_2$new(tlik2)

tlik2$updater$update(tlik2, task)

#lapply(param$clever_covariates(task, for_fitting = T), data.table)
est <- param$estimates(task, fold_number = "full")
IC_higher <- mean(est$IC)
psi_higher <- mean(est$psi)



tlik_old <- Targeted_Likelihood$new(likelihood, updater = list(cvtmle = F))
param_old <- Param_TSM_2$new(tlik_old)
tlik_old$updater$update(tlik_old, task)

est <- param_old$estimates(task, fold_number = "full")
IC_old <- mean(est$IC)
psi_old <- mean(est$psi)


print(c(psi_higher, psi_old))
return(c(psi_higher, psi_old))
}



