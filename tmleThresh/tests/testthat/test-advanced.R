
rmults <- function(n, size, prob) {
  size = size[1]
  prob = prob[1]
  apply(rmultinom(n = n, 1:size, rep(prob, size))==1, 2, which)

}
library(tmleThresh)
library(tmle3)
library(simcausal)
library(sl3)
options("sl3.save.training" = TRUE)
getOption("sl3.save.training")
n = 5000
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
data <- sim(setD, n = n)
data <- as.data.table(data)
node_list <- list("W" = c("Wa", "Wb", "Wc"), "A" = "S", "Y" = "Y")
npsem <- make_thresh_npsem(node_list)
task <- make_thresh_task(data, npsem)
xg <- Lrnr_xgboost$new()
hal <- Lrnr_hal9001$new(max_degree = 2)
lrnr_glm <- Lrnr_glm$new()
learner_list <- list("A" = xg, "Y" =hal)
A <- task$get_tmle_node("A")
cutoffs <- as.vector(quantile(A, seq(0.05, 0.95, length.out = 5)))
lik <- make_thresh_likelihood(task, learner_list = learner_list, cv = F, cutoffs= cutoffs)
lik$factor_list$Y
lik$factor_list
head(compute_thresh_estimate(lik, task))


cutoffs

for(cutoff in cutoffs) {
  keep <- data$S >=cutoff
  print(c(mean(data$pY[keep]), mean(data$Y[keep])))

}
colMeans(matrix(compute_thresh_estimate(lik, task), ncol = 5))
head(matrix(compute_thresh_estimate(lik, task), ncol = 5))








