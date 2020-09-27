context("Basic running of program")


library(data.table)
library(origami)
library(delayed)
library(sl3)
library(hal9001fast)
library(tmle3)
rmults <- function(n, size, prob) {
  size = size[1]
  prob = prob[1]
  apply(rmultinom(n = n, 1:size, rep(prob, size))==1, 2, which)

}

n = 15000
library(simcausal)
D <- DAG.empty()
D <- D +
  node("W1", distr = "rmults", size = 15 , prob = 1/15) +
  node("W", distr = "rconst", const = W1 + 1) +
  node("S", distr = "runif",  min = -5, max = 5) +
  node("Y", distr = "rbinom", size =1 , prob = plogis((S ))) +
  node("pY", distr = "rconst",const = plogis((S )))


setD <- set.DAG(D)
dat <- sim(setD, n = n)

D <- DAG.empty()
D <- D +
  node("W1", distr = "rmults", size = 15 , prob = 1/15) +
  node("W", distr = "rconst", const = W1 + 1) +
  node("S", distr = "runif",  min = -5, max = 5) +
  node("Y", distr = "rbinom", size =1 , prob = plogis((S ))) +
  node("pY", distr = "rconst",const = plogis((S )))


setD <- set.DAG(D)
truedata <- sim(setD, n = 50000)



cutoffs <-  quantile(dat$S, c(0.05, 0.1, 0.2, 0.4, 0.6))
results <- c()
for(cutoff in cutoffs) {
  res <- sum(truedata$Y * (truedata$S >= cutoff)) / sum(truedata$S >= cutoff)
  results <- c(results, res)
}



npsem <- list(define_node("W", "W", c()),
              define_node("A_learned", c("S"), c("A")),
              define_node("A", c("S"), c()),
              define_node("Y", "Y", c("A")))

task <- tmle3_Task$new(as.data.table(dat), npsem)
lrnr_bin <- make_learner(Lrnr_glm_fast)
cv_lrnr <- make_learner(Pipeline, Lrnr_cv$new(make_learner(Lrnr_glm), full_fit = T))

factor_list <- list(W = LF_emp$new("W"), A_learned = LF_fit$new("A_learned", cv_lrnr, type = "mean"))
likelihood <- Likelihood$new(factor_list)
likelihood <- likelihood$train(task)
generator_R <- learner_marker_task_generator(learned_marker_node = "A_learned", node = "Y", data_adaptive = F)

generator_A <- learner_marker_task_generator(learned_marker_node = "A_learned", node = "A", data_adaptive = F)
tmp <- (generator_R(task, likelihood)$revere_fold_task("full"))

lf_A <- LF_derived2$new("A", Lrnr_CDF$new(make_learner(Lrnr_xgboost),
                                         10,  cutoffs, cv = F), likelihood, generator_A, type = "mean")
lf_R <- LF_derived2$new("Y", Lrnr_thresh$new(Lrnr_cv$new(make_learner(Pipeline, lrnr_bin, make_learner(Lrnr_wrapper, 5))),
                                            strata_variable = "A",
                                            cutoffs = cutoffs), likelihood, generator_R, type = "mean")

lf_R <- LF_derived2$new("Y", Lrnr_thresh$new( make_learner(Lrnr_xgboost),
                                              strata_variable = "A", cv = F,
                                              cutoffs = cutoffs), likelihood, generator_R, type = "mean")



likelihood$add_factors(list(lf_A, lf_R))
cutoffs
cf_task <- task

cf_data <- data.table(rep(10 + max(cutoffs), cf_task$nrow))

setnames(cf_data, "A")
cf_data$id <- cf_task$id
cf_data$t <- cf_task$time
cf_task <- cf_task$generate_counterfactual_task(UUIDgenerate(), cf_data)

Yout <- as.data.table(matrix(likelihood$get_likelihood(task, "Y", fold_number = "full"), nrow = n))
Yout_cf <- as.data.table(matrix(likelihood$get_likelihood(cf_task, "Y", fold_number = "full"), nrow = n))
head(Yout)
head(Yout_cf)

Aout <- as.data.table(matrix(likelihood$get_likelihood(task, "A", fold_number = "full"), nrow = n))
print(unique(Aout))
# These should  be close



np_results <- c()
for(cutoff in cutoffs) {
  res <- sum(dat$Y * (dat$S >= cutoff)) / sum(dat$S >= cutoff)
  np_results <- c(np_results, res)
}
# These should all be close
out <- colMeans(as.matrix(Yout_cf))
out

np_results
assertthat::assert_that(max(abs(out - np_results)) < 0.01)

results <- c()
for(cutoff in cutoffs) {
  res <- mean(dat$S < cutoff)
  results <- c(results, res)
}
results
out <- colMeans(Aout)
out
assertthat::assert_that(max(abs(out - results)) < 0.01)


length(unlist(Aout))

