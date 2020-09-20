
make_revere <- function(task, likelihood, type = "gen") {
  genf <- make_generator(likelihood, type)
  return(sl3_revere_Task$new(genf, task))
}

derived_generator <- function(type = "univ") {
  gen <- function(tmle_task, likelihood) {
    return(make_revere(tmle_task, likelihood, type = type))
  }
  return(gen)
}


make_generator <- function(likelihood, type = "IPW") {

  gen_task_univ <- function(tmle_task, fold_number) {
    task <- tmle_task$get_regression_task("RR")
    X <- task$X
    g <- likelihood$get_likelihood(tmle_task, "A", fold_number)
    A <- tmle_task$get_tmle_node("A", format = T)$A
    Q <- likelihood$get_likelihood(tmle_task, "R", fold_number)
    cf_task1 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=1))
    cf_task0 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=0))
    ER1 <- likelihood$get_likelihood(cf_task1, "R", fold_number)
    ER0 <- likelihood$get_likelihood(cf_task0, "R", fold_number)

    data <- cbind(data.table(Y = task$Y, A = A, Q = Q, g = 1/g, ER1 = ER1, ER0 = ER0, weights = task$weights), X)
    covariates <- c(colnames(X), "A")
    outcome  <- "Y"
    new_task <- sl3_Task$new(data, covariates = covariates, outcome = outcome,  weights = "weights", folds = task$folds)
    return(new_task)
  }

  gen_task_IPW <- function(tmle_task, fold_number) {
    task <- tmle_task$get_regression_task("RR")
    X <- task$X
    R <- tmle_task$get_tmle_node("R", format = T)$R
    A <- tmle_task$get_tmle_node("A", format = T)$A
    g <- likelihood$get_likelihood(tmle_task, "A", fold_number)
    Q <- likelihood$get_likelihood(tmle_task, "R", fold_number)
    cf_task1 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=1))
    cf_task0 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=0))
    ER1 <- likelihood$get_likelihood(cf_task1, "R", fold_number)
    ER0 <- likelihood$get_likelihood(cf_task0, "R", fold_number)
    weights <- R/g
    covariates <- colnames(X)
    outcome <- "Y"
    outcome_type <- "binomial"
    data <- cbind(data.table(Y = A,  weightsold = task$weights, weights = weights, Q = Q, g = g, R = R, A = A, ER1 = ER1, ER0 = ER0), X)
    new_task <- sl3_Task$new(data, covariates = covariates, outcome = outcome, outcome_type = outcome_type, weights = "weightsold", folds = task$folds)
    return(new_task)
  }

  gen_task_plugin <- function(tmle_task, fold_number) {
    task <- tmle_task$get_regression_task("RR")
    X <- task$X
    R <- tmle_task$get_tmle_node("R", format = T)$R
    A <- tmle_task$get_tmle_node("A", format = T)$A
    g <- likelihood$get_likelihood(tmle_task, "A", fold_number)
    Q <- likelihood$get_likelihood(tmle_task, "R", fold_number)
    cf_task1 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=1))
    cf_task0 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=0))
    ER1 <- likelihood$get_likelihood(cf_task1, "R", fold_number)
    ER0 <- likelihood$get_likelihood(cf_task0, "R", fold_number)
    weights <- (ER1 + ER0)
    Y <- ER1/weights
    covariates <- colnames(X)
    outcome <- "Y"
    outcome_type <- "continuous"
    data <- cbind(data.table(Y = Y, weightsold = task$weights,  weights = weights, Q = Q, g = g, R = R, A = A, ER1 = ER1, ER0 = ER0), X)
    new_task <- sl3_Task$new(data, covariates = covariates, outcome = outcome, outcome_type = outcome_type, weights = "weightsold", folds = task$folds)
    return(new_task)
  }

  gen_task_general <- function(tmle_task, fold_number) {
    task <- tmle_task$get_regression_task("RR")
    X <- task$X
    R <- tmle_task$get_tmle_node("R", format = T)$R
    A <- tmle_task$get_tmle_node("A", format = T)$A
    g <- likelihood$get_likelihood(tmle_task, "A", fold_number)
    Q <- likelihood$get_likelihood(tmle_task, "R", fold_number)
    cf_task1 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=1))
    cf_task0 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=0))
    ER1 <- likelihood$get_likelihood(cf_task1, "R", fold_number)
    ER0 <- likelihood$get_likelihood(cf_task0, "R", fold_number)
    weightsIPW <- R/g * task$weights
    weightsplugin <- (ER1 + ER0)
    YIPW <- A
    Yplugin <- ER1/weightsplugin
    covariates <- colnames(X)
    outcome <- c("Yplugin", "YIPW")
    weights <- c(tmle_task$nodes$weights)
    data <- cbind(data.table(YIPW = YIPW, Yplugin = Yplugin, weightsplugin = weightsplugin, weightsIPW = weightsIPW, Q = Q, g = g, R = R, A = A, ER1 = ER1, ER0 = ER0), X)
    data <- cbind(task$get_data(,weights), data)
    new_task <- sl3_Task$new(data, covariates = covariates, outcome = outcome, weights = weights, folds = task$folds)
    return(new_task)
  }

  if(type == "IPW") {
    return(gen_task_IPW)
  } else if (type == "plugin") {
    return(gen_task_plugin)
  } else if(type == "univ"){
    return(gen_task_univ)
    }
  else {
    return(gen_task_general)
  }

}

make_eff_loss <- function(tmle_task, likelihood) {
  tmle_task <- tmle_task
  likelihood <- likelihood

  efficient_loss = function(preds, Y) {

    LRR <- preds

    R <- tmle_task$get_tmle_node("R", format = T)$R
    A <- tmle_task$get_tmle_node("A", format = T)$A

    g <- likelihood$get_likelihood(tmle_task, "A", "validation")
    lik <- likelihood

    cf_task1 <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A=1))
    cf_task0 <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A=0))
    ER1 <- lik$get_likelihood(cf_task1, "R", "validation")
    ER0 <- lik$get_likelihood(cf_task0, "R", "validation")
    ER <- lik$get_likelihood(tmle_task, "R", "validation")

    C1 <- A/g * (R - ER) + ER1
    C2 <- C1 + (1-A)/g * (R - ER) + ER0
    loss <- C1*-1*LRR + C2 * log(1 + exp(LRR))


    return(loss)
  }

}
