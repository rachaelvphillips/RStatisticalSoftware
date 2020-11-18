

get_task_list <- function(data, node_list, threshold_upper, threshold_lower = threshold_upper) {
  covariates <- node_list[["W"]]
  treatment <- node_list[["A"]]
  outcome <- node_list[["Y"]]
  weights <- node_list[["weights"]]
  data <- as.data.table(data)
  A <- data[[treatment]]
  if(!missing(threshold_lower)) {
    ind_lower <- as.numeric(A < threshold_lower)
  }
  if(!missing(threshold_upper)) {
    ind_upper <- as.numeric(A >= threshold_upper)
  }
  upper_var <- paste0(treatment, ">=", "u")
  lower_var <- paste0(treatment, "<", "l")
  data[[upper_var]] <- ind_upper
  data[[lower_var]] <- ind_lower
  print(data)
  cfu <- data.table::copy(data)
  cfl <- data.table::copy(data)
  cfu[[upper_var]] <- 1
  cfl[[lower_var]]  <- 1
  task_u_Y <- sl3_Task$new(data, covariates = c(covariates, upper_var), outcome = outcome, weights = weights )

  task_l_Y <- sl3_Task$new(data, covariates = c(covariates, lower_var), outcome = outcome, weights = weights )

  task_u_Y_cf <- sl3_Task$new(cfu, covariates = c(covariates, upper_var), outcome = outcome, weights = weights )

  task_l_Y_cf <- sl3_Task$new(cfl, covariates = c(covariates, lower_var), outcome = outcome, weights = weights )

  task_A_u <- sl3_Task$new(data, covariates = c(covariates), outcome = upper_var, weights = weights )
  task_A_l <- sl3_Task$new(data, covariates = c(covariates), outcome = lower_var, weights = weights )

  task_list <- list(data = data, Y = list(train_l = task_l_Y, train_u = task_u_Y, cfl = task_l_Y_cf, cfu = task_u_Y_cf), A = list(train_u = task_A_u, train_l = task_A_l))
  return(task_list)
}

get_preds <- function(task_list, lrnr_A = NULL, lrnr_Y = NULL) {
  if(is.null(lrnr_A))  {
    lrnr_A <- Lrnr_glmnet$new()
  }
  if(is.null(lrnr_Y))  {
    lrnr_Y <- Lrnr_glmnet$new()
  }

  lrnr_A_l <- lrnr_A$train(task_list[["A"]][["train_l"]])
  lrnr_A_u <- lrnr_A$train(task_list[["A"]][["train_u"]])
  lrnr_Y_l <- lrnr_Y$train(task_list[["Y"]][["train_l"]])
  lrnr_Y_u <- lrnr_Y$train(task_list[["Y"]][["train_u"]])

  g1_l <- lrnr_A_l$predict(task_list[["A"]][["train_l"]])
  g1_u <- lrnr_A_l$predict(task_list[["A"]][["train_u"]])
  Q_l <- lrnr_Y_l$predict(task_list[["Y"]][["train_l"]])
  Q_u <- lrnr_Y_u$predict(task_list[["Y"]][["train_u"]])
  Q1_l <- lrnr_Y_l$predict(task_list[["Y"]][["cfl"]])
  Q1_u <- lrnr_Y_u$predict(task_list[["Y"]][["cfu"]])
  data.table(g1_l= g1_l, g1_u = g1_u, Q_l =  Q_l, Q_u = Q_u, Q1_l = Q1_l, Q1_u = Q1_u)
}

do_update <- function(preds, task_list, node_list) {
  data <- task_list$data
  treatment <- node_list[["A"]]
  Y <- data[[node_list[["Y"]]]]

  if(!is.null(node_list[[weights]])) {
    weights <- data[[node_list[[weights]]]]
  } else{
    weights <- rep(1, nrow(data))
  }
  upper_var <- paste0(treatment, ">=", "u")
  lower_var <- paste0(treatment, "<", "l")
  H_l <- as.matrix(data[[lower_var]]/preds$g1_l)
  H_u <- as.matrix(data[[upper_var]]/preds$g1_u)
  eps_l <- suppressWarnings(coef(glm(Y ~ H_l - 1, offset =  qlogis(preds[["Q1_l"]]), data = list(H_l = H_l, Y = data[[node_list[["Y"]]]]), family = binomial(), start = rep(0, ncol(H_l)), weights = weights)))
  eps_u <- suppressWarnings(coef(glm(Y ~ H_u - 1, offset =  qlogis(preds[["Q1_u"]]), data = list(H_u = H_u, Y = data[[node_list[["Y"]]]]), family = binomial(), start = rep(0, ncol(H_l)), weights = weights)))

  Q1_l <- plogis( qlogis(preds[["Q1_l"]]) + eps_l/preds$g1_l )
  Q1_u <- plogis( qlogis(preds[["Q1_u"]]) + eps_u/preds$g1_u )
  Q_l <- as.vector(plogis( qlogis(preds[["Q_l"]]) + eps_l * H_l ))
  Q_u <- as.vector(plogis( qlogis(preds[["Q_u"]]) + eps_u * H_u))
  new_preds <- data.table(g1_l = preds$g1_l,  g1_u = preds$g1_u,Q_u = Q_u, Q_l = Q_l, Q1_l = Q1_l, Q1_u = Q1_u)
  IC_l <- H_l * (Y - new_preds$Q_l) * weights
  IC_u <- H_u * (Y - new_preds$Q_u) * weights
  psi_l <- weighted.mean(Q1_l, weights/sum(weights))
  psi_u <- weighted.mean(Q1_u, weights/sum(weights))

  estimates <- list(psi = c(psi_l, psi_u), IC = cbind(IC_l, IC_u) )

  return(estimates)
}



