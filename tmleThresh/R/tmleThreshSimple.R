
f_log_rr <- function(x, dx) {
  log(x[[2]]) - log(x[[1]])
}
df_log_rr <- function(x, dx) {
  dx[,2] / x[2] - dx[,1] / x[1]
}
rr_transform <- exp
delta_param_RR <- list(
  type = "RR",
  name = function(names) sprintf("RR(%s/%s)", names[[2]], names[[1]]),
  f = f_log_rr,
  df = df_log_rr,
  transform = rr_transform
)

RR_estimates <- function(tmle_ests) {
  lower <- tmle_ests$lower
  upper <- tmle_ests$upper
  for(i in seq_along(lower)) {
    est <- c(lower[[i]]$psi, upper[[i]]$psi)
    IC <- cbind(lower[[i]]$IC, upper[[i]]$IC)
    transf_est <- delta_param_RR$f(est, IC)
    transf_IC <- delta_param_RR$df(est, IC)
    radius <- 1.96 * apply(transf_IC,2, sd) / sqrt(nrow(transf_IC))
    CI <- delta_param_RR$transform(c(transf_est - radius, transf_est, transf_est + radius))
  }
}

thresholdTMLE <- function(data_full, node_list, thresholds_upper = NULL, thresholds_lower= NULL, biased_sampling_strata = NULL, biased_sampling_indicator = NULL, lrnr_A = Lrnr_glmnet$new(), lrnr_Y = Lrnr_glmnet$new()) {
  upper_list <- list()
  lower_list <- list()

  data_full <- as.data.table(data_full)
  data_full$id <- seq_len(nrow(data_full))
  if(!is.null(biased_sampling_indicator)) {
    data <- data_full[data_full[[biased_sampling_indicator]] ==1]
    data_full$grp <- data_full[[biased_sampling_strata]]
    biased_sampling <- TRUE
  } else {
    biased_sampling <- FALSE
    data <- data_full
  }


  for(i in seq_along(thresholds_upper)) {
    print(i)
    task_list <- get_task_list(data, node_list, thresholds_upper[i], NULL)

    preds <- get_preds(task_list, lrnr_A, lrnr_Y)

    ests <- do_update(preds, task_list, node_list)
    IC <- ests$IC

    if(biased_sampling) {
    IC_full <- matrix(0, nrow = nrow(data_full), ncol = ncol(IC))
    IC_full[data_full[[biased_sampling_indicator]] == 1,] <- IC * ests$weights
    proj_dat <- data.table(grp = data[[biased_sampling_strata]], IC = IC)
    IC_names <- setdiff(colnames(proj_dat), "grp")
    proj_dat <- proj_dat[, lapply(.SD, mean), by = "grp"]
    data_proj <- merge(data_full, proj_dat, by = "grp")
    data_proj <- data_proj[order(data_proj$id)]
    IC_proj <-  data_proj[,IC_names, with = F] * ( as.numeric(data_full[[biased_sampling_indicator]] == 1) * data_full[[node_list$weights]] - 1)
    ests$IC_IPCW <- as.matrix(IC_full)
     IC_full <- IC_full - IC_proj
    ests$IC <- as.matrix(IC_full)
    }
    upper_list[[i]] <- ests
  }

  for(i in seq_along(thresholds_lower)) {
    task_list <- get_task_list(data, node_list, NULL, thresholds_lower[i])
    preds <- get_preds(task_list, lrnr_A, lrnr_Y)
    ests <- do_update(preds, task_list, node_list)
    IC <- ests$IC

    if(biased_sampling) {
      IC_full <- matrix(0, nrow = nrow(data_full), ncol = ncol(IC))
      IC_full[data_full[[biased_sampling_indicator]] == 1,] <- IC * ests$weights
      proj_dat <- data.table(grp = data[[biased_sampling_strata]], IC = IC)
      IC_names <- setdiff(colnames(proj_dat), "grp")
      proj_dat <- proj_dat[, lapply(.SD, mean), by = "grp"]
      data_proj <- merge(data_full, proj_dat, by = "grp")
      data_proj <- data_proj[order(data_proj$id)]
      IC_proj <-  data_proj[,IC_names, with = F] * ( as.numeric(data_full[[biased_sampling_indicator]] == 1) * data_full[[node_list$weights]] - 1)
      ests$IC_IPCW <- as.matrix(IC_full)
      IC_full <- IC_full - IC_proj
      ests$IC <- as.matrix(IC_full)
    }
    upper_list[[i]] <- ests
  }
  return(list(upper = upper_list, lower = lower_list, thresholds_upper = thresholds_upper, thresholds_lower= thresholds_lower))
}


get_task_list <- function(data, node_list, threshold_upper, threshold_lower) {
  covariates <- node_list[["W"]]
  treatment <- node_list[["A"]]
  outcome <- node_list[["Y"]]
  weights <- node_list[["weights"]]
  data <- as.data.table(data)
  A <- data[[treatment]]
  if(!missing(threshold_lower) && !is.null(threshold_lower)) {
    ind_lower <- as.numeric(A < threshold_lower)
    lower_var <- paste0(treatment, "<", "l")
    data[[lower_var]] <- ind_lower
    cfl <- data.table::copy(data)
    cfl[[lower_var]]  <- 1
    task_l_Y <- sl3_Task$new(data, covariates = c(covariates, lower_var), outcome = outcome, weights = weights )
    task_l_Y_cf <- sl3_Task$new(cfl, covariates = c(covariates, lower_var), outcome = outcome, weights = weights )
    task_A_l <- sl3_Task$new(data, covariates = c(covariates), outcome = lower_var, weights = weights )

  }else {
    task_l_Y <- NULL
    task_l_Y_cf <- NULL
    task_A_l <- NULL
  }
  if( !missing(threshold_upper) && !is.null(threshold_upper)) {
    ind_upper <- as.numeric(A >= threshold_upper)
    upper_var <- paste0(treatment, ">=", "u")
    data[[upper_var]] <- ind_upper
    cfu <- data.table::copy(data)
    cfu[[upper_var]] <- 1
    task_u_Y <- sl3_Task$new(data, covariates = c(covariates, upper_var), outcome = outcome, weights = weights )
    task_u_Y_cf <- sl3_Task$new(cfu, covariates = c(covariates, upper_var), outcome = outcome, weights = weights )
    task_A_u <- sl3_Task$new(data, covariates = c(covariates), outcome = upper_var, weights = weights )

  } else {
    task_u_Y <- NULL
    task_u_Y_cf <- NULL
    task_A_u <- NULL
    }

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

  if(!is.null(task_list[["A"]][["train_l"]])) {
    lrnr_A_l <- lrnr_A$train(task_list[["A"]][["train_l"]])
    lrnr_Y_l <- lrnr_Y$train(task_list[["Y"]][["train_l"]])
    g1_l <- lrnr_A_l$predict(task_list[["A"]][["train_l"]])
    Q_l <- lrnr_Y_l$predict(task_list[["Y"]][["train_l"]])
    Q1_l <- lrnr_Y_l$predict(task_list[["Y"]][["cfl"]])

  }else {
    Q1_l <- NULL
    Q_l <- NULL
    g1_l <- NULL
  }
  if(!is.null(task_list[["A"]][["train_u"]])) {
    lrnr_A_u <- lrnr_A$train(task_list[["A"]][["train_u"]])
    lrnr_Y_u <- lrnr_Y$train(task_list[["Y"]][["train_u"]])
    g1_u <- lrnr_A_u$predict(task_list[["A"]][["train_u"]])
    Q_u <- lrnr_Y_u$predict(task_list[["Y"]][["train_u"]])
    Q1_u <- lrnr_Y_u$predict(task_list[["Y"]][["cfu"]])

  } else {
    Q1_u <- NULL
    Q_u <- NULL
    g1_u <- NULL
  }
  data.table(g1_l= g1_l, g1_u = g1_u, Q_l =  Q_l, Q_u = Q_u, Q1_l = Q1_l, Q1_u = Q1_u)
}


do_update <- function(preds, task_list, node_list) {
  data <- task_list$data
  treatment <- node_list[["A"]]
  Y <- data[[node_list[["Y"]]]]

  if(!is.null(node_list[["weights"]])) {
    weights <- data[[node_list[["weights"]]]]
  } else{
    weights <- rep(1, nrow(data))
  }
  if(!is.null(preds$g1_l)) {
    lower_var <- paste0(treatment, "<", "l")
    H_l <- as.matrix(data[[lower_var]]/preds$g1_l)
    lst <- as.data.frame(list(H_l = H_l, Y = data[[node_list[["Y"]]]]))
    eps_l <- (coef(glm(Y ~ H_l - 1, offset =  qlogis(preds[["Q_l"]]), data = lst, family = binomial(), start = rep(0, ncol(H_l)), weights = weights)))
    Q1_l <- plogis( qlogis(preds[["Q1_l"]]) + eps_l/preds$g1_l )
    Q_l <- as.vector(plogis( qlogis(preds[["Q_l"]]) + eps_l * H_l ))
    IC_l <- H_l * (Y - Q_l)
    psi_l <- weighted.mean(Q1_l, weights/sum(weights))

  } else {
    psi_l <- NULL
    IC_l <- NULL
  }

  if(!is.null(preds$g1_u)) {
    upper_var <- paste0(treatment, ">=", "u")
    H_u <- as.matrix(data[[upper_var]]/preds$g1_u)
    eps_u <- suppressWarnings(coef(glm(Y ~ H_u - 1, offset =  qlogis(preds[["Q_u"]]), data = list(H_u = H_u, Y = data[[node_list[["Y"]]]]), family = binomial(), start = rep(0, ncol(H_u)), weights = weights)))
    Q1_u <- plogis( qlogis(preds[["Q1_u"]]) + eps_u/preds$g1_u )
    Q_u <- as.vector(plogis( qlogis(preds[["Q_u"]]) + eps_u * H_u))
    IC_u <- H_u * (Y - Q_u)
    psi_u <- weighted.mean(Q1_u, weights/sum(weights))

  }  else {
    psi_u <- NULL
    IC_u <- NULL
  }


  estimates <- list(psi = c(psi_l, psi_u), IC = cbind(IC_l, IC_u), weights = weights)

  return(estimates)
}
