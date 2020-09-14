#' @import hal9001fast
#' @importFrom sl3 args_to_list
#' @importFrom uuid UUIDgenerate
#' @export
LRR_risk <- R6Class(
  classname = "LRR_risk",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(likelihood,  smoothness = 0, max_degree = 3, IPW_as_initial = T, type = c("efficient", "IPW", "grad"), grad_use_best = F, num_fits = 1,...){
      type <- match.arg(type)
      params <- sl3::args_to_list()
      params$type <- type
      private$.params <- params
      private$.useIPW <- IPW_as_initial
      private$.tmle_param <- Param_RR$new(likelihood, risk_function = self)
    },
    design_matrix = function(tmle_task, basis_list = NULL) {
       X <- tmle_task$get_tmle_node("W", format = T)
       task <- tmle_task$get_regression_task("RR")
        X <- task$X
       smoothness <- self$params$smoothness
       smoothness <- rep(c(smoothness)[[1]], ncol(X))
       if(is.null(basis_list)) {
         basis_list <- hal9001fast::enumerate_basis(as.matrix(X), max_degree = self$params$max_degree, order_map = c(smoothness), bins = rep(350, ncol(X)), include_zero_order = F, include_lower_order = F)
       }
       x_basis <- hal9001fast::make_design_matrix(as.matrix(X), basis_list)
       return(list(x_basis = x_basis, basis_list = basis_list))

    },
    IPW_risk = function(tmle_task, fold_number = "full", fit = NULL) {
      task <- tmle_task$get_regression_task("RR")
      R <- tmle_task$get_tmle_node("R", format = T)$R
      A <- tmle_task$get_tmle_node("A", format = T)$A
      g <- self$likelihood$get_likelihood(tmle_task, "A", fold_number)
      keep <- R==1
      weights <- (R/g)
      Y <- A
      f <- self$predict(tmle_task, type = "LRR", fit = fit)

      loss <- R * (-A*f + log(1 + exp(f)))/ g
      return(mean(loss))
    },
    plugin_risk = function(tmle_task, fold_number = "full") {
      lik <- self$likelihood
      cf_task1 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=1))
      cf_task0 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=0))
      ER1 <- lik$get_likelihood(cf_task1, "R", fold_number)
      ER0 <- lik$get_likelihood(cf_task0, "R", fold_number)
      ER <- lik$get_likelihood(tmle_task, "R", fold_number)
      LRR <- self$predict(tmle_task, type = "LRR")

      loss <- ER1 * -1 * LRR + (ER1 + ER0) * log(1 + exp(LRR))
      return(mean(loss))
    },
    efficient_risk = function(tmle_task, fold_number = "full", cache = F, fit = NULL) {
      #basis_list <- private$.basis_list
      #beta <- as.vector(private$.beta)
      #lst <- self$design_matrix(tmle_task, basis_list = basis_list)
      #x_basis <- lst$x_basis
      R <- tmle_task$get_tmle_node("R", format = T)$R
      A <- tmle_task$get_tmle_node("A", format = T)$A
      g <- self$likelihood$get_likelihood(tmle_task, "A", fold_number)
      lik <- self$likelihood
      cf_task1 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=1))
      cf_task0 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=0))
      ER1 <- lik$get_likelihood(cf_task1, "R", fold_number)
      ER0 <- lik$get_likelihood(cf_task0, "R", fold_number)
      ER <- lik$get_likelihood(tmle_task, "R", fold_number)
      LRR <- self$predict(tmle_task, type = "LRR", fit = fit)
      C1 <- A/g * (R - ER) + ER1
      C2 <- C1 + (1-A)/g * (R - ER) + ER0
      risk  <- function(beta) {
        #f <- as.vector(x_basis %*% beta[-1]) + beta[1]
        f <- LRR
        loss = C1*-1*f + C2 * log(1 + exp(f))

        return(mean(loss))
      }
      risk <- risk(beta)
      if(cache){
        private$.risk_history <- c(private$.risk_history, risk)
      }
      return(risk)

    },
    gradient_descent_update = function(tmle_task, fold_number = "full", search_eps = F) {
      eps <- private$.eps

      if(!is.null(eps)) {

        if(abs(eps) < 1e-7){
          return("converged")
        }

      }


       basis_list <- private$.basis_list
      beta <- private$.beta

      lst <- self$design_matrix(tmle_task, basis_list = basis_list)
      x_basis <- lst$x_basis
      R <- tmle_task$get_tmle_node("R", format = T)$R
      A <- tmle_task$get_tmle_node("A", format = T)$A
      g <- self$likelihood$get_likelihood(tmle_task, "A", fold_number)
      lik <- self$likelihood
      cf_task1 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=1))
      cf_task0 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=0))
      ER1 <- lik$get_likelihood(cf_task1, "R", fold_number)
      ER0 <- lik$get_likelihood(cf_task0, "R", fold_number)
      ER <- lik$get_likelihood(tmle_task, "R", fold_number)
      LRR <- self$predict(tmle_task, type = "LRR")

      C1 <- A/g * (R - ER) + ER1
      C2 <- C1 + (1-A)/g * (R - ER) + ER0
      Z <- -C1 + C2 * exp(LRR) / (1 + exp(LRR))
      risk  <- function(beta) {
        f <- as.vector(x_basis %*% beta[-1]) + beta[1]
        loss = C1*-1*f + C2 * log(1 + exp(f))
        return(mean(loss))
      }
      intercept <- beta[1]
      beta_free <- beta[-1]
      D <-  c(mean(Z) * intercept ,  as.vector(Matrix::crossprod(x_basis , Z/ length(Z)) * beta_free))

      D_star <- D - (sum(D * abs(beta))/ sum(beta^2))*abs(beta)
      norm <-   (norm(D_star, type = "2") / sqrt(length(D_star)))
      if(search_eps){
        print(norm)
      }
      if(norm < 1e-5) {
        return("converged")
      }
      D_star <- D_star /norm
      cur_risk <- risk(beta)
      eps <- private$.eps

      eps_max <- 1/max(abs(D_star))

      if(!is.null(eps)) {
        new_beta <- (1+eps*D_star) * beta
        new_risk <- risk(new_beta)
        better <- new_risk <= cur_risk
        if(!better){

        }
      } else {
        better <- T

      }
      if(is.null(eps) | search_eps | !better) {
        if(is.null(eps)) {
          eps <- eps_max/2
        }

        search_set <- -1 * c(exp(seq(log(2*abs(eps)), log(.000000001), length.out = 20)))
        #search_set <- c( -1 * search_set, search_set)
        search_set <- c(0.00005, search_set)
        risks <- c()
        for (ep in search_set) {
          new_beta <- (1+ep*D_star) * beta
          risks <- c(risks, risk(new_beta))
        }
        eps <- search_set[which.min(risks)]
        private$.eps <- eps
        new_risk <- min(risks)
        if(new_risk > cur_risk) {
          warning("fit not improved in gradient update")
          eps <- 0
        }
        new_beta <- (1+eps*D_star) * beta
      }
      private$.beta <- new_beta
      return("Not converged")

    },
    # Targets the R node in the likelihood such that it solves
    # specific score equations of the efficient risk
    update_likelihood = function(tmle_task, fold_number = "full", iter = 1) {
      updater <- self$updater
      for (i in 1:iter) {
        updater$set_estimates(tmle_task, fold_number)
        updater$update_step(self$likelihood, tmle_task, fold_number)
      }
    },
    # Initializes RR hal fit at optimal values
    initialize_best = function(tmle_task, fold_number = "full") {
      type  <- self$params$type
      self$type <- "none"
      private$.useIPW <- T
      self$updater$set_estimates(tmle_task, fold_number)
      self$updater$update_step(self$likelihood, tmle_task, fold_number)
      fit1 <- list(beta = private$.beta, basis_list = private$.basis_list)
      risk1 <- self$efficient_risk(tmle_task, fold_number)
      private$.useIPW <- F
      self$type <- "plugin"
      self$updater$set_estimates(tmle_task, fold_number)
      self$type <- type
      fit2 <- list(beta = private$.beta, basis_list = private$.basis_list)
      risk2 <- self$efficient_risk(tmle_task, fold_number)
      if(risk1 < risk2) {
        private$.basis_list <- fit1$basis_list
        private$.beta <- fit1$beta
      }

    },
    # Uses HAL to obtain an estimate of the LRR.
    # Either to obtain an initial fit or obtain an improved fit.
    train = function(tmle_task, fold_number = "full", type = NULL) {
      if(is.null(type)) {
        type <- self$params$type
      }

      if(type == "IPW" | private$.useIPW) {
        private$.useIPW <- F
        self$train_IPW(tmle_task, fold_number)
      } else if (type == "grad" & !is.null(private$.beta)) {
        first = T
        for(i in 1:100) {
          v = self$gradient_descent_update(tmle_task, fold_number, search_eps = first)
          if(v == "converged"){
            private$.eps <- NULL
            v = self$gradient_descent_update(tmle_task, fold_number, search_eps = T)
            if(v == "converged") {
              break
            }

          }
          first <- F
        }
      }
      else {
        self$train_plugin(tmle_task, fold_number)
      }
      return(self$efficient_risk(tmle_task, fold_number, cache = T))
    },
    # Estimates LRR via the IPW risk
    train_IPW = function(tmle_task, fold_number = "full") {

      task <- tmle_task$get_regression_task("RR")
      R <- tmle_task$get_tmle_node("R", format = T)$R
      A <- tmle_task$get_tmle_node("A", format = T)$A
      g <- self$likelihood$get_likelihood(tmle_task, "A", fold_number)
      keep <- R==1
      weights <- (R/g)
      Y <- A
      lst <- self$design_matrix(tmle_task)
      x_basis <- lst$x_basis
      basis_list <- lst$basis_list
      fitOnce <- function(i) {

        fit <- glmnet::cv.glmnet(x_basis[keep,], Y[keep], family = "binomial", weights = weights[keep], standardize = F, nlambda = 125)
        return(fit$lambda.min)
      }
      num_fits <- self$params$num_fits
      lambdas <- future.apply::future_lapply(1:num_fits, fitOnce)

      lambda <- mean(unlist(lambdas))

      fit <- glmnet::glmnet(x_basis[keep,], Y[keep],  family = "binomial", weights = weights[keep], standardize = F, lambda = lambda)

      private$.beta <-  as.vector(coef(fit))
      private$.basis_list <- basis_list
      beta <- private$.beta
      keep <- which(beta[-1]!=0)
      private$.basis_list <- basis_list[keep]
      private$.beta <- beta[c(1, keep+1)]
      private$.IPWfit <- list(basis_list = private$.basis_list, beta = private$.beta)

      return(invisible(list(basis_list = private$.basis_list, beta = private$.beta)))
    },
    # Estimates the LRR using the plug-in risk function with nuisance parameter E[R|A,W]
    # Primarily useful to obtain a hal fit for E[R|A=1,W] / E[R|A = 0,W] when E[R|A,W] might have been fit
    # with a different machine learning algorithm
    train_plugin = function(tmle_task, fold_number = "full"){

      task <- tmle_task$get_regression_task("RR")
      R <- tmle_task$get_tmle_node("R", format = T)$R
      lik <- self$likelihood
      cf_task1 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=1))
      cf_task0 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=0))
      ER1 <- lik$get_likelihood(cf_task1, "R", fold_number)
      ER0 <- lik$get_likelihood(cf_task0, "R", fold_number)
      weights <- ER1 + ER0
      Y <- ER1/weights
      if(all(Y == Y[1])) {
        print(as.data.table(Y))
        stop("Constants outcome. Glmnet will break down.")
      }
      lst <- self$design_matrix(tmle_task)

      x_basis <- lst$x_basis

      basis_list <- lst$basis_list
      fitOnce <- function(i, return_fit = F) {
        suppressWarnings(fit <- glmnet::cv.glmnet(x_basis, Y, family = binomial(), weights = weights, standardize = F, nlambda = 125))
        lambda.min <- fit$lambda.min
        if(isTRUE(all.equal(lambda.min,  min(fit$lambda)))) {
          # In case near border (early stopping)
          lambda = exp(seq(log(fit$lambda.min), log(.00005 * fit$lambda.min), length.out = 100))
          fit <- glmnet::cv.glmnet(x_basis, Y,  family = binomial(), weights = weights, standardize = F, lambda = lambda )
        }
        if(return_fit) return(fit)
        return(fit$lambda.min)
      }
      num_fits <- self$params$num_fits
      lambdas <- future.apply::future_lapply(1:num_fits, fitOnce)
      lambda <- mean(unlist(lambdas))
      suppressWarnings(fit <- glmnet::glmnet(x_basis, Y,  family = binomial(), weights = weights, standardize = F, lambda = lambda ))

      #suppressWarnings(fit <- glmnet::glmnet(x_basis, Y, family = binomial(), weights = weights, standardize = F, lambda = lambda.min))

      private$.beta <-  as.vector(coef(fit))
      private$.basis_list <- basis_list

      beta <- private$.beta
      keep <- which(beta[-1]!=0)
      private$.basis_list <- basis_list[keep]
      private$.beta <- beta[c(1, keep+1)]
      private$.pluginFit <- list(basis_list = private$.basis_list, beta = private$.beta)
      return(invisible(list(basis_list = private$.basis_list, beta = private$.beta)))

    },
    # This outputs a prediction of the RR by using E[R|A,W] as a plug-in estimator
    predict_derived = function(tmle_task, fold_number = "full") {
      lik <- self$likelihood
      cf_task1 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=1))
      cf_task0 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=0))
      ER1 <- lik$get_likelihood(cf_task1, "R", fold_number)
      ER0 <- lik$get_likelihood(cf_task0, "R", fold_number)
      return(ER1 / ER0)
    },
    # This outputs a prediction of the RR based on however HAL was fit.
    predict = function(tmle_task, type = c("LRR", "RR"), fit_object = NULL) {

      type <- match.arg(type)
      if(!is.null(fit_object)) {
        basis_list <- fit_object$basis_list
        beta <- fit_object$beta
      } else {
        basis_list <- private$.basis_list
        beta <- private$.beta
      }

      lst <- self$design_matrix(tmle_task, basis_list = basis_list)
      x_basis <- lst$x_basis
      predictions <-  as.vector(beta[1] + x_basis %*% beta[-1])
      if(type == "RR") {
        predictions <- exp(predictions)
      }
      return(predictions)
    },
    get_fit = function(type = c("cur", "IPW", "plugin")) {
      type <- match.arg(type)
      if(type == "cur") {
        return(self$fit_object)
      } else if (type == "IPW") {
        return(private$.IPWfit)
      } else {
        return(private$.pluginFit)
      }
    }


  ),
  private = list(
    .params = NULL,
    .beta = NULL,
    .basis_list = NULL,
    .useIPW = NULL,
    .eps = NULL,
    .risk_history = c(),
    .tmle_param = NULL,
    .IPWfit = NULL,
    .pluginFit = NULL
  ),
  active = list(
    params = function(){
      private$.params
    },
    type = function(type = NULL) {
      if(!is.null(type)){
        private$.params$type <- type
      }
      return(private$.params$type)
    },
    risk_history = function(){
      private$.risk_history
    },
    tmle_param = function(){
      private$.tmle_param
    },
    updater = function(){
      self$likelihood$updater
    },
    likelihood = function(){
      private$.params$likelihood
    },
    fit_object = function(){
      return(list(fit = private$.beta, basis_list = private$.basis_list))
    }
  )
)
