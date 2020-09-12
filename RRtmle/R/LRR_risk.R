#' @import hal9001fast
#' @importFrom sl3 args_to_list
#' @importFrom uuid UUIDgenerate
#' @export
LRR_risk <- R6Class(
  classname = "LRR_risk",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(likelihood, smoothness = 0, max_degree = 3, IPW_as_initial = T, type = c("efficient", "IPW"), num_fits = 1,...){
      type <- match.arg(type)
      params <- sl3::args_to_list()
      params$type <- type
      private$.params <- params
      private$.useIPW <- IPW_as_initial
    },
    design_matrix = function(tmle_task, basis_list = NULL) {
       X <- tmle_task$get_tmle_node("W", format = T)
       smoothness <- self$params$smoothness
       smoothness <- rep(c(smoothness)[[1]], ncol(X))
       if(is.null(basis_list)) {
         basis_list <- hal9001fast::enumerate_basis(as.matrix(X), max_degree = self$params$max_degree, order_map = c(smoothness), bins = rep(350, ncol(X)), include_zero_order = T, include_lower_order = T)
       }
       x_basis <- hal9001fast::make_design_matrix(as.matrix(X), basis_list)
       return(list(x_basis = x_basis, basis_list = basis_list))

    },
    train = function(tmle_task, fold_number = "full") {
      if(self$params$type == "IPW" | private$.useIPW) {
        private$.useIPW <- F
        return(self$train_IPW(tmle_task, fold_number))
      } else {
        return(self$train_plugin(tmle_task, fold_number))
      }
    },
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

        fit <- glmnet::cv.glmnet(x_basis[keep,], Y[keep], family = "binomial", weights = weights[keep], standardize = F, nlambda = 200)
        return(fit$lambda.min)
      }
      num_fits <- self$params$num_fits
      lambdas <- future.apply::future_lapply(1:num_fits, fitOnce)
      print(lambdas)
      lambda <- mean(unlist(lambdas))

      fit <- glmnet::glmnet(x_basis[keep,], Y[keep],  family = "binomial", weights = weights[keep], standardize = F, lambda = lambda)

      private$.fit <-  fit
      private$.basis_list <- basis_list
    },
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

      lst <- self$design_matrix(tmle_task)

      x_basis <- lst$x_basis

      basis_list <- lst$basis_list
      fitOnce <- function(i) {
        suppressWarnings(fit <- glmnet::cv.glmnet(x_basis, Y, family = binomial(), weights = weights, standardize = F, nlambda = 100))
        lambda.min <- fit$lambda.min
        if(isTRUE(all.equal(lambda.min,  min(fit$lambda)))) {
          # In case near border (early stopping)
          lambda = exp(seq(log(fit$lambda.min), .001*log(fit$lambda.min), length.out = 25))
          fit <- glmnet::cv.glmnet(x_basis, Y,  family = binomial(), weights = weights, standardize = F, lambda = lambda )
        }
        return(fit$lambda.min)
      }
      num_fits <- self$params$num_fits
      lambdas <- future.apply::future_lapply(1:num_fits, fitOnce)
      print(unlist(lambdas))
      lambda <- mean(unlist(lambdas))
      suppressWarnings(fit <- glmnet::glmnet(x_basis, Y,  family = binomial(), weights = weights, standardize = F, lambda = lambda ))


      #suppressWarnings(fit <- glmnet::glmnet(x_basis, Y, family = binomial(), weights = weights, standardize = F, lambda = lambda.min))

      private$.fit <-  fit
      private$.basis_list <- basis_list
    },
    predict = function(tmle_task, type = c("LRR", "RR")) {
      type <- match.arg(type)
      basis_list <- private$.basis_list
      fit <- private$.fit
      lst <- self$design_matrix(tmle_task, basis_list = basis_list)
      x_basis <- lst$x_basis
      predictions <- predict(fit, newx = x_basis, type = "link", s = "lambda.min")
      if(type == "RR") {
        predictions <- exp(predictions)
      }
      return(predictions)
    }

  ),
  private = list(
    .params = NULL,
    .fit = NULL,
    .basis_list = NULL,
    .useIPW = NULL
  ),
  active = list(
    params = function(){
      private$.params
    },
    likelihood = function(){
      private$.params$likelihood
    },
    fit_object = function(){
      return(list(fit = private$.fit, basis_list = private$.basis_list))
    }
  )
)
