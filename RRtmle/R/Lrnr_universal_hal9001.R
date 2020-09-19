
#' @export
Lrnr_universal_hal9001 <- R6Class(
  classname = "Lrnr_universal_hal9001", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 2, thresh = "sample_size",...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .train = function(task) {
      max_degree <- self$params$max_degree
      X <-as.matrix(task$X)
      Y <- task$Y
      outcome_type <- task$outcome_type$glm_family()
      basis_list <- hal9001fast::enumerate_basis(X, max_degree, bins = rep(350, ncol(X)))
      x_basis <- hal9001fast::make_design_matrix(X, basis_list)
      glmnet_fit <- glmnet::cv.glmnet(x_basis, Y, standardize = F, family = outcome_type)
      best_beta <- coef(glmnet_fit, s = "lambda.min")
      basis_to_check <- which(best_beta[-1] !=0)
      best_lambda <- glmnet_fit$lambda.min
      lambdas <- glmnet_fit$lambda

      lambdas <- lambdas[lambdas <= best_lambda]

      lambdas <- sort(lambdas, decreasing = T)

      preds <- predict(glmnet_fit, newx = x_basis, s = lambdas, type = "response")
      coefs <- as.matrix(coef(glmnet_fit, s = lambdas))[-1,]
      x_basis <- x_basis[,basis_to_check]
      new_lambda = NULL
      for(i in seq_along(lambdas)) {

        residual <- as.vector(Y - preds[,i])
        n <- length(residual)
        if(outcome_type == "binomial"){

          scores <- as.vector(Matrix::crossprod(x_basis, residual) / n)
          print(data.table(scores = scores))
          if(self$params$thresh == "sample_size") {
            passed <- scores <= sd(residual)/n
          } else if(is.numeric(self$params$thresh)) {

            passed <- scores <= sd(residual)/(n^(self$params$thresh))
          } else {
            passed <- scores <= sd(residual)/sqrt(n)/log(n)/2
          }
          if(sum(!passed) == 0) {
            print(sum(coefs[,i]!=0))
            new_lambda = lambdas[i]
            break
          }
        }
      }
      if(is.null(new_lambda)) {
        new_lambda <- min(lambdas)
      }

      fit_object = list(basis_list = basis_list, fit = glmnet_fit, lambda = new_lambda)
      return(fit_object)
    },
    .predict = function(task) {
      fit_object <- self$fit_object
      fit <- fit_object$fit
      X <- as.matrix(task$X)
      basis_list <- fit_object$basis_list
      x_basis <- hal9001fast::make_design_matrix(X, basis_list)
      return(predict(fit, newx = x_basis, s = fit_object$lambda, type = "response"))
    }
  )
)
