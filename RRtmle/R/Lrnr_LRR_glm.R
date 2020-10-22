Lrnr_LRR_glm <- R6Class(
  classname = "Lrnr_LRR_glm", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 3,
                          smoothness_degree = 0,
                          iter = 200,
                          method = NULL,
                          lasso = F,
                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("binomial"),

    .train = function(task) {
      print("here")
      method <- self$params$method
      X <- cbind(1,task$X)

      Y <- task$Y
      if(!is.null(method)) {
        weights <- task$get_data(,c("weightsIPW", "weightsplugin"))
        if(method == "IPW") {
          weights <- weights[, weightsIPW]
          Y <- Y[, YIPW]
        }
        else if(method == "plugin") {
          weights <- weights[, weightsplugin]
          Y <- Y[, Yplugin]
        }
        weights <- weights * task$weights

      } else {
        weights <- task$weights
      }
      if(self$params$lasso) {
        fit_object <- glmnet::cv.glmnet(as.matrix(X), Y, family = binomial(), weights = weights, intercept = F)
        coefs <- coef(fit_object, s = "lambda.min")
       } else {
         print("k")
        fit_object <- speedglm::speedglm.wfit(Y, as.matrix(X), family = binomial(), weights = weights, intercept = F)
        coefs <- fit_object$coef
       }


      fit_object <- list(coef = as.vector(coefs))
      return(fit_object)
    },
    .predict = function(task = NULL) {
      fit_obj <- self$fit_object
      coef <- as.vector(fit_obj$coef)
      X <- as.matrix(cbind(1,task$X))

      predictions <- X %*% coef
      return(predictions)
    },
    .required_packages = c("glmnet", "speedglm")
  )
)
