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
      method <- self$params$method
      X <- task$X_intercept

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
      } else {
        weights <- task$get_data(,"weights")
      }
      weights <- weights * task$weights
      if(self$params$lasso) {
        fit_object <- glmnet::cv.glmnet(as.matrix(X), Y, family = binomial(), weights = weights, intercept = F)
        coefs <- coef(fit_object, s = "lambda.min")
       } else {
        fit_object <- speedglm::speedglm.wfit(Y, as.matrix(X), family = binomial(), weights = weights, intercept = F)
        coefs <- fit_obj$coef
       }
      fit_object <- list(coef = coefs)
      return(fit_object)
    },
    .predict = function(task = NULL) {
      fit_obj <- self$fit_object
      coef <- fit_obj$coef
      X <- task$X_intercept
      predictions <- X %*% coef
      return(predictions)
    }
  )
)
