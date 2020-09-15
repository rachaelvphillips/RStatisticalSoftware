Lrnr_LRR_glm <- R6Class(
  classname = "Lrnr_LRR_glm", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 3,
                          smoothness_degree = 0,
                          iter = 200,
                          method = NULL,
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
      fit_object <- speedglm::speedglm.wfit(Y, as.matrix(X), family = binomial(), weights = weights, intercept = F)

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
