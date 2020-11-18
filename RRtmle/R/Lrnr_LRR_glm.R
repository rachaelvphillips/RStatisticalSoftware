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
      print("filter")
      keep <- task$weights!=0

      task <- task[task$weights!=0]
      print(task$nrow)
      method <- self$params$method
      X <- cbind(1,task$X)

      print(data.table(task$Y))
      Y <- (task$Y)
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
      Y <- as.vector(Y)
      X <- as.matrix(X)
      if(self$params$lasso) {
        print("lasso")

        family <-  "binomial"
        family <-   binomial()
        fit_object <- cv.glmnet(as.matrix(X), Y, family =  family, weights = weights, intercept = F)
        print(fit_object)
        print("Done")
        coefs <- coef(fit_object, s = "lambda.min")[-1]
        print(length(coefs))
        print("Done")
       } else {

        fit_object <- speedglm::speedglm.wfit(Y, as.matrix(X), family = binomial(), weights = weights, intercept = F)
        coefs <- fit_object$coef
       }


      fit_object <- list(coef = as.vector(coefs))
      return(fit_object)
    },
    .predict = function(task = NULL) {
      print("predict")
      fit_obj <- self$fit_object
      coef <- as.vector(fit_obj$coef)
      X <- as.matrix(cbind(1,task$X))

      predictions <- X %*% coef

      return(predictions)
    },
    .required_packages = c("speedglm")
  )
)
