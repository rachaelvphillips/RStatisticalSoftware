
Lrnr_glmnet_relaxed <- R6Class(
  classname = "Lrnr_glmnet_relaxed",
  inherit = Lrnr_base, portable = TRUE, class = TRUE,
  public = list(
    initialize = function(lambda = NULL, type.measure = "deviance", nfolds = 10,
                          alpha = 1, nlambda = 100, use_min = TRUE, relax_lambda = "best", ...) {
      super$initialize(params = args_to_list(), ...)
    }
  ),

  private = list(
    .properties = c("continuous", "binomial", "categorical", "weights"),

    .train = function(task) {
      args <- self$params
      outcome_type <- self$get_outcome_type(task)

      if (is.null(args$family)) {
        args$family <- outcome_type$glm_family(return_object = F)
      }

      if (args$family %in% "quasibinomial") {
        args$family <- "gaussian"
        warning(paste(
          "Lrnr_glmnet doesn't understand outcome_type =",
          "'quasibinomial'; fitting glmnet with family='gaussian'",
          "instead."
        ))
      }

      # specify data
      args$x <- as.matrix(task$X)
      args$y <- outcome_type$format(task$Y)

      if (task$has_node("weights")) {
        args$weights <- task$weights
      }

      if (task$has_node("offset")) {
        args$offset <- task$offset
      }


      fit_object <- sl3:::call_with_args(
        glmnet::cv.glmnet, args,
        names(formals(glmnet::glmnet))
      )
      fit_object$glmnet.fit$call <- NULL
      relax_lambda <- self$params$relax_lambda
      if(relax_lambda == "best") {
        relax_lambda = "lambda.min"
      } else{
        relax_lambda = min(fit_object$lambda)
      }

      coefs <- coef(fit_object, s = relax_lambda)
      print(length(coefs))
      keep <- coefs[-1] != 0
      print(sum(keep))

      x <- cbind(rep(1, nrow(args$x)), args$x[,keep, drop = F])
      args$family <- outcome_type$glm_family(return_object = T)

      glm_fit <- glm.fit(as.matrix(x), args$y, family = args$family, weights = args$weights, intercept = F)
      coefs <- glm_fit$coefficients

      fit_object <- list(keep = keep, coefs = coefs)
      return(fit_object)
    },

    .predict = function(task) {
      outcome_type <- private$.training_outcome_type
      if (self$params$use_min) {
        lambda <- "lambda.min"
      } else {
        lambda <- "lambda.1se"
      }
      keep <- self$fit_object$keep
      coefs <- as.vector(self$fit_object$coefs)
      coefs[is.na(coefs)] <- 0
      if(length(coefs) > 1){
        X <- as.matrix(task$X)[,keep, drop = F]

        link <- X %*% coefs[-1] + coefs[1]
      } else {
        X <- as.matrix(task$X)[,keep, drop = F]
        link <- rep(1, nrow(X)) * coefs[1]
      }

      if(outcome_type$type == "binomial") {
        return(plogis(link))
      } else{
        return(link)
      }

      return(predictions)
    },
    .required_packages = c("glmnet", "speedglm")
  )
)
