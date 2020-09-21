
#' @export
Lrnr_undersmooth <- R6Class(
  classname = "Lrnr_undersmooth", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(basis_generator = fourier_basis(), standardize = F, stopping_criterion = function(n, residual) {return(1/sqrt(n)/log(n))}, offset_learner = NULL, variable.control = list(), offset_as_covariate = FALSE, ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .train_sublearners = function(task) {
      lrnr <- self$params$offset_learner
      if(is.null(lrnr)) {
        return(NULL)
      }
      lrnr <- delayed_learner_train(lrnr, task)
      return(lrnr)
    },
    .train = function(task, offset_fit = NULL) {
      train <- function(task) {
        variable.control <- self$params$variable.control
        basis_generator <- self$params$basis_generator
        outcome_type <- task$outcome_type$glm_family()
        covariates <- task$nodes$covariates
        if(!is.null(variable.control$weights)) {
          weights <- task$get_data(,variable.control$weights) * task$weights
          covariates <- setdiff(covariates, variable.control$weights)
        }
        else if(task$has_node("weights")) {
          weights <- task$weights
        } else {
          weights <- NULL
        }
        if(!is.null(offset_fit)) {
          offset <- offset_fit$predict(task)
        }
        else if(!is.null(variable.control$offset)){
          offset <- task$get_node(variable.control$offset)
          covariates <- setdiff(covariates, variable.control$offset)
        }
        else if(task$has_node("offset")) {
          offset <- task$get_node("offset")
        } else {
          offset <- NULL
        }

        if(!is.null(variable.control$covariates)) {
          X <- as.matrix(task$get_data(,variable.control$covariates))
        } else {
          X <- as.matrix(task$X[,covariates, with = F])
        }
        if(self$params$offset_as_covariate) {
          X <- cbind(offset, X)
          offset <- NULL
        }
        design_generator <- basis_generator(X)
        x_basis <- design_generator(X)
        Y <- task$Y


        glmnet_fit <- glmnet::cv.glmnet(x_basis, Y, standardize = self$params$standardize, family = outcome_type, weights = weights, offset = offset)
        best_beta <- coef(glmnet_fit, s = "lambda.min")
        basis_to_check <- which(best_beta[-1] !=0)
        best_lambda <- glmnet_fit$lambda.min
        lambdas <- glmnet_fit$lambda

        lambdas <- lambdas[lambdas <= best_lambda]

        lambdas <- sort(lambdas, decreasing = T)

        preds <- predict(glmnet_fit, newx = x_basis, s = lambdas, type = "response", newoffset = offset)
        coefs <- as.matrix(coef(glmnet_fit, s = lambdas))[-1,]
        x_basis <- x_basis[,basis_to_check]
        new_lambda = NULL
        stopping_criterion <- self$params$stopping_criterion

        for(i in seq_along(lambdas)) {

          residual <- as.vector(Y - preds[,i])
          n <- length(residual)
          if(outcome_type == "binomial" |  outcome_type == "gaussian") {
            scores <- as.vector(Matrix::crossprod(x_basis, residual) / n)
            print(data.table(scores = scores))
            passed <- scores <= stopping_criterion(n, residual)
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
        fit_object = list(offset_fit = offset_fit, design_generator = design_generator, fit = glmnet_fit, lambda = new_lambda, has_offset = !is.null(offset))
      }
      stratify <- self$params$stratify_by
      if(!is.null(stratify)) {
        covariates <- setdiff(task$nodes$covariates, stratify)
        A <- as.factor(task$get_data(,stratify)[[1]])
        levels <- levels(A)
        fit_objects <- list()
        for(level in levels) {
          new_task <- task[which(A==level)]$next_in_chain(covariates = covariates)
          fit_objects[[level]] <- train(new_task)
        }
      } else {
        levels <- NULL
        fit_objects <- train(task)
      }

      return(list(fit_objects = fit_objects, levels = levels))
    },
    .predict = function(task) {
      fit_objects <- self$fit_object
      levels <- fit_objects$levels
      fit_objects <- fit_objects$fit_objects
      predict_once <- function(task, fit_object) {
        fit <- fit_object$fit
        offset_fit <- fit_object$offset
        design_generator <- fit_object$design_generator
        has_offset <- fit_object$has_offset
        variable.control <- self$params$variable.control
        if(!is.null(variable.control$covariates)) {
          X <- as.matrix(task$get_data(,variable.control$covariates))
        } else {
          X <- as.matrix(task$X)
        }


        if(has_offset | self$params$offset_as_covariate) {
          offset <- NULL
          if(!is.null(offset_fit)) {
            offset <- offset_fit$predict(task)
          } else if(task$has_node("offset")) {
            offset <- task$get_node("offset")
          }
          if(!is.null(offset) & self$params$offset_as_covariate) {
            X <- cbind(offset, X)
            offset <- NULL
          }
          x_basis  <- design_generator(X)
          preds <- predict(fit, newx = x_basis, s = fit_object$lambda, type = "response", newoffset = offset)
        } else {
          x_basis  <- design_generator(X)
          preds <- predict(fit, newx = x_basis, s = fit_object$lambda, type = "response")
        }
        return(as.vector(preds))
      }
      stratify <- self$params$stratify_by
      if(!is.null(stratify)) {
        covariates <- setdiff(task$nodes$covariates, stratify)
        A <- as.factor(task$get_data(,stratify)[[1]])
        predictions <- rep(0, task$nrow)
        for(level in levels) {
          fit_obj <- fit_objects[[level]]
          index <- which(A==level)
          new_task <- task[index]$next_in_chain(covariates = covariates)
          predictions[index] <- predict_once(new_task, fit_obj)
        }
      } else {
        predictions <- predict_once(task, fit_objects)
      }
      return(predictions)

    }))


