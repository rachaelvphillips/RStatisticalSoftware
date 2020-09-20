

#' @export
Lrnr_universal <- R6Class(
  classname = "Lrnr_universal", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(sieve_learner = make_learner(Lrnr_glmnet), max_degree = 2, through_weights = T, through_offset = T, nbasis = 100,...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .train_sublearners = function(task) {
      sieve_learner <- self$params$sieve_learner
      if(inherits(task, "delayed")){
        task <- task$compute()
      }
      max_degree <- self$params$max_degree
      g <- task$get_data(, "g")[[1]]
      Q <- task$get_data(, "Q")[[1]]
      A <- task$get_data(, "A")[[1]]
      X <- task$X
      Y <- task$Y

      if(self$params$through_weights) {
        if(self$params$through_offset) {
          covariates <- setdiff(task$nodes$covariates, "A")
          offset <- "Q"
        } else {
          covariates <- setdiff(c(task$nodes$covariates, "Q"), "A")
          offset <- NULL
        }
        weights <- task$weights / g
        column_names <- task$add_columns(data.table(weights = weights))
        next_task <- task$next_in_chain(covariates = covariates,offset = offset, column_names = column_names, weights = "weights")
      } else {
        if(self$params$through_offset) {
          covariates <- setdiff(c(task$nodes$covariates, "g"), "A")
          offset <- "Q"
        } else {
          covariates <- setdiff(c(task$nodes$covariates, "Q", "g"), "A")
          offset <- NULL
        }
        next_task <- task$next_in_chain(covariates = covariates, offset = offset)
      }
      X <- next_task$X

      covariates <- c()
      basis_list <- list()
      num_basis <- NULL
      var_names <- c()
      for(name in colnames(X)) {

        sd_val <- sd(X[[name]])
        if(sd_val < 1e-9) {
          next
        }
        var_names <- c(var_names, name)
        rangeval <- c(min(X[[name]]) - sd_val/2, max(X[[name]]) +sd_val/2)


        basis <- fda::create.fourier.basis(rangeval = rangeval, nbasis = self$params$nbasis)
        basis_list[[name]] <- basis

        new_data <- as.data.table(fda::eval.basis(X[[name]], basis))
        num_basis <- length(new_data)

        colnames(new_data) <- paste0(name, "_", colnames(new_data))
        covariates <- c(covariates, colnames(new_data))
        #covariates <- setdiff(covariates, grep("const", covariates, value = T))


        column_names <- next_task$add_columns(new_data)
        next_task <- next_task$next_in_chain(covariates = covariates, column_names = column_names)
      }

      private$.basis_list <- basis_list

      X <- next_task$X


      inters <- list()
      if(max_degree > 1){
        for(degree in c(max_degree)) {
          #inter_base <- setdiff(colnames(X), grep("const", colnames(X), value = T))
          #inter <- unlist(apply(combn(inter_base, degree),2, function(x){list(x)}), recursive = F)
          #inters[[as.character(degree)]] <- inter
          vars <- colnames(X)

          num_orig <- length(var_names)

          Xchange <- (matrix(as.vector(as.matrix(X)), ncol = num_orig))
          k <- num_basis

          form <- formula(paste0("~.^", max_degree))
          data <- model.frame(form, data = as.data.frame(Xchange))
          Xnew <- model.matrix(form, data = data)

          Xnew <- Xnew[,-1]
          Xnew <- do.call(cbind, unlist(apply(Xnew, 2, function(v) {list(matrix(v, ncol = k))}), recursive = F))
          Xnew <- data.table(Xnew)







          raw_data <- copy(next_task$internal_data$raw_data)
          suppressWarnings(raw_data <- alloc.col(raw_data, n = ncol(Xnew) + 1))
          shared_data <- Shared_Data$new(raw_data, force_copy = F)
          next_task <- sl3_Task$new(shared_data, nodes = next_task$nodes, column_names = next_task$column_names, folds = next_task$folds, outcome_type = next_task$outcome_type, row_index = next_task$row_index)
          new_columns <- next_task$add_columns(Xnew)
          next_task <- next_task$next_in_chain(covariates = colnames(Xnew), column_names = new_columns)
          #next_task <- next_task$add_interactions(inter)
        }
      }


      private$.inters <- inters
      private$.var_names <- var_names
      task0 <- next_task[A==0]
      task1 <- next_task[A==1]

      lrnr <- sieve_learner



      bundle <- bundle_delayed(list(inters = inters, var_names = var_names, basis_list = basis_list, fits = bundle_delayed(list("R0" = delayed_learner_train(lrnr, task0), "R1" = delayed_learner_train(lrnr, task1)))))

      return(bundle)



    },
    .process_task = function(task, basis_list, inters, var_names) {
      max_degree <- self$params$max_degree
      g <- task$get_data(, "g")[[1]]
      if(self$params$through_weights) {
        if(self$params$through_offset) {
          covariates <- setdiff(task$nodes$covariates, "A")
          offset <- "Q"
        } else {
          covariates <- setdiff(c(task$nodes$covariates, "Q"), "A")
          offset <- NULL
        }
        weights <- task$weights / g
        column_names <- task$add_columns(data.table(weights = weights))
        next_task <- task$next_in_chain(covariates = covariates,offset = offset, column_names = column_names, weights = "weights")
      } else {
        if(self$params$through_offset) {
          covariates <- setdiff(c(task$nodes$covariates, "g"), "A")
          offset <- "Q"
        } else {
          covariates <- setdiff(c(task$nodes$covariates, "Q", "g"), "A")
          offset <- NULL
        }
        next_task <- task$next_in_chain(covariates = covariates, offset = offset)
      }
      X <- next_task$X



      covariates <- c()
      num_basis <- NULL
      for(name in var_names) {
        basis <- basis_list[[name]]
        new_data <- as.data.table(fda::eval.basis(X[[name]], basis))
        num_basis <- ncol(new_data)
        colnames(new_data) <- paste0(name, "_", colnames(new_data))
        covariates <- c(covariates, colnames(new_data))
        #covariates <- setdiff(covariates, grep("const", covariates, value = T))
        column_names <- next_task$add_columns(new_data)
        next_task <- next_task$next_in_chain(covariates = covariates, column_names = column_names)
      }
      X <- next_task$X
      num_orig <-  length(var_names)


      if(max_degree > 1){
        for(degree in c(max_degree)) {
          #inter_base <- setdiff(colnames(X), grep("const", colnames(X), value = T))
          #inter <- unlist(apply(combn(inter_base, degree),2, function(x){list(x)}), recursive = F)
          #inters[[as.character(degree)]] <- inter
          vars <- colnames(X)


          Xchange <- (matrix(as.vector(as.matrix(X)), ncol = num_orig))
          k <- num_basis

          form <- formula(paste0("~.^", max_degree))

          data <- model.frame(form, data = as.data.frame(Xchange))
          Xnew <- model.matrix(form, data = data)
          Xnew <- Xnew[,-1]
          Xnew <- do.call(cbind, unlist(apply(Xnew, 2, function(v) {list(matrix(v, ncol = k))}), recursive = F))
          Xnew <- data.table(Xnew)

          raw_data <- copy(next_task$internal_data$raw_data)
          raw_data <- alloc.col(raw_data, n = ncol(Xnew) + 100)
          shared_data <- Shared_Data$new(raw_data, force_copy = F)
          next_task <- sl3_Task$new(shared_data, nodes = next_task$nodes, column_names = next_task$column_names, folds = next_task$folds, outcome_type = next_task$outcome_type, row_index = next_task$row_index)
          new_columns <- next_task$add_columns(Xnew)

          next_task <- next_task$next_in_chain(covariates = colnames(Xnew), column_names = new_columns)
          #next_task <- next_task$add_interactions(inter)
        }
      }
      return(next_task)

    },
    .train = function(task, sublearners) {
      fits <- sublearners$fits


      fit_object <- sublearners
      return(fit_object)
    },
    .predict = function(task) {
      A <- task$get_data(,"A")
      fit_object <- self$fit_object
      fits <- fit_object$fits

      basis_list <- fit_object$basis_list
      var_names = fit_object$var_names
      inters <- fit_object$inters

      new_task <- private$.process_task(task, basis_list, inters, var_names)


      task0 <- new_task[A==0]

      task1 <- new_task[A==1]
      fit0 <- fits[[1]]
      fit1 <- fits[[2]]
      preds0 <- fit0$predict(task0)
      preds1 <- fit1$predict(task1)





      preds <- c(preds0, preds1)
      indices <- c(which(A==0), which(A==1))
      preds <- preds[order(indices)]
      return(preds)

    },
    .inters = NULL,
    .basis_list = NULL,
    .var_names = NULL,
    .required_packages = c("glmnet", "fda")
  )
)






#' @export
Lrnr_LRR_glm <- R6Class(
  classname = "Lrnr_LRR_glm", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(
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
        coefs <- coef(fit_object, s = "lambda.min")[-1]
      } else {

        fit_object <- speedglm::speedglm.wfit(Y, as.matrix(X), family = binomial(), weights = weights, intercept = F)
        coefs <- fit_object$coef
      }
      fit_object <- list(coef = as.vector(coefs))
      return(fit_object)
    },
    .predict = function(task = NULL) {
      fit_obj <- self$fit_object
      coef <- fit_obj$coef
      X <- task$X_intercept

      predictions <- as.matrix(X) %*% coef
      return(predictions)
    },
    .chain = function(task = NULL) {
      fit_obj <- self$fit_object
      coef <- fit_obj$coef[-1]
      covariates <- task$nodes$covariates
      if(length(coef)!= length(covariates)) {
        stop("oops")
      }
      keep <- which(coef != 0)
      covariates <- covariates[keep]
      return(task$next_in_chain(covariates = covariates))

    }
  )
)
