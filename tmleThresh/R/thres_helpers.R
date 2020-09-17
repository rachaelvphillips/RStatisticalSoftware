
#' @export
threshold_likelihood <- function(tmle_task, learner_list, cutoffs, bins = 10, cv = T) {
  # covariates
  W_factor <- define_lf(LF_emp, "W")

  # treatment (bound likelihood away from 0 (and 1 if binary))
  A_type <- tmle_task$npsem[["A"]]$variable_type
  if (A_type$type == "continous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }


  A_factor <- define_lf(LF_derived, "A", learner = Lrnr_CDF$new(learner_list[["A"]], bins, cutoffs, cv = cv),likelihood, generator_R,  type = "mean", bound = A_bound)

  # outcome
  Y_factor <- LF_derived$new("Y", Lrnr_thresh$new(learner_list[["Y"]], tmle_task$npsem[["A"]]$variables, cutoffs =cutoffs, cv = cv ),likelihood, generator_R, type = "mean")


  # construct and train likelihood
  factor_list <- list(W_factor, A_factor, Y_factor)

  # add outcome censoring factor if necessary
  if (!is.null(tmle_task$npsem[["Y"]]$censoring_node)) {
    if (is.null(learner_list[["delta_Y"]])) {
      stop("Y is subject to censoring, but no learner was specified for censoring mechanism delta_Y")
    }
    #TODO
    delta_Y_factor <- define_lf(LF_fit, "delta_Y", learner = learner_list[["delta_Y"]], type = "mean", bound = c(0.025, 1))
    factor_list <- c(factor_list, delta_Y_factor)
  }

  if (!is.null(tmle_task$npsem[["A"]]$censoring_node)) {
    stop("A is subject to censoring, this isn't supported yet")
  }

  if (!is.null(tmle_task$npsem[["W"]]$censoring_node)) {
    stop("W is subject to censoring, this isn't supported yet")
  }

  likelihood_def <- Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  return(likelihood)
}



#Construct lrnr that chains the S column
#' @export
Lrnr_thresh <- R6::R6Class(
  classname = "Lrnr_thresh", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(lrnr = make_learner(Lrnr_glm), strata_variable, cutoffs,cv = T,
                          ...) {
      params <- args_to_list()

      lrnr <- make_learner(Pipeline, make_learner(Lrnr_chainer, cutoffs, strata_variable), lrnr_bin, make_learner(Lrnr_wrapper, length(cutoffs)))
      if(cv) {
        lrnr <- Lrnr_sl$new(lrnr, make_learner(Lrnr_solnp, metalearner_linear_multivariate, loss_squared_error_multivariate))
      }
      params$lrnr <- lrnr

      super$initialize(params = params, ...)
    },
    predict_fold = function(task, fold_number = "full") {
      return(private$.predict_fold(task, fold_number))
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "cv"),

    .train_sublearners = function(task) {

      lrnr <- self$params$lrnr

      #task <- private$.process_task(task, training = T )
      if(!("cv" %in% lrnr$properties)) {
        task <- task$revere_fold_task("full")
      }

      return(delayed_learner_train(lrnr, task))
    },
    .process_task = function(task, training = T) {

      args <- self$params
      cutoffs <- args$cutoffs

      strata_variable <- args$strata_variable

      if(inherits(task, "delayed")) {
        task <- task$compute()
      }
      if(inherits(task, "sl3_revere_Task")) {
        new_generator <- function(task, fold_number) {
          task <- task$revere_fold_task(fold_number)
          data <- task$data
          cutoffs <- args$cutoffs
          data_list <- list()

          for(cutoff in cutoffs) {
            Xcopy <- copy(data)
            Xcopy$bin <- cutoff
            Xcopy$Ind <- as.numeric(Xcopy[[strata_variable]] >= cutoff)
            Xcopy[[strata_variable]] <- NULL
            data_list[[as.character(cutoff)]] <- Xcopy
          }
          data <- rbindlist(data_list)

          nodes <- task$nodes
          nodes$id <- NULL
          nodes$t <- NULL
          nodes$covariates <- union(setdiff(task$nodes$covariates, strata_variable), c("Ind", "bin"))
          task <- sl3_Task$new(data, nodes = nodes)

          return(task)
        }
        task <- sl3_revere_Task$new(new_generator, task)


      } else {
        args <- self$params
        strata_variable <- args$strata_variable

        data <- task$data
        cutoffs <- args$cutoffs
        data_list <- list()

        for(cutoff in cutoffs) {
          Xcopy <- copy(data)
          Xcopy$bin <- cutoff
          Xcopy$Ind <- as.numeric(Xcopy[[strata_variable]] >= cutoff)
          Xcopy[[strata_variable]] <- NULL
          data_list[[as.character(cutoff)]] <- Xcopy
        }
        data <- rbindlist(data_list)

        nodes <- task$nodes
        nodes$id <- NULL
        nodes$t <- NULL
        nodes$covariates <- union(setdiff(task$nodes$covariates, strata_variable), c("Ind", "bin"))
        task <- sl3_Task$new(data, nodes = nodes)

      }
      return(task)
    },
    .train = function(task, fit) {

      return(list(lrnr = fit))

    },
    .predict_fold = function(task, fold_number) {
      args <- self$params
      cutoffs <- args$cutoffs

      strata_variable <- args$strata_variable
      #task <- private$.process_task(task, F)

      if(!("cv" %in% self$fit_object$lrnr$properties)) {
        task <- task$revere_fold_task("full")
        predictions <- self$fit_object$lrnr$predict(task)
      } else {
        predictions <- self$fit_object$lrnr$predict_fold(task, fold_number)
      }

      #predictions <- matrix(predictions, ncol = length(cutoffs))
      predictions <- sl3::unpack_predictions(predictions)

      predictions <- as.vector(predictions)
      return(predictions)
    },
    .predict = function(task = NULL) {
      args <- self$params
      cutoffs <- args$cutoffs
      print("here")
      strata_variable <- args$strata_variable
      #task <- private$.process_task(task, F)
      if(!("cv" %in% lrnr$properties)) {
        task <- task$revere_fold_task("full")
        predictions <- self$fit_object$lrnr$predict(task)
      } else {
        predictions <- self$fit_object$lrnr$predict_fold(task, "full")
      }
      print("here")
      print(self$fit_object$lrnr)

      predictions <- sl3::unpack_predictions(predictions)
      predictions <- as.vector(predictions)
      return(predictions)
    }
  )
)
#' @export
Lrnr_CDF <- R6::R6Class(
  classname = "Lrnr_CDF", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(lrnr, num_bins, threshs, type = "left-continuous", cv=T,
                          ...) {

      params <- args_to_list()

      super$initialize(params = params, ...)
    },
    predict_fold = function(task, fold_number = "full") {
      private$.predict_fold(task, fold_number)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial"),

    .train_sublearners = function(task) {
      args <- self$params
      lrnr <- args$lrnr
      num_bins <- args$num_bins

      out <- private$.process_task(task, training = T)
      folds <- out$folds
      task <- out$task
      if(args$cv) {
        lrnr <- Lrnr_sl$new(lrnr, sl3:::default_metalearner(list(type = "binomial")), folds = folds)
       # lrnr <- Lrnr_cv$new(lrnr, folds = folds)

      }
      if(!("cv" %in% lrnr$properties)) {
        task <- task$revere_fold_task("full")
      }

      cutoffs <- out$cutoffs


      lrnr <- delayed_learner_train(lrnr, task)

      private$cutoffs <- cutoffs
      return(lrnr)
    },
    .process_task = function(task, cutoffs = NULL, training = T) {

      args <- self$params
      num_bins <- args$num_bins

      if(inherits(task, "delayed")) {
        task <- task$compute()
      }
      if(inherits(task, "sl3_revere_Task")) {
        task1 <- task$revere_fold_task("validation")
        data <- task1$data
        Y <- task1$Y
        if(is.null(cutoffs)) {
          cutoffs <- as.vector(quantile(Y, seq(0, 1, length.out = num_bins)))
        }

        folds <- task1$folds
        orig_ids <- seq_len(task1$nrow)
        cluster_ids <- rep(orig_ids, length(cutoffs[-1]))
        new_folds <- id_folds_to_folds(folds, cluster_ids)




        new_generator <- function(task, fold_number) {
          old_task <- task

          task <- task$revere_fold_task(fold_number)
          folds <- task$folds
          orig_ids <- seq_len(task$nrow)

          data <- task$data
          Y <- task$Y

          Y <- (cutoffs[findInterval(Y, cutoffs, left.open = self$params$type != "left-continuous", all.inside
                                     = T)])
          #nested revere task


          data_list <- list()

          for(cutoff in (cutoffs[-1])) {
            Xcopy <- copy(data)
            Xcopy$bin <- cutoff
            Xcopy$in_bin <- as.numeric(Y==cutoff)
            Xcopy$Y <- Y
            data_list[[as.character(cutoff)]] <- Xcopy
          }
          data <- rbindlist(data_list)
          cluster_ids <- rep(orig_ids, length(cutoffs[-1]))
          new_folds <- id_folds_to_folds(folds, cluster_ids)
          if(training) {
            keep <- data$bin <= data$Y

            weights <- as.numeric(data$bin <= data$Y) * task$weights
            data$weights <- weights
            #data <- data[keep]

          } else {


            weights <-  rep(1, nrow(data)) * task$weights
            data$weights <- weights
          }

          nodes <- task$nodes
          nodes$covariates <- union(task$nodes$covariates, c("bin"))
          nodes$outcome <- "in_bin"
          nodes$weights <- "weights"
          pooled_task <- sl3_Task$new(data, nodes = nodes, weights = "weights")

          return(pooled_task)
        }
        pooled_task <- sl3_revere_Task$new(new_generator, task)
      } else {
        data <- task$data
        Y <- task$Y
        cutoffs <- as.vector(quantile(Y, seq(0, 1, length.out = num_bins)))

        Y <- (cutoffs[findInterval(Y, cutoffs, left.open = self$params$type != "left-continuous", all.inside
                                   = T)])
        #nested revere task


        data_list <- list()

        for(cutoff in (cutoffs[-1])) {
          Xcopy <- copy(data)
          Xcopy$bin <- cutoff
          Xcopy$in_bin <- as.numeric(Y==cutoff)
          Xcopy$Y <- Y
          data_list[[as.character(cutoff)]] <- Xcopy
        }
        data <- rbindlist(data_list)
        data <- data[data$bin <= data$Y]
        nodes <- task$nodes
        nodes$covariates <- union(task$nodes$covariates, c("bin"))
        nodes$outcome <- "in_bin"
        pooled_task <- sl3_Task$new(data, nodes = nodes)
        new_folds <- NULL
        cutoffs = list("full" = cutoffs)
      }
      return(list(task = pooled_task, cutoffs = cutoffs, folds = new_folds))
    },
    .train = function(task, fit) {

        return(list(lrnr = fit))
    },
    .predict_fold = function(task, fold_number) {

      args <- self$params
      cutoffs <- private$cutoffs
      #cutoffs <- cutoffs[[as.character(fold_number)]]
      orig_cutoffs <- cutoffs
      #cutoffs <- c(-Inf,cutoffs)
      threshs <- args$threshs
      lrnr <- self$fit_object$lrnr


      out <- private$.process_task(task, cutoffs, training = F)
      pooled_task <- out$task

      # This additional revere thing shouldnt be needed but iti s
      predictions <- matrix(lrnr$predict_fold(pooled_task$revere_fold_task(fold_number), fold_number), nrow = task$revere_fold_task("validation")$nrow)
      predictions <- cbind(rep(0, task$revere_fold_task("validation")$nrow), 1 - t(apply(1-predictions, 1, cumprod)))

      #Interpolate to get values at desired cutoffs
      predictions <- rbindlist(lapply(1:nrow(predictions), function(i){

        approx(orig_cutoffs, as.vector(predictions[i,]), xout = as.vector(threshs), rule = 2, yright = 1, yleft = 0 )
      }))

      predictions <- as.vector(predictions[,2][[1]])

      predictions <- (matrix(predictions, ncol = length(threshs), byrow = T))
      return(as.vector(predictions))
      #predictions <- sl3::pack_predictions(predictions)
      predictions <- sl3::pack_predictions(predictions)
    },
    .predict = function(task = NULL) {

      args <- self$params
      cutoffs <- private$cutoffs

      orig_cutoffs <- cutoffs
      #cutoffs <- c(-Inf,cutoffs)
      threshs <- args$threshs
      lrnr <- self$fit_object$lrnr

      pooled_task <- private$.process_task(task, cutoffs)$task

      predictions <- matrix(lrnr$predict_fold(pooled_task, "validation"), nrow = task$nrow)

      predictions <- cbind(rep(0, task$nrow), 1 - t(apply(1-predictions, 1, cumprod)))

      #Interpolate to get values at desired cutoffs
      predictions <- rbindlist(lapply(1:nrow(predictions), function(i){

        approx(orig_cutoffs, as.vector(predictions[i,]), xout = as.vector(threshs), rule = 2, yright = 1, yleft = 0 )
      }))

      predictions <- as.vector(predictions[,2][[1]])

      predictions <- (matrix(predictions, ncol = length(threshs), byrow = T))
      return(as.vector(predictions))
      #predictions <- sl3::pack_predictions(predictions)
      predictions <- sl3::pack_predictions(predictions)

      return(predictions)
    },
    cutoffs = NULL
  )
)



