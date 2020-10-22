
Lrnr_LRR_IPW_chainer <- R6Class(
  classname = "Lrnr_LRR_IPW_chainer", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(g1_var = "g1", treatment_var = "A", outcome_var = "R", sieve_learner = Lrnr_fourier$new(fourier_basis(2,1)), ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    },

    print = function() {
      print(self$name)
    }
  ),

  private = list(
    .properties = c("continuous", "binomial", "categorical", "weights", "offset"),

    .train = function(task) {
      trt <- self$params$treatment_var
      outcome <- self$params$outcome_var
      g1_var <- self$params$g1_var
      lrnr <- self$params$sieve_learner

      task_g <- task$next_in_chain(offset = g1_var, covariates = setdiff(task$nodes$covariates, trt), outcome = trt)
      lrnr <- lrnr$train(task_g)
      fit_object <- list()
      fit_object$lrnr <- lrnr
      return(fit_object)
    },

    .predict = function(task = NULL) {
      stop("Nothing to predict")
    },
    .chain = function(task = NULL) {
      trt <- self$params$treatment_var
      outcome <- self$params$outcome_var
      g1_var <- self$params$g1_var
      lrnr <- self$fit_object$lrnr
      task_g <- task$next_in_chain(offset = g1_var, covariates = setdiff(task$nodes$covariates, trt), outcome = trt)
      gnew <- lrnr$predict(task_g)
      covariates_IPW <- setdiff(task$covariates, trt)
      weights <- task$get_data(,outcome)[[1]]/gnew * task$weights

      new_data <- data.table(gnew, weights)
      colnames(new_data) <- c(g1_var, "weights")
      column_names <- task$add_columns(new_data)

      task <- task$next_in_chain(column_names = column_names, weights = "weights", outcome = trt, covariates = covariates_IPW)
      return(task)
    }
  )
)




Lrnr_LRR_plugin_chainer <- R6Class(
  classname = "Lrnr_LRR_plugin_chainer", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(Q_var = "Q", Q1_var = "ER1", Q0_var = "ER0", treatment_var = "A", outcome_var = "R", sieve_learner = Lrnr_fourier$new(fourier_basis(2,1)), ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    },

    print = function() {
      print(self$name)
    }
  ),

  private = list(
    .properties = c("continuous", "binomial", "categorical", "weights", "offset"),

    .train = function(task) {
      trt <- self$params$treatment_var
      outcome <- self$params$outcome_var
      Q_var <- self$params$Q_var
      lrnr <- self$params$sieve_learner

      task_q <- task$next_in_chain(offset = Q_var, covariates = union(task$nodes$covariates, trt), outcome = outcome)
      lrnr <- lrnr$train(task_q)
      fit_object <- list()
      fit_object$lrnr <- lrnr
      return(fit_object)
    },

    .predict = function(task = NULL) {
      stop("Nothing to predict")
    },
    .chain = function(task = NULL) {
      trt <- self$params$treatment_var
      outcome <- self$params$outcome_var
      Q_var <- self$params$Q_var
      Q0_var <- self$params$Q0_var
      Q1_var <- self$params$Q1_var

      lrnr <- self$fit_object$lrnr

      task_q <- task$next_in_chain(offset = Q_var, covariates = union(task$nodes$covariates, trt), outcome = outcome)
      cf_data1 <- data.table(rep(1, task_q$nrow))
      names(cf_data1) <- trt
      cf_data0 <- data.table(rep(0, task_q$nrow))
      names(cf_data0) <- trt
      column_map1 <- task_q$add_columns(cf_data1)
      column_map0 <- task_q$add_columns(cf_data0)
      task_q1 <- task_q$next_in_chain(column_names = column_map1, offset = Q1_var)
      task_q0 <- task_q$next_in_chain(column_names = column_map0, offset = Q0_var)

      Qnew <- lrnr$predict(task_q)
      Qnew0 <- lrnr$predict(task_q0)
      Qnew1 <- lrnr$predict(task_q1)
      weights <- (Qnew0 + Qnew1)*task$weights
      covariates_plugin <- setdiff(task$covariates, trt)
      outcome_val <- Qnew1 / (Qnew0 + Qnew1)
      new_data <- data.table(outcome_val, weights)
      colnames(new_data) <- c("outcome", "weights")
      column_names <- task$add_columns(new_data)

      task <- task$next_in_chain(column_names = column_names, weights = "weights", outcome = "outcome", covariates = covariates_plugin)
      return(task)
    }
  )
)











