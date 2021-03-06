
#' @export
Lrnr_hal9001_fixed <- R6::R6Class(
  classname = "Lrnr_hal9001_fixed", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 3,
                          fit_type = "glmnet",
                          n_folds = 10,
                          use_min = TRUE,
                          reduce_basis = NULL,
                          return_lasso = TRUE,
                          return_x_basis = FALSE,
                          basis_list = NULL,
                          cv_select = TRUE,
                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "weights", "ids"),

    .train = function(task) {
      args <- self$params

      outcome_type <- self$get_outcome_type(task)

      if (is.null(args$family)) {
        args$family <- args$family <- outcome_type$glm_family()
      }

      X <- as.matrix(task$X)
      # X <- apply(X, 2, function(v) {
      #   vals <- sort(unique(quantile(v, seq(0,1, length.out = 10))))
      #   v <- vals[findInterval(v, vals)]
      #   return(v)
      # })
      args$X <- X
      args$Y <- outcome_type$format(task$Y)
      args$yolo <- FALSE

      if (task$has_node("weights")) {
        args$weights <- task$weights
      }

      if (task$has_node("offset")) {
        args$offset <- task$offset
      }

      if (task$has_node("id")) {
        args$id <- task$id
      }

      fit_object <- do.call(hal9001::fit_hal, args)
      return(fit_object)
    },
    .predict = function(task = NULL) {
      predictions <- predict(self$fit_object, new_data = as.matrix(task$X))
      if (!is.na(safe_dim(predictions)[2])) {
        p <- ncol(predictions)
        colnames(predictions) <- sprintf("lambda_%0.3e", self$params$lambda)
      }
      return(predictions)
    },
    .required_packages = c("hal9001")
  )
)
