
#' @importFrom R6 R6Class
#'
#' @export

Lrnr_hal9001 <- R6Class(
  classname = "Lrnr_hal9001", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 2,
                          fit_type = "glmnet",
                          n_folds = 10,
                          use_min = TRUE,
                          reduce_basis = NULL,
                          return_lasso = TRUE,
                          return_x_basis = FALSE,
                          basis_list = NULL,
                          cv_select = TRUE,
                          lower.limits = -Inf,
                          upper.limits=Inf,
                          ...) {
      params <- args_to_list()

      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial"),

    .train = function(task) {
      args <- self$params

      outcome_type <- self$get_outcome_type(task)

      if (is.null(args$family)) {
        args$family <- args$family <- outcome_type$glm_family()
      }

      args$X <- as.matrix(task$X)
      args$Y <- outcome_type$format(task$Y)
      args$yolo <- FALSE
      args$'...' = list("lower.limits" = args$lower.limits, "upper.limits" = args$upper.limits)

      if (task$has_node("weights")) {
        args$weights <- task$weights
      }

      if (task$has_node("offset")) {
        args$offset <- task$offset
      }
      print(names(args))
      fit_object <- hal9001::fit_hal(args$X, args$Y, max_degree = args$max_degree, yolo = args$yolo, cv_select = args$cv_select, family = args$family, lower.limits = args$lower.limits, upper.limits = args$upper.limits, reduce_basis = args$reduce_basis)
      return(fit_object)
    },
    .predict = function(task = NULL) {
      predictions <- predict(self$fit_object, new_data = as.matrix(task$X))
      return(predictions)
    },
    .required_packages = c("hal9001")
  )
)
