Lrnr_LRR_subst <- R6Class(
  classname = "Lrnr_LRR_subst", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(
                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c( "RRtmle"),

    .train = function(task) {
      fit_object <- list()
      return(fit_object)
    },
    .predict = function(task = NULL) {
      fit_obj <- self$fit_object
      ER1 <- task$get_data(,"ER1")
      ER0 <- task$get_data(,"ER0")

      return(log(ER1/ER0))
    }
  )
)
