ltmle3_Node <- R6Class(
  classname = "ltmle3_Node",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Node,
  public = list(
    initialize = function(name, variables, parents = c(), time = NULL, summary_functions = NULL,
                          variable_type = NULL, censoring_node = NULL, scale = FALSE) {
      if(!is.list(summary_functions)){
        summary_functions <- list(summary_functions)
      }
      private$.ltmle_params = list(time = time, summary_functions = summary_functions)
      super$initialize(name, variables, parents,
                       variable_type, censoring_node = censoring_node, scale)
    },
    print = function() {
      node_class <- class(self)[1]
      cat(sprintf("%s: %s\n", node_class, self$name))
      cat(sprintf("\tVariables: %s\n", paste(self$variables, collapse = ", ")))
      cat(sprintf("\tParents: %s\n", paste(self$parents, collapse = ", ")))
      cat(sprintf("\tSummary Measures: %s\n", paste(sapply(self$summary_functions, function(f){f$name}), collapse = ", ")))

    }
  ),

  active = list(
    summary_functions = function(){
      private$.ltmle_params$summary_functions
    },
    time = function(){
      private$.ltmle_params$time
    }


  ),
  private = list(
    .ltmle_params = NULL
  )
)

define_lnode <- ltmle3_Node$new
