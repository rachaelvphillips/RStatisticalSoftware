
ltmle3_Node <- R6Class(
  classname = "ltmle3_Node",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Node,
  public = list(
    initialize = function(name, variables, parents = c(), time = NULL, summary_functions = NULL,
                          node_type = NULL,  at_risk_summary_function = NULL, include_competing_risks = F, times_to_pool = NULL, at_risk_node_name = NULL, variable_type = NULL, scale = FALSE) {
      if(!is.list(summary_functions)){
        summary_functions <- list(summary_functions)
      }
      if(time == "pooled" & is.null(times_to_pool)){
        stop("You specified that this node represents a pooled over time likelihood factor but did not supply the times_to_pool argument.")
      }
      # Handle counting process at risk summary functions
      if(!(is.null(node_type)) ) {
        # if(node_type == "counting_process"){
        #   # If counting process ensure that self is in risk set
        #   competing_set <- unique(variables)
        #   at_risk_summary_function$set_name(paste(paste(variables, collapse = "_"), "at_risk", sep = "_") )
        #   # Add last value of this counting process and at risk indicator to summary_functions
        #   # TODO check for duplicate summary functions
        #   last_val_summary <- make_summary_measure_last_value(variables, strict_past = T)
        #   last_val_summary$set_name(paste(paste(variables, collapse = "_"), "at_risk_last_val", sep = "_") )
        #   at_risk_competing_summary <- NULL
        #
        #
        #   # Store variables to be omitted from training set
        #   if(!is.null(at_risk_summary_function)) {
        #     at_risk_name <- at_risk_summary_function$name
        #   } else {
        #     at_risk_name <- NULL
        #   }
        #   if(include_competing_risks){
        #     at_risk_competing_summary <- make_summary_measure_competing_indicator(competing = competing_set, strict_past = T)
        #     at_risk_vars <- list(last_val = last_val_summary$name, at_risk_competing = at_risk_competing_summary$name, at_risk = at_risk_name)
        #     summary_functions <- c(summary_functions, last_val_summary, at_risk_competing_summary, at_risk_summary_function)
        #   } else {
        #
        #   at_risk_vars <- list(last_val = last_val_summary$name, at_risk_competing = NULL, at_risk = at_risk_name)
        #   summary_functions <- c(summary_functions, last_val_summary, at_risk_summary_function)
        # }
        #
        # }
      } else if (!is.null(at_risk_summary_function)){
        at_risk_summary_function$set_name(paste(paste(variables, collapse = "_"), "at_risk", sep = "_") )

        at_risk_name <- at_risk_summary_function$name
        last_val_summary <- make_summary_measure_last_value(variables)
        last_val_summary$set_name(paste(paste(variables, collapse = "_"), "at_risk_last_val", sep = "_") )

        at_risk_vars <- list(last_val = last_val_summary$name, at_risk = at_risk_name)
        summary_functions <- c(summary_functions, last_val_summary, at_risk_summary_function)
      } else{
        at_risk_vars <- NULL
      }
      private$.ltmle_params = list(at_risk_variable = at_risk_variable, time = time, times_to_pool = times_to_pool, summary_functions = summary_functions, node_type = node_type, at_risk_vars = at_risk_vars)
      super$initialize(name, variables, parents,
                       variable_type, censoring_node = NULL, scale)
    },
    print = function() {
      node_class <- class(self)[1]
      cat(sprintf("%s: %s\n", node_class, self$name))
      cat(sprintf("\tVariables: %s\n", paste(self$variables, collapse = ", ")))
      cat(sprintf("\tParents: %s\n", paste(self$parents, collapse = ", ")))
      cat(sprintf("\tSummary Measures: %s\n", paste(unlist(sapply(self$summary_functions, function(f){f$name})), collapse = ", ")))

    }
  ),

  active = list(
    summary_functions = function(){
      private$.ltmle_params$summary_functions
    },
    time = function(){
      private$.ltmle_params$time
    },
    node_type = function(){
      private$.ltmle_params$node_type
    },
    at_risk_vars = function(){
      private$.ltmle_params$at_risk_vars
    },
    times_to_pool = function(){
      private$.ltmle_params$times_to_pool
    },
    at_risk_node_name = function(){
      private$.ltmle_params$at_risk_node_name
    }



  ),
  private = list(
    .ltmle_params = NULL

  )
)

define_lnode <- ltmle3_Node$new
