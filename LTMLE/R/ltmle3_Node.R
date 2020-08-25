
ltmle3_Node <- R6Class(
  classname = "ltmle3_Node",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Node,
  public = list(
    initialize = function(name, variables, parents = c(), time = NULL, summary_functions = NULL,
                          at_risk_map = NULL, missing_row_implies_not_at_risk = T, times_to_pool = NULL, variable_type = NULL, scale = FALSE) {
      if(!is.list(summary_functions)){
        summary_functions <- list(summary_functions)
      }
      if(time == "pooled" & is.null(times_to_pool)){
        stop("You specified that this node represents a pooled over time likelihood factor but did not supply the times_to_pool argument.")
      }
      # Handle counting process at risk summary functions
      if (!is.null(at_risk_map)){
        if(!is.character(at_risk_map)){
          at_risk_map$set_name(paste(paste(variables, collapse = "_"), "at_risk", sep = "_") )
        }
      }
      private$.ltmle_params = list(missing_row_implies_not_at_risk = missing_row_implies_not_at_risk, at_risk_map = at_risk_map, time = time, times_to_pool = times_to_pool, summary_functions = summary_functions)
      super$initialize(name, variables, parents,
                       variable_type, censoring_node = NULL, scale)
    },
    print = function() {
      node_class <- class(self)[1]
      cat(sprintf("%s: %s\n", node_class, self$name))
      cat(sprintf("\tVariables: %s\n", paste(self$variables, collapse = ", ")))
      cat(sprintf("\tParents: %s\n", paste(self$parents, collapse = ", ")))
      cat(sprintf("\tSummary Measures: %s\n", paste(unlist(sapply(self$summary_functions, function(f){f$name})), collapse = ", ")))

    },
  risk_set = function(data, time){
    #Assumes data == data[t <= time,] and time is single number
    #Computes, for this node, the id's of those in data at risk of changing their value at this time
    at_risk_map <- self$at_risk_map
    missing_not_at_risk <- private$.ltmle_params$missing_row_implies_not_at_risk
    if(missing_not_at_risk){
      keep_id <- unique(data[data$t == time, id])
      data <- data[data$id %in% keep_id,]
    }
    if(is.null(at_risk_map)) risk_set <- unique(data$id)
    if(is.character(self$at_risk_map)) {
      if(missing_not_at_risk){
        #If those missing rows are not at risk
        #then only check for those with rows at this time
        risk_set <- data[t == time & data[,at_risk_map,with = F, drop = T]==1, "id", with = F, drop = T][[1]]
      } else{
        #Otherwise find the last value of risk indicator
        data <- data[!duplicated(data$id, fromLast = T), c("id", at_risk_map), with = F]
        #data <- data[, slice_tail(.SD), by = id, .SDcols = at_risk_map]
        risk_set <- data$id[data[[at_risk_map]] ==1]
      }
    } else{
      risk_set <-data[which(at_risk_map$summarize(data,time)[,at_risk_map$name, with = F, drop = T]==1), c("id"), with = F, drop = T][[1]]

    }
    return(unlist(risk_set, use.names = F))

  }
  ),
  active = list(
    summary_functions = function(){
      private$.ltmle_params$summary_functions
    },
    missing_not_at_risk = function(){
      private$.ltmle_params$missing_row_implies_not_at_risk
    },
    time = function(){
      private$.ltmle_params$time
    },
    at_risk_map = function(){
      private$.ltmle_params$at_risk_map
    },
    times_to_pool = function(){
      private$.ltmle_params$times_to_pool
    }



  ),
  private = list(
    .ltmle_params = NULL

  )
)

define_lnode <- ltmle3_Node$new
