

#' A likelihood that represents a collapsed likelihood,
#' meaning that some factors of the underlying likelihood are collapsed into one factor by
#' multiplication of probabilities.
#' @export
Collapsed_Likelihood <- R6Class(
  classname = "Collapsed_Likelihood",
  portable = TRUE,
  class = TRUE,
  inherit = Likelihood,
  active = list(
    nodes = function() {
      return(names(private$.collapse_by))
    },
    hidden_likelihood = function(){
      return(private$.hidden_likelihood)
    }
  ),
  public = list(

    initialize = function(hidden_likelihood, collapse_by, cache = NULL, ...) {
      params <- args_to_list()
      factor_list <- list()
      for(node in names(collapse_by)){
        hidden_nodes <- collapse_by[[node]]
        factors <- lapply(hidden_nodes, function(hidden) {
          hidden_likelihood$factor_list[[hidden]]
          })
        factor_list <- c(factor_list, list(LF_collapsed$new(node, factors)))
      }
      #names(factor_list) <- names(collapse_by)
      params$factor_list <- factor_list
      private$.hidden_likelihood <- hidden_likelihood

      private$.collapse_by <- collapse_by

      super$initialize(factor_list, cache = cache, ...)
    },
    get_hidden_nodes = function(node){
      return(private$.get_hidden_nodes(node))
    }

  # get_likelihood = function(tmle_task, node, fold_number = "full") {
  #   hidden_nodes <- private$.get_hidden_nodes(node)
  #   if(length(hidden_nodes)==1){
  #     return(self$hidden_likelihood$get_likelihood(tmle_task, hidden_nodes, fold_number))
  #   }
  #
  #   hidden_likelihoods <- do.call(cbind, lapply(hidden_nodes, self$hidden_likelihood$get_likelihood, tmle_task = tmle_task, fold_number = fold_number))
  #   likelihood <- unlist(apply(hidden_likelihoods, 1, prod))
  #   likelihood <- data.table(likelihood)
  #   setnames(likelihood, node)
  #   return(likelihood)
  # }
  ),
  private = list(
    .hidden_likelihood = NULL,
    .collapse_by = NULL,
    .get_hidden_nodes = function(node){
      return(private$.collapse_by[[node]])
    }
  )

)
