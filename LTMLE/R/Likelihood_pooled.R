

#' Likelihood object that supports a list of pooled LF objects (e.g. LF_fit_pooled)
#' All output behavior is identical to the standard Likelihood object.
#' @export
Likelihood_pooled <- R6Class(
  classname = "Likelihood_pooled",
  portable = TRUE,
  class = TRUE,
  inherit = Likelihood,
  active = list(
    factor_list_pooled = function(){
      self$params$factor_list_pooled
    }
  ),
  public = list(
    initialize = function(factor_list, cache = NULL, ...) {
      params <- args_to_list()
      if (inherits(factor_list, "LF_base")) {
        factor_list <- list(factor_list)
      }

      factor_names <- lapply(factor_list, `[[`, "name")
      factor_list_unpooled <- list()
      for(i in seq_along(factor_names)){
        names <- c(factor_names[[i]])
        factor <- factor_list[[i]]
        for(name in names){
          factor_list_unpooled$name <- factor
        }
      }
      names(factor_list) <-  sapply(factor_names, function(name) paste(name, collapse = "%"))
      #factor_list is used for training
      #factor_list_unpooled is used for prediction
      #Assumes that training is done via mutation

      #factor_list appears just like a standard factor_list except it contains duplicate LF_factors
      #corresponding to the pooled nodes.
      params$factor_list <- factor_list_unpooled
      #factor_list_pooled only contains the unique factors
      params$factor_list_pooled <- factor_list
      if (is.null(cache)) {
        cache <- Likelihood_cache$new()
      }
      private$.cache <- cache

      super$initialize(params)
    },
    validate_task = function(tmle_task) {
      assert_that(is(tmle_task, "tmle3_Task"))

      factor_list <- self$factor_list
      factor_names <- names(factor_list)
      task_nodes <- names(tmle_task$npsem)
      if (!all(factor_names %in% task_nodes)) {
        stop("factor_list and task$npsem must have matching names")
      }
    },
    get_likelihood = function(tmle_task, node, fold_number = "full") {
      likelihood_factor <- self$factor_list[[node]]
      # first check for cached values for this task

      likelihood_values <- self$cache$get_values(likelihood_factor, tmle_task, fold_number, node)
      # note the above contains all values from the pooled task, not just this node.

      if (is.null(likelihood_values)) {
        # if not, generate new ones
        likelihood_values <- likelihood_factor$get_likelihood(tmle_task, fold_number)
        nodes <- names(likelihood_values)
        # Cache all likelihood values for all nodes in likelihood_values
        for(node in nodes) {
          self$cache$set_values(likelihood_factor, tmle_task, 0, fold_number, likelihood_values[, node, with = F], node)
        }
      }
      #Subset to only likelihood values of this node
      likelihood_values <- likelihood_values[, node, with = F]

      return(likelihood_values)
    },
    get_likelihoods = function(tmle_task, nodes = NULL, fold_number = "full") {
      if (is.null(nodes)) {
        nodes <- self$nodes
      }

      if (length(nodes) > 1) {
        all_likelihoods <- lapply(nodes, function(node) {
          self$get_likelihood(tmle_task, node, fold_number)
        })
        likelihood_dt <- as.data.table(all_likelihoods)
        setnames(likelihood_dt, nodes)
        return(likelihood_dt)
      } else {
        return(self$get_likelihood(tmle_task, nodes[[1]], fold_number))
      }
    }
  )
)
