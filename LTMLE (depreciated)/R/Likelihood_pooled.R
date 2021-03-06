

#' Likelihood object that supports a list of pooled LF objects (e.g. LF_fit_pooled)
#' All output behavior is identical to the standard Likelihood object.
#' @export
Likelihood_pooled <- R6Class(
  classname = "Likelihood_pooled",
  portable = TRUE,
  class = TRUE,
  inherit = Lrnr_base,
  active = list(
    factor_list_pooled = function(){
      return(self$params$factor_list_pooled)
    },
    factor_list = function() {
      return(self$params$factor_list)
    },
    nodes = function() {
      return(names(self$factor_list))
    },
    cache = function() {
      return(private$.cache)
    },
    censoring_nodes = function() {
      return(private$.censoring_nodes)
    }
  ),
  private = list(
    .train_sublearners = function(tmle_task) {
      factor_fits <- lapply(self$factor_list_pooled, function(factor) {
        factor$delayed_train(tmle_task) })
      result <- bundle_delayed(factor_fits)
      return(result)
    },
    .train = function(tmle_task, factor_fits) {
      factor_list <- self$factor_list_pooled
      for (i in seq_along(factor_list)) {
        factor_list[[i]]$train(tmle_task, factor_fits[[i]])
      }
      # TODO: mutating factor list of Lrnr_object instead of returning a fit
      #       which is not what sl3 Lrnrs usually do

      censoring_nodes <- lapply(tmle_task$npsem, function(node) {
        node$censoring_node$name
      })

      names(censoring_nodes) <- names(tmle_task$npsem)
      private$.censoring_nodes <- censoring_nodes
      return("trained")
    },
    .predict = function(tmle_task) {
      stop("predict method doesn't work for Likelihood. See Likelihood$get_likelihoods for analogous method")
    },
    .chain = function(tmle_task) {
      stop("chain method doesn't work for Likelihood. Currently, no analogous functionality")
    },
    .cache = NULL,
    .censoring_nodes = NULL
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
          factor_list_unpooled[[name]] <- factor
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

    print = function() {
      lapply(self$factor_list, print)
      invisible(NULL)
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
    get_likelihood = function(tmle_task, node, fold_number = "full", drop_id = F, drop_time = F, to_wide = F, expand = expand) {
      # The nonexpanded likelihood is simply for fast montecarlo simulation when you don't need the degenerate probabilities
      # As a result, these likelihoods will not be cached.
      likelihood_factor <- self$factor_list[[node]]
      # first check for cached values for this task

      likelihood_values <- self$cache$get_values(likelihood_factor, tmle_task, fold_number, node ="")
      if(!expand & !is.null(likelihood_values)) {
        #Only store the full likelihood
        #Regression task should be cached so this is cheap
        keep <- tmle_task$get_regression_task(target_node, expand = T)$data$at_risk == 1
        likelihood_values <- likelihood_values[keep,]
      }

      # note the above contains all values from the pooled task, not just this node.

      if (is.null(likelihood_values)) {
        # if not, generate new ones

        likelihood_values <- likelihood_factor$get_likelihood(tmle_task, fold_number, expand = expand)

        #nodes <- setdiff(names(likelihood_values), c("id", "t"))
        #names_of <- colnames(likelihood_values)
        #keep_cols <- intersect(c("t", "id", grep(node,names_of , value = T)), names_of)
        # Cache all likelihood values for all nodes in likelihood_values.
        #for(node in nodes) {
        if(expand){
          self$cache$set_values(likelihood_factor, tmle_task, 0, fold_number, likelihood_values, node = "")
        }
        #}
      }
      #Subset to only likelihood values of this node
      #In case node is pooled, then grab all times
      names_of <- colnames(likelihood_values)
      keep_cols <- intersect(c("t", "id", grep(node, names_of, value = T)), names_of)


      likelihood_values <- likelihood_values[,  keep_cols, with = F]

      if(to_wide & length(unique(likelihood_values$t))==1){
        likelihood_values$t <- NULL
      }
      else if(to_wide){
        likelihood_values <- reshape(likelihood_values, idvar = "id", timevar = "t", direction = "wide")

      }
      if(drop_id & "id" %in% colnames(likelihood_values)) likelihood_values$id <- NULL
      if(drop_time & "t" %in% colnames(likelihood_values)) likelihood_values$t <- NULL
      return(likelihood_values)
    },
    get_likelihoods = function(tmle_task, nodes = NULL, fold_number = "full", drop_id = F, drop_time = F, to_wide = F, expand = T) {
      if (is.null(nodes)) {
        nodes <- self$nodes
      }

      if (length(nodes) > 1) {
        all_likelihoods <- lapply(nodes, function(node) {
          self$get_likelihood(tmle_task, node, fold_number, to_wide = to_wide, expand =  expand)
        })
        contains_t <- all(unlist(lapply(all_likelihoods, function(lik){
          "t" %in% colnames(lik)
        })))
        if(contains_t){
          likelihood_dt <- all_likelihoods %>% reduce(full_join, c("id", "t"))#as.data.table(all_likelihoods)
        } else{
          all_likelihoods <- lapply(all_likelihoods, function(lik){
            if(!(all(c("id", "t") %in% colnames(lik)))){
              if("t" %in% colnames(lik)) lik$t <- NULL
              return(lik)
            }
            if(length(unique(lik$t))==1){
              if("t" %in% colnames(lik)) lik$t <- NULL
              return(lik)
            }
            reshape(lik, idvar = "id", timevar = "t", direction = "wide")
          })

          likelihood_dt <- all_likelihoods %>% reduce(full_join, c("id"))#as.data.table(all_likelihoods)
          if("t" %in% colnames(likelihood_dt)) likelihood_dt$t <- NULL
        }
        #setnames(likelihood_dt, nodes)
        if(drop_id) likelihood_dt$id <- NULL
        if(drop_time & t %in% colnames(likelihood_dt)) likelihood_dt$t <- NULL
        return(likelihood_dt)
      } else {
        likelihood <- self$get_likelihood(tmle_task, nodes[[1]], fold_number, to_wide = to_wide, drop_id = drop_id, drop_time = drop_time, drop = T)
        return(likelihood)
      }
    },
    get_possible_counterfactuals = function(nodes = NULL) {
      # get factors for nodes
      factor_list <- self$factor_list
      if (!is.null(nodes)) {
        factor_list <- factor_list[nodes]
      }

      all_levels <- lapply(factor_list, function(likelihood_factor) {
        likelihood_factor$variable_type$levels
      })
      all_levels <- all_levels[!(sapply(all_levels, is.null))]
      level_grid <- expand.grid(all_levels)
      return(level_grid)
    },
    base_train = function(task, pretrain) {
      self$validate_task(task)
      fit_object <- private$.train(task, pretrain)
      new_object <- self$clone() # copy parameters, and whatever else
      new_object$set_train(fit_object, task)
      return(new_object)
    },
    add_factors = function(factor_list) {
      if (inherits(factor_list, "LF_base")) {
        factor_list <- list(factor_list)
      }

      factor_names <- sapply(factor_list, `[[`, "name")

      # train factors if necessary
      factor_list <- lapply(factor_list, train_lf, self$training_task)

      # add factors to list of factors
      private$.params$factor_list[factor_names] <- factor_list
    },
    sample = function(tmle_task = NULL, sample_lib = NULL) {
      # for now assume nodes are in order
      # TODO: order nodes based on dependencies
      if (is.NULL(sample_lib = NULL)) {
        nodes <- names(self$factor_list)
        sample_lib <- rep(list(NULL), length(nodes))
        names(sample_lib) <- nodes
      }

      for (node in names(self$factor_list)) {
        tmle_task <- factor_list$node$sample(tmle_task, sample_lib$node)
      }

      return(tmle_task)
    }
  )
)
