#' @importFrom purrr reduce
#' @importFrom dplyr full_join
#' @export
ltmle3_Task <- R6Class(
  classname = "ltmle3_Task",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Task,
  active = list(
    data = function() {
      all_variables <- unlist(lapply(self$npsem, `[[`, "variables"))
      self$get_data(columns = unique(c("id" , "t", all_variables)))
    }
  ),
  public = list(
    initialize = function(data, npsem, id = "id", ...) {
      super$initialize(data, npsem, id = "id", ...)
    },
    get_tmle_node = function(node_name, format = FALSE, include_time = F) {
      # node as dt vs node as column
      # scaling
      # caching that accounts for these
      # keep defaults the same
      # use this for get regession task
      # format variables (using format Y ) when
      # categorical should be formatted as factors
      # what does the ate and tsm spec do here
      cache_key <- sprintf("%s_%s_%s", node_name, format, include_time)

      cached_data <- get0(cache_key, private$.node_cache, inherits = FALSE)
      if (!is.null(cached_data)) {
        return(cached_data)
      }
      tmle_node <- self$npsem[[node_name]]
      node_var <- tmle_node$variables
      if (is.null(node_var)) {
        return(data.table(NULL))
      }
      data <- self$get_data(self$row_index, c("t", node_var))

      time <- private$.npsem[[node_name]]$time
      data <- data[data$t == time,]
      t <- data$t
      data$t <- NULL



      if ((ncol(data) == 1)) {
        data <- unlist(data, use.names = FALSE)
      }

      if (format == TRUE) {

        var_type <- tmle_node$variable_type

        data <- var_type$format(data)
        data <- self$scale(data, node_name)
        data <- data.table(data)
        setnames(data, node_var)


      }

      if(include_time){
        data <- data.table(cbind(t ,data))
        setnames(data, c("t", node_var))
      }



      assign(cache_key, data, private$.node_cache)

      return(data)
    },
    get_regression_task = function(target_node, scale = FALSE, drop_censored = FALSE, is_time_variant = FALSE) {
      npsem <- self$npsem
      target_node_object <- npsem[[target_node]]
      outcome <- target_node_object$variables

      parent_names <- target_node_object$parents
      parent_nodes <- npsem[parent_names]
      #times <- as.vector(sapply(parent_nodes, function(node) node$time))
      outcome_data <- self$get_tmle_node(target_node, format = TRUE)

      time <- target_node_object$time
      past_data <- self$get_data()
      IsTreatmentNode <- stringr::str_detect(target_node, "[AY]")
      skip <- F
      if(IsTreatmentNode){
        past_data <- past_data[past_data$t <= time,]

      }
      else {
        if(time ==1){
          past_data <- data.table(NULL)
          covariates <- c()
          skip = T
        }
        else{
          past_data <- past_data[past_data$t <= time - 1,]

        }
      }


      all_covariate_data <- past_data
      if(!skip){
        summary_measures <- target_node_object$summary_functions

        all_covariate_data <- lapply(summary_measures, function(fun){
          if(IsTreatmentNode){
            #If this summary measure depends on treatment then set past to t-1
            #This ensures that no summary measures use the outcome
            #And allows summary measures to depend on the most recent L
            #Which happens at the same time as A.
            if(length(intersect(fun$column_names, outcome))>0){
              all_covariate_data <- all_covariate_data[all_covariate_data$t <= time - 1,]
            }
          }
          return(fun$summarize(all_covariate_data))}
        )
        all_covariate_data <- all_covariate_data %>% purrr::reduce(dplyr::full_join, by = "id")
        print(all_covariate_data)

        covariates <- setdiff(colnames(all_covariate_data), "id")
        t_col <- data.table("t" = rep(time, nrow(all_covariate_data)))
        all_covariate_data <- cbind(t_col, all_covariate_data)
      }

      nodes <- self$nodes

      node_data <- self$get_data(, unlist(nodes))
      # Keep only node_data for each individual once
      # Assumes that these nodes are constant in time
      node_data <- node_data[!duplicated(node_data$id),]
      nodes$outcome <- outcome
      nodes$covariates <- covariates
      regression_data <- do.call(cbind, list(all_covariate_data, outcome_data, node_data))
      print(regression_data)
      regression_task <- sl3_Task$new(
        regression_data,
        nodes = nodes,
        outcome_type = target_node_object$variable_type,
        folds = self$folds
      )

      return(regression_task)
    },
    generate_counterfactual_task = function(uuid, new_data) {
      # for current_factor, generate counterfactual values
      node_names <- names(new_data)

      node_variables <- sapply(
        node_names,
        function(node_name) {
          self$npsem[[node_name]]$variables
        }
      )

      node_times <- sapply(
        node_names,
        function(node_name) {
          self$npsem[[node_name]]$time
        }
      )
      node_index <- lapply(
        node_times,
        function(time) {
          which(self$get_data()$t==time)
        }
      )

      old_data <- data.table::copy(self$get_data()[, unique(node_variables), with = F])

      lapply(seq_along(node_index), function(i){
        index <- node_index[[i]]
        var <- node_variables[[i]]
        set(old_data, index, var, new_data[,node_names[[i]],with=F])
      })

      new_data <- old_data

      #setnames(new_data, node_names, node_variables)

      new_task <- self$clone()
      new_column_names <- new_task$add_columns(new_data, uuid)

      new_task$initialize(
        self$internal_data, self$npsem,
        column_names = new_column_names,
        folds = self$folds,
        row_index = self$row_index
      )

      return(new_task)
    }
  )
)
