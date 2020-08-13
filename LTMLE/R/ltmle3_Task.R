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
    },
    weights = function() {
      # Assumes weights are constant in time
      super$weights[!duplicated(self$data$id)]
    },
    offset = function() {
      # Assumes offset are constant in time
      super$offset[!duplicated(self$data$id)]
    }
  ),
  public = list(
    initialize = function(data, npsem, id = "id", time = "t",  ...) {

      #If time and id are not named as "t" and "id" then
      #add columns
      if(id!="id"){
        data[,"id", with = F] <- data[,id, with = F]
        id = "id"
      }
      if(time!="t"){
        data[, "t", with = F] <- data[,time, with = F]
        time = "t"
      }

      super$initialize(data, npsem, id = "id", time = "t", add_censoring_indicators = F,  ...)
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
    get_regression_task = function(target_node, scale = FALSE, drop_censored = FALSE, is_time_variant = F) {

      # If target_node specifies multiple nodes
      # then return the pooled regression task obtained from each node-specific regression task, if possible.
      if(length(target_node)>1){
        all_tasks <- lapply(target_node, self$get_regression_task, scale, drop_censored , is_time_variant)
        all_nodes <- lapply(all_tasks, function(task) task$nodes)
        time_is_node <- sapply(all_nodes, function(node) !is.null(node$time))
        assertthat::assert_that(all(time_is_node))
        assertthat::assert_that(all.equal(unique(unlist(all_nodes, use.names = F)), unlist(all_nodes[[1]], use.names = F)))
        #If they contain the same nodes we can merge
        pooled_data <- do.call(rbind, lapply(all_tasks, function(task) task$data))
        nodes <- all_nodes[1]
        # Make sure time is included as covariate
        nodes$covariates <- union("t", nodes$covariates)


        pooled_regression_task <- sl3_Task$new(
          pooled_data,
          nodes = nodes,
          outcome_type = self$npsem[[target_node[1]]]$variable_type,
          folds = self$folds
        )
        return(pooled_regression_task)
      }


      npsem <- self$npsem
      target_node_object <- npsem[[target_node]]
      outcome <- target_node_object$variables
      #get t-1 data and then get all t varialbes of past
      time <- target_node_object$time
      past_data <- self$get_data()
      past_data <- past_data[,]

      if(drop_censored & !is.null(target_node_object$censoring_node)){
        # Censoring node should represent a n*t binary vector which is 1 for individual n
        # at time t if the individual is censored/dies at that time.
        # The time of failure/censoring is extracted from this vector.
        censoring_node_var <- target_node_object$censoring_node$variables
        assertthat::assert_that(all(past_data[, censoring_node_var, with = F][[1]] %in% c(0,1)), msg = "Error: drop_censored is T for non-binary censoring variable ")
        getT <- function(X){
          event_index <- which(X[,censoring_node_var,with=F]==1)
          event_index <- min(event_index)
          if(length(event_index)==0){
            return(Inf)
          }
          event_t <- X[event_index, "t", with = F][[1]]

          return(event_t)
        }
        failure_times <- past_data[, getT(.SD), by = id, .SDcols = c("t", censoring_node_var)][[2]]

        dont_keep <- which(failure_times < time)
      } else {
        dont_keep <- c()
      }

      parent_names <- target_node_object$parents
      parent_nodes <- npsem[parent_names]


      times <- as.vector(sapply(parent_nodes, function(node) node$time))
      parent_covariates <- as.vector(sapply(parent_nodes, function(node) node$variables))
      add_to_past <- unique(parent_covariates[times == time])

      outcome_data <- self$get_tmle_node(target_node, format = TRUE)
      past_data <- past_data[past_data$t <= time,]
      ignore_cols <- setdiff(colnames(past_data), c("t", "id", add_to_past))
      # Safeguard to ensure no future data leakage
      # Handling down the line should ensure that this step isn't ever needed
      # But better safe than sorry
      set(past_data,  which(past_data$t == time), ignore_cols, NA)
      set(past_data, , setdiff(colnames(past_data), c("t", "id", unique(parent_covariates))) , NA)

      skip <- F
      if(length(parent_covariates) == 0){
        past_data <- data.table(NULL)
        covariates <- c()
        skip <- T
      }




      all_covariate_data <- past_data
      if(!skip){
        summary_measures <- target_node_object$summary_functions

        all_covariate_data <- lapply(summary_measures, function(fun){

          #If this summary measure depends on treatment then set past to t-1
          #This ensures that no summary measures use the outcome
          #And allows summary measures to depend on the most recent L
          #Which happens at the same time as A.

          assertthat::assert_that(all(fun$column_names %in% c("id" , "t", parent_covariates)), msg = sprintf("Error: Summary_meaure depends on covariates not found in parent nodes. Summary measure: %s & Parent covariates: %s", fun$column_names, parent_covariates))
          if(length(intersect(fun$column_names, ignore_cols)) ==  length((fun$column_names))){
            all_covariate_data <- all_covariate_data[all_covariate_data$t <= time - 1,]
          }

          return(fun$summarize(all_covariate_data))}
        )
        all_covariate_data <- all_covariate_data %>% purrr::reduce(dplyr::full_join, by = "id")


        covariates <- setdiff(colnames(all_covariate_data), "id")
        t_col <- data.table(rep(time, nrow(all_covariate_data)))
        colnames(t_col) <- "t"
        if("t" %in% colnames(all_covariate_data))  all_covariate_data$t <- NULL

        all_covariate_data <- cbind(t_col, all_covariate_data)
      }
      if(is_time_variant){
        # TODO Not sure if this is the correct use of is_time_variant
        covariates <- union("t", covariates)
      }

      nodes <- self$nodes

      node_data <- self$get_data(, unlist(nodes))
      # Keep only node_data for each individual at the time of this tmle node
      node_data <- node_data[self$data$t == time,]
      nodes$outcome <- outcome
      nodes$covariates <- covariates
      nodes$time <- "t"
      regression_data <- do.call(cbind, list(all_covariate_data, outcome_data, node_data))
      if(drop_censored & !is.null(target_node_object$censoring_node) & length(dont_keep)>0){
        regression_data <- regression_data[-dont_keep,, with = F]
      }
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
