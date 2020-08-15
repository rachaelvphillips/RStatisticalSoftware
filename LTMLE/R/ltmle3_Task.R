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
      keep <- !duplicated(self$data$id)
      id <- self$data$id[keep]
      #Return weights in order of id
      super$weights[keep][order(id)]
    },
    offset = function() {
      keep <- !duplicated(self$data$id)
      id <- self$data$id[keep]
      # Assumes offset are constant in time
      # return offset
      super$offset[keep][order(id)]
    }
  ),
  public = list(
    initialize = function(data, npsem, id = "id", time = "t",  ...) {

      #If time and id are not named as "t" and "id" then
      # add columns
      if(id!="id"){
        data[,"id", with = F] <- data[,id, with = F]
        id = "id"
      }
      if(time!="t"){
        data[, "t", with = F] <- data[,time, with = F]
        time = "t"
      }
      # Missing values may not be censoring, just information that the value hasn't changed since last time point
      # Censoring should not be handled automatically.
      super$initialize(data, npsem, id = "id", time = "t", add_censoring_indicators = F,  ...)
    },
    get_tmle_node = function(node_name, format = FALSE, include_time = F, include_id = T) {
     # returns node value
     # note if row is missing for person at time specified by node
      # then it is assumed the individual was not monitored and was not subject to change
      # the last known value of the node (at previous time) will be used instead

       # node as dt vs node as column
      # scaling
      # caching that accounts for these
      # keep defaults the same
      # use this for get regession task
      # format variables (using format Y ) when
      # categorical should be formatted as factors
      # what does the ate and tsm spec do here
      cache_key <- sprintf("%s_%s_%s_%s", node_name, format, include_time, include_id)

      cached_data <- get0(cache_key, private$.node_cache, inherits = FALSE)
      if (!is.null(cached_data)) {
        return(cached_data)
      }
      tmle_node <- self$npsem[[node_name]]
      node_var <- tmle_node$variables

      if (is.null(node_var)) {
        return(data.table(NULL))
      }

      data <- self$get_data(self$row_index, c("t", "id", node_var))
      if(format == TRUE & !is.null(tmle_node$node_type) ){
        if(tmle_node$node_type == "counting_process") {
          # Convert counting process format to hazard outcome format

          to_hazard <- function(v) {
            first_jump <- min(which(v==1))
            v[] = 0
            v[first_jump] =1
            return(v)
          }
          data_as_haz <- data[, lapply(.SD, to_hazard), by = id, .SDcols = node_var]
          data_as_haz$t <- data$t
          data <- data_as_haz

          time <- private$.npsem[[node_name]]$time
          all_ids <- unique(data$id)
          orig_data <- data

          last_obs_value <- function(X){
            max_time_index <- which.max(X$t)
            return(X[max_time_index,])
          }
          # Handle that those whose outcome was not subject to change do not have a row in dataset
          # Therefore, extract last observed value of node_var up until this time.
          data <- data[data$t <= time, ]
          # get last observed value for each person (either at this time if there is a row or not)
          data <- data[, last_obs_value(.SD), by = id, .SDcols = c("t", node_var)]
          # If person was not monitored at this time then there hazard is 0.
          data[which(data$t!= time), node_var] <- 0
          # set time to current time
          data$t = time

        }
      } else {
        time <- private$.npsem[[node_name]]$time
        all_ids <- unique(data$id)
        orig_data <- data

        last_obs_value <- function(X){
          max_time_index <- which.max(X$t)
          return(X[max_time_index,])
        }
        # Handle that those whose outcome was not subject to change do not have a row in dataset
        # Therefore, extract last observed value of node_var up until this time.
        data <- data[data$t <= time, ]
        # get last observed value for each person (either at this time if there is a row or not)
        data <- data[, last_obs_value(.SD), by = id, .SDcols = c("t", node_var)]
        # set time to current time
        data$t = time
      }





      # # Check those who were not subject to change/monitoring (e.g. missing row at this time)
      # not_monitored <- setdiff(all_ids, id)
      # if(length(not_monitored)>0){
      #   # Get last observed value of outcome, could happen at different times for each person
      #   data_not_monitored <- orig_data[which(orig_data$id %in% not_monitored), ]
      #   last_obs_value <- function(X){
      #     max_time_index <- which.max(X$t)
      #     return(X[max_time_index, ,])
      #   }
      #   data_not_monitored <- orig_data[, last_obs_value(.SD), by = id, .SDcols = c("t", node_var)]
      #   # Match up column order
      #   data_not_monitored <- data_not_monitored[, colnames(data),with = F]
      #   #merge data by rows
      #   data <- rbind(data, data_not_monitored)
      #
      # }
      #Order by id
      data <- data[order(data$id),]

      id <- data$id
      data$id <- NULL

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
        if(!is.data.table(data)){
          data <- data.table(data)
          setnames(data, node_var)
        }
        data$t <- t


      }
      if(include_id){
        if(!is.data.table(data)){
          data <- data.table(data)
          setnames(data, node_var)
        }
        data$id <- id
      }




      assign(cache_key, data, private$.node_cache)

      return(data)
    },
    get_regression_task = function(target_node, scale = FALSE, drop_censored = FALSE, is_time_variant = F) {
      # TODO Not sure what is_time_variant is for.

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
      time <- target_node_object$time
      at_risk_vars <- target_node_object$at_risk_vars
      past_data <- self$get_data()
      past_data <- past_data[,]


      parent_names <- target_node_object$parents
      parent_nodes <- npsem[parent_names]
      all_vars <- unique(unlist(lapply(npsem, `[[`, "variables")))

      times <- as.vector(sapply(parent_nodes, function(node) node$time))
      parent_covariates <- as.vector(sapply(parent_nodes, function(node) node$variables))
      past_same_time_vars <- unique(parent_covariates[times == time])

      # Note that those with missing rows will be included in outcome_data.
      # There value will be set to last measured value.
      outcome_data <- self$get_tmle_node(target_node, format = TRUE, include_id = T)
      #Subset to data up until current time
      past_data <- past_data[past_data$t <= time,]
      # Find ids with observation row at this time. E.g. those who were monitored at this time
      # These people are implicitely in the at-risk set.
      risk_set <- c()
      monitored_ids <- past_data[, any(.SD==time), by = id, .SDcols = c("t")]
      monitored_ids <- monitored_ids[which(monitored_ids[[2]]), id]

      risk_set <- unlist(c(risk_set, monitored_ids))
      # Extract those who were not monitored at this time to be dropped.
      #not_monitored_ids <- monitored_ids[-which(monitored_ids[[2]]), id][[1]]

      # Ensure that variables explicitely mentioned as nodes in npsem are not included in data
      # For example W -> A where W and A both happen at time t. Then the node W should not observe
      # ... the value of A at time t. Setting to NA ensures not future data leakage for summary measures.
      strict_past_vars <- setdiff(all_vars, c("t", "id", past_same_time_vars))
      # Safeguard to ensure no future data leakage
      # Needed for summary functions that depend on both current time values and past time values
      set(past_data,  which(past_data$t == time), strict_past_vars, NA)

      skip <- F
      if(length(parent_covariates) == 0){
        past_data <- (NULL)
        covariates <- c()
        skip <- T

      }



      #This contains all people who have had some observation in the past  up until now
      all_covariate_data <- past_data
      if(!skip){
        summary_measures <- target_node_object$summary_functions

        all_covariate_data <- lapply(summary_measures, function(fun){

          #If this summary measure depends on treatment then set past to t-1
          #This ensures that no summary measures use the outcome
          #And allows summary measures to depend on the most recent L
          #Which happens at the same time as A.
          subset_time <- time

          if(!is.null(strict_past_vars) & any(fun$column_names %in% strict_past_vars)){
            if(!fun$strict_past){
              #
              #warning("Summary measure is not based on strict past and does not respect the time ordering of npsem. Manually handling this.")
              #subset_time <- subset_time - 1
            }
          }

          return(fun$summarize(all_covariate_data, subset_time))}
        )

        all_covariate_data <- all_covariate_data %>% purrr::reduce(dplyr::full_join, by = "id")
        #all_covariate_data$id <- NULL
        covariates <- setdiff(colnames(all_covariate_data), "id")
        # TODO generate names of at_risk_vars with uuid so that no colissions happen with real data
        covariates <- setdiff(covariates, at_risk_vars)

        # If risk set is specified for node then handle this.
        if(!is.null(at_risk_vars)){
          # if any of risk_indicator_cols is 0 then person is no longer at risk
          risk_indicator_cols <-  c(at_risk_vars$at_risk_competing, at_risk_vars$at_risk)
          # Those at risk should have a row of only 1's.
          keep <- which(rowSums(all_covariate_data[, risk_indicator_cols, with = F])==length(risk_indicator_cols))
          still_at_risk_ids <- all_covariate_data[keep, id]

          # The risk_set shrinks
          risk_set <- intersect(risk_set, still_at_risk_ids)

          # Indicator whether person is still at risk
          #all_covariate_data$at_risk <- as.numeric(all_covariate_data$id %in% still_at_risk_ids)
          # last measured value of outcome for each person
          # Set outcome of those who are not at risk to last observed value
          last_val <- all_covariate_data[,at_risk_vars$last_val, with = F]
          if(!is.null(target_node_object$node_type)){
            if(target_node_object$node_type == "counting_process") {
              # If counting process then last value  should be set to zero (for hazard)
              # technically only true for those not in risk set (but only needed for those people so its fine)
              last_val[] <- 0
            }
          }
          # remove all at_risk based covariates which aren't part of the machine learning
          setDT(all_covariate_data)
          set(all_covariate_data, ,"last_val", last_val)
          set(all_covariate_data, , at_risk_vars$last_val,  NULL)
          set(all_covariate_data, , at_risk_vars$at_risk_competing,  NULL)
          set(all_covariate_data, , at_risk_vars$at_risk,  NULL)



        } else {
          setDT(all_covariate_data)

          all_covariate_data[, last_val := NA]
        }
        # Add column specifying those at risk
        #all_covariate_data$at_risk <- as.numeric(all_covariate_data$id %in% risk_set)



        #t_col <- data.table(rep(time, nrow(all_covariate_data)))
        #colnames(t_col) <- "t"
        if("t" %in% colnames(all_covariate_data))  all_covariate_data$t <- NULL

        #all_covariate_data <- cbind(t_col, all_covariate_data)
      }
      if(is_time_variant){
        # TODO Not sure if this is the correct use of is_time_variant
        #covariates <- union("t", covariates)
      }

      nodes <- self$nodes

      node_data <- self$get_data(, unlist(nodes))
      # Keep only node_data for each individual at the time of this tmle node
      node_data <- node_data[self$get_data()$t == time,]
      nodes$outcome <- outcome
      nodes$covariates <- covariates
      nodes$time <- "t"
      # Merge all data by id. Full join so all rows are kept and NA's added.
      # NA rows will necessarily be dropped if drop_censored = T.
      # Otherwise only the outcome variables should be NA for unobserved/unmonitored outcomes.
      # If one is interested in getting predictions for people not observed at this time
      # then one must add a row with the desired outcome.
      if(skip) {
        regression_data <-  list(outcome_data, node_data) %>% reduce(full_join, "id")

      } else {
        regression_data <-  list(all_covariate_data, outcome_data, node_data) %>% reduce(full_join, "id")

      }
      regression_data$at_risk <- as.numeric(regression_data$id %in% risk_set)
      # For those who had missing rows, set their last_val to the current value of node_var.
      set_last_val_to_current<- intersect(which(is.na(regression_data$last_val)), which(!(regression_data$id %in% risk_set)))
      if(length(set_last_val_to_current) > 0 ){
        set(regression_data,set_last_val_to_current, "last_val", regression_data[set_last_val_to_current, outcome, with = F] )

      }
      # Only keep rows of those in risk set
      if(drop_censored){
        regression_data <- regression_data[which(regression_data$id %in% risk_set),]
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
