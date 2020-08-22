#' @importFrom purrr reduce
#' @importFrom dplyr full_join
#' @importFrom tmle3 tmle3_Task
#' @export
ltmle3_Task <- R6Class(
  classname = "ltmle3_Task",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Task,
  private = list(
    .force_at_risk = F,
    .non_data_columns = NULL
  ),
  active = list(
    data = function() {
      all_variables <- unlist(lapply(self$npsem, `[[`, "variables"))
      self$get_data(columns = unique(c(private$.non_data_columns, "id" , "t", all_variables)))
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
    },
    extra_summary_measure_columns = function(){
      private$.non_data_columns
    },
    force_at_risk = function(){
      private$.force_at_risk
    }
  ),
  public = list(
    initialize = function(data, npsem, id = "id", time = "t", extra_summary_measure_columns = "at_risk", force_at_risk = F,  ...) {

      #If time and id are not named as "t" and "id" then
      # add columns
      private$.force_at_risk = force_at_risk
      private$.non_data_columns = extra_summary_measure_columns

      if(!inherits(data, "Shared_Data")){
        if(id!="id"){
          data[,"id", with = F] <- data[,id, with = F]
          id = "id"
        }
        if(time!="t"){
          data[, "t", with = F] <- data[,time, with = F]
          time = "t"
        }
        data <- setorder(data, id)
        data <- setorder(data, t)
        shared_data <- data
      } else{
        shared_data <- data
      }

      # Missing values may not be censoring, just information that the value hasn't changed since last time point
      # Censoring should not be handled automatically.

      super$initialize(shared_data, npsem, id = "id", time = "t", add_censoring_indicators = F,  ...)
      private$.uuid <- digest::digest(self$data)
    },
    get_tmle_node = function(node_name, format = FALSE, include_time = F, include_id = T, force_time_value = "No") {
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
      cache_key <- sprintf("%s_%s_%s_%s_%s", node_name, format, include_time, include_id, force_time_value)

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

      if(is.numeric(force_time_value)){
        time <- force_time_value
      } else {
        time <- tmle_node$time
      }
      if(time == "pooled"){
        #If pooled then get value of variable at all times. Use the force_time argument.
        #all_times <- unique(data$t)
        #times <- seq(min(all_times), max(all_times), 1)
        times <- tmle_node$times_to_pool

        nodes_vals_for_all_times <- lapply(times, self$get_tmle_node, node_name = node_name, format=format,  include_time = T, include_id = T )
        return(do.call(rbind, nodes_vals_for_all_times))
      }


      if(FALSE & format == TRUE & !is.null(tmle_node$node_type) ){
        if(tmle_node$node_type == "counting_process") {

        }
      } else {
        all_ids <- unique(data$id)
        orig_data <- data

        last_obs_value <- function(X){
          max_time_index <- nrow(X)

          return(X[max_time_index,])
        }
        # Handle that those whose outcome was not subject to change and do not have a row in dataset
        # Therefore, extract last observed value of node_var up until this time.
        data <- data[data$t <= time, ]
        # get last observed value for each person (either at this time if there is a row or not)
        # TODO If people are censored then they will end up having data drawn here (last observed value)
        # TODO Should we treat censoring and at_risk sets differently?
        data <- data[, last_obs_value(.SD), by = id, .SDcols = c("t", node_var)]
        # set time to current time
        data$t = time
      }

      #Order by id
      #data <- data[order(data$id),]

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
    get_regression_task = function(target_node, scale = FALSE, drop_censored = FALSE, is_time_variant = F, force_time_value = "No") {
      # TODO Not sure what is_time_variant is for.




      if(!is.numeric(force_time_value)){

        cache_key <- sprintf("%s_%s_%s", target_node, scale, is_time_variant)
        cached_data <- get0(cache_key, private$.node_cache, inherits = FALSE)
        if (!is.null(cached_data)) {
          if(drop_censored){
            cached_data <- cached_data[which(cached_data$get_data()$at_risk==1),]
          }
          return(cached_data)
        }
      }


      # If target_node specifies multiple nodes
      # then return the pooled regression task obtained from each node-specific regression task, if possible.
      # E.g. pool over time.
      if(length(target_node)>1){

        all_tasks <- lapply(target_node, self$get_regression_task, scale, drop_censored , is_time_variant)
        all_nodes <- lapply(all_tasks, function(task) task$nodes)
        time_is_node <- sapply(all_nodes, function(node) !is.null(node$time))
        assertthat::assert_that(all(time_is_node))
        assertthat::assert_that(all.equal(unique(unlist(all_nodes, use.names = F)), unlist(all_nodes[[1]], use.names = F)))
        #If they contain the same nodes we can merge
        pooled_data <- do.call(rbind, lapply(all_tasks, function(task) task$get_data()))
        pooled_data <- pooled_data[order(pooled_data$id)]
        pooled_data <- pooled_data[order(pooled_data$t)]
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
      if(is.numeric(force_time_value)){
        time <- force_time_value
      } else {
        time <- target_node_object$time

      }

      at_risk_vars <- target_node_object$at_risk_vars
      past_data <- self$data
      past_data <- past_data[,]
      if(time == "pooled"){
        # TODO summary measures ae expensive to compute. The task cache helps.

        # If node is pooled across time then get pooled regression task
        times <- target_node_object$times_to_pool

        #times <- seq(min(times), max(times), 1)
        all_tasks <- lapply(times, self$get_regression_task, target_node = target_node, scale = scale, drop_censored = F, is_time_variant = is_time_variant )

        all_nodes <- lapply(all_tasks, function(task) task$nodes)
        time_is_node <- sapply(all_nodes, function(node) !is.null(node$time))
        assertthat::assert_that(all(time_is_node))
        assertthat::assert_that(all.equal(unique(unlist(all_nodes, use.names = F)), unlist(all_nodes[[1]], use.names = F)))
        #If they contain the same nodes we can merge
        pooled_data <- do.call(rbind, lapply(all_tasks, function(task) task$get_data()))
        nodes <- all_nodes[[1]]
        # Make sure time is included as covariate
        nodes$covariates <- union("t", nodes$covariates)
        pooled_data <- pooled_data[order(pooled_data$id)]
        pooled_data <- pooled_data[order(pooled_data$t)]

        pooled_regression_task <- sl3_Task$new(
          pooled_data,
          nodes = nodes,
          outcome_type = self$npsem[[target_node[1]]]$variable_type,
          folds = self$folds
        )

        if(!is.numeric(force_time_value)){
          #Store tasks
          assign(cache_key, pooled_regression_task, private$.node_cache)
        }
        if(drop_censored){
          pooled_regression_task <- pooled_regression_task[which(pooled_data$at_risk==1),]
        }
        return(pooled_regression_task)
        }

      parent_names <- target_node_object$parents
      parent_nodes <- npsem[parent_names]

      if(F&is.null(unlist(target_node_object$summary_functions))){
        # No summary functions so simply stack node values of parents
        parent_data <- do.call(cbind, lapply(parent_names, self$get_tmle_node, include_id = F, include_time = F, format = T))
        outcome_data <- self$get_tmle_node(target_node, format = TRUE, include_id = T, include_time = (time == "pooled"), force_time_value = force_time_value)

        data <- cbind(parent_data, outcome_data)
        data$at_risk <- self$get_tmle_node(target_node_object$at_risk_node_name, format = TRUE, include_id = T, include_time = (time == "pooled"), force_time_value = force_time_value)
        # TODO last value
        #data[, paste0("last_val",  setdiff(colnames(outcome_data), c("id", "t"))), with = F] <-
        # TODO custom nodes not handled
        regression_task <- sl3_Task$new(
          Shared_Data$new(data, force_copy = F),
          id = "id",
          time = "t",
          covariates = colnames(parent_data),
          outcome = setdiff(colnames(outcome_data), c("id", "t")),
          outcome_type = target_node_object$variable_type,
          folds = self$folds
        )
        return(regression_task)
      }


      all_vars <- unique(unlist(lapply(npsem, `[[`, "variables")))

      times <- as.vector(sapply(parent_nodes, function(node) node$time))
      parent_covariates <- as.vector(sapply(parent_nodes, function(node) node$variables))
      past_same_time_vars <- unique(parent_covariates[times == time])

      # Note that those with missing rows will be included in outcome_data.
      # There value will be set to last measured value.
      outcome_data <- self$get_tmle_node(target_node, format = TRUE, include_id = T, include_time = (time == "pooled"), force_time_value = force_time_value)

      past_data <- past_data[past_data$t <= time,]

      risk_set <- c()

      risk_set <- c(risk_set, unique(past_data$id))

      skip <- F
      if(length(parent_covariates) == 0 & is.null(unlist(target_node_object$summary_functions))){
        past_data <- data.table(rep(NA, nrow(outcome_data)))
        setnames(past_data, paste0("last_val_", setdiff(colnames(outcome_data), c("t", "id"))))
        past_data$id <- outcome_data$id

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
          rowsums <- rowSums(all_covariate_data[, risk_indicator_cols, with = F])
          keep <- which(is.na(rowsums) | rowsums == length(risk_indicator_cols))
          still_at_risk_ids <- all_covariate_data[keep, id]
          # The risk_set shrinks
          risk_set <- intersect(risk_set, still_at_risk_ids)

          # Indicator whether person is still at risk
          #all_covariate_data$at_risk <- as.numeric(all_covariate_data$id %in% still_at_risk_ids)
          # last measured value of outcome for each person
          # Set outcome of those who are not at risk to last observed value
          last_val <- all_covariate_data[,at_risk_vars$last_val, with = F]
          if(F &!is.null(target_node_object$node_type)){
            if(target_node_object$node_type == "counting_process") {
              # If counting process then last value  should be set to zero (for hazard)
              # technically only true for those not in risk set (but only needed for those people so its fine)
              last_val[] <- 0
            }
          }
          # remove all at_risk based covariates which aren't part of the machine learning
          setDT(all_covariate_data)

          names_last_val <- paste("last_val", outcome, sep  = "_")
          set(all_covariate_data, , names_last_val, last_val)
          set(all_covariate_data, , at_risk_vars$last_val,  NULL)
          if(!is.null( at_risk_vars$at_risk_competing)) set(all_covariate_data, , at_risk_vars$at_risk_competing,  NULL)
          set(all_covariate_data, , at_risk_vars$at_risk,  NULL)
          if(drop_censored){
            all_covariate_data <- all_covariate_data[which(all_covariate_data$id %in% risk_set),]
          }

        } else {
          setDT(all_covariate_data)
          names_last_val <- paste("last_val", outcome, sep  = "_")
          all_covariate_data[, (names_last_val) := NA]
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
      #TODO since number of time rows vary per person, only time-indepdent nodes make sense
      # Keep only node_data for each individual at the time of this tmle node
      node_data <- node_data[!duplicated(node_data$id),]
      node_data$t = time
      nodes$outcome <- outcome
      nodes$covariates <- covariates
      nodes$time <- "t"
      # Merge all data by id. Full join so all rows are kept and NA's added.
      # NA rows will necessarily be dropped if drop_censored = T.
      # Otherwise only the outcome variables should be NA for unobserved/unmonitored outcomes.
      # If one is interested in getting predictions for people not observed at this time
      # then one must add a row with the desired outcome.
      if(skip) {
        regression_data <-  list(all_covariate_data, outcome_data, node_data) %>% reduce(full_join, "id")

      } else {
        regression_data <-  list(all_covariate_data, outcome_data, node_data) %>% reduce(dplyr::inner_join, "id")

      }
      if(private$.force_at_risk){
        regression_data$at_risk <- 1
      } else{
        regression_data$at_risk <- as.numeric(regression_data$id %in% risk_set)

      }
      # For those who had missing rows, set their last_val to the current value of node_var.
      set_last_val_to_current<- intersect(which(is.na(regression_data$last_val)), which(!(regression_data$id %in% risk_set)))
      if(length(set_last_val_to_current) > 0 ){
        set(regression_data,set_last_val_to_current, "last_val", regression_data[set_last_val_to_current, outcome, with = F] )

      }
      # Only keep rows of those in risk set


      regression_data$t = time
      regression_data <- regression_data[order(regression_data$id), ]
      regression_data <- regression_data[order(regression_data$t), ]
      regression_data <- Shared_Data$new(regression_data, force_copy = F)

      regression_task <- sl3_Task$new(
        regression_data,
        nodes = nodes,
        outcome_type = target_node_object$variable_type,
        folds = self$folds
      )
      if(!is.numeric(force_time_value)){
        #Store tasks
        assign(cache_key, regression_task, private$.node_cache)
      }

      if(drop_censored){
        #regression_data <- regression_data[which(regression_data$id %in% risk_set),]

        regression_task <- regression_task[which(regression_task$data$id %in% risk_set),]
      }

      return(regression_task)
    },
    generate_counterfactual_task = function(uuid, new_data, force_at_risk = NULL, through_data=F) {
      if(!through_data){
        if(!("t" %in% colnames(new_data))){

          #If not in id/t format. Then convert to id/t format and handle node names
          node <-  setdiff(colnames(new_data), c("id", "t"))
          dat <- self$get_tmle_node(node, include_time = T, include_id = T)
          node_vars <- sapply(
            node,
            function(node_name) {
              self$npsem[[node_name]]$variables
            }
          )

          set(dat, , node, new_data[,node,with=F])
          setnames(dat, node, node_vars)
          new_data <- dat
        }

        data <- data.table::copy(self$data)

        if(all(c("id", "t") %in% colnames(new_data))){
          #This is needed
          node <-  setdiff(colnames(new_data), c("id", "t"))
          helper <- function(x){
            id <- x$id
            t <- x$t
            index <- which(data$id == id &  data$t==t)

            if(length(index)>0){
              #If t and id in data then replace
              if(data[index,node, with = F][[1]]==x[, node, with = F][[1]]){
                #If same value do nothing
                return()
              }
              set(data, index, node, x[, node, with = F])
              return(1)
            } else {
              #otherwise add new row
              #first extract most recent row for values of other variables
              indices<-which(data$id == id & data$t <= t)
              index <- indices[which.max(data[indices, "t", with = F][[1]])]
              new_row <- data[index,]
              if(new_row[,node, with = F][[1]]==x[, node, with = F][[1]]){
                #If same value do nothing
                return()
              }
              #set new node value and time value
              set(new_row, , node, x[, node, with = F])
              set(new_row, , "t", x[, "t", with = F])
              data <<- rbind(data, new_row)
              return(0)
            }

          }
          #This will update the data object with new values of new_data
          new_data[, helper(.SD), by = c("t", "id"), .SDcols = colnames(new_data)]

          #Shouldn't be needed but clean duplicate rows.
          data <- data[!duplicated(data[, -c("t")])]

          setorder(data, id)
          setorder(data, t)

          #shared_data <- Shared_Data$new(data, force_copy = F)

          new_task <- self$clone()


          new_task$initialize(
            data, self$npsem,
            folds = self$folds,
            row_index = self$row_index,
            t = "t",
            id = "id",
            nodes = self$nodes,
            force_at_risk = ifelse(is.null(force_at_risk), self$force_at_risk, force_at_risk),
            extra_summary_measure_columns = private$.non_data_columns
          )
          return(new_task)
        }


      }

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
      is_pooled <- "pooled" %in% unlist(node_times)

      node_index <- lapply(
        node_times,
        function(time) {
          if(is_pooled) {

            return(1:nrow(self$data))
          }
          which(self$data$t==time)
        }
      )

      old_data <- data.table::copy(self$data[, unique(node_variables), with = F])

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
        row_index = self$row_index,
        force_at_risk = force_at_risk
      )

      return(new_task)
    }

  )
)
