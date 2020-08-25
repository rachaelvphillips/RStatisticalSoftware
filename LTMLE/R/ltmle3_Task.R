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
    .non_data_columns = NULL,
    .thin = NULL
  ),
  active = list(
    data = function() {
      all_variables <- unlist(lapply(self$npsem, `[[`, "variables"))

      if(!is.null(private$.non_data_columns)){
        columns <- unique(c(private$.non_data_columns, "id" , "t", all_variables))
      } else{
        columns <- unique(c( "id" , "t", all_variables))

      }
      self$get_data(columns =columns)
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
    thin = function(){
      private$.thin
    },
    extra_summary_measure_columns = function(){
      private$.non_data_columns
    },
    force_at_risk = function(){
      private$.force_at_risk
    }
  ),
  public = list(
    initialize = function(data, npsem, id = "id", time = "t", extra_summary_measure_columns = "at_risk", force_at_risk = F, thin = T, ...) {

      #If time and id are not named as "t" and "id" then
      # add columns
      private$.force_at_risk = force_at_risk
      private$.non_data_columns = extra_summary_measure_columns
      private$.thin = thin

      if(!inherits(data, "Shared_Data")){
        if(id!="id"){
          data[,"id", with = F] <- data[,id, with = F]
          id = "id"
        }
        if(time!="t"){
          data[, "t", with = F] <- data[,time, with = F]
          time = "t"
        }
        data <- setkey(data, id, t)
        shared_data <- data
      } else{
        shared_data <- data
      }

      # Missing values may not be censoring, just information that the value hasn't changed since last time point
      # Censoring should not be handled automatically.

      super$initialize(shared_data, npsem, id = "id", time = "t", add_censoring_indicators = F,  ...)
      private$.uuid <- digest::digest(self$data)
    },
    get_tmle_node = function(node_name, format = FALSE, include_time = F, include_id = T, force_time_value = "No", expand = F, compute_risk_set = T) {
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
      cache_key <- sprintf("%s_%s_%s_%s", node_name, format, force_time_value, expand)

      cached_data <- get0(cache_key, private$.node_cache, inherits = FALSE)
      if (!is.null(cached_data)) {
        if(!include_time){
          cached_data$t <- NULL
        }
        if(!include_id){
          cached_data$id <- NULL
        }
        return(cached_data)
      }
      tmle_node <- self$npsem[[node_name]]
      node_var <- tmle_node$variables

      if (is.null(node_var)) {
        return(data.table(NULL))
      }



      if(is.numeric(force_time_value)){
        time <- force_time_value
      } else {
        time <- tmle_node$time
      }

      if(time == "pooled"){
        at_risk_map <- tmle_node$at_risk_map
        time <- tmle_node$times_to_pool
        if(expand  | !is.null(tmle_node$at_risk_map) | !tmle_node$missing_not_at_risk){
          data <- print(lapply(time, self$get_tmle_node, node_name= node_name, format = format, include_time = T, include_id = T, expand = expand))
          #setkey(data, id , t)
          return(data)
        }
        else {
          data <- self$data
          data <- data[t %in% time , c("id", "t", node_var), with = F ]
        }

      } else {

        data <- self$data[t <= time, ]

        if(compute_risk_set){
          risk_set <- tmle_node$risk_set(data, time)
        }
        data <- data[, c("id", "t", node_var), with = F]

        if(expand){
          # Get most recent value for all
          data <- data[, last(.SD), by = id]
          if(compute_risk_set){
            data$at_risk <- as.numeric(data$id %in% risk_set)
            set(data, ,"last_val", unlist(self$data[t <= time - 1, last(.SD), by = id, .SDcols = c(node_var)][[2]], use.names = F))
          }

        } else {
          # Get most recent value for all those at risk
          if(compute_risk_set){
            data <- data[id %in% risk_set, last(.SD), by = id]
          } else {
            data <- data[, last(.SD), by = id]
          }

        }

      }
      set(data,, "t", time)
      if (format == TRUE) {
        data_node <- data[, node_var, with = F]
        if ((ncol(data_node) == 1)) {
          data_node <- unlist(data_node, use.names = FALSE)
        }
        var_type <- tmle_node$variable_type
        data_node <- var_type$format(data_node)
        data_node <- self$scale(data_node, node_name)
        set(data,, node_var, data_node)
      }

      assign(cache_key, data, private$.node_cache)

      if(!include_time){
        data$t <- NULL
      }
      if(!include_id){
        data$id <- NULL
      }



      return(data)
    },
    get_regression_task = function(target_node, scale = FALSE, drop_censored = FALSE, is_time_variant = F, force_time_value = "No", expand = F) {
      # TODO Not sure what is_time_variant is for.

      if(!is.numeric(force_time_value)){
        cache_key <- sprintf("%s_%s_%s_%s", target_node, scale, is_time_variant, expand)
        cached_data <- get0(cache_key, private$.node_cache, inherits = FALSE)
        if (!is.null(cached_data)) {
          return(cached_data)
        }
      }

      # If target_node specifies multiple nodes
      # then return the pooled regression task obtained from each node-specific regression task, if possible.
      # E.g. pool over time.
      if(length(target_node)>1){
        all_tasks <- lapply(target_node, self$get_regression_task, scale, drop_censored , is_time_variant, expand = expand)
        all_nodes <- lapply(all_tasks, function(task) task$nodes)
        time_is_node <- sapply(all_nodes, function(node) !is.null(node$time))
        pooled_data <- do.call(rbind, lapply(all_tasks, function(task) task$get_data()))
        setkey(pooled_data, id, t)
        nodes <- all_nodes[[1]]
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

      past_data <- self$data

      if(time == "pooled"){
        # TODO summary measures ae expensive to compute. The task cache helps.

        # If node is pooled across time then get pooled regression task
        times <- target_node_object$times_to_pool
        all_tasks <- lapply(times, self$get_regression_task, target_node = target_node, scale = scale, drop_censored = F, is_time_variant = is_time_variant, expand = expand )
        all_nodes <- lapply(all_tasks, function(task) task$nodes)
        time_is_node <- sapply(all_nodes, function(node) !is.null(node$time))
        assertthat::assert_that(all(time_is_node))
        assertthat::assert_that(all.equal(unique(unlist(all_nodes, use.names = F)), unlist(all_nodes[[1]], use.names = F)))
        #If they contain the same nodes we can merge
        pooled_data <- do.call(rbind, lapply(all_tasks, function(task) task$get_data()))
        nodes <- all_nodes[[1]]
        # Make sure time is included as covariate
        nodes$covariates <- union("t", nodes$covariates)
        setkey(pooled_data, id, t)
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
        return(pooled_regression_task)
      }

      parent_names <- target_node_object$parents
      parent_nodes <- npsem[parent_names]

      if(is.null(unlist(target_node_object$summary_functions))){
        # No summary functions so simply stack node values of parents
        parent_data <- do.call(cbind, lapply(parent_names, self$get_tmle_node, include_id = T, include_time = F, format = T, expand = T, compute_risk_set = F))
        outcome_data <- self$get_tmle_node(target_node, format = TRUE, include_id = T, include_time = T, force_time_value = force_time_value, expand = expand)
        data <- merge(outcome_data, parent_data, by = id, all.x = T)
        covariates <- colnames(parent_data)
        outcome = setdiff(colnames(outcome_data), c("id", "t", "last_val", "at_risk"))
        if((time == "pooled")){
          covariates <- c(covariates, "t")
        }
        outcome_index <-  match(outcome, colnames(outcome_data))
        if(length(parent_data)>0){
          cov_index <- ncol(outcome_data) + match(covariates, colnames(parent_data))
        } else {
         cov_index <- c()
        }
        #Due to time indexing, we do not have unique column names.
        #In order to support pooling across time, we shouldn't use node names as column names
        #important that outcome variable name doesnt change
        setnames(data, make.unique(colnames(data)))
        covariates <- colnames(data)[cov_index]
        outcome <- colnames(data)[outcome_index]

        regression_task <- sl3_Task$new(
          Shared_Data$new(data, force_copy = F),
          id = "id",
          time = "t",
          covariates = covariates,
          outcome = outcome,
          outcome_type = target_node_object$variable_type,
          folds = self$folds
        )

        return(regression_task)
      }


      all_vars <- unique(unlist(lapply(npsem, `[[`, "variables")))

      times <- as.vector(sapply(parent_nodes, function(node) node$time))
      parent_covariates <- as.vector(sapply(parent_nodes, function(node) node$variables))

      # Note that those with missing rows will be included in outcome_data.
      # There value will be set to last measured value.
      outcome_data <- self$get_tmle_node(target_node, format = TRUE, include_id = T, include_time = (time == "pooled"), force_time_value = force_time_value, expand = expand)

      past_data <- past_data[t <= time & id %in% outcome_data$id,]



      if(length(parent_covariates) != 0 | !is.null(unlist(target_node_object$summary_functions))){
        summary_measures <- target_node_object$summary_functions
        all_covariate_data <- lapply(summary_measures, function(fun){
          return(fun$summarize(past_data, time))}
        )
        all_covariate_data <- all_covariate_data %>% purrr::reduce(merge, by = "id")
        covariates <- setdiff(colnames(all_covariate_data), "id")
        if("t" %in% colnames(all_covariate_data))  all_covariate_data$t <- NULL

      } else {
        past_data <- data.table(rep(NA, nrow(outcome_data)))
        past_data$id <- outcome_data$id
        covariates <- c()
      }

      nodes <- self$nodes
      node_data <- self$get_data(, unlist(nodes))

      #TODO since number of time rows vary per person, only time-indepdent nodes make sense
      # Keep only node_data for each individual at the time of this tmle node
      node_data <- node_data[node_data$id %in% outcome_data$id, last(.SD), by = id]

      node_data$t = time
      nodes$outcome <- outcome
      nodes$covariates <- covariates
      nodes$time <- "t"
      # Merge all data by id. Full join so all rows are kept and NA's added.
      # NA rows will necessarily be dropped if drop_censored = T.
      # Otherwise only the outcome variables should be NA for unobserved/unmonitored outcomes.
      # If one is interested in getting predictions for people not observed at this time
      # then one must add a row with the desired outcome.

      regression_data <-  list(all_covariate_data, outcome_data, node_data) %>% reduce(merge, "id")
      regression_data$t = time

      regression_data <- Shared_Data$new(regression_data, force_copy = F)

      regression_task <- sl3_Task$new(
        regression_data,
        nodes = nodes,
        outcome_type = target_node_object$variable_type,
        folds = self$folds
      )
      if(!is.numeric(force_time_value)){
        assign(cache_key, regression_task, private$.node_cache)
      }

      return(regression_task)
    },
    generate_counterfactual_task = function(uuid, new_data, force_at_risk = NULL, through_data=F) {
      #If single value then
      if(!through_data){

        if(!("t" %in% colnames(new_data))){

          #If not in id/t format. Then convert to id/t format and handle node names
          node <-  setdiff(colnames(new_data), c("id", "t"))
          dat <- self$get_tmle_node(node, include_time = T, include_id = T)
          node <-  setdiff(colnames(new_data), c("id", "t"))
          node_vars <- sapply(
            node,
            function(node_name) {
              self$npsem[[node_name]]$variables
            }
          )

          set(dat, , unlist(node_vars, use.names = F), new_data[,node,with=F])

          #setnames(dat, node, node_vars)
          new_data <- dat
        }


        data <- data.table::copy(self$get_data())
        node <-  setdiff(colnames(new_data), c("id", "t"))
        has_row <- which(unlist(data[.(new_data$id, new_data$t), !is.na(node[[1]]), with = F], use.names = F))
        append_row_data <- new_data[-has_row]
        alter_row_data <- new_data[has_row]
        data[.(alter_row_data$id, alter_row_data$t), .(node) :=  alter_row_data[, .(node)]]
        if(nrow(append_row_data) > 0){
          data <- rbind(data, append_row_data, fill = T)
          setkey(data, id, t)
          setnafill(data, "locf")

        }

          new_task <- self$clone()


          new_task$initialize(
            data, self$npsem,
            folds = self$folds,
            row_index = self$row_index,
            t = "t",
            id = "id",
            nodes = self$nodes,
            force_at_risk = ifelse(is.null(force_at_risk), self$force_at_risk, force_at_risk),
            extra_summary_measure_columns = private$.non_data_columns,
            thin = self$thin
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
        force_at_risk = force_at_risk,
        thin = self$thin
      )

      return(new_task)
    }

  )
)
