
# Generated data-adaptive covariate "A"
learner_marker_task_generator <- function(learned_marker_node = "A", learned_marker_var = "A", marker_node = "A", node = "Y", data_adaptive = T, is_time_variant = F) {
  outA <- function(tmle_task, likelihood, for_training = F) {
    generator <- function(task, fold_number) {
      if(!data_adaptive) {
        preds <- tmle_task$get_tmle_node("A")
      } else {
        preds <- likelihood$get_likelihood(tmle_task,  learned_marker_node, fold_number)
      }

      new_data <- as.data.table(preds)
      setnames(new_data, learned_marker_node)

      column_names <- task$add_columns(new_data)

      new_task <- task$next_in_chain(outcome = learned_marker_node, column_names = column_names)
      data <- new_task$data

      return(new_task)
    }

    sl3_revere_Task$new(generator, tmle_task$get_regression_task(node))
  }

  outR <- function(tmle_task, likelihood, for_training = F) {
    generator <- function(task, fold_number) {
      if(!data_adaptive) {
        preds <- tmle_task$get_tmle_node("A")
        id <- tmle_task$get_tmle_node("A", include_id = T)$id
      } else {
        id <- tmle_task$get_tmle_node("A", include_id = T)$id
        preds <- likelihood$get_likelihood(tmle_task,  learned_marker_node, fold_number)
      }
      new_data <- as.data.table(preds)

      setnames(new_data, learned_marker_var)
      new_data$id <- id
      id_task <- data.table(id =task$data$id)
      new_data <- merge(id_task, new_data, by = "id", all.x = T)
      new_data$id <- NULL

      column_names <- task$add_columns(new_data)
      # Removes the covariates corresponding with A
      covariates <- sort(union(setdiff(task$nodes$covariates, tmle_task$npsem[[marker_node]]$variables), learned_marker_var))
      #Assumes that no extra covariates have been added.
      new_task <- task$next_in_chain(covariates = covariates, column_names = column_names)

      return(new_task)
    }

    sl3_revere_Task$new(generator, tmle_task$get_regression_task(node, drop_censored = for_training, expand = !for_training, is_time_variant = is_time_variant))
  }
  if(node != "A") {
    return(outR)
  } else if(node == "A") {
    return(outA)
  }
}

