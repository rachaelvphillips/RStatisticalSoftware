
#' @export
Lrnr_subsemble_hal <- R6::R6Class(
  classname = "Lrnr_subsemble_hal",
  portable = TRUE, class = TRUE, inherit = Lrnr_subsemble,
  public = list(



  ),
  active = list(

  ),
  private = list(
    .chain = function(task) {
      predictions <- self$predict(task)


      predictions <- as.data.table(predictions)
      # Add predictions as new columns
      task <- task$revere_fold_task("full")
      new_col_names <- task$add_columns(predictions, self$fit_uuid)
      # new_covariates = union(names(predictions),task$nodes$covariates)
      next_task = task$next_in_chain(
        covariates = names(predictions),
        column_names = new_col_names
      )
      fit_objects = lapply(self$fit_object$learner_fits, function(item){ tmp = item$fit_object$fit
      if(is.null(tmp)){
        tmp = item$fit_object
      }
      return(tmp)})
      next_task = next_task$add_to_bucket(fit_objects, "hal_fits")
      next_task = next_task$add_to_bucket(task, "task")
      return(next_task)

    }
  )
)

