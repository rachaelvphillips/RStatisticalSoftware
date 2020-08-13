
#' @export
Lrnr_subsemble_hal_metalearner <- R6::R6Class(
  classname = "Lrnr_subsemble_hal_metalearner", inherit = Lrnr_subsemble_metalearner,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(num_groups, num_preds,...) {
      meta_lrnr = make_learner(Lrnr_hal9001_subsemble_collapser)
      params <- args_to_list()
      params$num_groups <- num_groups
      params$num_preds <- num_preds
      params$meta_lrnr <- meta_lrnr
      do.call(super$initialize, params)
      private$.name="Subsemble_metalearner_hal"

    }
  ),
  private = list(
    .train_sublearners = function(task) {
      len_pred <- self$params$num_preds
      num_groups = self$params$num_groups
      # generate training subtasks
      # X <- task$X
      # ncol_X = ncol(task$X)
      ncol_X = num_groups*len_pred



      ids <- rep(1:len_pred, num_groups)
      levels <- split(1:ncol_X, ids)

      meta_lrnr <- self$params$meta_lrnr
      num_groups <- self$params$num_groups

      # Create set of meta learners
      learners <- replicate(len_pred, meta_lrnr$clone())
      learner_names <- sapply(1:len_pred, function(i){

        lrnr = learners[[i]]
        return(paste0("meta_", i, "->", lrnr[["name"]]))
      })

      private$.learner_names <- learner_names
      private$.learners <- learners


      subtasks <- lapply(1:len_pred, function(i) {
        learner = learners[[i]]
        level = levels[[i]]
        subsemble_meta_task <- function(task){
          assert_that(ncol(task$X)==ncol_X)
          return(task$next_in_chain(covariates = colnames(task$X)[level])$add_to_bucket(i, "lambda_index"))
        }

        delayed_task <- delayed_fun(subsemble_meta_task)(task)
        delayed_task$name <- "subsemble_meta_task"
        if (learner$is_trained) {
          return(learner)
        } else {
          delayed_learner <- delayed_learner_train(learner, delayed_task)
          delayed_learner$expect_error <- TRUE
          return(delayed_learner)
        }
      })
      return(bundle_delayed(subtasks))
    }
  )
)
