

#' @import delayed
#' @import sl3
#' @import origami
#' @export
Lrnr_subsemble_metalearner <- R6::R6Class(
  classname = "Lrnr_subsemble_metalearner", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(meta_lrnr = make_learner(Lrnr_nnls), num_groups,...) {

      params <- args_to_list()
      params$num_groups <- num_groups
      private$.name="Subsemble_metalearner"
      super$initialize(params = params, ...)
    },
    print = function() {
      print(self$name)
      print(self$fits)
    }
  ),

  active = list(
    fits = function() {
      fit_object <- private$.fit_object
      if (!is.null(fit_object)) {
        data.table::data.table(lrnrs = fit_object$lrnrs, weights = fit_object$x)
      } else {
        data.table::data.table(lrnrs = character(0), weights = numeric(0))
      }
    }
  ),

  private = list(
    .properties = c("subsemble"),
    .learner_names = NULL,
    .learners = NULL,
    .train_sublearners = function(task) {

      num_groups = self$params$num_groups


      delayed_levels <- function(task, num_groups) {
        if(inherits(task, "delayed")) {
          task <- task$compute()
        }
        ncol_X <- ncol(task$X)
        len_pred <- ncol_X/num_groups
        if(len_pred%%1 != 0) {
          stop("Number of predictions is not a multiple of the number of groups.")
        }
        ids <- rep(1:len_pred, num_groups)
        levels <- split(1:ncol_X, ids)
        return(levels)
      }

      levels <- delayed_fun(delayed_levels)(task, num_groups)

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

        subsemble_meta_task <- function(task, levels){
          if(inherits(levels, "delayed")) {
            levels <- levels$compute()
          }
          level <- levels[[i]]
          assert_that(ncol(task$X)==ncol_X)
          return(task$next_in_chain(covariates = colnames(task$X)[level]))
        }

        delayed_task <- delayed_fun(subsemble_meta_task, levels)(task)
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
    },

    .train = function(task, trained_sublearners) {

      is_error <- sapply(trained_sublearners, function(result) {
        inherits(result, "error") || inherits(result, "try-error")
      })
      learner_errors <- trained_sublearners[is_error]
      errored_learners <- self$params$learners[is_error]

      for (i in seq_along(errored_learners)) {
        message <- learner_errors[[i]]
        learner <- errored_learners[[i]]
        warning(sprintf(
          "%s failed with message: %s. It will be removed from the stack",
          learner$name, message
        ))
      }
      if (all(is_error)) {
        stop("All learners in stack have failed")
      }
      learner_names <- private$.learner_names
      learner_names <- sapply(1:length(learner_names), function(i){
        lrnr = trained_sublearners[[i]]
        return(paste0(learner_names[[i]],"_", "base->",stringr::str_extract(colnames(task$X)[i], "Lrnr_[^_]*"),"_", stringr::str_remove(colnames(task$X)[i], ".*_strata_[0-9]+_" )))
      })
      private$.learner_names <- learner_names

      learner_names <- learner_names[!is_error]
      names(trained_sublearners) <- learner_names

      fit_object <- list(
        learner_fits = trained_sublearners,
        learner_errors = learner_errors, is_error = is_error
      )
      return(fit_object)
    },

    .predict = function(task = NULL) {



      is_error <- private$.fit_object$is_error
      learner_fits <- private$.fit_object$learner_fits[!is_error]
      learners <- self$params$learners[!is_error]
      learner_names <- private$.learner_names[!is_error]
      n_to_pred <- task$nrow
      n_learners <- length(learner_names)


      current_fit <- learner_fits[[1]]
      current_preds <- current_fit$base_predict(task)
      current_names <- learner_names[1]
      n_to_pred <- safe_dim(current_preds)[2]
      if(is.na(n_to_pred)) n_to_pred <- 1
      num_groups <- self$params$num_groups
      pred_ids <- rep(1:n_to_pred, num_groups)


      ## Cannot use := to add columns to a null data.table (no columns),
      ## hence we have to first seed an initial column, then delete it later
      learner_preds <- data.table::data.table(
        current_preds = current_preds
      )

      if (!is.na(safe_dim(current_preds)[2]) &&
          safe_dim(current_preds)[2] > 1) {
        current_names <- paste0(current_names, "_", names(current_preds))
        stopifnot(length(current_names) == safe_dim(current_preds)[2])
      }

      setnames(learner_preds, names(learner_preds), current_names)

      for (i in seq_along(learner_fits)) {
        current_fit <- learner_fits[[i]]
        if (i > 1) {
          current_preds <- current_fit$base_predict(task)
        }

        current_names <- learner_names[i]
        if (!is.na(safe_dim(current_preds)[2]) &&
            safe_dim(current_preds)[2] > 1) {
          current_names <- paste0(learner_names[i], "_", names(current_preds))
          stopifnot(length(current_names) == safe_dim(current_preds)[2])
        }
        set(learner_preds, j = current_names, value = current_preds)
        invisible(NULL)
      }



      return(learner_preds)
    },
    .chain = function(task) {
      predictions <- self$predict(task)
      predictions <- as.data.table(predictions)
      # Add predictions as new columns
      task <- task$revere_fold_task("full")
      new_col_names <- task$add_columns(predictions, self$fit_uuid)
      # new_covariates = union(names(predictions),task$nodes$covariates)
      return(task$next_in_chain(
        covariates = names(predictions),
        column_names = new_col_names
      ))
    },
    .required_packages = c("nnls")
  ),
)
