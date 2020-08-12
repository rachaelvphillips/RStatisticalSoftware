
make_learner_subsemble <- function(lrnr, task, k = 5, meta_lrnr = make_learner(Lrnr_nnls)){
  n <- task$nrow
  folds <- compute_subsemble_folds(n,k)
  strata_ids <- attr(folds , 'strata_ids')
  lrnr_sub <- make_learner(Lrnr_subsemble, lrnr = lrnr, strata_ids = strata_ids)
  meta_lrnr <- make_learner(Lrnr_subsemble_meta, meta_lrnr = meta_lrnr, strata_ids = strata_ids)
  final_lrnr <- make_learner(Pipeline, lrnr_sub, meta_lrnr)
  new_task <- task$next_in_chain(folds = folds)
  return(list(subsemble_lrnr = final_lrnr, cv_task = new_task))
}
compute_subsemble_folds <- function(n,k){
  strata = (sample(1:n))%%k +1
  folds = origami::make_folds(n, strata_ids =strata )
  attr(folds , 'strata_ids') = strata
  return(folds)
}
Lrnr_subsemble <- R6::R6Class(
  classname = "Lrnr_subsemble",
  inherit = Lrnr_base,
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(lrnr, strata_ids) {

      levels = split(1:length(strata_ids), strata_ids)
      num_groups = length(levels)
      learners <- replicate(num_groups, lrnr$clone())

      # catch learner names and make unique if there's repetition
      learner_names <- sapply(1:num_groups, function(grp){
        lrnr = learners[[grp]]
        paste0(lrnr[["name"]], "_strata_", grp)
        })
      if (any(duplicated(learner_names))) {
        learner_names <- make.unique(learner_names, sep = "_")
      }
      private$.learner_names <- learner_names
      params <- list(learners = learners, strata_ids = strata_ids, levels = levels, num_groups = num_groups )


      super$initialize(params = params)
    },
    print = function() {
      if (is.null(private$.fit_object)) {
        print(private$.learner_names)
        # lapply(self$params$learners, print)
      } else {
        print(private$.learner_names)
        # lapply(private$.fit_object, print)
      }
    },
    update_errors = function(is_error) {
      private$.fit_object$is_error <- is_error
    }
  ),

  active = list(
    name = function() {
      # learners = self$params$learners
      # learner_names = sapply(learners, function(learner) learner$name)
      # name = paste(learner_names, collapse="x")
      name <- "Stack"
      return(name)
    },
    learner_fits = function() {
      result <- self$fit_object$learner_fits
      return(result)
    }
  ),

  private = list(
    # modified names of learners
    .learner_names = NULL,

    .train_sublearners = function(task) {
      # generate training subtasks
      learners <- self$params$learners
      levels <- self$params$levels
      subtasks <- lapply(1:length(learners), function(i) {
        learner = learners[[i]]
        level = levels[[i]]
        if (learner$is_trained) {
          return(learner)
        } else {
          delayed_learner <- delayed_learner_train(learner, task$subset_task(row_index = level))
          delayed_learner$expect_error <- TRUE
          return(delayed_learner)
        }
      })
      return(bundle_delayed(subtasks))
    },
    .train = function(task, trained_sublearners) {
      # check fits for errors
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

      learner_names <- private$.learner_names[!is_error]
      names(trained_sublearners) <- learner_names

      fit_object <- list(
        learner_fits = trained_sublearners,
        learner_errors = learner_errors, is_error = is_error
      )
      return(fit_object)
    },
    .predict = function(task) {
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
    }
  )
)
