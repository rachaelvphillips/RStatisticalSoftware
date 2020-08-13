
#' @import delayed
#' @import sl3
#' @import origami
#' @import data.table
#'
#' @export
make_learner_subsemble <- function(lrnr, task, k = 5, num_preds, meta_lrnr = make_learner(Lrnr_nnls)){
  n <- task$nrow
  folds <- compute_subsemble_folds(n,k)
  strata_ids <- attr(folds , 'strata_ids')
  lrnr_sub <- make_learner(Lrnr_subsemble, lrnr = lrnr, strata_ids = strata_ids)
  meta_lrnr <- make_learner(Lrnr_subsemble_metalearner, meta_lrnr = meta_lrnr, num_groups = k, num_preds = num_preds)
  final_lrnr <- make_learner(Pipeline, lrnr_sub, meta_lrnr)
  new_task <- task$next_in_chain(folds = folds)
  return(list(subsemble_lrnr = final_lrnr, cv_task = new_task))
}
#' @export
compute_subsemble_folds <- function(n,k){
  strata = (sample(1:n))%%k +1
  folds = origami::make_folds(n, strata_ids =strata )
  attr(folds , 'strata_ids') <- strata
  return(folds)
}
#' @export
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
      name <- "Subsemble"
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
    .properties = c("subsemble", "subsemble_metalearner"),
    .train_sublearners = function(task) {
      # generate training subtasks
      learners <- self$params$learners
      levels <- self$params$levels
      subtasks <- lapply(1:length(learners), function(i) {
        learner = learners[[i]]
        level = levels[[i]]
        subsemble_task <- function(task){
          if(!is.null(task$row_index)){
            sub_rows = as.vector(na.omit(match(level, task$row_index)))
          }
          else{
            sub_rows = level
          }

          return(task$subset_task(row_index = sub_rows))
        }

        delayed_task <- delayed_fun(subsemble_task)(task)
        delayed_task$name <- "subsemble_task"
        if (learner$is_trained) {
          return(learner)
        } else {
          delayed_learner <- delayed_learner_train(learner,delayed_task)
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




      ## Cannot use := to add columns to a null data.table (no columns),
      ## hence we have to first seed an initial column, then delete it later
      learner_preds <- data.table::data.table(
        ..delete = rep(NA, n_to_pred)
      )



      for (i in seq_along(learner_fits)) {
        current_fit <- learner_fits[[i]]

        current_preds <- rep(NA, n_to_pred)
        try({
          current_preds <- current_fit$base_predict(task)
        })


        pred_names <- private$.name_preds(learner_names[i], current_preds)


        set(learner_preds, j = pred_names, value = as.data.table(current_preds))
        invisible(NULL)
      }
      learner_preds$..delete <- NULL


      return(learner_preds)
    },
    .name_preds = function(learner_name, preds) {
      current_names <- learner_name
      if (!is.na(safe_dim(preds)[2]) &&
          safe_dim(preds)[2] > 1) {
        prednames <- colnames(preds)
        if (is.null(prednames)) {
          prednames <- sprintf("col_%03d", seq_len(ncol(preds)))
        }
        current_names <- paste0(current_names, "_", prednames)
        stopifnot(length(current_names) == safe_dim(preds)[2])
      }

      return(current_names)
    }
  )
)
