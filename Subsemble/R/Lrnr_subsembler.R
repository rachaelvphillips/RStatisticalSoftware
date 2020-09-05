


Lrnr_subsemble <-



Lrnr_subsemble_base <- R6::R6Class(
  classname = "Lrnr_subsemble_base",
  inherit = Lrnr_base,
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(learners, k) {
      if(length(learners) == k){
        learners <- learners
      } else {
        learners <- lapply(1:k, function(i) learners[[1]]$clone())
      }

      # catch learner names and make unique if there's repetition
      learner_names <- sapply(1:k, function(grp){
        lrnr <- learners[[grp]]
        paste0(lrnr[["name"]], "_strata_", grp)
      })
      if (any(duplicated(learner_names))) {
        learner_names <- make.unique(learner_names, sep = "_")
      }
      private$.learner_names <- learner_names
      params <- list(learners = learners, k = k )


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
      k <- self$params$k
      compute_levels <- function(task){
        ids <- (task$id)
        uniq_ids <- unique(ids)
        strata_ids <- split(uniq_ids, sample(1:length(uniq_ids))%%k +1)
        row_index_list <- lapply(strata_ids, function(grp) {
          which(ids %in% grp)
        })
        return(row_index_list)
      }
      row_index_list <- delayed_fun(compute_levels)(task)
      subtasks <- lapply(1:k, function(i) {
        learner = learners[[i]]

        subsemble_task <- function(task, row_index_list){
          if (inherits(row_index_list, "Delayed")) {
            row_index_list <- row_index_list$compute()
          }

          return(task$subset_task(row_index = row_index_list[[i]]))
        }

        delayed_task <- delayed_fun(subsemble_task)(task, row_index_list)
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




Lrnr_subsemble_metalearner <- R6::R6Class(
  classname = "Lrnr_subsemble_metalearner", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(meta_lrnr = make_learner(Lrnr_nnls), num_groups, num_preds = 1,...) {

      params <- args_to_list()
      params$num_groups <- num_groups
      params$num_preds <- num_preds
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
          return(task$next_in_chain(covariates = colnames(task$X)[level]))
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
