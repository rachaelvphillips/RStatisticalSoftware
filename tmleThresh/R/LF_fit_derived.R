#' Derived Likelihood Factor Estimated from Data + Other Likelihood values, using sl3.
#'
#' Uses an \code{sl3} learner to estimate a likelihood factor from data.
#' Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Likelihood objects
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_lf(LF_fit, name, learner, ..., type = "density")}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name in the nodes specified by \code{\link{tmle3_Task}$npsem}
#'     }
#'     \item{\code{learner}}{An sl3 learner to be used to estimate the factor
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{type}}{character, either "density", for conditional density or, "mean" for conditional mean
#'     }
#'     }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{learner}}{The learner or learner fit object}
#'     }
#'
#' @export
LF_fit_derived <- R6::R6Class(
  classname = "LF_fit_derived",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, learner, base_likelihood, task_generator, cast_to_long = T, byrow = T, ..., type = "density") {
      super$initialize(name, ..., type = type)
      private$.learner <- learner
      private$.base_likelihood <- base_likelihood
      private$.task_generator <- task_generator
      private$.cast_to_long <- cast_to_long
      private$.byrow <- byrow
    },
    delayed_train = function(tmle_task) {
      # call task generator
      derived_task <- self$task_generator(tmle_task, self$base_likelihood)

      # just return prefit learner if that's what we have
      # otherwise, make a delayed fit and return that
      if (self$learner$is_trained) {
        return(self$learner)
      }

      outcome_node <- self$name
      learner_fit <- delayed_learner_train(self$learner, derived_task)
      return(learner_fit)
    },
    train = function(tmle_task, learner_fit) {
      super$train(tmle_task)
      private$.learner <- learner_fit
    },
    get_mean = function(tmle_task, fold_number, ...) {
      derived_task <- self$task_generator(tmle_task, self$base_likelihood)
      learner <- self$learner
      if(self$cast_to_long) {
        preds <- as.data.table(as.vector(learner$predict_fold(derived_task, fold_number)))
      } else {
        preds <- as.data.table(learner$predict_fold(derived_task, fold_number))

      }

      setnames(preds, self$name)

      return(preds)
    },
    get_density = function(tmle_task, fold_number, ...) {
      derived_task <- self$task_generator(tmle_task, self$base_likelihood)
      learner <- self$learner
      if(self$cast_to_long) {
        preds <- learner$predict_fold(derived_task, fold_number)

        if(self$byrow) {
          preds <- as.vector(t(preds))
        } else {
          preds <- as.vector(preds)
        }
        preds <- as.data.table(preds)
      } else {
        preds <- as.data.table(learner$predict_fold(derived_task, fold_number))

      }
      setnames(preds, self$name)

      # todo: think about derived task with other outcome types (this assumes continuous)
      return(preds)
    }
  ),
  active = list(
    learner = function() {
      return(private$.learner)
    },
    task_generator = function() {
      return(private$.task_generator)
    },
    base_likelihood = function() {
      return(private$.base_likelihood)
    },
    cast_to_long = function(){
      return(private$.cast_to_long)
    },
    byrow = function(){
      return(private$.byrow)
    }
  ),
  private = list(
    .name = NULL,
    .learner = NULL,
    .base_likelihood = NULL,
    .task_generator = NULL,
    .cast_to_long = NULL,
    .byrow = NULL
  )
)
