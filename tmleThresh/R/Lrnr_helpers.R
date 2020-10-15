# Some helper



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
    initialize = function(name, learner, base_likelihood, task_generator, cast_to_long = T, byrow = F, ..., type = "density") {
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
      preds <- sl3::unpack_predictions(learner$predict_fold(derived_task, fold_number))
      if(self$cast_to_long) {
        if(self$byrow) {
          preds <- as.vector(t(preds))
        } else {
          preds <- as.vector(preds)
        }
      }
      preds <- as.data.table(preds)
      setnames(preds, self$name)

      return(preds)
    },
    get_density = function(tmle_task, fold_number, ...) {
      derived_task <- self$task_generator(tmle_task, self$base_likelihood)
      learner <- self$learner
      if(self$cast_to_long) {
        preds <- learner$predict_fold(derived_task, fold_number)
        preds <- sl3::unpack_predictions(preds)
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



Lrnr_wrapper <- R6Class(
  classname = "Lrnr_wrapper",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(ncol, pack = T, ...) {

      params <- list(ncol = ncol, pack = pack, ...)
      super$initialize(params = params, ...)
    }
  ),

  active = list(
    name = function() {

      name <- "Wrapper"
      return(name)
    }
  ),
  private = list(

    .train = function(task) {
      fit_object = list()
      return(fit_object)
    },

    .predict = function(task) {
      ncol <- self$params$ncol

      X <- ((task$X))

      out <- as.data.table(apply(X, 2, function(v) {
        if(self$params$pack){
          v <- matrix(v, ncol = ncol)
          predictions <- pack_predictions(v)
        } else {
          predictions <- v
        }
        return((predictions))
      }))


      return(out)
    },
    .required_packages = NULL
  )
)

Lrnr_chainer <- R6Class(
  classname = "Lrnr_chainer",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(cutoffs, strata_variable, ...) {

      params <- list(strata_variable = strata_variable, cutoffs = cutoffs,...)
      super$initialize(params = params, ...)
    }
  ),

  active = list(
    name = function() {

      name <- "Wrapper"
      return(name)
    }
  ),
  private = list(

    .train = function(task) {
      fit_object = list()
      return(fit_object)
    },

    .chain = function(task) {

      args <- self$params
      cutoffs <- args$cutoffs
      strata_variable <- args$strata_variable

      if(inherits(task, "delayed")) {
        task <- task$compute()
      }
      if(inherits(task, "sl3_revere_Task")) {
        print(stop("this shouldnt be"))
        new_generator <- function(task, fold_number) {
          task <- task$revere_fold_task(fold_number)
          data <- task$data
          cutoffs <- args$cutoffs
          data_list <- list()

          for(cutoff in cutoffs) {
            Xcopy <- copy(data)
            Xcopy$bin <- cutoff
            Xcopy$Ind <- as.numeric(Xcopy[[strata_variable]] >= cutoff)
            Xcopy[[strata_variable]] <- NULL
            data_list[[as.character(cutoff)]] <- Xcopy
          }
          data <- rbindlist(data_list)

          nodes <- task$nodes

          nodes$covariates <- union(setdiff(task$nodes$covariates, strata_variable), c("Ind", "bin"))
          task <- sl3_Task$new(data, nodes = nodes)

          return(task)
        }
        task <- sl3_revere_Task$new(new_generator, task)


      } else {
        args <- self$params
        strata_variable <- args$strata_variable

        data <- task$data
        cutoffs <- args$cutoffs
        data_list <- list()

        for(cutoff in cutoffs) {
          Xcopy <- copy(data)
          Xcopy$bin <- cutoff
          Xcopy$Ind <- as.numeric(Xcopy[[strata_variable]] >= cutoff)
          Xcopy[[strata_variable]] <- NULL
          data_list[[as.character(cutoff)]] <- Xcopy
        }

        data <- rbindlist(data_list)

        nodes <- task$nodes

        nodes$covariates <- union(setdiff(task$nodes$covariates, strata_variable), c("Ind", "bin"))
        task <- sl3_Task$new(data, nodes = nodes)
        print(task)
        print(table(task$Y))
      }

      return(task)
    },
    .required_packages = NULL
  )
)

