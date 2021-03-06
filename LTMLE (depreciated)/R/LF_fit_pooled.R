#' Likelihood Factor Estimated from Data using sl3.
#'
#' Uses an \code{sl3} learner to estimate a likelihood factor from data.
#' Similar to LF_fit object but allows vector-valued names (corresponding to multiple nodes sharing the same regression_task variables but at different times)
#' This learner will pool the regression tasks across time. This depends on ltmle3_Task which has the pooled task functionality.
#' Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @importFrom AR AR.Sim
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
#'TODO If outcome is multivariate then this assumes learner will handle the multi outcome and return predictions as a matrix

#' @export

LF_fit_pooled <- R6Class(
  classname = "LF_fit_pooled",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, learner, is_time_variant = FALSE, ..., type = "density") {
      super$initialize(name, ..., type = type)
      private$.learner <- learner
      # TODO: add parameter is_time_variant
      private$.is_time_variant <- is_time_variant
    },
    delayed_train = function(tmle_task) {
      # just return prefit learner if that's what we have
      # otherwise, make a delayed fit and return that


      if (self$learner$is_trained) {
        return(self$learner)
      }

      outcome_node <- self$name

      # fit scaled task for bounded continuous
      learner_task <- tmle_task$get_regression_task(outcome_node, scale = TRUE, drop_censored=TRUE,
                                                    is_time_variant = self$is_time_variant)

      learner_fit <- delayed_learner_train(self$learner, learner_task)

      return(learner_fit)
    },
    train = function(tmle_task, learner_fit) {
      #super$train(tmle_task)
      #get first node for outcome type

      tmle_nodes <- lapply(self$names, function(node) tmle_task$npsem[[node]])
      private$.variable_type <- lapply(tmle_nodes, function(node) node$variable_type)
      private$.training_task <- tmle_task
      private$.learner <- learner_fit
    },
    shape_predictions = function(learner_task, preds){
      # Reshapes predictions if this is a pooled_task
      # returns the predictions as n x len(self$names) data.table.
      if(length(self$name)==1){
        preds <- data.table(preds)
        setnames(preds, self$name)
        return(preds)
      }
      ids <- learner_task$get_node("id")
      preds_list  <- lapply(unique(ids), function(id){
        # Returns the vector of predictions for a single person
        preds[which(ids == id)]
      })
      # Stack preds so each row corresponds to the prediction for a single person
      if(length(preds_list) == 1){
        preds_stacked <- data.table(preds_list[[1]])
      } else {
        preds_stacked <- data.table(do.call(rbind, preds_list))
      }

      setnames(preds_stacked, self$name)
      return(preds_stacked)
    },
    get_mean = function(tmle_task, fold_number, expand = T) {
      # TODO: prediction is made on all data, so is_time_variant is set to TRUE
      learner_task <- tmle_task$get_regression_task(self$name, is_time_variant = TRUE, expand = T)
      learner <- self$learner

      preds <- learner$predict_fold(learner_task, fold_number)

      # unscale preds (to handle bounded continuous)
      preds_unscaled <- tmle_task$unscale(preds, self$name)
      preds <- data.table(preds_unscaled)
      preds$id <- learner_task$data$id
      preds$t <- learner_task$data$t
      setnames(preds, c(learner_task$nodes$outcome, "id", "t"))
      if(to_wide){
        preds <- reshape(preds, idvar = "id", timevar = "t", direction = "wide")
        setnames(preds, c("id", self$name))
      }
      if(drop_id & "id" %in% colnames(preds)) preds$id <- NULL
      if(drop_time & "t" %in% colnames(preds)) preds$t <- NULL
      if(drop & ncol(preds) == 1) preds <- unlist(preds, use.names = F)
      return(preds)
    },
    get_density = function(tmle_task, fold_number, check_at_risk = T, to_wide = T, drop_id = F, drop_time = F, drop = T, expand = T) {
      # TODO: prediction is made on all data, so is_time_variant is set to TRUE

      learner_task <- tmle_task$get_regression_task(self$name, is_time_variant = TRUE, expand = expand)

      learner <- self$learner

      preds <- learner$predict_fold(learner_task, fold_number)

      outcome_type <- self$learner$training_task$outcome_type
      Y <- learner_task$Y


      observed <- outcome_type$format(learner_task$Y)
      data <-  learner_task$get_data()


      if (outcome_type$type == "binomial") {
        likelihood <- ifelse(observed == 1, preds, 1 - preds)
      } else if (outcome_type$type == "categorical") {
        unpacked <- sl3::unpack_predictions(as.vector(preds))
        index_mat <- cbind(seq_along(observed), observed)
        likelihood <- unpacked[index_mat]
      } else if (outcome_type$type == "continuous") {
        likelihood <- preds
      } else {
        stop(sprintf("unsupported outcome_type: %s", outcome_type$type))
      }

      if(check_at_risk & "at_risk" %in% colnames(data) ) {
        # By setting check_at_risk = F then one can obtain the counterfactual predictions
        # conditioned on everyone being at risk.
        names_last_val <- paste0("last_val_", learner_task$nodes$outcome)
        # TODO Support multivariate outcome
        assertthat::assert_that(all(names_last_val %in% colnames(data)), msg = "If at_risk is a column then last_val must be as well.")
        not_at_risk <- which(data$at_risk == 0)
        if(length(not_at_risk)>0){
          lst_val <- data[not_at_risk, names_last_val, with = F]
          #If not at risk then equal to last value with prob 1
          likelihood[not_at_risk] <- as.numeric(observed[not_at_risk] == lst_val)
        }

      }
      likelihood <- data.table(likelihood)
      likelihood$id <- learner_task$data$id
      likelihood$t <- learner_task$data$t

      setnames(likelihood, c(learner_task$nodes$outcome, "id", "t"))
      if(to_wide){
        likelihood <- reshape(likelihood, idvar = "id", timevar = "t", direction = "wide")
        setnames(likelihood, c("id", self$name))
      }
      if(drop_id & "id" %in% colnames(likelihood)) likelihood$id <- NULL
      if(drop_time & "t" %in% colnames(likelihood)) likelihood$t <- NULL
      if(drop & ncol(likelihood) == 1) likelihood <- unlist(likelihood, use.names = F)


      return(likelihood)
    },
    get_likelihood = function(tmle_task, fold_number = "full", expand = T) {
      if (self$type == "mean") {
        values <- self$get_mean(tmle_task, fold_number, expand = expand)
      } else {
        values <- self$get_density(tmle_task, fold_number, expand = expand)
      }
      if (!is.null(self$bound)) {
        values <- bound(values, self$bound)
      }

      return(values)
    },

    sample = function(tmle_task, n_samples = NULL, fold_number = "full") {
      warning("This function has not been implemented for pooled tasks.")
      # TODO: fold
      # TODO: pooled tasks
      if (is.null(tmle_task)) {
        tmle_task <- self$training_task
      }
      if (is.null(n_samples)) {
        return(tmle_task)
      }

      learner_task <- tmle_task$get_regression_task(self$name)
      learner <- self$learner

      outcome_type <- learner$training_task$outcome_type

      if (outcome_type$type == "binomial") {
        # TODO: option to return task
        # TODO: think about how folds should be structured on resample
        # need to keep ids the same
        # probably also predict using training set fits
        preds <- learner$predict_fold(learner_task, "full")
        values <- sapply(preds, function(p) rbinom(n_samples, 1, p))
      } else if (outcome_type$type == "categorical") {
        preds <- learner$predict_fold(learner_task, "full")
        unpacked <- sl3::unpack_predictions(as.vector(preds))
        values <- apply(
          unpacked, 1,
          function(probs) {
            apply(
              rmultinom(n_samples, 1, probs) == 1, 2,
              function(onehots) outcome_type$levels[which(onehots)]
            )
          }
        )
      } else if (outcome_type$type == "continuous") {
        if ("sampling" %in% learner$properties) {
          values <- learner$sample(learner_task, n_samples, "full")
        } else {
          values <- matrix(nrow = n_samples, ncol = tmle_task$nrow)
          for (i in 1:tmle_task$nrow) {
            subject <- tmle_task[i]
            f_X <- function(a) {
              cf_data <- data.table(a)
              setnames(cf_data, names(cf_data), self$name)
              subject_a <- subject$generate_counterfactual_task(UUIDgenerate(), cf_data)

              pred <- learner$predict_fold(subject_a$get_regression_task(self$name), "full")
              likelihood <- unlist(pred)

              return(likelihood)
            }
            samples <- AR::AR.Sim(n_samples, f_X,
                              xlim = c(min(learner$training_task$Y), max(learner$training_task$Y))
            )
            values[, i] <- samples
          }
        }
      } else {
        stop(sprintf("unsupported outcome_type: %s", outcome_type$type))
      }
      values <- t(values)
      return(values)
    }
  ),
  active = list(
    learner = function() {
      return(private$.learner)
    },
    is_time_variant = function() {
      return(private$.is_time_variant)
    }
  ),
  private = list(
    .name = NULL,
    .learner = NULL,
    .is_time_variant = NULL
  )
)
