#' Defines a TML Estimator (except for the data)
#'
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_Threshold <- R6Class(
  classname = "tmle3_Spec_Threshold",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(method = c("Psi_W", "cond_mean"), threshold_function = function(A) {as.vector(quantile(A, seq(0.1, 0.9, length.out = 10)))}
                          , cdf_bins = 10, data_adaptive = F, cv = data_adaptive, ...) {
      method <- match.arg(method)
      super$initialize(
        method = method, weights = weights, threshold_function = threshold_function,
        cdf_bins = cdf_bins, data_adaptive = data_adaptive, cv = cv, ...
      )
    },
    make_tmle_task = function(data, node_list) {
      npsem <- make_thresh_npsem(node_list, data_adaptive = self$options$data_adaptive)
      tmle_task <- make_thresh_task(data, npsem, weights = node_list$weights)
      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      # produce trained likelihood when likelihood_def provided
      cdf_bins <- self$options$cdf_bins
      threshold_function <- self$options$threshold_function

      if (!is.null(self$options$likelihood_override)) {
        likelihood <- self$options$likelihood_override$train(tmle_task)
      } else {
        if(self$options$method == "Psi_W") {
          likelihood <- make_thresh_likelihood(tmle_task, learner_list, threshold_function, cdf_bins, cv = self$options$cv, marker_learner = learner_list[["A_learned"]])

        } else {
          likelihood <- make_thresh_likelihood_eff(tmle_task, learner_list)
        }
      }

      return(likelihood)
    },
    make_targeted_likelihood = function(likelihood, updater) {
      if(self$options$method == "cond_mean") {
        submodel_type_by_node = list("Y" = "logistic", "A" = "EIC")
      } else {
        submodel_type_by_node <- "logistic"
      }
      print(submodel_type_by_node)
      targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater, submodel_type_by_node = submodel_type_by_node)
      return(targeted_likelihood)
    },
    make_updater = function(lower_bound = 0, one_dimensional = F, constrain_step = F, delta_epsilon_Y = 1e-2, delta_epsilon_A = 0.05,...) {
      if(self$options$method == "cond_mean") {
        constrain_step = T
        one_dimensional = T
        delta_epsilon = list("Y" = delta_epsilon_Y, "A" = function(x) {
          res = 0.1/max(abs(x))

          res <- min(res, delta_epsilon_A)
          return(res)
        })
      } else {

      }
      updater <- tmle3_Update$new(lower_bound = lower_bound, one_dimensional = one_dimensional, constrain_step = constrain_step, delta_epsilon = delta_epsilon, ...)
    },
    make_params = function(tmle_task, likelihood, ...) {
      thresholds <- self$options$threshold_function
      thresholds <- thresholds(tmle_task$get_tmle_node("A"))

      if(self$options$method == "Psi_W") {
        thres <- Param_thresh$new(likelihood, ...)

      } else {
        thres <- Param_thresh_eff$new(likelihood, thresholds = thresholds, ...)
      }
      tmle_params <- list(thres)
      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

#' All Treatment Specific Means
#'
#' O=(W,A,Y)
#' W=Covariates
#' A=Treatment (binary or categorical)
#' Y=Outcome (binary or bounded continuous)
#' @importFrom sl3 make_learner Lrnr_mean
#' @param treatment_level the level of A that corresponds to treatment
#' @param control_level the level of A that corresponds to a control or reference level
#' @export
tmle_Threshold <- function(threshold_function = function(A) {as.vector(quantile(A, seq(0.05, 0.95, length.out = 10)))}
                           , cdf_bins = 10, data_adaptive = F, cv = data_adaptive, ...) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_Threshold$new(threshold_values, num_thresholds, cdf_bins, ...)
}
