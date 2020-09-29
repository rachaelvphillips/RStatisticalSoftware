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
    initialize = function(threshold_function = function(A) {as.vector(quantile(A, seq(0.05, 0.95, length.out = 10)))}
                          , cdf_bins = 10, data_adaptive = F, cv = data_adaptive, ...) {
      super$initialize(
        threshold_function = threshold_function,
        cdf_bins = cdf_bins, data_adaptive = data_adaptive, cv = cv, ...
      )
    },
    make_tmle_task = function(data, node_list, weights = NULL) {
      npsem <- make_thresh_npsem(node_list, data_adaptive = self$options$data_adaptive)
      tmle_task <- make_thresh_task(data, npsem, weights = weights)
      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      # produce trained likelihood when likelihood_def provided
      cdf_bins <- self$options$cdf_bins
      threshold_function <- self$options$threshold_function

      if (!is.null(self$options$likelihood_override)) {
        likelihood <- self$options$likelihood_override$train(tmle_task)
      } else {
        likelihood <- make_thresh_likelihood(tmle_task, learner_list, threshold_function, cdf_bins, cv = self$options$cv, marker_learner = learner_list[["A_learned"]])
      }

      return(likelihood)
    },
    make_params = function(tmle_task, likelihood) {
      thres <- Param_thresh$new(likelihood)
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
