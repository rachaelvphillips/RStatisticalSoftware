#' Fitting Intercept Models
#'
#' This learner provides fitting procedures for intercept models. Such models
#' predict the outcome variable simply as the mean of the outcome vector.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom assertthat assert_that is.count is.flag
#'
#' @export
#'
#' @keywords data
#'
#' @return Learner object with methods for training and prediction. See
#'  \code{\link{Lrnr_base}} for documentation on learners.
#'
#' @format \code{\link{R6Class}} object.
#'
#' @family Learners
#'
#' @section Parameters:
#' \describe{
#'   \item{\code{...}}{Not used.}
#' }
#'
#' @template common_parameters
#
Lrnr_mean <- R6Class(
  classname = "Lrnr_mean", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(marker_node, thresholds...) {
      params <- list(...)
      super$initialize(params = params, ...)
    },

    print = function() {
      print(self$name)
    }
  ),

  private = list(
    .properties = c("continuous", "binomial", "categorical", "weights", "offset"),

    .train = function(task) {



      return(fit_object)
    },

    .predict = function(task = NULL) {
      predictions <- rep(private$.fit_object$mean, task$nrow)

      if (self$fit_object$training_offset) {
        offset <- task$offset_transformed(NULL, for_prediction = TRUE)
        predictions <- predictions + offset
      }

      predictions <- as.matrix(predictions, ncol = 1)
      return(predictions)
    }
  )
)
