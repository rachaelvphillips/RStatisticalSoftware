

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
#' @export
LF_collapsed <- R6Class(
  classname = "LF_collapsed",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, factor_list , ..., type = "density") {
      super$initialize(name, ..., type = type)
      private$.factor_list<- factor_list

    },
    get_mean = function(tmle_task, fold_number) {
      stop("Error: Conditional means cannot be collapsed.")
      learner_task <- tmle_task$get_regression_task(self$name, scale = FALSE)
      preds <- self$mean_fun(learner_task)

      return(preds)
    },
    get_density = function(tmle_task, fold_number) {
      factors <- self$factor_list
      preds <- lapply(factors, function(factor) {
        factor$get_density(tmle_task, fold_number)
        })
      preds <- do.call(cbind, preds)
      if(is.data.table(preds)){
      preds <- unlist(apply(preds, 1, prod))
      }
      preds <- data.table(preds)
      setnames(preds, self$name)
      return(preds)
    }
  ),
  active = list(
    factor_list = function() {
      return(private$.factor_list)
    }
  ),
  private = list(
    .name = NULL,
    .factor_list = NULL
  )
)
