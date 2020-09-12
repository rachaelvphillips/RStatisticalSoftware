
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @export
Param_RR <- R6Class(
  classname = "Param_RR",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, outcome_node = "R", ...) {
      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)
      cf_task <- observed_likelihood$training_task
      # cf_data where everyone is the maximum level of of the node so that they are above threshhold in every group

      cf_task1 <- cf_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=1))
      cf_task0 <- cf_task$generate_counterfactual_task(UUIDgenerate(), data.table(A=0))


      observed_likelihood$get_likelihood(cf_task1, "R")
      observed_likelihood$get_likelihood(cf_task0, "R")

      private$.risk_function <- LRR_risk$new(observed_likelihood, ...)
      # get initial LRR
      #risk <- self$risk_function
      #risk$train(tmle_task, fold_number)


    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", for_fitting = F, refit = F) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      risk <- self$risk_function
      if(refit) {
        risk$train(tmle_task, fold_number)
      }
      LRR <- risk$predict(tmle_task)
      A <- tmle_task$get_tmle_node("A")
      g <- self$observed_likelihood$get_likelihood(tmle_task, "A")
      g <- bound(g, c(0.005, 1))
      H1 <- A/g * LRR
      H2 <- (1/g) * log(1 + exp(LRR))
      HA <- cbind(H1, H2)

      return(list(R = HA))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      lik <- self$observed_likelihood
      ER <- lik$get_likelihood(tmle_task, "R", fold_number)
      R <- tmle_task$get_tmle_node("R")
      #Assumes estimates is called before clever covariates every time
      HA <- (self$clever_covariates(tmle_task, fold_number, refit = T)$R)
      H <- rowSums(HA)
      IC <- H*(R - ER)
      risk <- self$risk_function
      psi <- risk$predict(tmle_task, type = "RR")
      # Returns the value of the RR function and the IC equation that needed to be solved.
      result <- list(psi = psi, IC  = as.matrix(HA*as.vector(R - ER)))
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- "RR"
      return(param_form)
    },
    risk_function = function(){
      private$.risk_function
    },
    update_nodes = function() {
      return(c(self$outcome_node))
    }
  ),
  private = list(
    .type = "RR",
    .supports_outcome_censoring = FALSE,
    .submodel_type_supported = c("logistic"),
    .risk_function = NULL

  )
)
