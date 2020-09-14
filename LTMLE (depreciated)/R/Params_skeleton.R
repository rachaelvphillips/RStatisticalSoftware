#' Base Class for Defining Parameters
#'
#' A parameter is a function of the likelihood. Once given a \code{\link{Likelihood}} object, a parameter will a value.
#' These objects also contain information about the efficient influence function (EIF) of a parameter, as well as its clever covariate(s).
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'

#' @export
Param_base <- R6Class(
  classname = "Param_base",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(observed_likelihood, ..., nodes_to_target = "Y") {
      # What is nodes to target for?
      nodes_to_target = c("A%B%C", "D")
      expanded_nodes_map <- c("A" = c("A%B%C",1), "B" = c("A%B%C",2), "C" = c("A%B%C",3), "D" = c("D",1))
      private$.expanded_nodes_map <- expanded_nodes_map
      private$.observed_likelihood <- observed_likelihood
      private$.nodes_to_target <- nodes_to_target

      # if (!is.null(observed_likelihood$censoring_nodes[[outcome_node]])) {
      #   if (!self$supports_outcome_censoring) {
      #     stop(sprintf(
      #       "%s has censoring mechanism, but %s does not yet support outcome node censoring",
      #       outcome_node, class(self)[1]
      #     ))
      #   }
      # }

      if (inherits(observed_likelihood, "Targeted_Likelihood")) {
        # register parameter with updater
        observed_likelihood$updater$register_param(self)
      } else if (inherits(observed_likelihood, "Likelihood")) {
        warning("Parameter was passed a non-Targeted Likelihood object so estimates cannot be updated from initial")
      } else {
        stop("Invalid Likelihood class: ", class(observed_likelihood))
      }
    },
    clever_covariates = function(tmle_task, fold_number = "full"){
      # If three nodes need to be targetted together then their key is "node1%node2%node3"
      # The element of list should be clever covariates (a list) and the EIC_comp for the bundle (a vector) ...
      # the EIC comp is the convergence criterion for the bundle and is also used for collapsing in one step
      result <- list("A%B%C" = list(H= H, D = EIC_comp), "D" = list(H = H, D = EIC_comp) )

    },
    submodel_info = function(){
      # For each bundle specify the loss and submodel that needs to be used
      return(list("A%B%C" = list(loss, submodel), "D" = list(loss, submodel)))
    },
    loss_function = function(estimate_list, observed_list){
      # Example of loss function which takes list
      loss <- lapply(seq_along(estimate_list), function(i){
        get_loss(estimate_list[[i]], observed_list[[i]])
      })
      return(rowSums(do.call(cbind, loss)))
    },
    submodel = function(epsilon, initial_list, H_list){
      # TODO Also handle when only one element is passed, and not list
      #  Example of submodel that updates multiple nodes with single epsilon
      submodel_helper = function(epsilon, initial, H) {
        plogis(qlogis(initial) + H %*% epsilon)
      }
      new_likelihoods <- lapply(seq_along(initial_list), function(i){
        submodel_helper(epsilon, initial_list[[i]], H_list[[i]] )
      })
      # returns list
      return(new_likelihoods)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      # Compute full estimate and EIC
      result <- c(Estimates, EIC)
      return(result)
    },
    print = function() {
      cat(sprintf("%s: %s\n", class(self)[1], self$name))
    }

  ),
  active = list(
    name = function() {
      return(private$.type)
    },
    type = function() {
      return(private$.type)
    },
    expanded_nodes_map = function(){
      private$.expanded_nodes_map
    },

    observed_likelihood = function() {
      return(private$.observed_likelihood)
    },
    nodes_to_target = function() {
      return(private$.nodes_to_target)
    },
    # TODO not sure what this is
    supports_outcome_censoring = function() {
      return(private$.supports_outcome_censoring)
    },
    targeted = function(){
      return(rep(private$.targeted, length(self$nodes_to_target)))
    }
  ),
  private = list(
    .type = "undefined",
    .observed_likelihood = NULL,
    .nodes_to_target = NULL,
    .expanded_nodes_map = NULL,
    .targeted = TRUE,
    .supports_outcome_censoring = FALSE
  )
)


