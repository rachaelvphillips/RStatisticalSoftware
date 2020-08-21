#' Defines an update procedure (submodel+loss function)
#'
#' Current Limitations:
#' loss function and submodel are hard-coded (need to accept arguments for these)
#' @section Constructor:
#'   \code{define_param(maxit, cvtmle, one_dimensional, constrain_step, delta_epsilon, verbose)}
#'
#'   \describe{
#'     \item{\code{maxit}}{The maximum number of update iterations
#'     }
#'     \item{\code{cvtmle}}{If \code{TRUE}, use CV-likelihood values when
#'        calculating updates.
#'     }
#'     \item{\code{one_dimensional}}{If \code{TRUE}, collapse clever covariates
#'        into a one-dimensional clever covariate scaled by the mean of their
#'        EIFs.
#'     }
#'     \item{\code{constrain_step}}{If \code{TRUE}, step size is at most
#'        \code{delta_epsilon} (it can be smaller if a smaller step decreases
#'        the loss more).
#'     }
#'     \item{\code{delta_epsilon}}{The maximum step size allowed if
#'        \code{constrain_step} is \code{TRUE}.
#'     }
#'     \item{\code{convergence_type}}{The convergence criterion to use: (1)
#'        \code{"scaled_var"} corresponds to sqrt(Var(D)/n)/logn (the default)
#'        while (2) \code{"sample_size"} corresponds to 1/n.
#'     }
#'     \item{\code{fluctuation_type}}{Whether to include the auxiliary covariate
#'        for the fluctuation model as a covariate or to treat it as a weight.
#'        Note that the option \code{"weighted"} is incompatible with a
#'        multi-epsilon submodel (\code{one_dimensional = FALSE}).
#'     }
#'     \item{\code{use_best}}{If \code{TRUE}, the final updated likelihood is set to the
#'        likelihood that minimizes the ED instead of the likelihood at the last update
#'        step.
#'     }
#'     \item{\code{verbose}}{If \code{TRUE}, diagnostic output is generated
#'        about the updating procedure.
#'     }
#'     }
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Update <- R6Class(
  classname = "tmle3_Update",
  portable = TRUE,
  class = TRUE,
  public = list(
    # TODO: change maxit for test
    initialize = function(maxit = 100, cvtmle = TRUE, one_dimensional = FALSE,
                          constrain_step = FALSE, delta_epsilon = 1e-4,
                          convergence_type = c("scaled_var", "sample_size"),
                          fluctuation_type = c("standard", "weighted"),
                          optim_delta_epsilon = TRUE,
                          use_best = FALSE,
                          verbose = FALSE) {
      private$.maxit <- maxit
      private$.cvtmle <- cvtmle
      private$.one_dimensional <- one_dimensional
      private$.constrain_step <- constrain_step
      private$.delta_epsilon <- delta_epsilon
      private$.convergence_type <- match.arg(convergence_type)
      private$.fluctuation_type <- match.arg(fluctuation_type)
      private$.optim_delta_epsilon <- optim_delta_epsilon
      private$.use_best = use_best
      private$.verbose <- verbose
    },
    collapse_covariates = function(ED, clever_covariates) {
      #ED <- ED_from_estimates(estimates)
      EDnormed <- ED / norm(ED, type = "2")
      collapsed_covariate <- clever_covariates %*% EDnormed
      collapsed_covariate <-matrix(collapsed_covariate,ncol=1)

      return(collapsed_covariate)
    },
    update_step = function(likelihood, tmle_task, fold_number = "full", update_nodes = NULL) {
      # Allow only updating subset of nodes per step
      if(is.null(update_nodes)){
        update_nodes <- self$update_nodes
      }
      unkeyed_nodes <- unlist(lapply(update_nodes, function(node) self$key_to_node_bundle(node)))
      current_step <- self$step_number + 1

      # initialize epsilons for this step
      na_epsilons <- as.list(rep(NA, length(unkeyed_nodes)))
      names(na_epsilons) <- unkeyed_nodes
      private$.epsilons[[current_step]] <- na_epsilons

      for (update_node in update_nodes) {
        # update_node coulde keyed multinode
        # get new submodel fit

        expanded_node <- self$key_to_node_bundle(update_node)
        submodel_data <- self$generate_submodel_data(
          likelihood, tmle_task,
          fold_number, update_node,
          drop_censored = TRUE, for_fitting = TRUE
        )
        #This will be a named epsilon vector
        new_epsilon <- self$fit_submodel(submodel_data)

        # update likelihoods
        if(length(expanded_node)==1){
          likelihood$update(as.vector(unlist(new_epsilon, use.names = F)), current_step, fold_number, update_node)
        } else {
          # If factors have shared epsilon
          for(node in expanded_node ){
            # TODO apply_update will fail
            likelihood$update(new_epsilon[[node]], current_step, fold_number, node)
          }
        }

        if (fold_number != "full") {
          # update full fit likelihoods if we haven't already
          if(length(expanded_node)==1){
            likelihood$update(as.vector(unlist(new_epsilon, use.names = F)), self$step_number, "full", update_node)
          } else {
            # If factors have shared epsilon
            for(node in expanded_node ){
              likelihood$update(new_epsilon[[node]],self$step_number, "full", node)
            }
          }
        }
        # Update shared epsilons
        if(length(expanded_node)==1){
          private$.epsilons[[current_step]][[update_node]] <- as.vector(unlist(new_epsilon, use.names = F))
        } else {
          for(node in expanded_node ){
            private$.epsilons[[current_step]][[node]] <- new_epsilon[[node]]
          }
        }
      }

      # update step number
      private$.step_number <- current_step
    },
    key_to_node_bundle = function(key){
      unlist(stringr::str_split(key, "%"))
    },
    node_bundle_to_key = function(node_bundle){
      paste(node_bundle, sep = "%")
    },


    generate_submodel_data = function(likelihood, tmle_task,
                                      fold_number = "full",
                                      update_node_key = "Y",
                                      drop_censored = FALSE, for_fitting = F) {

      # Generates a sub_model data in list form if needed
      update_node <- self$key_to_node_bundle(update_node_key)

      clever_covariates_info <- lapply(self$tmle_params, function(tmle_param) {
        tmle_param$clever_covariates(tmle_task, fold_number, for_fitting = for_fitting)[[update_node_key]]
      })
      clever_covariates <- lapply(clever_covariates_info, `[[`, "H")
      # Note that clever covariates are not collapsed in this function.
      # The EIC component is absorbed into the epsilon vector

      submodel_info <- lapply(self$tmle_params, function(tmle_param) {
        tmle_param$submodel_info()
      })

      reduce_f <- function(x,y){
        out = list()
        names_out <- self$key_to_node_bundle(update_node_key)
        for(i in seq_along(x)){
          out[[names_out[i]]] <- cbind(x[[i]], y[[i]])
        }
        return(out)
      }
      # Merge clever covariates for each node across parameters
      covariates_dt_list <- clever_covariates %>% reduce(reduce_f)
      if(!is.list(covariates_dt_list)) covariates_dt_list <- list(covariates_dt_list)
      #covariates_dt <- do.call(cbind, node_covariates)
      # Stack EIC node-specific components for each target parameter
      #EIC_comps_stacked <- unlist(EIC_comps)
      node_submodel_info <- lapply(submodel_info, `[[`, update_node_key)
      # TODO check that parameters share the same loss functions and submodels
      node_submodel_info <- node_submodel_info[[1]]
      loss <- node_submodel_info$loss
      submodel <- node_submodel_info$submodel
      family <- node_submodel_info$submodel_family
      offset_transform <- node_submodel_info$offset_transform


      observed_list <- lapply(update_node,function(node){
        observed <- unlist(tmle_task$get_tmle_node(node)[,node, with = F], use.names = F)
        observed <- tmle_task$scale(observed, node)
      })
      id_list = list()
      t_list = list()
      initial_list <- lapply(update_node , function(node, task, fold){

        initial <- likelihood$get_likelihood(node, tmle_task=task, fold_number = fold, verify_in_sync = F )
        id_list <<- c(id_list, list(initial$id ))
        t_list <<- c(t_list, list(initial$t ))
        initial <- initial[,node,with=F]
        initial <- unlist(initial[,node,with=F], use.names=F)
        initial <- tmle_task$scale(initial, node)
        initial <- bound(initial, 0.005)
        }, task = tmle_task, fold = fold_number)



      names(observed_list) <- update_node
      names(covariates_dt_list) <- update_node
      names(initial_list) <- update_node
      names(t_list) <- update_node
      names(id_list) <- update_node
      at_risk_list<- lapply(update_node, function(node) which(unlist(tmle_task$get_regression_task(node)$get_data(,"at_risk"), use.names=F)==1))

      if(drop_censored){
        observed_list <- lapply(seq_along(observed_list), function(i) {observed_list[[i]][at_risk_list[[i]]]})
        initial_list <- lapply(seq_along(initial_list), function(i) {initial_list[[i]][at_risk_list[[i]]]})
        t_list <- lapply(seq_along(t_list), function(i) {t_list[[i]][at_risk_list[[i]]]})
        id_list <- lapply(seq_along(id_list), function(i) {id_list[[i]][at_risk_list[[i]]]})
        covariates_dt_list <- lapply(seq_along(covariates_dt_list), function(i) {covariates_dt_list[[i]][at_risk_list[[i]],,drop = F]})
        at_risk_list <-  lapply(seq_along(at_risk_list), function(i) {at_risk_list[[i]][at_risk_list[[i]]]})
      }
      names(observed_list) <- update_node
      names(covariates_dt_list) <- update_node
      names(initial_list) <- update_node
      names(t_list) <- update_node
      names(id_list) <- update_node
      names(at_risk_list) <- update_node
      if(for_fitting) {
        training_task <- likelihood$training_task
      } else {
        training_task <- NULL
      }
      submodel_data <- list(
        observed = observed_list,
        H = covariates_dt_list,
        initial = initial_list,
        loss = loss,
        submodel = submodel,
        t_list = t_list,
        id_list = id_list,
        at_risk_list = at_risk_list,
        offset_transform = offset_transform,
        family = family,
        key = update_node_key,
        node = update_node,
        training_task = training_task
      )


      return(submodel_data)
    },
    fit_submodel = function(submodel_data) {
      update_node <- submodel_data$node

      if (self$one_dimensional) {
        #Collapsing is part of fitting
        observed_task <- submodel_data$training_task

        EIC_comps <- lapply(update_node, function(node){
          unlist(lapply(self$tmle_params, function(tmle_param) {
            tmle_param$get_EIC_component(observed_task, node)
          }))})


        names(EIC_comps) <- submodel_data$node
        new_H <- lapply(seq_along(submodel_data$H),  function(i){
          self$collapse_covariates(EIC_comps[[i]], submodel_data$H[[i]])
        } )
        names(new_H) <-  submodel_data$node
        submodel_data$H <- new_H
        #covariates_dt <- self$collapse_covariates(EIC_comps_stacked, covariates_dt)
      }

      observed <- unlist(submodel_data$observed, use.names=F )
      initial <- unlist(submodel_data$initial, use.names=F )
      H <- do.call(rbind, submodel_data$H)
      loss_function <- submodel_data$loss
      submodel <- submodel_data$submodel
      # family and offset are part of submodel and are passed through PARAM
      family <- submodel_data$family
      offset <- submodel_data$offset_transform(initial, observed)
      submodel_data_pooled = list(H = H, observed = observed,  initial = initial)



      if (self$constrain_step) {
        ncol_H <- ncol(H)
        if (!(is.null(ncol_H) || (ncol_H == 1))) {
          stop(
            "Updater detected `constrain_step=TRUE` but multi-epsilon submodel.\n",
            "Consider setting `collapse_covariates=TRUE`"
          )
        }


        risk <- function(epsilon) {
          submodel_estimate <- self$apply_submodel(submodel_data, epsilon, submodel)
          loss <- loss_function(submodel_estimate, submodel_data$observed)
          mean(loss)
        }


        if (self$optim_delta_epsilon) {
          optim_fit <- optim(
            par = list(epsilon = self$delta_epsilon), fn = risk,
            lower = 0, upper = self$delta_epsilon,
            method = "Brent"
          )
          epsilon <- optim_fit$par
        } else {
          epsilon <- self$delta_epsilon
        }

        risk_val <- risk(epsilon)
        risk_zero <- risk(0)

        # # TODO: consider if we should do this
        # if(risk_zero<risk_val){
        #   epsilon <- 0
        # }

        if (self$verbose) {
          cat(sprintf("risk_change: %e ", risk_val - risk_zero))
        }
      } else {
        if (self$fluctuation_type == "standard") {
          suppressWarnings({
            submodel_fit <- glm(observed ~ H - 1, submodel_data_pooled,
                                offset = offset,
                                family = family,
                                start = rep(0, ncol(H))
            )
          })
        } else if (self$fluctuation_type == "weighted") {
          if (self$one_dimensional) {
            suppressWarnings({
              submodel_fit <- glm(observed ~ -1, submodel_data_pooled,
                                  offset = offset,
                                  family = family,
                                  weights = as.numeric(H),
                                  start = rep(0, ncol(H))
              )
            })
          } else {
            warning(
              "Updater detected `fluctuation_type='weighted'` but multi-epsilon submodel.\n",
              "This is incompatible. Proceeding with `fluctuation_type='standard'`."
            )
            suppressWarnings({
              submodel_fit <- glm(observed ~ H - 1, submodel_data_pooled,
                                  offset = offset,
                                  family = family,
                                  start = rep(0, ncol(H))
              )
            })
          }
        }
        epsilon <- coef(submodel_fit)

        # NOTE: this protects against collinear covariates
        # (which we don't care about, we just want an update)
        epsilon[is.na(epsilon)] <- 0
      }

      if (self$verbose) {
        max_eps <- epsilon[which.max(abs(epsilon))]
        cat(sprintf("(max) epsilon: %e ", max_eps))
      }
      if (self$one_dimensional) {
        epsilon <- lapply(EIC_comps, function(comp) comp*epsilon)
        names(epsilon) <- update_node
      } else{
        epsilon <- as.list(replicate(list(epsilon), length(update_node)))
        names(epsilon) <- update_node
      }

      return(epsilon)
    },

    apply_submodel = function(submodel_data, epsilon, submodel) {
      # TODO
      submodel(epsilon, submodel_data$initial, submodel_data$H, submodel_data$observed)
    },
    debundle_submodel = function(bundle, node){
      submodel_data <- list(
      observed = bundle$observed[node],
      H = bundle$H[node],
      initial = bundle$initial[node],
      loss = bundle$loss,
      submodel = bundle$submodel,
      t = bundle$t_list[node],
      id = bundle$id_list[node],
      at_risk_list = bundle$at_risk_list[node],
      offset_transform = bundle$offset_transform,
      family = bundle$family
      )
      return(submodel_data)
    },
    apply_update = function(tmle_task, likelihood, fold_number, new_epsilon, update_node) {
      # Update_node will always be a single node (never a bundled nodes key)
      # TODO However, submodel_data won't work for nodes contained in bundled nodes

      submodel_data <- self$generate_submodel_data(
        likelihood, tmle_task,
        fold_number, update_node, drop_censored = FALSE, for_fitting = F
      )

      # TODO make sure that submodel data/clever covariates are cached.
      # Currently the above is computing the submodel data for all nodes in bundle...
      # Alternatively we could do apply_update in bundles, but then changed to Targeted_lik
      # would need to be made.

      submodel_data <- self$debundle_submodel(submodel_data, update_node)


      updated_likelihood <- self$apply_submodel(submodel_data, new_epsilon, submodel_data$submodel)

      # This is absorbed into clever covariates
      #updated_likelihood[-submodel_data$at_risk_list[[update_node]]] <- submodel_data$initial[[node]][-submodel_data$at_risk_list[[update_node]]]
      t <- unlist(submodel_data$t)
      id <- unlist(submodel_data$id)

      if(any(!is.finite(updated_likelihood[[update_node]]))){
        stop("Likelihood was updated to contain non-finite values.\n
             This is likely a result of unbounded likelihood factors")
      }
      # un-scale to handle bounded continuous
      updated_likelihood <- tmle_task$unscale(
        updated_likelihood,
        update_node
      )
      updated_likelihood <- data.table(t = t, id = id, unlist(updated_likelihood))
      setnames(updated_likelihood, c("t", "id", update_node))
      return(updated_likelihood)
    },
    check_convergence = function(tmle_task, fold_number = "full") {
      # TODO For each update_node, get the eic components (from clever cov result)
      # check each component separately
      # TODO change everything

      #estimates <- self$current_estimates

      n <- length(unique(tmle_task$id))
      if ( self$convergence_type == "scaled_var") {
        # TODO need to compute EIC variance for each parameter once
        # NOTE: the point of this criterion is to avoid targeting in an overly
        #       aggressive manner, as we simply need check that the following
        #       condition is met |P_n D*| / SE(D*) =< max(1/log(n), 1/10)
        estimates <- private$.EIC_sd

        norm_weights <- estimates
        ED_threshold <- 1 / sqrt(n)/log(n)

      } else if (self$convergence_type == "sample_size") {


        ED_threshold <- 1 / sqrt(n)/log(n)
      }
      #ED_threshold <- 1 / sqrt(n)/log(n)
      # get |P_n D*| of any number of parameter estimates
      list_of_EIC_norms <- lapply(self$update_nodes, function(node_key) {
        nodes <- self$key_to_node_bundle(node_key)
        comps <- lapply(nodes, function(node){
          unlist(lapply(self$tmle_params, function(param) param$get_EIC_component(tmle_task, node)))
        })
        comp <- Reduce(`+`, comps)
        if (self$convergence_type == "sample_size"){
          wghts <- 1/sqrt(length(comp) )
        } else {
          wghts <- 1/norm_weights[[node_key]]
          wghts <- wghts/ sqrt(sum(wghts^2))
        }

        return(norm(comp*(wghts), type = "2"))
      })



      ED_criterions <- as.vector(unlist(list_of_EIC_norms))
      print(data.table(ED_criterions))

      if (self$verbose) {
        cat(sprintf("max(abs(ED)): %e\n", ED_criterions))
      }
      passed <- self$update_nodes[ED_criterions <= ED_threshold]
      #Returns nodes that converged

      return(passed)

    },
    update_best = function(likelihood){
      current_step <- self$step_number
      ED <- private$.EDs[[current_step]]
      ED_2_norm <- sqrt(sum(ED^2))
      if(ED_2_norm<private$.best_ED){
        likelihood$cache$update_best()
        private$.best_ED <- ED_2_norm
      }
    },
    update = function(likelihood, tmle_task, update_spec = NULL) {
      # update_spec is a list of length(maxit) where each element specifies which nodes
      # to be updates at the iteration. This allows user to specify specific targeting strategies
      # TODO handle early convergence
      update_fold <- self$update_fold
      vars <- lapply(self$tmle_params, function(param) param$get_EIC_var(tmle_task, update_fold))

      list_of_EIC_vars <- lapply(self$update_nodes, function(node_key) {

        sd_eic <- sqrt(unlist(lapply(vars, `[[`, node_key)))


        return(sd_eic)
      })
      names(list_of_EIC_vars) <- self$update_nodes
      private$.EIC_sd <- list_of_EIC_vars
      if(is.null(update_spec)) {
        maxit <- private$.maxit
        update_spec <- as.list(rep(list(self$update_nodes), maxit))
      }
      # seed current estimates
      #private$.current_estimates <- lapply(self$tmle_params, function(tmle_param) {
       # tmle_param$estimates(tmle_task, update_fold)
      #})
      # TODO make sure EIC variances are computed during first iteration and stored
      # Maybe in param object
      converged_nodes <- c()
      count = 1
      for (subset_nodes in update_spec) {
        # Only update nonconverged nodes
        subset_nodes <- setdiff(subset_nodes, converged_nodes)

        count = count+1
        self$update_step(likelihood, tmle_task, update_fold, subset_nodes)

        # update estimates based on updated likelihood
        # TODO current estimates not needed. we only need EIC comp from clever cov
        #private$.current_estimates <- lapply(self$tmle_params, function(tmle_param) {
         # tmle_param$estimates(tmle_task, update_fold)
        #})
        #TODO check convergence per factor.
        # Maybe subset update_spec so that only necessary nodes are updated.
        converged_nodes <- self$check_convergence(tmle_task, update_fold)
        if(length(converged_nodes) == length(self$update_nodes)){
          # If all nodes converged then stop
          break
        }


        if(self$use_best){
          self$update_best(likelihood)
        }
      }

      if(self$use_best){
        self$update_best(likelihood)
        likelihood$cache$set_best()
      }

    },
    register_param = function(new_params) {
      if (inherits(new_params, "Param_base")) {
        new_params <- list(new_params)
      }
      private$.tmle_params <- c(private$.tmle_params, new_params)
      private$.targeted_components <- unlist(lapply(private$.tmle_params, `[[`, "targeted"))
      new_update_nodes <- unlist(lapply(new_params, `[[`, "update_nodes"))
      private$.update_nodes <- unique(c(
        private$.update_nodes,
        new_update_nodes
      ))
    }
  ),
  active = list(
    epsilons = function() {
      return(private$.epsilons)
    },
    EDs = function() {
      return(private$.EDs)
    },
    tmle_params = function(new_params = NULL) {
      if (!is.null(new_params)) {
        if (inherits(new_params, "Param_base")) {
          new_params <- list(new_params)
        }
        private$.tmle_params <- new_params
        private$.update_nodes <- unique(unlist(lapply(
          new_params, `[[`,
          "update_nodes"
        )))
      }
      return(private$.tmle_params)
    },
    update_nodes = function() {
      return(private$.update_nodes)
    },
    update_fold = function() {
      if (self$cvtmle) {
        # use training predictions on validation sets
        update_fold <- "validation"
      } else {
        # use predictions from full fit
        update_fold <- "full"
      }
    },
    step_number = function() {
      return(private$.step_number)
    },
    maxit = function() {
      return(private$.maxit)
    },
    cvtmle = function() {
      return(private$.cvtmle)
    },
    one_dimensional = function() {
      return(private$.one_dimensional)
    },
    constrain_step = function() {
      return(private$.constrain_step)
    },
    delta_epsilon = function() {
      return(private$.delta_epsilon)
    },
    convergence_type = function() {
      return(private$.convergence_type)
    },
    fluctuation_type = function() {
      return(private$.fluctuation_type)
    },
    optim_delta_epsilon = function() {
      return(private$.optim_delta_epsilon)
    },
    use_best = function() {
      return(private$.use_best)
    },
    verbose = function() {
      return(private$.verbose)
    },
    current_estimates = function(){
      return(private$.current_estimates)
    }
  ),
  private = list(
    .epsilons = list(),
    .EDs = list(),
    .best_ED = Inf,
    .tmle_params = NULL,
    .update_nodes = NULL,
    .step_number = 0,
    # TODO: change maxit for test
    .maxit = 100,
    .cvtmle = NULL,
    .one_dimensional = NULL,
    .constrain_step = NULL,
    .delta_epsilon = NULL,
    .optim_delta_epsilon = NULL,
    .convergence_type = NULL,
    .fluctuation_type = NULL,
    .use_best = NULL,
    .verbose = FALSE,
    .targeted_components = NULL,
    .current_estimates = NULL,
    .EIC_sd = NULL
  )
)
