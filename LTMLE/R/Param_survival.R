#' Survival Curve
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
#' @section Constructor:
#'   \code{define_param(Param_survival, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood}}{the counterfactual likelihood for this treatment
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention
#'     }
#' }
#' @export
Param_survival <- R6Class(
  classname = "Param_survival",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list, ..., outcome_node, target_times = NULL) {
      # TODO: check outcome_node, current I(T<=t, delta=1), need I(T=t, delta=1)
      # W A processN processA
      private$.cf_likelihood <- make_CF_Likelihood(observed_likelihood, intervention_list)
      times <- setdiff(sort(unique(observed_likelihood$training_task$time)),0)
      private$.times <- times
      if(is.null(target_times)){
        target_times <- times
        private$.targeted <- rep(TRUE, length(target_times))
      } else {
        private$.targeted <-  rep(TRUE, length(target_times))
      }
      private$.target_times <- target_times
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
    },
    long_to_mat = function(x,id, time){
      dt <- data.table(id=id,time=time,x=as.vector(x))
      wide <- dcast(dt, id~time, value.var="x")
      mat <- as.matrix(wide[,-1,with=FALSE])
      return(mat)
    },
    hm_to_sm = function(hm){
      # TODO: check
      sm <- t(apply(1-hm,1,cumprod))
      # sm <- cbind(1,sm[,-ncol(sm)])
      return(sm)
    },
    clever_covariates_internal = function(tmle_task = NULL, fold_number = "full", subset_times = FALSE, for_fitting = T) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      # Intervene on censoring/riskset so that hazard probabilities are returned
      # In particular, set the counting processes to 1 and specify that risk set should not be checked


      cf_data <- tmle_task$data[, c("processN", "processA")]
      cf_data[, "processN"] <- 1
      cf_data[, "processA"] <- 1
      cf_task = task$generate_counterfactual_task(tmle_task$uuid, new_data = cf_data, T)

      intervention_nodes <- names(self$intervention_list)
      intervention_nodes = "A"
      g_A <- self$observed_likelihood$get_likelihoods(cf_task, intervention_nodes, fold_number, drop_id = T)

      g_A <- unlist(g_A[,intervention_nodes, with = F])

      # I(A=1)
      #TODO this
      #cf_pA <- self$cf_likelihood$get_likelihoods(cf_task, intervention_nodes, fold_number, drop_id = T)
      cf_pA <- tmle_task$get_tmle_node("A")[,A]

      #Should already be matrix. Assumed to be in order of time
      Q <- self$observed_likelihood$get_likelihoods(cf_task, c("processN"), fold_number, drop_id = T, to_wide = T)
      # TODO: make bound configurable


      Q <- as.matrix(Q)
      Q <- bound(Q, 0.005)


      G <- self$observed_likelihood$get_likelihoods(cf_task, c("processA"), fold_number, drop_id = T, to_wide = T)
      G <- as.matrix(G)


      Q_surv <- as.matrix(self$hm_to_sm(Q))
      G_surv <- as.matrix(self$hm_to_sm(G))

      # fix t-1
      G_surv <- cbind(rep(1, nrow(G_surv)),G_surv[,-ncol(G_surv)])

      time <- cf_task$data$t

      times <- setdiff(sort(unique(time)),0)
      #t_mat = matrix(rep(ks,nrow(cf_pA)), nrow = nrow(cf_pA), ncol = length(ks), byrow =T)

      # Compute clever covariates for one  tgt time in long format
      # stacked by vertically in chunks of time

      hk_at_tgt <- function(tgt_time) {
        index_set_to_zero <- which(times > tgt_time)


        # Get survival at target time
        Q_surv_tgt <- Q_surv[,which(tgt_time ==times)]
        # This is the clever covariate matrix at a single target time


        clever_cov = -1 *(cf_pA*Q_surv_tgt/g_A)/ (G_surv*Q_surv)
        if(length(index_set_to_zero) >0){
          clever_cov[,index_set_to_zero] <- 0
        }

        return(as.vector(clever_cov))
      }
      target_times <- self$target_times
      # bind the clever covariates for each tgt time by columns
      HA <- as.matrix(do.call(cbind, lapply(target_times, hk_at_tgt)))

      # get observed counting process
      # Contains columns t and id

      observed_N <- tmle_task$get_tmle_node("processN", include_time = T, include_id = T)
      at_risk <- tmle_task$get_regression_task("processN", drop_censored = F)$get_data(,"at_risk")

      # Only those at risk have likelihood updated.


      if(for_fitting){
      observed_N_wide <- reshape(observed_N, idvar = "id", timevar = "t", direction = "wide")
      to_dNt <- function(v){

        dt = data.table(matrix(c(0,diff(unlist(v))), nrow = 1))
        colnames(dt) <- paste0("d", colnames(v))
        return(dt)
      }
      # Converts to dNt format
      observed_dN_wide <- observed_N_wide[, to_dNt(.SD), by = id]
      observed_dN_wide$id <- NULL



      jump_time <- apply(observed_N_wide, 1, function(v){
        time <- which(v==1)
        if(length(time) ==0){
          return(Inf)
        }
        return(min(time))
      })
      t_mat <- matrix(1:ncol(observed_dN_wide), ncol = ncol(observed_dN_wide), nrow = nrow(observed_dN_wide), byrow = T )
      # TODO dont recompute this every time
      ind = apply(((t_mat <= jump_time)),2,as.numeric)


      # compute scaled and zeroed residuals vector


      residuals = as.vector((as.matrix(observed_dN_wide) - Q)*ind/nrow(observed_N_wide))


      #  Compute EIC component for one step/convergence criterion
      D1 = colSums(HA*residuals)
      new_comps <- list(list(processN = D1 ))
      names(new_comps) <- tmle_task$uuid
      new_comps
      private$.D_cache <- c(private$.D_cache, new_comps)

      }
      zero_rows <- which(at_risk == 0)
      # TODO dont compute HA for all people
      HA[zero_rows,] <- 0
      #Returns clever covariates and EIC component
      return(list(processN = list(H = HA)))
    },
    submodel_info = function(){
      #Returns list of submodel info
      # This includes submodel function, loss function, family object (for glm), offset_transform function (for glm)
      return(list(processN = list(submodel_family = self$submodel_family, offset_transform = self$offset_transform, loss = self$loss_function, submodel = self$submodel )))
    },
    offset_transform = function(initial, observed){
      # convert initial likelihood to to hazard format
      return(qlogis(ifelse(observed == 1, initial, 1 - initial)))
    },
    submodel_family = function(){
      return(binomial())
    },
    submodel = function(epsilon, initial, H, observed) {
      # Logistic submodel requires p/1-p
      observed <- observed$processN
      initial <- initial$processN
      initial <- ifelse(observed ==1, initial, 1 - initial)
      H <- H$processN
      result <- plogis(qlogis(initial) + H %*% epsilon)
      # convert hazard format to likelihood format
      result <- ifelse(observed ==1, result, 1 - result)
      list(processN = result)
    },
    loss_function = function(estimate, observed) {
      # Assumes estimate and observed are single (long) vectors
      -log(estimate$processN)
    },
    clever_covariates = function(tmle_task, fold_number = "full", for_fitting = T){
      self$clever_covariates_internal(tmle_task, fold_number, subset_times = TRUE, for_fitting = for_fitting)
    },
    get_EIC_component = function(task, node) {
      return(private$.D_cache[[task$uuid]][[node]])
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      #TODO everything
      result <- list(psi = psi, IC = IC_long)
      return(result)
    }
  ),
  active = list(
    # TODO: modify
    name = function() {
      param_form <- sprintf("E[P(T > %s|%s, W)]", self$times, self$cf_likelihood$name)
      return(param_form)
    },
    cf_likelihood = function() {
      return(private$.cf_likelihood)
    },
    intervention_list = function() {
      return(self$cf_likelihood$intervention_list)
    },
    update_nodes = function() {
      return(self$outcome_node)
    },
    times = function() {
      return(private$.times)
    },
    target_times = function() {
      return(private$.target_times)
    }
  ),
  private = list(
    .D_cache = list(),
    .type = "survival",
    .cf_likelihood = NULL,
    .supports_outcome_censoring = TRUE,
    .times = NULL,
    .target_times = NULL
  )
)
