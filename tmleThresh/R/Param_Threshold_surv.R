#' Estimates the threshold function E_W E[Y | A >=c, W] for a range of given values threshold values c
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @importFrom dplyr near
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_ATT, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention.
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood_treatment}}{the counterfactual likelihood for the treatment
#'     }
#'     \item{\code{cf_likelihood_control}}{the counterfactual likelihood for the control
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention
#'     }
#' }
#' @export
Param_thresh_survival <- R6Class(
  classname = "Param_thresh_survival",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, thresh_node = "A", outcome_node = c("Nt"), type = 1, target_times = NULL) {
      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)

      #compute_thresh_estimate(likelihood = observed_likelihood, fold_number = "full", type = type)
      #compute_thresh_estimate(likelihood = observed_likelihood, fold_number = "validation", type = type)


      private$.censoring_node <-"At"
      private$.thresh_node <- thresh_node
      cutoffs <- observed_likelihood$factor_list$Nt$learner$cutoffs

      private$.cutoffs <- cutoffs
      private$.strict_threshold <- F
      private$.target_times <- target_times
      #rivate$.cf_task <- cf_task
      private$.type_est <- type

    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", for_fitting = F) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      thresh_node <- private$.thresh_node
      censoring_node <- private$.censoring_node
      cutoffs <- private$.cutoffs

      cdfS <- as.vector(self$observed_likelihood$get_likelihood(tmle_task, thresh_node, fold_number))

      cdfS <- bound(cdfS, c(0.0005, .9995))
      cdfS_wide <- matrix(cdfS, ncol = length(cutoffs))
      if("A_learned" %in% names(tmle_task$npsem)) {
        S <- self$observed_likelihood$get_likelihood(tmle_task, "A_learned", fold_number)
      } else {
        S <- tmle_task$get_tmle_node(thresh_node)
      }


      if(private$.type_est == 1) {
        indS <- as.vector(unlist(lapply(cutoffs, function(cutoff) {as.numeric(S >= cutoff)})))
      } else if(private$.type_est == 0) {
        indS <- as.vector(unlist(lapply(cutoffs, function(cutoff) {as.numeric(S < cutoff)})))
      }
      indS_wide <- matrix(indS, ncol = length(cutoffs))

      #Uses
      times <- sort(tmle_task$npsem[["Nt"]]$time)

      if(!tmle_task$force_at_risk) {
        cf_task <- tmle_task$clone()
        cf_task$force_at_risk <- T
      } else {
        cf_task <- tmle_task
      }

      Q <- matrix(self$observed_likelihood$get_likelihood(cf_task, "Nt", fold_number), ncol = length(cutoffs))
      # Survival of N

      Q <- bound(Q, 0.0005)
      Q_surv <- apply(Q, 2, function(v) {
        #Survival stacked by person
        as.vector((apply( 1 -matrix(v, ncol = length(times), byrow = T), 1, cumprod)))
      })
      #Survival left continuous of A
      G <- self$observed_likelihood$get_likelihood(cf_task, "At", fold_number)
      G <- bound(G, 0.0005)
      G_surv <- apply( 1 -matrix(G, ncol = length(times), byrow = T), 1, cumprod)
      G_surv <- as.vector(rbind(rep(1, ncol(G_surv)), G_surv[-1,]))


      if(nrow(Q) != nrow(Q_surv)) {
        stop("Row numbers of surv and hazard dont match.")
      }



      # Gets stacked clever covariates

      hk_at_tgt <- function(tgt_time) {
        return(as.vector(do.call(cbind, lapply(1:ncol(Q_surv), function(j) {
          cdfS <- 1 - cdfS_wide[,j]

          indS <- indS_wide[,j]
          Q_surv <- Q_surv[,j]
          G_surv <- G_surv
          Q_surv <- matrix(Q_surv, ncol = length(times), byrow = F)
          G_surv <- matrix(G_surv, ncol = length(times), byrow = F)
          index_set_to_zero <- which(times > tgt_time)

          # Get survival at target time
          Q_surv_tgt <- Q_surv[,which(tgt_time ==times)]
          # This is the clever covariate matrix at a single target time


          clever_cov = -1 *(indS*Q_surv_tgt/cdfS)/ (G_surv*Q_surv)
          if(length(index_set_to_zero) >0){
            clever_cov[,index_set_to_zero] <- 0
          }
          return(as.vector(t(clever_cov)))
        }))))
      }
      target_times <- self$target_times
      if(is.null(target_times)) {
        target_times <- times
        private$.target_times <- times
      }
      HA <- as.matrix(do.call(cbind, lapply(target_times, hk_at_tgt)))
      at_risk <- unlist(tmle_task$get_regression_task("Nt", drop_censored = F, expand = T)$get_data(,"at_risk"))

      zero_rows <- as.numeric(at_risk == 1)
      HA <- HA * zero_rows
      observed_N <- tmle_task$get_tmle_node("Nt", include_time = T, include_id = T, expand= T)

      IC_N <- NULL
      if(for_fitting) {
        observed_N_wide <- reshape(observed_N, idvar = "id", timevar = "t", direction = "wide")
        observed_N_wide$id <- NULL
        # TODO this is crazy slow, code better and dont recompute
        to_dNt <- function(v){
          dt = c(0,diff((v)))
          return(dt)
        }
        # Converts to dNt format
        observed_dN_wide <- t(apply(observed_N_wide, 1, to_dNt))




        jump_time <- nrow(observed_N_wide) - apply(observed_N_wide, 1, sum) + 1
        t_mat <- matrix(1:ncol(observed_dN_wide), ncol = ncol(observed_dN_wide), nrow = nrow(observed_dN_wide), byrow = T )
        # TODO dont recompute this every time
        ind = t_mat <= jump_time
        residuals = as.vector((as.vector(t(observed_dN_wide)) - Q)*as.vector(t(ind)))
        IC_N <- HA*residuals
        IC_N <- do.call(cbind, unlist(apply(IC_N, 2, function(v) {
          list(matrix(v, ncol = length(cutoffs)))
        }), recursive = F))
      }

      n = nrow(observed_N)
      k = length(cutoffs)
      H <- matrix(0, nrow = nrow(HA), ncol = k*ncol(HA))
      for(i in 1:k){
        first <- (i-1)*n + 1
        last <- (i)*n
        j = 1:ncol(HA) + (i-1)*ncol(HA)
        H[first:last,j] <- HA[first:last,]
      }

      return(list(Nt = H, IC = list(Nt = IC_N)))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      cf_task <- tmle_task$clone()
      cf_task$force_at_risk <- TRUE
      thresh_node <- private$.thresh_node
      censoring_node <- private$.censoring_node
      cutoffs <- private$.cutoffs

      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))

      # clever_covariates happen here (for this param) only, but this is repeated computation
      IC_N <- self$clever_covariates(tmle_task, fold_number)$IC[[self$outcome_node]]

      Y <- matrix(tmle_task$get_tmle_node(self$outcome_node, impute_censoring = TRUE), nrow = tmle_task$nrow)

      #get E[Y|A>=1, W]
      EY <- matrix(self$observed_likelihood$get_likelihood(tmle_task, self$outcome_node, fold_number), nrow = tmle_task$nrow)

      #get E[Y|A>=1, W]
      EY1 <- matrix(compute_thresh_estimate(self$observed_likelihood, tmle_task, type = private$.type_est, fold_number = fold_number, return_estimate = F), nrow = tmle_task$nrow)
      #EY1 <- matrix(self$observed_likelihood$get_likelihood(cf_task, self$outcome_node, fold_number), nrow = tmle_task$nrow)

      psi <- self$empirical_mean(tmle_task, EY1)


      IC <- IC_N  + t((t(EY1)  - psi))
      #weights <- tmle_task$get_regression_task(self$outcome_node)$weights
      result <- list(psi = psi, IC = IC )
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("ATE[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
      return(param_form)
    },
    cf_task = function() {
      return(private$.cf_task)
    },

    update_nodes = function() {
      return(c(self$outcome_node))
    },
    target_times = function() {
      private$.target_times
    }
  ),
  private = list(
    .type = "Threshold",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .supports_outcome_censoring = TRUE,
    .submodel_type_supported = c("logistic"),
    .supports_weights = T,
    .censoring_node = NULL,
    .thresh_node = NULL,
    .cutoffs = NULL,
    .strict_threshold = NULL,
    .cf_task = NULL,
    .type_est = NULL,
    .target_times = NULL
  )
)
