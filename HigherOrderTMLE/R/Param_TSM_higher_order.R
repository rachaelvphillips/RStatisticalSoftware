#' Average Treatment Effect
#'
#' Parameter definition for the Average Treatment Effect (ATE).
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
Param_TSM_higher_order <- R6Class(
  classname = "Param_TSM_higher_order",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, likelihood_tilde, outcome_node = "Y") {
      super$initialize(observed_likelihood, list(), outcome_node)
      private$.likelihood_tilde <- likelihood_tilde


    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", for_fitting = F) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      #cf_task0 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A = 0))

      cache_key <- tmle_task$uuid
      if(cache_key %in% private$.cf_tasks_uuid) {
        cf_task1 <- tmle_task
      } else {
        cf_task1 <- get0(cache_key, private$.task_cache, inherits = FALSE)
        if(is.null(cf_task1)) {
          cf_task1 <-  tmle_task$generate_counterfactual_task(paste0(tmle_task$uuid, "1"), data.table(A = 1))
          assign(cache_key, cf_task1, private$.task_cache)
          private$.cf_tasks_uuid <- c(private$.cf_tasks_uuid, cf_task1$uuid)
        }
      }

      A <- tmle_task$get_tmle_node("A")
      Y <- tmle_task$get_tmle_node("Y")
      Y <-  tmle_task$scale(Y, "Y")
      Q_tilde <- tmle_task$scale(self$likelihood_tilde$get_likelihood(tmle_task, "Y", fold_number), "Y")
      Q_tilde1  <- tmle_task$scale(self$likelihood_tilde$get_likelihood(cf_task1, "Y", fold_number), "Y")
      g_tilde <- self$likelihood_tilde$get_likelihood(tmle_task, "A", fold_number)
      g <- self$observed_likelihood$get_likelihood(tmle_task, "A", fold_number)
      g1 <- self$observed_likelihood$get_likelihood(cf_task1, "A", fold_number)
      g_tilde1 <- self$likelihood_tilde$get_likelihood(cf_task1, "A", fold_number)

      Q <- tmle_task$scale(self$observed_likelihood$get_likelihoods(tmle_task, "Y", fold_number), "Y")
      Q1 <- tmle_task$scale(self$observed_likelihood$get_likelihoods(cf_task1, "Y", fold_number), "Y")

      bounds <- 0.00005
      if(any(Q < bounds)) {
        Q <- bound(Q, bounds)
        warning("Q bounded")
      }
      if(any(g < bounds)) {
        g <- bound(g, bounds)
        warning("g bounded")
      }
      if(any(g1 < bounds)) {
        g1 <- bound(g1, bounds)
        warning("g1 bounded")
      }
      if(any(Q1 < bounds)) {
        Q1 <- bound(Q1, bounds)
        warning("Q1 bounded")
      }


      # First order TMLE
      #####
      if(for_fitting) {
        observed <- Y
        observed <-  Q_tilde1
        offset <- stats::qlogis(Q1)

        weights <- g_tilde1/g1
        if(any(weights > 10)){
          warning("weights bounded")
          weights <- bound(weights, c(0,10))
        }

        H <- as.matrix(1/g1)
        if(any(abs(H) > 50)) {
          warning("bound H")
          H <- bound(H, c(-50,50))
        }
        data <- list(H = H, observed = observed)





        suppressWarnings(epsilon <- coef(glm(observed ~ H -1, data, offset = offset, start = rep(0, ncol(data$H)), family = binomial(), weights = tmle_task$get_regression_task("Y")$weights * weights)))
        if(is.na(epsilon)) {
          epsilon <- 0
        }
        if(is.infinite(epsilon)) {
          warning("infinite eps")
        }
        private$.first_epsilon <- epsilon

      } else {
        epsilon <- private$.first_epsilon
      }

      Q_tmle <- stats::plogis(stats::qlogis(Q) + as.matrix(A/g1) %*% epsilon)
      Q_tmle1 <- stats::plogis(stats::qlogis(Q1) + as.matrix(1/g1) %*% epsilon)
      if(all(Q_tmle==0)) {
        print(head(Q_tmle))
        print(head(g1))
        print(epsilon)
        stop("tmle 0")
      }
      if(all(Q_tmle1==0)) {
        print(quantile(Q_tmle))
        print(quantile(g1))
        print(quantile(Q_tilde1))
        print(epsilon)
        print(quantile(1/g1))
        print(quantile(H))
        print(quantile(weights))
        print(quantile(offset))
        stop("tmle1 0")
      }
      if(any(is.na(Q_tmle1))) {
        stop("NA qtmle1")
      }
      if(any(is.na(Q_tmle))) {
        stop("NA qtmle1")
      }
      Q_tmle[is.na(Q_tmle)] <- 0
      Q_tmle1[is.na(Q_tmle1)] <- 0
      #####


      # Second order TMLE
      #####
      Q1Q_1 <- Q1 * (1 - Q1)
      Q1Q_1_tmle <- Q_tmle1 * (1 - Q_tmle1)

      inv_term <- 1/mean(g_tilde1 / g1^2 *Q1Q_1_tmle)
      prod_term <- mean( Q1Q_1_tmle / g1)
      C_tilde_1 <- inv_term * prod_term
      C_g <- A/g1
      #####
      # C_y <- A/g1 * Q1Q_1_tmle /Q1Q_1
      # - C_tilde_1 * C_g * (g_tilde1/g1) * (Q1Q_1_tmle)/Q1Q_1

      C_y <- A/g1 * Q1Q_1_tmle /Q1Q_1 * (1 - C_tilde_1 * g_tilde1/g1)
      C_y_tilde <- g_tilde1/g1 * Q1Q_1_tmle /Q1Q_1 * (1 - C_tilde_1 * g_tilde1/g1)
      if(FALSE) {
        C_y <- cbind( C_y,A/g1)
      }

      C_y <- bound(C_y, c(-50,50))
      # C_a <- -epsilon * Q1Q_1_tmle / g1^2 -
      #   (C_tilde_1 * g_tilde1/g1^2) * (Q_tilde1 - Q_tmle1) +
      #   (C_tilde_1 * g_tilde1/g1^3) * epsilon * Q1Q_1_tmle

      C_a <- -epsilon * Q1Q_1_tmle / g1^2 * (1 - C_tilde_1 * g_tilde1/g1) -
        C_tilde_1 * g_tilde1/g1^2 * (Q_tilde1 - Q_tmle1)

      C_a <- bound(C_a, c(-50,50))

      ICY <- as.matrix(C_y) * as.vector(Y - Q)
      ICA <- C_a * (A - g)

      ICA_tilde <- C_a * (g_tilde1 - g1)
      ICY_tilde <- C_y_tilde * (Q_tilde - Q)
      if(any(is.infinite(ICY))) {
        stop("unboundedyinf")
      }
      if(any(is.infinite(ICA))) {
        stop("unboundedainf")
      }

      if(any(is.na(ICY))) {
        stop("unboundedyna")
      }
      if(any(is.na(ICA))) {
        stop("unboundedana")
      }
      #####
      return(list(Y = as.matrix(C_y), A = as.matrix(C_a),
                 IC= list(A = ICA, Y = ICY), new_Q = list(Q_tmle, Q_tmle1),
                 ICtilde = list(
                   A = ICA_tilde,
                   Y = ICY_tilde
                 )))
    },
    update_last = function(tmle_task, fold_number = "full"){
      HA <- self$clever_covariates(tmle_task, fold_number, for_fitting = T)
      likelihood_factor <- self$observed_likelihood$factor_list[["Y"]]
      cf_task1 <- get0(tmle_task$uuid, private$.task_cache, inherits = FALSE)

      self$observed_likelihood$cache$set_values(likelihood_factor, tmle_task, self$observed_likelihood$updater$step_number, fold_number, HA$new_Q[[1]], node = "Y")
      self$observed_likelihood$cache$set_values(likelihood_factor, cf_task1, self$observed_likelihood$updater$step_number, fold_number, HA$new_Q[[2]], node = "Y")
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {

      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      cache_key <- tmle_task$uuid

      cf_task1 <- get0(cache_key, private$.task_cache, inherits = FALSE)
        if(is.null(cf_task1)) {
          cf_task1 <-  tmle_task$generate_counterfactual_task(paste0(tmle_task$uuid, "1"), data.table(A = 1))
          assign(cache_key, cf_task1, private$.task_cache)
          private$.cf_tasks_uuid <- c(private$.cf_tasks_uuid, cf_task1$uuid)
        }


      # clever_covariates happen here (for this param) only, but this is repeated computation
      HA <- self$clever_covariates(tmle_task, fold_number, for_fitting = T)


      # todo: make sure we support updating these params

      Y <- tmle_task$get_tmle_node("Y", impute_censoring = TRUE)
      A <- tmle_task$get_tmle_node("A")
      Q <- self$observed_likelihood$get_likelihood(tmle_task, "Y", fold_number)
      g <- self$observed_likelihood$get_likelihood(tmle_task, "A", fold_number)

      EY1 <- self$observed_likelihood$get_likelihood(cf_task1, "Y", fold_number)
     # EY0 <- self$observed_likelihood$get_likelihood(cf_task0, self$outcome_node, fold_number)
      psi <- self$empirical_mean(tmle_task, EY1)
      #psi <- mean(EY1 - EY0)

      IC <- HA$Y[,1] * (Y - Q) +  HA$A * (A - g) + EY1 - psi

      if(any(is.infinite(IC))) {
        stop("unbounded")
      }
      if(any(is.infinite(IC))) {
        stop("unbounded")
      }


      result <- list(psi = psi, IC = as.matrix(IC), ICtilde = cbind(HA$ICtilde$A, HA$ICtilde$Y))
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("ATE[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
      return(param_form)
    },
    likelihood_tilde = function() {
      return(private$.likelihood_tilde)
    },
    update_nodes = function() {
      return(c("Y", "A"))
    }
  ),
  private = list(
    .type = "ATE",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .supports_outcome_censoring = TRUE,
    .likelihood_tilde = NULL,
    .first_epsilon = NULL,
    .task_cache = new.env(),
    .cf_tasks_uuid = c()
  )
)
