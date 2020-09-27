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
Param_thresh_eff <- R6Class(
  classname = "Param_thresh_eff",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, thresh_node = "A", outcome_node = "Y", type = 1, thresholds = NULL, num_bins = 100 ) {
      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)



      range <- range(observed_likelihood$training_task$get_tmle_node("A"))
      A_grid <- seq(from = range[1], to = range[2], length.out = num_bins)
      private$.A_grid <- A_grid
      private$.censoring_node <- (observed_likelihood$censoring_nodes[[outcome_node]])
      private$.thresh_node <- thresh_node
      if(is.null(thresholds)) {
        thresholds <- A_grid[seq(2, length(A_grid) - 5, 10)]
      }
      private$.thresholds <- thresholds
      private$.strict_threshold <- F
      private$.type_est <- type
      private$.cache <- new.env()

    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", for_fitting = F) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      t = proc.time()
      thresholds <- private$.thresholds
      likelihood <- self$observed_likelihood
      A_grid <- private$.A_grid
      thresh_node <- private$.thresh_node
      censoring_node <- private$.censoring_node
      thresholds <- private$.thresholds
      #g <- likelihood$get_likelihood(tmle_task, "A")
      Q  <- likelihood$get_likelihood(tmle_task, "Y")
      Y <- tmle_task$get_tmle_node("Y")
      A <- tmle_task$get_tmle_node("A")
      variable_type <- tmle_task$npsem$A$variable_type$type

      long_task <- self$make_long_task(tmle_task)
      is_long <- F
      if(long_task$uuid == tmle_task$uuid) {
        is_long <- T
      }
      g_long <- likelihood$get_likelihood(long_task, "A")
      Q_long <- likelihood$get_likelihood(long_task, "Y")

      long_g_id <- data.table(A = long_task$get_tmle_node("A", include_id = F), g = g_long, id = NA)
      set(long_g_id, , "id", long_task$get_data(,"trueid"))

      Q_long_id <- data.table(Q = Q_long, id = NA)
      set(Q_long_id,,  "id", long_task$get_data(,"trueid"))


      if(variable_type == "continuous") {
        A_diff <- long_g_id[,  diff(A), by = id]
        g_int <- na.omit(long_g_id[, frollsum(.SD, 2, na.rm=T, align = "right"), by = id,.SDcols = "g"])[[2]] /2
        Q_int <- na.omit(Q_long_id[, frollsum(.SD, 2, na.rm=T, align = "right"), by = id, .SDcols = "Q"])[[2]] /2
        A_int <- long_g_id[, A, by = id][[2]]
        integral <- data.table(A_diff,  g = g_int, Q = Q_int)
        setnames(integral, c("id", "A_diff", "g", "Q"))
        cdf <- integral[, c(0,cumsum(A_diff*g)), by = id]
        psi_W <- integral[,  c(rev(cumsum(rev(A_diff * g * Q ))),0), by = id]
        psi_W$A <- A_int
        setnames(psi_W, c("id", "int", "A"))


      } else {
        A_diff <- 1
        g_int <- long_g_id[, g, by = id][[2]]
        Q_int <- Q_long_id[, Q, by = id][[2]]
        A_int <- long_g_id[, A, by = id][[2]]
        integral <- data.table(A_diff,  g = g_int, Q = Q_int)
        setnames(integral, c("id", "A_diff", "g", "Q"))
        cdf <- integral[, c(0, cumsum(g)[-length(g)]), by = id]
        psi_W <- integral[,  rev(cumsum(rev( g * Q ))), by = id]
        psi_W$A <- A_int
        setnames(psi_W, c("id", "int", "A"))
      }

      trueid <- Q_long_id$id

      #For each id, CDF values at A_grid

      indices <- sapply(thresholds, function(a) {
        which.min( abs(A_grid - a))
      })
      a_mat <- do.call(cbind, lapply(thresholds, function(a) {
        rep(a, length(unique(tmle_task$id)))
      }))
      a_mat <- a_mat <=A
      a_long <- as.vector(a_mat)
      cdfb <- cdf[[2]]
      psi_Wb <- psi_W[, int]
      cdfb <- as.vector(matrix(cdfb, ncol = length(A_grid), byrow = T)[,indices])
      psi_Wb <- as.vector(matrix(psi_Wb, ncol = length(A_grid), byrow = T)[,indices])
      #print(data.table(matrix(psi_Wb, ncol = length(thresholds))))
      survlong <-  1 - cdfb
      psi_Wlong <- psi_Wb / survlong
      if(is_long) {
        Qwide <- matrix(Q, nrow = length(unique(trueid)))
        amat <- a_mat

        IC_A <- do.call(rbind, lapply(1:ncol(Qwide), function(i) {
          Q <- Qwide[,i]
          result <- matrix((Q - psi_Wlong )  / survlong, ncol =length(thresholds))
          return(result)
        })) * a_mat

      } else {
        IC_long <- (Q - psi_Wlong ) * a_long / survlong
        IC_A <- matrix(IC_long, ncol = length(thresholds))
      }


      psi_W <-  matrix(psi_Wlong, ncol = length(thresholds))
      #print(data.table(matrix(survlong, ncol = length(thresholds))))
      # results <- lapply(thresholds, function(a) {
      #   index <- which.min( abs(A_grid - a))
      #   surv <- 1 - cdf[, .SD[index] , by = id][[2]]
      #
      #   psi_W <- psi_W[, int[index], by = id][[2]] / surv
      #
      #
      #   IC <- (Q - psi_W ) * (A>=a) / surv
      #
      #   return(list(IC = IC, psi_W = psi_W))
      # })
      # IC_A <-  as.matrix(do.call(cbind, lapply(results, `[[`, "IC")))

      #phi_W <-  as.matrix(do.call(cbind, lapply(results, `[[`, "psi_W")))
      if(is_long) {
        survwide <- matrix(survlong, ncol = length(thresholds))
        H_Y <- (do.call(cbind, lapply(1:ncol(a_mat), function(i) {
          a <- a_mat[,i]
          surv <- survwide[,i]
          return(a/surv)
        })))
      } else {
        Hlong <- a_long / survlong
        H_Y <- matrix(Hlong, ncol = length(thresholds))

      }



      # H_Y <- as.matrix(do.call(cbind, lapply(thresholds, function(thresh) {
      #   ind <- as.numeric(A >= thresh)
      #   index <- which.min( abs(A_grid - thresh))
      #
      #   surv <- 1 - cdf[, .SD[index] , by = id][[2]]
      #   H_Y <- ind/surv
      # })))


      EIC_Y <- NULL
      EIC_A <- NULL
      if(for_fitting) {
        EIC_Y <- H_Y * (Y - Q)
        EIC_A <- IC_A
      }
      return(list(Y = H_Y, A = IC_A, psi_W = psi_W,
                  IC = list(Y = EIC_Y, A = EIC_A)))
    },
    make_long_task = function(tmle_task) {
      if(tmle_task$uuid %in% private$.long_task_uuids) {
        return(tmle_task)
      }
      key <- tmle_task$uuid
      cached_task <- get0(key, private$.cache, inherits = FALSE)
      if(!is.null(cached_task)){
        return(cached_task)
      }

      data <- tmle_task$data
      nodes <- tmle_task$nodes
      A_grid <- private$.A_grid
      var <- tmle_task$npsem[["A"]]$variables
      new_data <- rbindlist(lapply(A_grid, function(a) {
        data <- copy(data)
        set(data, , var, a)
        return(data)
      }))
      new_data$trueid <- new_data$id
      new_data$id <- sort(as.factor(seq_len(nrow(new_data))))
      long_task <- tmle3_Task$new(new_data, tmle_task$npsem, id = "id", time = "t", nodes = nodes)
      assign(key, long_task, private$.cache)
      private$.long_task_uuids <- c(private$.long_task_uuids, long_task$uuid)
      return(long_task)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      thresh_node <- private$.thresh_node
      censoring_node <- private$.censoring_node
      thresholds <- private$.thresholds

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
      return(c("Y", "A"))
    }
  ),
  private = list(
    .type = "Threshold",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .supports_outcome_censoring = TRUE,
    .submodel_type_supported = list("Y" = "logistic", "A" = "EIC"),
    .supports_weights = T,
    .censoring_node = NULL,
    .thresh_node = NULL,
    .thresholds = NULL,
    .strict_threshold = NULL,
    .cf_task = NULL,
    .type_est = NULL,
    .A_grid = NULL,
    .cache = NULL,
    .long_task_uuids = c()

  )
)
