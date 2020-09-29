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
    initialize = function(observed_likelihood, thresh_node = "A", outcome_node = "Y", type = 1, thresholds = NULL, num_bins = 200, discretize_type = c("mix", "equal_mass", "equal_range"), discretize_g = T) {

      discretize_type <- match.arg(discretize_type)
      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)



      range <- (observed_likelihood$training_task$get_tmle_node("A"))
      variable_type <- observed_likelihood$training_task$npsem$A$variable_type$type
      if(variable_type!="continuous") {
        A_grid <- sort(unique(range))
        thresholds <- A_grid
      }
      else {
      #A_grid <- sl3:::make_bins(range, "equal_mass", num_bins)
        if(discretize_type == "mix") {
          A_grid1 <- sl3:::make_bins(range, "equal_mass", round(num_bins * 0.85))
          min_diff1 <- min(abs(diff(A_grid1)))
          A_grid2 <- sl3:::make_bins(range, "equal_range", round(num_bins * 0.3))
          min_diff2 <- min(abs(diff(A_grid2)))
          A_grid <- sort(union(A_grid1, A_grid2))
          min_diff <- min(min_diff2, min_diff1)
          print(length(A_grid))
          remove <- setdiff(which(diff(A_grid) <min_diff) +1, c(length(A_grid)) )
          remove <- intersect(remove, seq(2, length(A_grid), 2))
          print(paste0("Removed: ", length(remove)))

          A_grid <- A_grid[-remove]
        } else {
          A_grid <- sl3:::make_bins(range, discretize_type, num_bins)
        }
      #A_grid <- seq(0, 1, length.out = num_bins)
      #A_grid <- sort(unique(quantile(range, A_grid, type = 3)))
      private$.censoring_node <- (observed_likelihood$censoring_nodes[[outcome_node]])
      private$.thresh_node <- thresh_node
      if(is.null(thresholds)) {
        if(length(A_grid) > 10) {
          thresholds <- quantile(range, seq(0.05, 0.95, length.out = 10))
          thresholds <- A_grid[findInterval(thresholds, A_grid)]
          #thresholds <- c(A_grid[190])

          } else {
          thresholds <- A_grid
        }
      }

      #print(thresholds %in% A_grid)

      A_grid <- sort(unique(union(A_grid, thresholds)))

      thresholds <- sort(unique(thresholds))
}

      private$.A_grid <- A_grid
      private$.thresholds <- thresholds
      private$.strict_threshold <- F
      private$.type_est <- type
      private$.cache <- new.env()
      private$.discretize_g <- discretize_g
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", for_fitting = F, node = NULL) {
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
      Q  <- as.vector(likelihood$get_likelihood(tmle_task, "Y"))
      g  <- as.vector(likelihood$get_likelihood(tmle_task, "A"))
      Y <- tmle_task$get_tmle_node("Y")
      A <- tmle_task$get_tmle_node("A")
      variable_type <- tmle_task$npsem$A$variable_type$type

      long_task <- self$make_long_task(tmle_task)

      is_long <- F
      if(long_task$uuid == tmle_task$uuid) {
        is_long <- T
      }

      g_long <- as.vector(likelihood$get_likelihood(long_task, "A"))
      if(any(g_long <0)) {
        stop("Negative density epsilon to big")
      }
      Q_long <- as.vector(likelihood$get_likelihood(long_task, "Y"))

      long_g_id <- data.table(A = long_task$get_tmle_node("A", include_id = F), g = g_long, id = NA)
      set(long_g_id, , "id", long_task$get_data(,"trueid"))
      setnames(long_g_id, c("A", "g", "id"))

      Q_long_id <- data.table(Q = unlist(Q_long), id = NA)

      set(Q_long_id,,  "id", long_task$get_data(,"trueid"))
      setnames(Q_long_id, c("Q", "id"))


      # Store common information based on long_task (only recompute when fitting)
      if(for_fitting){
        private$.cache_shared <- list()
      }
      cdf <- private$.cache_shared$cdf
      psi_W <- private$.cache_shared$psi_W
      integral <- private$.cache_shared$integral
      if(!is.null(cdf)) {
        #print(data.table(matrix(cdf, nrow = 1000, byrow = F)))

      }
      if(is.null(cdf)) {
        if(variable_type == "continuous") {
          A_diff <- long_g_id[,  diff(A), by = id]
          #\g_int <- na.omit(long_g_id[, frollsum(.SD, 2, na.rm=T, align = "right"), by = id,.SDcols = "g"])[[2]] /2
          g_int <- (long_g_id[, .SD[-nrow(.SD)], by = id, .SDcols = "g"])[[2]]

          #TODO do we want our density/mean and integrals to be exact?
          if(private$.discretize_g == T && for_fitting){
            g_discrete <- matrix(g_int, nrow = length(unique(A_diff$id)), byrow = T)
            match_index <- findInterval(A, A_grid, all.inside = T)
            index_mat <- cbind(1:length(A), match_index)
            g_discrete <- g_discrete[index_mat]
            print("risk")
            print(mean(-log(g_discrete)))
            print(mean(-log(g)))
            g <- g_discrete
          }
          #g_int <- na.omit(long_g_id[, .SD[-.N], by = id,.SDcols = "g"])[[2]]
          #Q_int <- na.omit(Q_long_id[, frollsum(.SD, 2, na.rm=T, align = "right"), by = id, .SDcols = "Q"])[[2]] /2
          Q_int <- (Q_long_id[, .SD[-nrow(.SD)], by = id, .SDcols = "Q"])[[2]]

          A_int <- long_g_id[, A, by = id][[2]]
          integral <- data.table(A_diff,  g = g_int, Q = Q_int)
          setnames(integral, c("id", "A_diff", "g", "Q"))
          integral[, A_g := A_diff*g]

          #cdf <- integral[, c(0,cumsum(A_diff*g)), by = id]
          cdf <- integral[, c(0,cumsum(A_g)), by = id]
          #cdf <- integral[,  c(rev(cumsum(rev(A_g ))),0), by = id]

          #old_cdf <- (integral[,  c(0,cumsum(A_g )), by = id])
          #setnames(cdf, c("id", "int"))
          #cdf <- cdf[,  int[1] - int, by = id]

          norm_C <- cdf[, last(.SD), by = id][[2]]

          if(any(abs(norm_C - 1) > 1e-6 )) {
            print(sum(norm_C > 1+1e-9))
            print(quantile(norm_C))

            warning("Density no longer integrates to 1")
          }

          #norm_C <- 1
          cdf <- cdf[[2]]
          cdf <- as.vector(matrix(cdf, ncol = length(A_grid), byrow = T))
          if(T & for_fitting & any(norm_C > 1)) {
            long_g_id$g <- long_g_id$g/norm_C
            g <- g/norm_C
            integral$g <- as.vector(t(matrix(integral$g,  ncol = length(A_grid)-1, byrow = T)/norm_C))
            cdf <- cdf /norm_C
          }

          psi_W <- integral[,  c(rev(cumsum(rev(A_g * Q ))),0), by = id]
          psi_W_other <- integral[,  c(0, cumsum(Q*A_g  )), by = id]
          setnames(psi_W_other, c("id", "int"))

          psi_W_other <- psi_W_other[,  int[.N] - int, by = id]
          psi_W <- psi_W_other
          psi_W$A <- A_int
          setnames(psi_W, c("id", "int", "A"))

          private$.cache_shared$cdf <- cdf
          private$.cache_shared$psi_W <- psi_W
          private$.cache_shared$integral <- integral
          likelihood_factor <- likelihood$factor_list[["A"]]
          step_number <- likelihood$cache$get_update_step(likelihood_factor, task, fold_number, node = "A")
          step_number_long <- likelihood$cache$get_update_step(likelihood_factor, long_task, fold_number, node = "A")

          if(any(norm_C > 1) ) {
            likelihood$cache$set_values(likelihood_factor, task, step_number, fold_number, g, node = "A")
            likelihood$cache$set_values(likelihood_factor, long_task, step_number_long, fold_number, long_g_id$g, node = "A")
          }


        } else {
          #A_int <- long_g_id[, A, by = id][[2]]
          A_diff <- 1
          integral <- data.table(long_g_id$id, A = A_diff,  g = long_g_id$g, Q = Q_long_id$Q)
          setnames(integral, c("id", "A_diff", "g", "Q"))
          integral[, A_g := g]
          cdf <- integral[, c(0, cumsum(g)[-length(g)]), by = id]
          norm_C <- integral[,  sum(g), by = id][[2]]
          print(quantile(norm_C))
          cdf <- cdf[[2]]
          cdf <- as.vector(matrix(cdf, ncol = length(A_grid), byrow = T))

          psi_W <- integral[,  rev(cumsum(rev( g * Q ))), by = id]
          #psi_W$A <- A_int
          setnames(psi_W, c("id", "int"))
          integral <- integral[,.SD, by = id]
        }

      }



      trueid <- Q_long_id$id

      #For each id, CDF values at A_grid

      cache_calc <- private$.cache_calcs
      if(is_long) {
        indices <- private$.cache_long$indices
        a_mat <- private$.cache_long$a_mat
      } else {
        indices <- private$.cache_short$indices
        a_mat <- private$.cache_short$a_mat
      }
      if(is.null(indices)) {
        indices <- sapply(thresholds, function(a) {
          which.min( abs(A_grid - a))
        })
        a_mat <- do.call(cbind, lapply(thresholds, function(a) {
          rep(a, length(unique(tmle_task$id)))
        }))
        if(is_long) {
          private$.cache_long$indices <- indices
          private$.cache_long$a_mat <- a_mat
        } else {
          private$.cache_short$indices <- indices
          private$.cache_short$a_mat <- a_mat
        }

      }


      a_mat <- a_mat <=A
      a_long <- as.vector(a_mat)

      #print(data.table(matrix(cdfb, ncol = length(A_grid), byrow = T)))
      psi_Wb <- psi_W[, int]
      cdfb <- as.vector(matrix(cdf, ncol = length(A_grid), byrow = F)[,indices])
      psi_Wb <- as.vector(matrix(psi_Wb, ncol = length(A_grid), byrow = T)[,indices])
      #print(data.table(matrix(psi_Wb, ncol = length(thresholds))))
      survlong <-  1 - cdfb
      if(min(survlong) < 0.0005) {
        warning("very small survival prob")
      }
      survlong <- bound(survlong, c(0.0005,10))

      psi_Wlong <- psi_Wb / survlong

      if(is.null(node) || (!is.null(node) & node == "A")) {
      if(is_long) {
        Qwide <- matrix(Q, nrow = length(unique(trueid)))
        amat <- a_mat

        IC_A <- do.call(rbind, lapply(1:ncol(Qwide), function(i) {
          Q <- Qwide[,i]
          result <- matrix((Q - psi_Wlong )  / survlong, ncol =length(thresholds))
          return(result)
        }))

        IC_A <- IC_A * a_mat
       # print(max(abs(IC_A)))
      #  print("IC")

        for(i in 1:ncol(IC_A)) {
        IC_id <- data.table(IC = IC_A[,i], id = long_g_id$id)
        IC_new <- IC_id[, IC[-.N], by = id]
        #IC_new <- IC_id[, IC, by = id]
        setnames(IC_new, c("id", "IC"))
        IC_new$A_g <- integral$A_g

        center <-  IC_new[, sum(IC* A_g), by = id][[2]]
       # print(max(center))
        #IC_A[,i] <- IC_A[,i] - center
        }


      } else {

        IC_long <- (as.vector(Q) - as.vector(psi_Wlong) ) * a_long / survlong
        IC_A <- matrix(IC_long, ncol = length(thresholds))

      }

      } else {
        IC_A <- NULL
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
      if(is.null(node) || (!is.null(node) & node == "Y")) {
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
      } else {
        H_Y <- NULL
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
        if(is.null(node) | node == "Y") {
          EIC_Y <- H_Y * (Y - Q)
        }
        if(is.null(node) | node == "A") {
          EIC_A <- IC_A
        }

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
      H <- self$clever_covariates(task, fold_number, for_fitting = T)
      IC <- H$IC
      IC_A <- IC$A
      IC_Y <- IC$Y
      psi_W <- H$psi_W
      psi <- colMeans(psi_W)
      IC_W <- t(t(psi_W) - psi)
      IC <- IC_W + IC_A + IC_Y
      result <- list(psi = psi, IC = IC)
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
    .cache_short = list(),
    .cache_long = list(),
    .cache_shared = list(),
    .long_task_uuids = c(),
    .discretize_g = NULL

  )
)
