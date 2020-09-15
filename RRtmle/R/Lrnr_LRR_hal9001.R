Lrnr_LRR_hal9001 <- R6Class(
  classname = "Lrnr_LRR_hal9001", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 3,
                          smoothness_degree = 0,
                          iter = 250,
                          method = NULL,
                          grad_desc = F,
                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("binomial"),

    .train = function(task) {
      method <- self$params$method
      X <- task$X

      Y <- task$Y
      if(!is.null(method)) {
        weights <- task$get_data(,c("weightsIPW", "weightsplugin"))
        if(method == "IPW") {
          weights <- weights[, weightsIPW]
          Y <- Y[, YIPW]
        }
        else if(method == "plugin") {
          weights <- weights[, weightsplugin]
          Y <- Y[, Yplugin]
        }
      } else {
        weights <- task$get_data(,"weights")
      }
      weights <- weights * task$weights
      smoothness <- self$params$smoothness_degree
      smoothness <- rep(c(smoothness)[[1]], ncol(X))
      basis_list <- hal9001fast::enumerate_basis(as.matrix(X), max_degree = self$params$max_degree, order_map = c(smoothness), bins = rep(350, ncol(X)), include_zero_order = F, include_lower_order = F)
      x_basis <- hal9001fast::make_design_matrix(as.matrix(X), basis_list)


      fit <- glmnet::cv.glmnet(x_basis, Y, family = binomial(), weights = weights, standardize = F, nlambda = 125)
      fit_object <- list(basis_list = basis_list, beta = as.vector(coef(fit)))
      keep <- which(fit_object$beta[-1]!=0)
      basis_list <- basis_list[keep]
      beta <- fit_object$beta[c(1, keep+1)]
      fit_object$basis_list <- basis_list
      fit_object$beta <- beta
      private$.beta <- beta
      private$.basis_list <- basis_list
      x_basis <- x_basis[,keep]
      if(!self$params$grad_desc){
        return(fit_object)
      }
      search_eps = T
      for(i in 1:self$params$iter) {

        out <- private$.gradient_descent_update(task, x_basis, search_eps)
        search_eps = F

        if(out == "converged") {
          break
        }
      }
      fit_object$beta <- private$.beta

      return(fit_object)
    },
    .predict = function(task = NULL) {
      fit_obj <- self$fit_object
      basis_list <- fit_obj$basis_list
      beta <- fit_obj$beta
      X <- task$X
      x_basis <- hal9001fast::make_design_matrix(as.matrix(X), basis_list)
      predictions <- x_basis %*% beta[-1] + beta[1]
      return(predictions)
    },
    .gradient_descent_update = function(task, x_basis, search_eps = F) {
      eps <- private$.eps

      if(!is.null(eps)) {

        if(abs(eps) < 1e-7){

          return("converged")
        }

      }



      beta <- private$.beta

      g <- task$get_data(,"g")[[1]]
      ER <- task$get_data(,"Q")[[1]]
      R <- task$get_data(,"R")[[1]]
      A <- task$get_data(,"A")[[1]]
      ER0 <- task$get_data(,"ER0")[[1]]
      ER1 <- task$get_data(,"ER1")[[1]]

      LRR <- as.vector(x_basis %*% beta[-1]) + beta[1]
      C1 <- A/g * (R - ER) + ER1
      C2 <- C1 + (1-A)/g * (R - ER) + ER0
      Z <- -C1 + C2 * exp(LRR) / (1 + exp(LRR))
      risk  <- function(beta) {
        f <- as.vector(x_basis %*% beta[-1]) + beta[1]
        loss = C1*-1*f + C2 * log(1 + exp(f))
        return(mean(loss))
      }

      intercept <- beta[1]
      beta_free <- beta[-1]
      D <-  c(mean(Z) * intercept ,  as.vector(Matrix::crossprod(x_basis , Z/ length(Z)) * beta_free))

      D_star <- D - (sum(D * abs(beta))/ sum(beta^2))*abs(beta)
      norm <-   (norm(D_star, type = "2") / sqrt(length(D_star)))

      if(norm < 1e-5) {

        return("converged")
      }
      D_star <- D_star / norm
      cur_risk <- risk(beta)
      eps <- private$.eps

      eps_max <- 1/max(abs(D_star))

      if(!is.null(eps)) {
        new_beta <- (1+eps*D_star) * beta
        new_risk <- risk(new_beta)
        better <- new_risk <= cur_risk
        if(!better){

        }
      } else {
        better <- T

      }
      if(is.null(eps) | search_eps | !better) {
        if(is.null(eps)) {
          eps <- eps_max/2
        }

        search_set <- -1 * c(exp(seq(log(2*abs(eps)), log(.000000001), length.out = 20)))
        #search_set <- c( -1 * search_set, search_set)
        search_set <- c(0.00005, search_set)
        risks <- c()
        for (ep in search_set) {
          new_beta <- (1+ep*D_star) * beta
          risks <- c(risks, risk(new_beta))
        }
        eps <- search_set[which.min(risks)]
        private$.eps <- eps
        new_risk <- min(risks)
        print(new_risk)
        if(new_risk > cur_risk) {
          warning("fit not improved in gradient update")
          eps <- 0
        }
        new_beta <- (1+eps*D_star) * beta
      }
      private$.beta <- new_beta
      return("Not converged")

    },
    .eps = NULL,
    .basis_list = NULL,
    .beta = NULL,
    .required_packages = c("hal9001fast")
  )
)
