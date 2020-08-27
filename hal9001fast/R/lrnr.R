

#' @importFrom R6 R6Class
#' @import sl3
Lrnr_hal9001fast <- R6::R6Class(
  classname = "Lrnr_hal9001fast", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 3,
                          formula = NULL,
                          reduce_basis = function(n){1/sqrt(n)},
                          keep_min = 3,
                          dcor_pval = 0.25,
                          cor_pval = 0.1,
                          smoothness_orders = 0,
                          bins = function(n){max(50,min(n/10,500))},
                          screen_basis_main_terms = F,
                          screen_basis_interactions = F,
                          lambda = NULL,
                          cv_select=F,
                          keep_cov = NULL,
                          max_num_two_way = 75000,
                          max_total_basis = 250000,
                          link_as_pred = F,
                          verbose = F,


                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial"),

    .train = function(task) {
      params <- self$params
      X <- as.matrix(task$X)
      Y <- as.vector(task$Y)
      formula <- params$formula
      smoothness_orders <- params$smoothness_orders
      if(length(smoothness_orders)>1){
        smoothness_orders <- smoothness_orders[match(colnames(X), names(smoothness_orders))]

      }
      bins <- params$bins
      if(is.function(bins)){
        bins <- round(bins(nrow(X)))
      }

      reduce_basis <- params$reduce_basis
      cv_select<- self$params$cv_select
      if(is.function(reduce_basis)){
        reduce_basis <- reduce_basis(nrow(task$X))
      }

      # Perform nonparametric variable screening by independence
      if(!is.null(params$dcor_pval)){
        if(nrow(task$X)<=200){
          params$dcor_pval = 0.5
        }
        else if(nrow(task$X)<=400){
            params$dcor_pval = 0.2
          }
        pvals <- apply(X,2, function(x,y){energy::dcorT.test(x,y)$p.value}, y = Y)

        remove <- which(pvals > params$dcor_pval)
        if(ncol(X)-length(remove)<params$keep_min){
          remove = which(rank(-pvals) <= (ncol(X) - params$keep_min))
        }
        if(!is.null(params$keep_cov)){
          remove = setdiff(remove, match(params$keep_cov,colnames(X)))
        }

        print(paste0("Number of covariates removed: ", length(remove)))
        #Effectively remove columns
        if(length(smoothness_orders)>1){
          smoothness_orders = smoothness_orders[-remove]

        }
        if(length(remove)>0){
          X = X[,-remove,drop=F]

        }
      }
      else{
        remove = NULL
      }
      print(dim(X))
      # perform hal fit
      if(!is.null(formula)){
        #Use formula spec if available
        data <- cbind(X,Y)
        colnames(data) = c(colnames(X), task$nodes$outcome)
        formula <- formula_hal(formula$formula, data, smoothness_orders, bins, custom_group = formula$custom_group)

        fit <- fit_halfast(formula = formula, screen_cor_pval = params$cor_pval, family = task$outcome_type$glm_family(),  lambda = params$lambda,
                           cv_select = cv_select,reduce_basis = reduce_basis,
                           max_num_two_way = params$max_num_two_way,
                           max_total_basis = params$max_total_basis, verbose = self$params$verbose)
      }
      else{

        fit <- fit_halfast(X = X, Y=Y, max_degree = params$max_degree, smoothness_orders=smoothness_orders, screen_cor_pval = params$cor_pval, family = task$outcome_type$glm_family(),
                           screen_basis_main_terms = params$screen_basis_main_terms,
                           screen_basis_interactions = params$screen_basis_interactions,
                           lambda = params$lambda, reduce_basis = reduce_basis, cv_select = cv_select, num_bins = bins,
                           max_num_two_way = params$max_num_two_way,
                           max_total_basis = params$max_total_basis,verbose = self$params$verbose)

      }





      fit_object <- list()
      fit_object$remove = remove
      fit_object$fit = fit
      return(fit_object)
    },
    .predict = function(task = NULL) {
      X = as.matrix(task$X)
      fit <- self$fit_object$fit
      remove <- self$fit_object$remove

      if(!is.null(remove) & length(remove)>0){
        X = as.matrix(X)
        X = X[,-remove]
      }
      if(self$params$link_as_pred){
        predictions = predict(fit$glmnet_lasso, newx = make_design_matrix(as.matrix(X), fit$basis_list), type = "link")
      } else {
        predictions = predict(fit, new_data = as.matrix(X))
      }

      return(predictions)
    },
    .required_packages = c("hal9001fast")
  )
)
