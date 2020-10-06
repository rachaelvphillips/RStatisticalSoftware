
make_learner_subsemble_hal <- function(lrnr, task, k = 5, num_preds){
  assertthat::assert_that("sl3_hal_Task" %in% class(task))
  meta_lrnr = make_learner(Lrnr_hal9001_subsemble_collapser)

  n <- task$nrow
  folds <- compute_subsemble_folds(n,k)
  strata_ids <- attr(folds , 'strata_ids')
  lrnr_sub <- make_learner(Lrnr_subsemble_hal, lrnr = lrnr, strata_ids = strata_ids)
  meta_lrnr <- make_learner(Lrnr_subsemble_hal_metalearner, meta_lrnr = meta_lrnr, num_groups = k, num_preds = num_preds)
  final_lrnr <- make_learner(Pipeline, lrnr_sub, meta_lrnr)
  new_task <- task$next_in_chain(folds = folds)
  return(list(subsemble_lrnr = final_lrnr, cv_task = new_task))
}

#' @importFrom R6 R6Class
#' @import sl3
#' @export
sl3_hal_Task <- R6::R6Class(
  classname = "sl3_hal_Task",
  portable = TRUE,
  class = TRUE,
  inherit = sl3_Task,
  public = list(
  add_to_bucket = function(item, name){
    private$.bucket[[name]] <- item
    return(invisible(self))
  }
  ),
  active = list(
    bucket = function(){
      private$.bucket
    }
  ),
  private = list(
    .bucket = list()
  )
)
#' Custom subsemble meta learner for hal
#' @importFrom R6 R6Class
#' @import sl3
#' @export
Lrnr_hal9001_subsemble_collapser <- R6::R6Class(
  classname = "Lrnr_hal9001_collapser", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(convex = T,...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial"),

    .train = function(task) {
      orig_X = as.matrix(task$bucket[["task"]]$X)

    args <- self$params

      x <- as.matrix(task$X)
      y <- task$Y
      if(task$outcome_type$glm_family() == "gaussian"){
        family = gaussian()
        nnls_fit <- nnls::nnls(as.matrix(x), y)
        if (args$convex == TRUE) {
          init_coef <- coefficients(nnls_fit)
          init_coef[is.na(init_coef)] <- 0
          if (sum(init_coef) > 0) {
            coef <- init_coef / sum(init_coef)
          } else {
            warning("All algorithms have zero weight", call. = FALSE)
            coef <- init_coef
          }

        } else{
          coef <- coefficients(nnls_fit)
        }
      } else if(task$outcome_type$glm_family() == "binomial"){
        family = binomial()
        coef <- coef(cv.glmnet(as.matrix(x), y, lower.limits = 0, family = "binomial"), s = "lambda.min")
        if (args$convex == TRUE) {
          if(sum(coef)>0){
            coef <- coef / sum(coef)
          }
        }
        }
      coef_meta <- as.vector(coef)
      print(coef_meta)
      throw_away = which(coef_meta==0)




      hal_fits = lapply(task$bucket[["hal_fits"]], function(item){
        item
      })
      if(length(throw_away)>0){
        print("heress")
        hal_fits = hal_fits[-throw_away]
        coef_meta = coef_meta[-throw_away]
      }


      index = task$bucket[["lambda_index"]]
      basis_betas <- lapply(hal_fits, function(fit){
       blist = fit$basis_list
       coefs = as.vector(fit$coefs[,index])
       a0 = coefs[1]
       betas = coefs[-1]

       keep = which(betas !=0)
       return(list(basis_list = blist[keep], a0 = a0, coefs = betas[keep]))
      })

      # Combined basis
      basis_list = unlist(lapply(basis_betas, function(item){
        return(item$basis_list)
      }), recursive = F)
      # Stacked beta vector
      coefs_stacked = unlist(lapply(1:length(basis_betas), function(i){
        item = basis_betas[[i]]
        coefs = item$coefs * coef_meta[i]
        return(coefs)
      }))

      a0s = sapply(1:length(basis_betas), function(i){
        item = basis_betas[[i]]
        item$a0 * coef_meta[i]
      })
      intercept = sum(a0s)
      x_basis = make_design_matrix(orig_X, basis_list)
      # Collapse duplicates
      copy_map = make_copy_map(x_basis)
      uniq_cols = as.numeric(names(copy_map))
      x_basis = x_basis[,uniq_cols]
      basis_list = basis_list[uniq_cols]
      beta = sapply(1:length(copy_map), function(i){
        dupes = copy_map[[i]]
        beta = sum(coefs_stacked[dupes])
        return(beta)
      })
      #x_basis = cbind(rep(1, nrow(x_basis)), x_basis)
      #beta = c(intercept, beta)
    #fit <- glmnet(x_basis,y,  lambda = hal_fits[[1]]$lambda_star[index], standardize= F, intercept = F, penalty.factor = c(0,rep(1,ncol(x_basis)-1)))
      # Get gradient descent updated beta (logistic not implemented yet)
    #new_beta =  as.vector(fit$beta)

      #new_beta = run_descent_gaussian(x_basis, y, beta, max_iter = 500, verbose = T)
      # print(data.frame(cbind(beta, new_beta)))
      # print(data.frame(cbind(sum(abs(beta)), sum(abs(new_beta)))))
      # print(data.frame(cbind(mean((y-as.vector(x_basis%*%beta))^2), mean((y-as.vector(x_basis%*%new_beta))^2))))
      fit_object = list(hal_fits = hal_fits, intercept = intercept, beta = beta, basis_list = basis_list)

      return(fit_object)
    },
    .predict = function(task = NULL) {
      fit_obj = self$fit_object
      beta = fit_obj$beta
      intercept = fit_obj$intercept
      basis_list = fit_obj$basis_list
      orig_X = as.matrix(task$bucket[["task"]]$X)

      x_basis = make_design_matrix(as.matrix(orig_X), basis_list)
      link = intercept + as.vector(x_basis %*% beta)
      if(task$outcome_type$glm_family() == "gaussian"){
        return(link)
      }
      else if(task$outcome_type$glm_family() == "binomial"){
        return(stats::plogis(link))
      }
      return(link)
    },
    .required_packages = c("hal9001")
  )
)
