
#' @import  purrr
#' @import  delayed
#' @import  future
#' @import  SIS
#' @import  delayed
#' @import  future
#' @import  sl3
#' @import  glmnet
#' @import  hal9001
#' @import  biglasso
#' @import  bigstatsr
#' @import  bigmemory
#' @import  qvalue
#' @import  sl3
#' @export
#'
#'


HALnet <- R6Class("HALnet", private = list(
  .params = NULL,
  .learners_spec = NULL,



  outcome_type = function(){
    return(self$learners_spec$trainingTask$outcome_type$type)
  },
  list_to_stack = function(lst){
    if(is.null(lst)){
      return()
    }
    lst=c(lst)
    if(length(lst)==1){
      return(lst[[1]])
    }
    return(do.call(function(...){
      return(make_learner(Stack, ...))
    }, lst))

  },
  set_hal_lrnrs = function(bins = 150){
    enforce_meta_monotone = self$params$enforce_meta_monotone
    if(enforce_meta_monotone){
      lower.limits = 0
    }
    else{
      lower.limits = -Inf
    }
    if(!self$params$include_meta_learner){

      print("Not including meta learners")
      private$.learners_spec$meta_learners = NULL


      return()
    }
    else if(toupper(self$params$which_meta_learner) %in% c("glm", "GLM")){
      private$.learners_spec$meta_learners = list(make_learner(Lrnr_glm))
      return()
    }



    lambda_path = self$learners_spec$lambda_path
    hal_matrix_learner = make_learner(Lrnr_hal9001, max_degree = 2, lambda = self$learners_spec$lambda_path, yolo=F,cv_select=F,lower.limits=lower.limits)
    hal_matrix_learner2 = make_learner(Lrnr_hal9001, max_degree = 2, lambda = self$learners_spec$lambda_path, yolo=F,cv_select=F)
    hal_matrix_learner3 = make_learner(Lrnr_hal9001, max_degree = 1, lambda = self$learners_spec$lambda_path, yolo=F,cv_select=F,lower.limits=lower.limits)

    # make_lambda_selector = function(index){
    #   make_learner(Lrnr_column_selector, column_index=index)
    # }

    # lambda_selectors_stack = private$list_to_stack(lapply(1:length(lambda_path), make_lambda_selector))
    # hal_learner = make_learner(Pipeline, hal_matrix_learner, lambda_selectors_stack)
      #The above isnt needed if we use Stackfixed.R
     hal_learner = c(hal_matrix_learner, hal_matrix_learner3)
    cv_hal_learner =  make_learner(Lrnr_hal9001, max_degree = 2, yolo=F,cv_select=T,lower.limits=lower.limits)
    hal_learners = make_learner(Pipeline, make_learner(Lrnr_discretizer, bins = bins), private$list_to_stack(c(hal_learner, cv_hal_learner)))


    Halgrad_lrnr = make_learner(Lrnr_HALgrad)
    private$.learners_spec$meta_learners = list(hal_learners, Halgrad_lrnr,make_learner(Lrnr_glm))


  },

  set_baseline_learners_enet = function(){
    learners = list()
    for(alpha in self$params$alpha){
      learners = c(learners, make_learner(BigLrnr_biglasso, alpha = alpha))
    }

    private$.learners_spec$base_learners= learners

  },
  set_screeners = function(){
    screeners_to_include = c(self$params$screeners_to_include)
    screeners = list()
    if("All" %in% screeners_to_include){
      screeners = c(screeners, "All" = "All")
    }

    if("screen.biglasso" %in% screeners_to_include){
      print("hi")
      for(alpha in c(self$params$screen.biglasso_alpha)){
        print("by")
        screeners = c(screeners, "screen.biglasso" =make_learner(BigLrnr_screener_biglasso, alpha = alpha, dfmax = self$params$screen.max_selected))
        print(screeners)
      }


    }
    if("screen.limma" %in% screeners_to_include){

      for(cutoff in c(self$params$screen.pval.fdr_cutoff)){
        screeners=c(screeners, "screen.limma" = make_learner(BigLrnr_screener_limma, cutoff = cutoff, max_num = self$params$screen.max_selected))
      }

    }
    if("screen.corP" %in% screeners_to_include){
      for(cutoff in c(self$params$screen.pval.fdr_cutoff)){
        screeners = c(screeners, "screen.corP" = make_learner(BigLrnr_screener_corP, minPvalue_fdr = cutoff, maxscreen = self$params$screen.max_selected))
      }

    }

    if("screen.randomPCA" %in% screeners_to_include){
      screeners = c(screeners, "screen.randomPCA" = make_learner(BigLrnr_randomPCA, K = self$params$screen.pca_dim))
    }


    private$.learners_spec$screeners = screeners
  },


  set_group_subsetters = function(grps, group_PCA_dim){
    num_groups = length(grps)
    name_groups = names(grps)
    group_subsetters = list()
    for(grp in grps){
      lrnr_grp = make_learner(BigLrnr_subset_covariates, bm_columns = as.data.frame(grp))
      if(is.numeric(group_PCA_dim)){
        lrnr_pca = make_learner(BigLrnr_randomPCA, group_dim = group_PCA_dim)
        lrnr_grp = make_learner(Pipeline,lrnr_grp ,lrnr_pca)
      }

      group_subsetters = c(group_subsetters, lrnr_grp)
    }
    private$.learners_spec$group_subsetters = group_subsetters
  },
  buildLibrary = function(){
    screeners = self$learners_spec$screeners
    library=list()
    for(name in names(screeners)){
      scrn = screeners[[name]]
      if(name=="All"){
        #Add learners with no screening
        library = c(library, private$getElement(NULL))

      }

      else{

        library = c(library, private$getElement(scrn) )
      }

    }

    library_stack = private$list_to_stack(library)
    return(library_stack)
  },
  buildSuperLearner = function(){

    sl3_library = private$buildLibrary()
    sl3_superlearner = make_learner(Pipeline, Lrnr_cv$new(sl3_library), make_learner(Lrnr_cv_selector))
    if(!is.null(self$learners_spec$subset_learner)){
      sl3_superlearner = make_learner(Pipeline, self$learners_spec$subset_learner, sl3_superlearner)
      sl3_library = make_learner(Pipeline, self$learners_spec$subset_learner, sl3_library)
    }
    private$.learners_spec$sl3_superlearner = sl3_superlearner
    private$.learners_spec$sl3_library = sl3_library

  },
  getElement = function(screener){
    meta_learner = private$list_to_stack(self$learners_spec$meta_learners)
    baselearners = self$learners_spec$base_learners
    group_subsetters = self$learners_spec$group_subsetters
    base = list()


    for(lrnr in baselearners){
      tmp = list()
      print(group_subsetters)
      if(length(group_subsetters)==0){
        tmp = lrnr
      }
      else{
        for(grp in group_subsetters){
          tmp = c(tmp, make_learner(Pipeline, grp, lrnr))
        }
      }


      tmp = private$list_to_stack(c(tmp))

      if(!is.null(meta_learner)){
        base = c(base, make_learner(Pipeline,tmp , meta_learner))
      }
      else{
        base = c(base, tmp)
      }

    }

    base = private$list_to_stack(base)

    if(is.null(screener)){

      return(base)
    }
    else{
      element = make_learner(Pipeline, screener, base)
    }

    return(element)

  },

  set_delayed_learner = function(nworkers){
    plan(multiprocess, workers = nworkers)
    private$.learners_spec$sl3_delayed_learner = delayed_learner_train(self$learners_spec$sl3_superlearner, self$learners_spec$trainingTask)
  },
  setTasks = function(){
    hal_task = self$params$methylTask

    private$.learners_spec$testTask = hal_task$testTask

    private$.learners_spec$trainingTask = hal_task$trainTask
  },
  setLambdaPath = function(lambda_max, eps = 0.01, K = 100){
    if(is.null(lambda_max)){
      lambda_max=1
    }

    private$.learners_spec$lambda_path = round(exp(seq(log(lambda_max)+3, log(lambda_max*eps),
                                                   length.out = K)), digits = 10)

  }




),
active = list(
  params = function(){
    return(private$.params)
  },
  groups = function(){
    return(self$params$methylTask$groups)
  },
  lambda_path = function(){
    return(self$learners_spec$lambda_path)
  },
  learners_spec = function(){
    return(private$.learners_spec)
  },
  superlearner = function(){
    return(self$learners_spec$sl3_superlearner)
  },
  library = function(){
    return(self$learners_spec$sl3_library)
  },
  testTask = function(){
    return(self$learners_spec$testTask)
  },
  trainingTask = function(){
    return(self$learners_spec$trainingTask)
  }
),
public = list(
  initialize =function(methylTask = NULL, screeners_to_include = c("screen.corP", "screen.biglasso"),   baseline_learner_choice = c("glmnet"), include_meta_learner = T, which_meta_learner = "HAL",
                       screen.pval.fdr_cutoff = c(0.1), screen.biglasso_alpha = c(0.5), screen.max_selected = 300000, screen.pca_dim = 1000,
                       include_groups=T, group_PCA_dim = NULL,
                     baseline_glmnet_lambda_choice = c("lambda.min"), alpha = c(1),
                     enforce_meta_monotone = T, eps = 1e-5, K = 200,
                     nworkers = 1){


    screeners_to_include = c(screeners_to_include)
    private$.params = args_to_list()
    private$.learners_spec = list()

    private$setLambdaPath(methylTask$lambda_max, eps, K)


    private$set_screeners()
    private$set_hal_lrnrs()
    private$set_baseline_learners_enet()
    private$buildSuperLearner()

    private$setTasks()
    private$set_delayed_learner(nworkers)




  },
  visualize = function(){
    return(plot(self$learners_spec$sl3_delayed_learner))
  },
  train_cv = function(parallel = T, also_full_model = T, verbose = T){

    if(parallel){
      nworkers = self$params$nworkers
      plan(multiprocess, workers = nworkers)
      delayed_lrnr= self$learners_spec$sl3_delayed_learner
      sched <- Scheduler$new(delayed_lrnr, FutureJob, nworkers = nworkers, verbose = verbose)
      cv_fit <- sched$compute()
      print("no")
      private$.learners_spec$sl3_superlearner = cv_fit
    }
    else{
      cv_fit = self$learners_spec$sl3_superlearner$train(self$learners_spec$trainingTask)

      private$.learners_spec$sl3_superlearner = cv_fit
    }
    if(also_full_model){

      self$train_full(parallel, verbose)

    }
    return(cv_fit)
  },
  train_full = function(parallel = T, verbose = F){
    if(parallel){
      lrnr = delayed_learner_train(self$learners_spec$sl3_library, self$learners_spec$trainingTask)
      sched <- Scheduler$new(lrnr, FutureJob, nworkers = self$params$nworkers, verbose = verbose)
      fit <- sched$compute()
      private$.learners_spec$sl3_library = fit
    }
    else{
      fit = self$learners_spec$sl3_library$train(self$learners_spec$trainingTask)
      private$.learners_spec$sl3_library = fit
    }

    return(fit)
  },
  predict = function(newTask = NULL, CV_pred = F){


    if(!CV_pred){
      lrnr =  self$learners_spec$sl3_library
    }
    else{
      lrnr =  self$learners_spec$sl3_superlearner
    }

    if(!CV_pred){

      coefs = which(self$learners_spec$sl3_superlearner$learner_fits$Lrnr_cv_selector$coefficients==1)
      if(is.null(newTask)){
        predmat = data.frame(lrnr$predict())
        return(list(predmat,predmat[,coefs, drop = F]))
      }
      predmat = data.frame(lrnr$predict(newTask))
      return(list(predmat,predmat[,coefs, drop=F]))
    }

    if(is.null(newTask)){
      return(lrnr$predict())
    }
    return(lrnr$predict(newTask))
  },


  #Replaces a subset of the fields in learners_spec with those contained in spec.
  #Can be used to specify custom base learners and meta learners.
  construct_super_learner_from_spec = function(spec){
    valid_names = intersect(names(spec), names(self$learners_spec))
    if(length(valid_names)==0){
      stop("Error: No valid names in spec.")
    }
    for(name in valid_names){
      item = spec[[name]]
      private$.learners_spec[[name]] = item
    }
    private$build_super_learner()
    private$set_delayed_learner(self$params$nworkers)



  },
  cv_risk = function(){
    self$learners_spec$sl3_superlearner$learner_fits$`CV_Pipeline(BigLrnr_screener_corP_0.1_5_3e+05->Pipeline(BigLrnr_biglasso_1->Stack))`$cv_risk(loss_squared_error)

  },
  validate = function(task = NULL, CV_pred = F){
    if(is.null(task)){
      task = private$methylTask$getTestTask()
    }

    truth = task$get_node("outcome")
    mat = self$predict(task, CV_pred = CV_pred, at_optimal = F )
    out = data.frame(t(apply(mat, MARGIN = 2, function(x){c(cor(x,truth),mean((truth-x)^2))})))
    colnames(out) = c("R", "MSE")
    out1=out[order(-out$R),]
    out2 = out[self$get_best_model_name(),]
    return(list(out1, out2))

  }


)
)
