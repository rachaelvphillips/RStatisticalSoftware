#' @import  bigmemory
#' @import  data.table
#' @export

BigLrnr_biglasso <- R6Class(
  classname = "BigLrnr_biglasso",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(alpha = 1, ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "bigmatrix"),

    .train = function(task) {
      if(!inherits(task,"sl3_BigTask")){
        stop("Task must be a BigTask.")
      }

      verbose <- getOption("sl3.verbose")
      params <- self$params
      alpha = params$alpha

      outcome_type = self$get_outcome_type(task)$type
      fbm =task$fbm
      row_indices = task$row_index_safe

      col_index = task$bm_indices

      Y=task$Y

      if(outcome_type=="continuous"){

        fit = bigstatsr::big_spLinReg(fbm,Y, ind.train = row_indices, ind.col = col_index, alphas = alpha)
      }
      else{
        fit = bigstatsr::big_spLogReg(fbm,Y, ind.train = row_indices, ind.col = col_index, alphas = alpha)
      }
      betas = as.matrix(do.call(cbind, lapply(fit[[1]], function(elem){return(as.vector(elem$beta))})))

      fit_object = list(fit= fit, bm_columns=task$bm_columns,col_index = task$bm_indices, beta = get_beta(betas))

      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      fit = self$fit_object$fit
      bm_columns = self$fit_object$bm_columns
      col_index=self$fit_object$col_index
      row_indices = task$row_index_safe
      fbm =task$fbm
      pred = predict(fit, fbm, ind.row = row_indices, ind.col =col_index[match(bm_columns,task$bm_columns)])

      return(pred)
    },
    .chain = function(task) {

      fit = self$fit_object$fit
      bm_columns = self$fit_object$bm_columns
      col_index=self$fit_object$col_index
      row_indices = task$row_index_safe
      fbm =task$fbm
      pred = as.data.table(as.vector(predict(fit, fbm, ind.row = row_indices, ind.col =col_index[match(bm_columns,task$bm_columns)])))


      col_map = task$add_columns(pred)

      nexttask = task$next_in_chain(covariates = c(task$nodes$covariates, colnames(pred)), column_names= col_map)
      return(nexttask)
    }






  )
)







#' @import  bigmemory
#' @import  data.table
#' @export
BigLrnr_screener_corP <- R6Class(
  classname = "BigLrnr_screener_corP",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(minPvalue_fdr=0.1, minscreen=5, maxscreen=Inf, ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "bigmatrix", "screener"),

    .train = function(task) {
      if(!inherits(task,"sl3_BigTask")){
        stop("Task must be a BigTask.")
      }
      verbose <- getOption("sl3.verbose")
      params <- self$params
      minPvalue <- params$minPvalue_fdr
      minscreen <- params$minscreen
      maxscreen <- params$maxscreen

      outcome_type = self$get_outcome_type(task)$type
      fbm =task$fbm
      row_indices = task$row_index_safe
      col_index = task$bm_indices
      col_names = task$bm_columns
      Y=task$Y

      if(outcome_type=="continuous"){

        fit_biglm = bigstatsr::big_univLinReg(fbm,Y,ind.train = row_indices,ind.col=col_index)
      }
      else{
        fit_biglm = bigstatsr::big_univLogReg(fbm,Y,ind.train = row_indices,ind.col=col_index)
      }

      listPvalue = as.vector(predict(fit_biglm, log10 = F))

      if(length(listPvalue)>=2000){
        qs = qvalue::qvalue(listPvalue)
        listPvalue = as.vector(qs$qvalues)
      }

      selected <- which(listPvalue <= minPvalue)

      if (length(selected) < minscreen) {
        selected =  order(listPvalue,decreasing=F)[1:minscreen]
      }
      else if(length(selected) > maxscreen){
        selected =  order(listPvalue,decreasing=F)[1:maxscreen]
      }

      selected_index = col_index[selected]
      selected_names = col_names[selected]
      fit_object = list(selected_index=selected_index, selected_names=selected_names)

      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      selected_names = self$fit_object$selected_names
      bm=task$bm
      return(as.data.table(bm[,selected_names]))
    },
    .chain = function(task) {
      selected_names = self$fit_object$selected_names
      newtask=task$next_in_chain(bm_columns = selected_names)
      return(newtask)
    }






  )
)








#' @import  bigmemory
#' @import  data.table
#' @export
BigLrnr_randomPCA <- R6Class(
  classname = "BigLrnr_randomPCA",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(center = T, scale = T, K = 10, ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "bigmatrix", "preprocessing"),

    .train = function(task) {
      if(!inherits(task,"sl3_BigTask")){
        stop("Task must be a BigTask.")
      }
      verbose <- getOption("sl3.verbose")
      params <- self$params
      center = params$center
      scale = params$scale
      K = params$K

      outcome_type = self$get_outcome_type(task)$type
      fbm =task$fbm
      row_indices = task$row_index_safe
      col_index = task$bm_indices
      col_names = task$bm_columns
      bigPCA = bigstatsr::big_randomSVD(fbm,
                     fun.scaling = big_scale(center = center, scale = scale), row_indices,
                     ind.col = col_index,
                     k = K )
      fit_object = list(bigPCA = bigPCA,col_index=col_index,col_names=col_names )
      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      row_indices = task$row_index_safe
      tsk_bm_columns = task$bm_columns

      col_index=self$fit_object$col_index
      col_names=self$fit_object$col_names
      task_col_indices = col_index[match(tsk_bm_columns,col_names)]
      bigPCA = self$fit_object$bigPCA
      fbm=task$fbm
      pred = as.data.table(predict(bigPCA, fbm, ind.row=row_indices,ind.col=task_col_indices))
      return(pred)
    },
    .chain = function(task) {
     pred = private$.predict(task)
      col_map=task$add_columns(pred)
      newtask = task$next_in_chain(covariates = c(task$nodes$covariates, colnames(pred)), column_names= col_map)
      return(newtask)
    }






  )
)





#' @import  bigmemory
#' @import  data.table
#' @export
BigLrnr_subset_covariates <- R6Class(
  classname = "BigLrnr_subset_covariates",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(data_columns = NULL, bm_columns = NULL,covariates = NULL, ...) {

      params <- args_to_list()
      params$bm_columns = unlist(bm_columns)
      params$covariates = unlist(covariates)
      params$data_columns = unlist(data_columns)
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "bigmatrix"),

    .train = function(task) {
      if(!inherits(task,"sl3_BigTask")){
        stop("Task must be a BigTask.")
      }

      verbose <- getOption("sl3.verbose")
      params <- self$params
      data_columns = params$data_columns
      bm_columns = params$bm_columns
      covariates = params$covariates
      task_bm_columns = task$bm_columns
      task_data_columns = task$nodes$covariates


      data_columns = intersect(task_data_columns,c(covariates,data_columns))
      bm_columns = intersect(task_bm_columns,c(covariates,bm_columns))
      if(length(data_columns) ==0){
        data_columns = NULL
      }
      if(length(bm_columns) ==0){
        bm_columns = NULL
      }



      fit_object = list(data_columns = data_columns, bm_columns = bm_columns)
      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      data_columns = self$fit_object$data_columns
      bm_columns = self$fit_object$bm_columns
      bm = task$bm
      row_index = task$row_index_safe
      pred = cbind(task$X[,data_columns, with = F], as.data.table(bm[row_index,bm_columns]))
      return(pred)
    },
    .chain = function(task) {
      data_columns = self$fit_object$data_columns
      bm_columns = self$fit_object$bm_columns
      newtask = task$next_in_chain(covariates = data_columns, bm_columns = bm_columns)
      return(newtask)
    }






  )
)




#' @import  bigmemory
#' @import  data.table
#' @export
BigLrnr_screener_limma <- R6Class(
  classname = "BigLrnr_screener_limma",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(cutoff = 0.1, max_num = Inf,...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c(
      "continuous", "binomial", "screener"
      , "bigmatrix"),

    .train = function(task) {
      if(!inherits(task,"sl3_BigTask")){
        stop("Task must be a BigTask.")
      }
      verbose <- getOption("sl3.verbose")
      params <- self$params
      max_num <- params$max_num
      cutoff <- params$cutoff
      y = task$get_node("outcome")
      bm = task$bm
      X = bm[task$row_index, task$bm_columns]
      col_names = colnames(bm)
      colnames(X) = col_names
      X=t(X)

      design <- as.matrix(cbind(rep(1, times = length(y)), y))

      mod_fit <- limma::lmFit(object = as.matrix(X), design = design)
      rm(design)
      rm(X)
      mod_fit <- limma::eBayes(mod_fit)

      # extract indices of relevant CpG sites
      tt_fit <- limma::topTable(mod_fit, coef = 2, num = max_num, sort.by = "p", p.value=cutoff)

      indices_pass <- as.numeric(which(col_names %in% rownames(tt_fit)))

      if(length(indices_pass) <=5 ){
        tt_fit <- limma::topTable(mod_fit, coef = 2, num = 10, sort.by = "p", p.value=1)
        indices_pass <- which(col_names %in% rownames(tt_fit))
      }

      screen_ind <- as.numeric(indices_pass)
      col_names = col_names[screen_ind]

      fit_object = list()
      fit_object$col_names = col_names
      fit_object$tt_fit = tt_fit
      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      bm = task$bm
      col_names= self$fit_object$col_names
      predictions <- as.matrix(bm[task$row_index,col_names])
      return(predictions)
    },

    .chain = function(task = NULL){
      nexttask = task$next_in_chain(bm_columns = col_names)
      return(nexttask)
    }
  )
)


#' @import  bigmemory
#' @import  data.table
#' @export
BigLrnr_screener_biglasso <- R6Class(
  classname = "BigLrnr_screener_biglasso",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(alpha = 1, df_max= 50000, add_to_data = T, ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "bigmatrix", "screener"),

    .train = function(task) {
      if(!inherits(task,"sl3_BigTask")){
        stop("Task must be a BigTask.")
      }

      verbose <- getOption("sl3.verbose")
      params <- self$params
      alpha = params$alpha

      outcome_type = self$get_outcome_type(task)$type
      fbm =task$fbm
      row_indices = task$row_index_safe
      df_max = params$df_max
      col_index = task$bm_indices
      Y=task$Y

      if(outcome_type=="continuous"){

        fit = bigstatsr::big_spLinReg(fbm,Y, ind.train = row_indices, ind.col = col_index, alphas = alpha, dfmax=df_max)
      }
      else{
        fit = bigstatsr::big_spLogReg(fbm,Y, ind.train = row_indices, ind.col = col_index, alphas = alpha,dfmax = df_max)
      }

      betas = bigstatsr::get_beta(as.matrix(do.call(cbind, lapply(fit[[1]], function(elem){return(as.vector(elem$beta))}))))
      keep_index = which(abs(betas) > 1e-9)
      keep_columns = task$bm_columns[keep_index]
      fit_object = list(fit=fit, selected = keep_columns)

      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")

      selected = self$fit_object$selected

      row_indices = task$row_index_safe
      bm =task$bm
      pred = as.data.table(bm[row_indices,selected])

      return(pred)
    },
    .chain = function(task) {
      selected = self$fit_object$selected

      if(self$params$add_to_data){


        col_map=task$add_columns(task$bm[task$row_index_safe,selected])
        newtask = task$next_in_chain(covariates = c(task$nodes$covariates, selected), column_names= col_map)
        return(newtask)

      }
      else{
        nexttask = task$next_in_chain(bm_columns = selected)

      }
      return(nexttask)
    }






  )
)




#' @import  bigmemory
#' @import  data.table
#' @export
Lrnr_screener_dist_cor <- R6Class(
  classname = "Lrnr_screener_dist_cor",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(method = NULL, algorithm = 50, ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "bigmatrix"),

    .train = function(task) {
      if(!inherits(task,"sl3_BigTask")){
        stop("Task must be a BigTask.")
      }

      verbose <- getOption("sl3.verbose")
      params <- self$params
      method = params$method
      max_num = params$max_num
      outcome_type = self$get_outcome_type(task)$type
      if(is.null(method)){
        method = ifelse(outcome_type == "continuous", "DC-SIS", "MV-SIS")
      }
      fbm =task$fbm
      row_indices = task$row_index_safe

      col_index = task$bm_indices


      Y=task$Y
      screenIID = function(mat, ind, Y, method){
        obj = VariableScreening::screenIID(mat[row_indices,ind], Y, method)
        return(obj$measurement)
      }

      importance_scores = bigstatsr::big_apply(fbm, a.FUN = screenIID, a.combine = c, ind = col_index,  Y = Y, method = method)
      ranks = rank(importance_scores)
      top_ind =which(ranks <= max_num)
      keep_ind = col_index[top_ind]

      fit_object = list(selected = task$bm_columns[keep_ind], ranks = importance_scores[top_ind])

      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      selected = self$fit_object$selected
      row_indices = task$row_index_safe


      return(task$bm[row_indices,selected])
    },
    .chain = function(task) {

      selected = self$fit_object$selected

      nexttask = task$next_in_chain(bm_columns = selected)
      return(nexttask)
    }
  )
)




