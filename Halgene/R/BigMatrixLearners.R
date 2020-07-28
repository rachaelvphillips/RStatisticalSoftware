




#Pairs correlation screening with biglasso
Lrnr_big_baseline <- R6Class(
  classname = "Lrnr_big_baseline",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(mat, group = NULL, dobiglasso= F,alpha = c(1, 0.5,0.05), screen_cor = F, minPvalue_fdr = 0.1, minscreen= 100, maxscreen = 600000, method = 'pearson', ...) {

      alpha[alpha < 1e-3] = 1e-3
      mat = bigmemory::describe(mat)
      params <- args_to_list()
      params$mat= mat
      params$group= group
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("screener", "continuous", "binomial"),

    .train = function(task) {
      print("hi")
      args <- self$params
      params = self$params
      screen_cor = params$screen_cor
      outcome_type <- self$get_outcome_type(task)
      method <- args$method
      alpha = params$alpha
      group = params$group
      group = as.vector(unlist(group))
      minPvalue <- args$minPvalue_fdr
      minscreen <- args$minscreen
      maxscreen <- args$maxscreen
      dobiglasso = args$dobiglasso
      X <- params$mat
      X = bigmemory::attach.big.matrix(X)
      Y <- outcome_type$format(task$Y)

      index = task$row_index
      if(is.null(index)){
        index = 1:nrow(X)
      }
      if(!bigmemory::is.filebacked(X)){
        #Z=bigmemory::deepcopy(X,  backingfile = "tmpback10101.bk")
      }
      else{
        Z=X
      }
      Z=X


      ######
      ####### Begin screening
      BM2FBM <- function(bm) {
        bigstatsr::FBM(nrow = bigmemory::nrow(bm), ncol = bigmemory::ncol(bm), type = "double",
                       backingfile = file.path(bigmemory::dir.name(bm), bigstatsr::sub_bk(bigmemory::file.name(bm))),
                       create_bk = FALSE)
      }
      FBM_X=BM2FBM(Z)
      if(screen_cor){


        if(outcome_type$type=="continuous"){

          fit_biglm = bigstatsr::big_univLinReg(FBM_X,Y,ind.train = index)
        }
        else{
          fit_biglm = bigstatsr::big_univLogReg(FBM_X,Y,ind.train = index)
        }

        listPvalue = as.vector(predict(fit_biglm, log10 = F))

        if(length(listPvalue)>=2000){
          qs = (qvalue::qvalue(listPvalue))

          listPvalue = as.vector(qs$qvalues)
        }

        selected <- which(listPvalue <= minPvalue)

        if (length(selected) < minscreen) {
          selected =  order(listPvalue,decreasing=F)[1:100]
        }
        if(length(selected) > maxscreen){
          selected =  order(listPvalue,decreasing=F)[1:maxscreen]
        }

      }
      else{
        selected = 1:ncol(X)

      }

      #######
      row_indices = index
      col_indices = selected


      if(!is.null(group)){
        group_indices = which(colnames(X) %in% group)
      }
      else{
        group_indices = col_indices
      }

      print(length(row_indices))
      print(length(intersect(col_indices, group_indices)))

      if(dobiglasso){
        print("biglasso")
        bm = FBM_X$bm()

        bm2 = bigmemory::deepcopy(bm, rows = row_indices, cols = intersect(col_indices, group_indices),backingfile = "tmpfile101.bk")
        if(outcome_type$type=="continuous"){
          print("biglasso2")
          fit = biglasso::cv.biglasso(bm2,Y, alpha = alpha,family = "gaussian", penalty = "enet")
        }
        else{
          fit = biglasso::cv.biglasso(bm2,Y, alpha = alpha,family = "binomial", penalty = "enet")
        }
      }
      else{


      if(outcome_type$type=="continuous"){

        fit = bigstatsr::big_spLinReg(FBM_X,Y, ind.train = row_indices, ind.col = intersect(col_indices, group_indices), alphas = alpha)
      }
      else{
        fit = bigstatsr::big_spLogReg(FBM_X,Y, ind.train = row_indices, ind.col = intersect(col_indices, group_indices), alphas = alpha)
      }
      }


      fit_object <- list(fit = fit, FBM_X = FBM_X, col_indices = intersect(col_indices, group_indices) )

      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      fit_obj = self$fit_object
      fit = fit_obj$fit
      FBM_X = fit_obj$FBM_X

      col_indices = fit_obj$col_indices
      row_indices = task$row_index
      if(is.null(row_indices)){
        row_indices = 1:FBM_X$nrow
      }
      if(self$params$dobiglasso){
        print("biglasso prdds")
        bm = fit_obj$bm
      Xs = bigmemory::deepcopy(bm, cols = col_indices, rows = row_indices,backingfile = "tmpfile102.bk")

      pred = predict(fit, Xs, type = "response")
      }
      else{
        pred = predict(fit, FBM_X, ind.row = row_indices, ind.col = col_indices)
      }
      pred = as.vector(pred)

      return(pred)

    }

  )
)
