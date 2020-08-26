


Lrnr_bigcorP_screener <- R6Class(
  classname = "Lrnr_bigcorP_screener",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(mat,minPvalue_fdr = 0.1, minscreen= 100, maxscreen = 10000, method = 'pearson', ...) {

      mat = bigmemory::describe(mat)
      params <- args_to_list()
      params$mat= mat
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("screener", "continuous", "binomial"),

    .train = function(task) {

      args <- self$params
      params = self$params
      outcome_type <- self$get_outcome_type(task)
      method <- args$method
      minPvalue <- args$minPvalue_fdr
      minscreen <- args$minscreen
      maxscreen <- args$maxscreen
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

      BM2FBM <- function(bm) {
        bigstatsr::FBM(nrow = bigmemory::nrow(bm), ncol = bigmemory::ncol(bm), type = "double",
                       backingfile = file.path(bigmemory::dir.name(bm), bigstatsr::sub_bk(bigmemory::file.name(bm))),
                       create_bk = FALSE)
      }
      FBM_X=BM2FBM(Z)
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


      # listPvalue <- as.vector(biganalytics::apply(X, 2, function(x, Y, method) {
      #   ifelse(var(x) <= 0, 1, cor.test(x[index], y = Y, method = method)$p.value)
      # }, Y = Y, method = method))

      # nextup = big_apply(tmpp,  function(X, ind, Y, method) {
      #   apply(X[,ind],MARGIN=2, FUN= function(x,Y, method){ ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method)$p.value)
      #   },Y = Y, method = method)}, Y = Y, method = method, a.combine = c)

      selected <- which(listPvalue <= minPvalue)

      if (length(selected) < minscreen) {
        selected[rank(listPvalue) <= minscreen] <- TRUE
      }
      if(length(selected) > maxscreen){
        selected =  order(listPvalue,decreasing=F)[1:maxscreen]
      }

      selected_names <- colnames(X)[selected]


      fit_object <- list(selected = selected_names)

      return(fit_object)
    },

    .predict = function(task = NULL) {

      verbose <- getOption("sl3.verbose")
      indices = task$row_index
      params <- self$params
      X = params$mat
      X = bigmemory::attach.big.matrix(X)
      if(is.null(indices)){
        indices = 1:nrow(X)
      }

      pred = as.data.table(X[indices,self$fit_object$selected])

      return(pred)

    },
    .chain = function(task) {

      keep = self$fit_object$selected
      mat = self$params$mat
      mat = bigmemory::attach.big.matrix(mat)
      indices = task$row_index
      if(is.null(indices)){
        indices = 1:nrow(mat)
      }
      meth_dat = as.data.table(mat[,keep])
      orig_data = task$get_data(row = 1:nrow(mat))

      dat = set(meth_dat,j=colnames(orig_data), value = orig_data)

      shared_data = Shared_Data$new(dat, force_copy = F)
      newtask = make_sl3_Task(data = shared_data, covariates = colnames(meth_dat), outcome = task$nodes$outcome, outcome_type= self$get_outcome_type(task)$type, row_index = indices)


      return(newtask)

    }

  )
)

Lrnr_biglasso <- R6Class(
  classname = "Lrnr_biglasso",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(alpha = 1,penalty= "enet", ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial"),

    .train = function(task) {

      verbose <- getOption("sl3.verbose")
      params <- self$params
      alpha = params$alpha
      penalty = params$penalty
      X = task$X
      tmp = bigmemory::as.big.matrix(as.matrix(X))

      if(alpha == 1){
        penalty = "lasso"
      }
      if(alpha == 0){
        penalty = "ridge"
      }
      y = as.vector(task$Y)

      outcome_type = self$get_outcome_type(task)$type
      fit_cv = biglasso::cv.biglasso(tmp, y,family = ifelse(outcome_type == "continuous", "gaussian", "binomial" ), output.time = T, alpha = alpha, penalty = penalty,nfolds=10)
      fit_object = list()
      fit_object$fit = fit_cv

      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      fit = self$fit_object$fit

      pred = as.data.table(as.vector(predict(fit, bigmemory::as.big.matrix(as.matrix(task$X)), type ="response")))

      return(pred)
    }





  )
)



Lrnr_biglasso_screener <- R6Class(
  classname = "Lrnr_biglasso_screener",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(mat = NULL, alpha = 1, penalty = "enet",...) {
      # mat = bigmemory::deepcopy(mat)
      #bigmemory::write.big.matrix(mat, "tempfile",col.names=T)
      mat = bigmemory::describe(mat)
      #mat = "tempfile"
      # if(!is.null(mat) & !bigmemory::is.big.matrix(mat)){
      #   mat = bigmemory::as.big.matrix(as.matrix(mat))
      # }
      params <- args_to_list()
      params$mat= mat
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("screener", "continuous", "binomial"),

    .train = function(task) {
      verbose <- getOption("sl3.verbose")
      params <- self$params
      mat = params$mat
      alpha = params$alpha
      penalty = params$penalty

      #mat = bigmemory::read.big.matrix(filename = mat, type = "double", header=T)
      mat = bigmemory::attach.big.matrix(mat)
      if(alpha == 1){
        penalty = "lasso"
      }
      if(alpha == 0){
        penalty = "ridge"
      }
      if(!is.null(mat)){

        tmp = mat
        if(!bigmemory::is.big.matrix(tmp)){

          tmp = bigmemory::as.big.matrix(as.matrix(mat))
        }
        index = as.vector(task$row_index)
        if(is.null(index)){
          index = 1:nrow(tmp)
        }
        tmp = bigmemory::deepcopy(tmp, rows = index)

      }
      else{
        X = task$X
        tmp = bigmemory::as.big.matrix(as.matrix(X))
        index = 1:nrow(tmp)
      }

      y = task$Y


      outcome_type = self$get_outcome_type(task)$type


      fit_cv = biglasso::cv.biglasso(tmp, as.vector(y),family = ifelse(outcome_type == "continuous", "gaussian", "binomial" ), output.time = T, alpha = alpha, penalty = penalty,ncores =1)
      index = fit_cv$min
      fit = fit_cv$fit
      coefs = fit$beta[,index]
      coefs = coefs[-1]
      keep = names(which(abs(coefs)>1e-15))
      if(length(keep) == 0){
        keep = colnames(tmp)[1:10]
      }
      fit_object = list()
      fit_object$keep = keep

      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      keep = self$fit_object$keep
      mat = self$params$mat
      mat = bigmemory::attach.big.matrix(mat)
      index = task$row_index
      if(is.null(index)){
        index = 1:nrow(mat)
      }
      if(!is.null(mat)){

        res = as.data.table(mat[index,keep])
        colnames(res) = keep
        return(res)
      }
      pred = task$X[, keep, with = FALSE, drop = FALSE]
      return(pred)
    },


    .chain = function(task) {
      keep = self$fit_object$keep
      mat = self$params$mat
      #mat=bigmemory::read.big.matrix(filename = mat,type="double",header=T)
      mat = bigmemory::attach.big.matrix(mat)
      index = task$row_index
      if(is.null(index)){
        index = 1:nrow(mat)
      }
      if(!is.null(mat)){

        keep = self$fit_object$keep


        indices = task$row_index
        if(is.null(indices)){
          indices = 1:nrow(mat)
        }
        meth_dat = as.data.table(mat[,keep])
        orig_data = task$get_data(row = 1:nrow(mat))

        dat = set(meth_dat,j=colnames(orig_data), value = orig_data)

        shared_data = Shared_Data$new(dat, force_copy = F)
        newtask = make_sl3_Task(data = shared_data, covariates = colnames(meth_dat), outcome = task$nodes$outcome, outcome_type= self$get_outcome_type(task)$type, row_index = indices)

        return(newtask)
      }
      return(task$next_in_chain(covariates = keep))
    }
  )
)





Lrnr_constant <- R6Class(
  classname = "Lrnr_constant",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(constant_mat,...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }

  ),
  private = list(
    .properties = c(
      "continuous", "binomial", "categorical", "weights"
    ),

    .train = function(task) {
      verbose <- getOption("sl3.verbose")
      params <- self$params

      constant_mat <- params$constant_mat

      fit_object = list()
      fit_object$constant_mat = constant_mat



      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      indices = task$row_index
      constant_mat= self$fit_object$constant_mat
      if(is.null(indices)){
        indices = 1:nrow(constant_mat)
      }
      predictions <- as.data.table(constant_mat[indices,])
      return(predictions)
    }
  )
)




Lrnr_screener_limma <- R6Class(
  classname = "Lrnr_screener_limma",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(cutoff = 0.1, max_num = 1000,...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c(
      "continuous", "binomial", "categorical", "weights"
    ),

    .train = function(task) {
      verbose <- getOption("sl3.verbose")
      params <- self$params
      max_num <- params$max_num
      cutoff <- params$cutoff
      y = task$get_node("outcome")
      X=task$X
      cpgNames = task$nodes$covariates
      X=t(X)

      design <- as.matrix(cbind(rep(1, times = length(y)), y))


      mod_fit <- limma::lmFit(object = as.matrix(X), design = design)
      rm(design)
      rm(X)
      mod_fit <- limma::eBayes(mod_fit)

      # extract indices of relevant CpG sites
      tt_fit <- limma::topTable(mod_fit, coef = 2, num = max_num, sort.by = "p", p.value=cutoff)

      indices_pass <- as.numeric(which(cpgNames %in% rownames(tt_fit)))

      if(length(indices_pass) <=5 ){
        print("none selected")
        tt_fit <- limma::topTable(mod_fit, coef = 2, num = 10, sort.by = "p", p.value=1)
        indices_pass <- which(cpgNames %in% rownames(tt_fit))
      }

      screen_ind <- as.numeric(indices_pass)
      cpgNames = cpgNames[screen_ind]

      fit_object = list()
      fit_object$cpgNames = cpgNames

      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      X <- as.matrix(task$X)

      predictions <- as.matrix(X[,self$fit_object$cpgNames])
      return(predictions)
    }
  )
)



Lrnr_discretizer <- R6Class(
  classname = "Lrnr_discretizer",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(bins = 50,...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c(
      "continuous", "binomial", "categorical", "weights"
    ),

    .train = function(task) {
      verbose <- getOption("sl3.verbose")
      params <- self$params
      column_index <- params$column_index
      bins <- params$bins

      # specify data
      X <- as.matrix(task$X)

      convertColumn = function(x){
        quants = seq(0,1,1/bins)
        q=quantile(x,quants)
        (q)
        nearest <- findInterval(x, q)
        x <- q[nearest]
        return(x)
      }
      quantizer = function(X){as.matrix(apply(X, MARGIN = 2, FUN =convertColumn))}
      fit_object <- list()
      fit_object$quantizer = quantizer


      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      X <- as.matrix(task$X)
      quantizer = self$fit_object$quantizer
      predictions <- as.matrix(quantizer(X))
      return(predictions)
    }
  )
)



##############
Lrnr_column_selector <- R6Class(
  classname = "Lrnr_column_selector",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(column_index,...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c(
      "continuous", "binomial", "categorical", "weights"
    ),

    .train = function(task) {
      verbose <- getOption("sl3.verbose")
      params <- self$params
      column_index <- params$column_index


      # specify data
      X <- as.matrix(task$X)


      fit_object <- list()

      fit_object$pred <- as.vector(X[,column_index])

      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      X <- as.matrix(task$X)

      predictions <- as.vector(X[, self$params$column_index])
      return(predictions)
    }
  )
)






#' Coefficient Magnitude Screener
#'
#' This learner provides screening of covariates based on the magnitude of
#' their estimated coefficients in a (possibly regularized) GLM.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#'
#' @export
#'
#' @keywords data
#'
#' @return Learner object with methods for training and prediction. See
#'  \code{\link{Lrnr_base}} for documentation on learners.
#'
#' @format \code{\link{R6Class}} object.
#'
#' @family Learners
#'
#' @section Parameters:
#' \describe{
#'   \item{\code{learner}}{An instantiated learner to use for estimating
#'     coefficients used in screening.}
#'   \item{\code{threshold = 1e-3}}{Minimum size of coefficients to be kept.
#'   If NULL, will select max_retain variables.}
#'   \item{\code{min_retain = 2}}{Minimum number of variables to be kept.}
#'   \item{\code{max_retain = NULL}}{Maximum number of variables to be kept.
#'   Default selects all with absolute coefficient above threshold.}
#'   \item{\code{verbose = FALSE}}{Print selected variables.}
#'   \item{\code{...}}{Other parameters passed to \code{learner}.}
#' }
#'
#'
#'
Lrnr_screener_selecter <- R6Class(
  classname = "Lrnr_screener_coefs",
  inherit = Lrnr_base, portable = TRUE, class = TRUE,
  public = list(
    initialize = function(nvar = NULL, ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("screener"),

    .train = function(task) {
      nvar = self$params$nvar
      covariates = task$nodes$covariates
      if(is.null(nvar) | nvar >  length(covariates)){
        nvar = length(covariates)
      }
      fit_object <- list(selected = covariates[1:nvar])
      return(fit_object)

    },

    .predict = function(task) {
      task$X[, private$.fit_object$selected, with = FALSE, drop = FALSE]
    },

    .chain = function(task) {
      return(task$next_in_chain(covariates = private$.fit_object$selected))
    },
    .required_packages = c()
  )
)




Lrnr_screener_randomForest_ordered <- R6Class(
  classname = "Lrnr_screener_randomForest",
  inherit = Lrnr_base, portable = TRUE, class = TRUE,
  public = list(
    initialize = function(nVar = 10,
                          ntree = 1000,
                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),

  private = list(
    .properties = c("binomial", "continuous", "categorical"),

    .train = function(task) {
      call_with_args <- function(fun, args, other_valid = list(), keep_all = FALSE,
                                 silent = FALSE, ignore = c()) {

        # drop ignore args
        args <- args[!(names(args)%in%ignore)]
        if (!keep_all) {
          # catch arguments to be kept
          formal_args <- names(formals(fun))
          all_valid <- c(formal_args, other_valid)

          # find invalid arguments based on combination of formals and other_valid
          invalid <- names(args)[which(!(names(args) %in% all_valid))]

          # subset arguments to pass
          args <- args[which(names(args) %in% all_valid)]

          # return warnings when dropping arguments
          if (!silent & length(invalid) > 0) {
            message(sprintf(
              "Learner called function %s with unknown args: %s. These will be dropped.\nCheck the params supported by this learner.",
              as.character(substitute(fun)), paste(invalid, collapse = ", ")
            ))
          }
        }
        do.call(fun, args)
      }
      args <- self$params
      outcome_type <- self$get_outcome_type(task)
      args$x <- task$X
      args$y <- outcome_type$format(task$Y)

      if (is.null(args$mtry)) {
        args$mtry <- floor(ncol(args$x))
      }
      if (outcome_type$type == "binomial") {
        args$y <- factor(args$y, levels = c(0, 1))
      }
      rf_fun <- getS3method(
        "randomForest", "default",
        envir = getNamespace("randomForest")
      )
      rf_object <- call_with_args(rf_fun, args)
      covariates <- task$nodes$covariates
      orders = order(rank(-rf_object$importance))
      covariates = covariates[orders]

      fit_object <- list(selected = covariates)
      return(fit_object)
    },

    .predict = function(task) {
      task$X[, private$.fit_object$selected, with = FALSE, drop = FALSE]
    },

    .chain = function(task) {
      return(task$next_in_chain(covariates = private$.fit_object$selected))
    },

    .required_packages = c("randomForest")
  )
)



Lrnr_screener_coefs <- R6Class(
  classname = "Lrnr_screener_coefs",
  inherit = Lrnr_base, portable = TRUE, class = TRUE,
  public = list(
    initialize = function(learner, threshold = 1e-3, min_retain = 2,
                          max_retain = NULL, verbose = FALSE, ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("screener"),

    .train = function(task) {
      print("training")
      # check arguments
      threshold <- self$params$threshold
      min_retain <- self$params$min_retain
      max_retain <- self$params$max_retain
      covs <- task$nodes$covariates
      print("step 1")
      if(is.null(threshold) & is.null(max_retain)){
        stop("threshold, max_retain or both must be provided")
      }

      if(is.null(min_retain)){
        stop("min_retain must be provided")
      }

      if(length(min_retain) >= length(covs)){
        stop("modify min_retain to be less than the number of covariates")
      }

      if(length(max_retain) > length(covs)){
        max_retain <- NULL
        warning("max_retain is greater than the number of covariates, modifying max_retain to NULL")
      }

      learner <- self$params$learner
      fit <- learner$train(task)
      print("Step 2")
      rf_object = fit$fit_object
      selected <- (rank(-rf_object$importance) <= max_retain)
      selected_names <- names(task$X)[selected]
      covariates <- task$nodes$covariates
      covariate_selected <- sapply(covariates, function(covariate) {
        any(grep(covariate, selected_names))
      })
      fit_object <- list(selected = covariates[covariate_selected])
      return(fit_object)

    },

    .predict = function(task) {
      task$X[, private$.fit_object$selected, with = FALSE, drop = FALSE]
    },

    .chain = function(task) {
      return(task$next_in_chain(covariates = private$.fit_object$selected))
    },
    .required_packages = c()
  )
)



##########

#' Truncates predictions to ensure loss function is bounded.
#'
#' @param preds A \code{numeric} vector of predictions to to be bounded.
#' @param bounds Either a \code{numeric} vector of length two, giving the
#'  closed interval (lower, upper), or just a lower bound. In the latter case,
#'  the upper bound is computed as (1 - lower). The default is 0.001.
#'
#' @return Truncated predictions.
#'
#' @keywords internal
bound <- function(preds, bounds = 0.001) {
  lower <- bounds[[1]]
  if (length(bounds) > 1) {
    upper <- bounds[[2]]
  } else {
    upper <- 1 - lower
  }
  preds_bounded <- pmin(pmax(preds, lower), upper)
  return(preds_bounded)
}

################################################################################

# if warning is in ignoreWarningList, ignore it; otherwise post it as usual
SuppressGivenWarnings <- function(expr, warningsToIgnore) {
  h <- function(w) {
    if (w$message %in% warningsToIgnore) invokeRestart("muffleWarning")
  }
  withCallingHandlers(expr, warning = h)
}

################################################################################

GetWarningsToSuppress <- function(update.step = FALSE) {
  warnings.to.suppress <- c(
    paste(
      "glm.fit: fitted probabilities numerically 0",
      "or 1 occurred"
    ),
    paste(
      "prediction from a rank-deficient fit may be",
      "misleading"
    ),
    "non-integer #successes in a binomial glm!",
    "the matrix is either rank-deficient or indefinite",
    "glm.fit: algorithm did not converge"
  )
  return(warnings.to.suppress)
}

################################################################################

#' Streamline Function Arguments
#'
#' Reduce a list of function argsuments by taking a function body and returning
#' only the arguments that belong to the function signature.
#'
#' @param Args A \code{list} of function arguments to be streamlined.
#' @param fun A \code{function} whose signature will be used to reduce the
#'  arguments passed in.
#'
#' @keywords internal
keep_only_fun_args <- function(Args, fun) {
  keepArgs <- intersect(names(Args), names(formals(fun)))
  # captures optional arguments given by user
  if (length(keepArgs) > 0) {
    Args <- Args[keepArgs]
  } else {
    Args <- NULL
  }
  return(Args)
}

################################################################################

#' Call with filtered argument list
#'
#' Call a function with a list of arguments, eliminating any that aren't
#' matched in the function prototype
#'
#' @param fun A \code{function} whose signature will be used to reduce the
#' @param args A \code{list} of function arguments to use.
#' @param other_valid A \code{list} of function arguments names that are valid,
#'   but not formals of \code{fun}.
#' @param keep_all A \code{logical} don't drop arguments, even if they aren't
#'   matched in either the function prototype or other_valid.
#' @param silent A \code{logical} indicating whether to pass \code{message}s
#'  when arguments not found in \code{formals} are passed to \code{fun}.
#' @param ignore A \code{character} vector indicating which arguments should be dropped
#'
#' @keywords internal
call_with_args <- function(fun, args, other_valid = list(), keep_all = FALSE,
                           silent = FALSE, ignore = c()) {

  # drop ignore args
  args <- args[!(names(args)%in%ignore)]
  if (!keep_all) {
    # catch arguments to be kept
    formal_args <- names(formals(fun))
    all_valid <- c(formal_args, other_valid)

    # find invalid arguments based on combination of formals and other_valid
    invalid <- names(args)[which(!(names(args) %in% all_valid))]

    # subset arguments to pass
    args <- args[which(names(args) %in% all_valid)]

    # return warnings when dropping arguments
    if (!silent & length(invalid) > 0) {
      message(sprintf(
        "Learner called function %s with unknown args: %s. These will be dropped.\nCheck the params supported by this learner.",
        as.character(substitute(fun)), paste(invalid, collapse = ", ")
      ))
    }
  }
  do.call(fun, args)
}

