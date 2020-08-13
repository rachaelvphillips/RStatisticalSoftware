
#' @import  bigmemory
#' @import  data.table
#' @import  R6
#' @export
HALnet_Task <- R6Class("HALnet_Task", private = list(

  ##########
  #Variables
  ##########
  .params = NULL,
  .split_indices = NULL,
  .trainTask = NULL,
  .testTask = NULL,
  .lambda_max = NULL,
  .max_cor_dim = 25000,


  ##########
  #FUNCTIONS
  ##########

  splitData = function(){

    num_rows = self$nrow
    split = self$params$split
    if(is.null(split) | split == 1){
      private$.split_indices = NULL
    }
    print("Splitting data into training and validation sets...")
    dt = sort(sample(num_rows, num_rows*split))
    private$.split_indices = dt
  },
  #Compute groups using heirarchical clustering.
  #If too big then clustering occurs at meta-level, after clustering via kmeans to smaller size.
  computeGroups = function(num_groups, select_optimal_number, selection_method,kmeans_pre_cluster){

    print("Computing groups...")
    fbm = self$fbm
    bm = self$bm

    bm_names = self$base_columns
    base_column_index = match(bm_names, colnames(bm))
    if(kmeans_pre_cluster | ncol(bm) > private$.max_cor_dim){
      if(is.null(self$params$base_predictors)){
        fbm_new = fbm
      }
      else{
        fbm_new = bigstatsr::big_copy(fbm, ind.col = base_column_index)
      }

      fbm_transp = bigstatsr::big_transpose(fbm_new)

      bm_transp = fbm_transp$bm()

      bigkmeans = biganalytics::bigkmeans(bm_transp, centers = min(private$.max_cor_dim, ncol(bm)/100), dist = "cosine")
      members = bigkmeans$cluster
      kmeans_centers = bigkmeans$centers
      corMat = cor(kmeans_centers)
      df = kmeans_centers


    }
    else{
      df = fbm[,base_column_index]
      corMat  = bigstatsr::big_cor(fbm, ind.col = base_column_index)
      corMat = corMat[,]

      members = NULL

    }
    dist= as.dist(abs(1-corMat))
    #Generate heir-cluster groups (possibly at meta level)
    if(select_optimal_number){
      hclust_choice = factoextra::fviz_nbclust(
        t(df),
        FUNcluster = factoextra::hcut,
        method = selection_method[1],
        diss = dist,
        k.max = num_groups,
        print.summary = F)
      num_groups = which.max(hclust_choice$data[,2])


    }


   print("clustering...")
    hclust_fit = hclust(dist, members = members, method = "complete")
    group_index = cutree(hclust_fit, num_groups)
    private$.params$hclust_fit  = hclust_fit

    print(group_index)
    groups = split(bm_names, group_index)
    private$.params$groups = groups
    print(paste0("We have generated ", length(groups), " groups."))
  },


  setLambdaMax = function(){
    lrnr = make_learner(BigLrnr_biglasso)

    lrnr_trained = lrnr$train(self$trainTask)
    pred = as.matrix(lrnr_trained$predict())

    Y = self$trainTask$Y
    bins = 15
    convertColumn = function(x){
      quants = seq(0,1,1/bins)
      q=quantile(x,quants)

      nearest <- findInterval(x, q)
      x <- q[nearest]
      return(x)
    }
    quantizer = function(X){as.matrix(apply(X, MARGIN = 2, FUN =convertColumn))}
    x = as.matrix(quantizer(pred))

    hal_fit = fit_hal(X=as.matrix(x), Y=as.vector(Y), family =ifelse(self$params$outcome_type == "continuous", "gaussian", "binomial" ),yolo=F, return_lasso = T, n_folds=10, max_degree=2)
    lambda_max = max(hal_fit$lambda)
    private$.lambda_max = lambda_max

  },
  setTasks = function(){
    bm = self$bm

    meta_data = self$meta_data

    split_indices = self$split_indices

    shared_dat = Shared_BigData$new(data=meta_data, bm=bm,force_copy = F)


    private$.trainTask = make_sl3_BigTask(data = shared_dat, bm_columns = self$base_columns, covariates = self$meta_columns, outcome = self$outcome_name, outcome_type= self$params$outcome_type, row_index = split_indices)

    if(is.null(split_indices)){
      private$.testTask = NULL
      return()
    }

    private$.testTask = make_sl3_BigTask(data = shared_dat, bm_columns = self$base_columns, covariates = self$meta_columns, outcome = self$outcome_name, outcome_type = self$params$outcome_type, row_index = setdiff(1:self$nrow, split_indices))

  }


),

active = list(
  params = function(){
    return(private$.params)
  },
  nrow = function(){
    return(nrow(self$meta_data))
  },
  meta_data = function(){
    return(self$params$meta_data)
  },
  bm = function(){
    return(attach.big.matrix(self$params$big_data))
  },
  fbm =function(){
    return(private$.trainTask$fbm)
  },
  meta_columns = function(){
    return(self$params$meta_predictors)
  },
  base_columns = function(){
    cols = self$params$base_predictors
    if(is.null(cols)){
      return(colnames(self$bm))
    }
    return(cols)
  },
  outcome_name = function(){
    return(self$params$outcome)
  },
  trainTask = function(){
    return(private$.trainTask)
  },
  testTask = function(){
    return(private$.testTask)
  },
  groups = function(){
    return(self$params$groups)
  },
  split_indices = function(){
    return(private$.split_indices)
  },
  lambda_max = function(){
    return(private$.lambda_max)
  }


),
public = list(
  initialize = function(meta_data, big_data, outcome = "Age.years", outcome_type = "continuous", split = 0.8, split_data = T, groups = NULL, compute_groups = F, num_groups = 20, select_optimal_number = F, selection_method = c("silhouette" , "wss", "gap_stat"),
                       kmeans_pre_cluster = F, base_predictors = NULL, meta_predictors = NULL, force_copy=T, big_data_backingfile= tempfile()) {
    big_data_backingfile_name = paste0(basename(big_data_backingfile), ".bk")
    big_data_backingfile_path = dirname(big_data_backingfile)

    #Turn data into data.table
    if(inherits(big_data, "big.matrix")){
      if(force_copy){
        big_data = bigmemory::deepcopy(big_data, backingfile = big_data_backingfile_name, backingpath = big_data_backingfile_path)
      }
    }
    else if (is.matrix(big_data) | inherits(big_data, "data.frame")) {

      big_data = as.big.matrix(as.matrix(big_data), backingfile = big_data_backingfile_name,backingpath = big_data_backingfile_path)

    }


    meta_data = as.data.table(meta_data)
    meta_data = meta_data[,c(base_predictors, outcome), with = F]


    params <- args_to_list()
    params$big_data = describe(big_data)
    private$.params = params
    if(split_data){
      private$splitData()
    }

    private$setTasks()
    private$setLambdaMax()
    if(compute_groups){
      private$computeGroups(num_groups,select_optimal_number, c(selection_method),kmeans_pre_cluster)
    }
  }



)

)
