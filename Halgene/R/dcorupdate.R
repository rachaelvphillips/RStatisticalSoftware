





#' @import  purrr
#' @import  delayed
#' @import  future
#' @import  SIS
#' @import  energy
#' @import  future
#' @import  sl3
#' @import  glmnet
#' @import  hal9001
#' @export
Lrnr_screener_dist_cor <- R6Class(
  classname = "Lrnr_screener_dist_cor",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(method = "backward_selection", num_to_select =30, cutoff = 0.01, verbose = F,...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial"),

    .train = function(task) {


      verbose <- getOption("sl3.verbose")
      params <- self$params
      method = params$method
      algorithm = method
      num_to_select = params$num_to_select
      outcome_type = self$get_outcome_type(task)$type
      verb = params$verbose
      fast_analysis = F
      cutoff = params$cutoff

      if(method == "forward_selection"){


        output = dcor_screen(as.data.frame(task$X), task$Y, num_to_select, search_set = 1:ncol(task$X), fast_analysis = F, perform_batch_forward_selection_based_filtering = F, perform_clustered_backward_selection = F, exact_algorithm = F, filter_search_set_to = 10000,  compute_pdcor_by_n_forward = 50, compute_pdcor_by_n_backward = 50, algorithm = method,  max_iterations = ifelse(algorithm=="forward_selection",  num_to_select, 25), 1,  dcor_retain_top_n = 1, pdcor_retain_top_n=1, rejection_dcor_cutoff = ifelse(ncol(task$X) >= 2*num_to_select, cutoff, 0), rejection_pdcor_cutoff =  ifelse(ncol(task$X) >= 2*num_to_select, cutoff, 0), increase_cutoff_by = 0, min_proportion_reduce = ifelse(ncol(task$X)>200,0.05,0), backward_selection_maximum_cluster_size= ifelse(fast_analysis,5,3), forward_selection_maximum_number_clusters = round(num_to_select/2), reduce_by_n_clusters = 1, cluster_by_dcor = F, recompute_clusters  = T, verbose = verb)
        }
        else if(method == "backward_selection"){

          output  = dcor_screen(as.data.frame(task$X), task$Y, num_to_select, search_set = 1:ncol(task$X), fast_analysis = F, perform_batch_forward_selection_based_filtering = T, perform_clustered_backward_selection = T, exact_algorithm = F, filter_search_set_to = 10000,  compute_pdcor_by_n_forward = 50, compute_pdcor_by_n_backward = 50, algorithm = "backward_selection", max_iterations = ifelse(algorithm=="forward_selection",  num_to_select, 25), selection_batch_size = 1,  dcor_retain_top_n = 5, pdcor_retain_top_n=1, rejection_dcor_cutoff = cutoff, rejection_pdcor_cutoff = cutoff, increase_cutoff_by = 0, min_proportion_reduce = ifelse(fast_analysis,0.05, ifelse(ncol(task$X) >= 2*num_to_select, 0.05, 0)), backward_selection_maximum_cluster_size= ifelse(fast_analysis,5,3), forward_selection_maximum_number_clusters = round(num_to_select/2), reduce_by_n_clusters = 1, cluster_by_dcor = F, recompute_clusters  = T, verbose = verb)
          }
      selected = output$final_selected_names
      fit_object = list()
      fit_object$selected = selected

      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      selected = self$fit_object$selected



      return(task$X[,selected, with = F])
    },
    .chain = function(task) {

      selected = self$fit_object$selected

      nexttask = task$next_in_chain(covariates = selected)
      return(nexttask)
    }
  )
)






#Selects num_to_select optimal variables for nonparametrically predicting y from the variables specified by search_set.
dcor_screen = function(X, y, num_to_select, search_set = 1:ncol(X), fast_analysis = F, perform_batch_forward_selection_based_filtering = T, perform_clustered_backward_selection = T, exact_algorithm = F, filter_search_set_to = 10000,  compute_pdcor_by_n_forward = ifelse(fast_analysis, ifelse(length(search_set)>= 10000, 50, 100), 250), compute_pdcor_by_n_backward = ifelse(fast_analysis, 150, 300), algorithm = c("backward_selection", "forward_selection", "batch_forward_selection"), initial_max_selected = ifelse(algorithm=="forward_selection", num_to_select, 7*num_to_select), max_iterations = ifelse(algorithm=="forward_selection",  num_to_select, 25), selection_batch_size = ifelse(algorithm=="forward_selection", 1,round(initial_max_selected/max_iterations)),  dcor_retain_top_n = 5, pdcor_retain_top_n=1, rejection_dcor_cutoff = 0.05, rejection_pdcor_cutoff = 0.05, increase_cutoff_by = 0, min_proportion_reduce = ifelse(fast_analysis,0.05,0.03), backward_selection_maximum_cluster_size= ifelse(fast_analysis,5,3), forward_selection_maximum_number_clusters = round(num_to_select/2), reduce_by_n_clusters = 1, cluster_by_dcor = F, recompute_clusters  = T, verbose = T){
  algorithm = c(algorithm)[1]
  #Obtain initial selection set
  final_selection_above_cutoff=NULL
  forward_selection = (algorithm == "batch_forward_selection")| (algorithm == "forward_selection")
  batch_forward_selection =  (algorithm == "batch_forward_selection")
  backward_selection = (algorithm == "backward_selection")

  if(forward_selection){
    if(algorithm == ("forward_selection")){
      selection_batch_size = 1

    }
    initial_max_selected = num_to_select
    max_iterations = ceiling(num_to_select/selection_batch_size)


  }

  if(exact_algorithm){
    min_proportion_reduce = -1
    rejection_dcor_cutoff=-1
    rejection_pdcor_cutoff = -1
    perform_clustered_backward_selection=F
  }
  if(verbose){
    print(paste0("Initial search set contains ", length(search_set), " variables"))
    print("Selecting initial set of variables and reducing search set...")
  }


  if(perform_batch_forward_selection_based_filtering | forward_selection ){
    first_out = batch_forward_selection(X, y, num_to_select, search_set, filter_search_set_to,  compute_pdcor_by_n_forward, selection_batch_size, dcor_retain_top_n, pdcor_retain_top_n, initial_max_selected, max_iterations, rejection_dcor_cutoff, rejection_pdcor_cutoff, increase_cutoff_by, min_proportion_reduce, verbose)

    selected= first_out$selected
    retain_for_later = first_out$retain_for_later
    if(batch_forward_selection){
      #After performing batch forward_selection as filtering do standard forward selection on the selected set.
      first_out = batch_forward_selection(X, y, num_to_select, selected, filter_search_set_to = filter_search_set_to,  compute_pdcor_by_n_forward, selection_batch_size = 1, dcor_retain_top_n, pdcor_retain_top_n, initial_max_selected = num_to_select, max_iterations = num_to_select, rejection_dcor_cutoff = -1, rejection_pdcor_cutoff = -1, increase_cutoff_by = -1, min_proportion_reduce = -1, verbose)
      selected= first_out$selected
      retain_for_later = first_out$retain_for_later
    }
  }
  else{
    selected = search_set
    retain_for_later = c()
  }
  first_selection = selected

  reduced_selected = selected

  if(length(reduced_selected) == num_to_select){
    results = list()
    results$final_selected_index = as.numeric(reduced_selected)
    results$final_selected_names = colnames(X)[as.numeric(reduced_selected)]
    results$final_selected_conditional_scores = as.vector(sapply(as.numeric(reduced_selected), function(j){sqrt(abs(pdcor(X[,j],y, X[,setdiff(as.numeric(reduced_selected), c(j))])))}))
    results$final_selected_marginal_scores = sapply(as.numeric(reduced_selected), function(j){sqrt(dcor2d(X[,j],y))})
    return(results)
  }
  if(verbose){
    print(paste0("Reducing selection set by clusters of (average) size ", backward_selection_maximum_cluster_size, "..." ))
  }

   if(perform_clustered_backward_selection){
     second_out = clustered_backward_selection(reduced_selected, X, y, num_to_select,  compute_pdcor_by_n_backward , rejection_dcor_cutoff ,rejection_pdcor_cutoff, backward_selection_maximum_cluster_size, reduce_by_n_clusters, cluster_by_dcor , recompute_clusters  , verbose,fast_analysis, forward_selection_maximum_number_clusters )
     reduced_selected = second_out$reduced_selected
     final_selection_above_cutoff = second_out$final_selection_above_cutoff
   }

  cluster_selected = reduced_selected


  if(verbose & perform_clustered_backward_selection){
    print(paste0("Cluster-based reduction process retained ", length(unique(reduced_selected)), " variables."))
    print("Now reducing selection set by one variable at a time until desired number of variables reached.")
  }
  if(verbose & perform_clustered_backward_selection){
    print(paste0("Adding ", length(retain_for_later), " of the variables retained from each iteration of the first selection process."))
  }

  reduced_selected = unique(c(retain_for_later, reduced_selected))



  entered = length(reduced_selected) !=0
  third_output= backward_selection(entered, X,y, reduced_selected, num_to_select, compute_pdcor_by_n_backward,rejection_dcor_cutoff, verbose,fast_analysis)

  if(is.null(final_selection_above_cutoff)){
    final_selection_above_cutoff = third_output$final_selection_above_cutoff
  }

  reduced_selected = third_output$reduced_selected


  if(verbose){
    print("Selection process finished.")
    print(paste0("We selected ", length(reduced_selected), " variables."))
    print("Computing final results...")
  }

  results = list()
  results$final_selected_index = as.numeric(reduced_selected)
  results$final_selected_names = colnames(X)[as.numeric(reduced_selected)]
  results$final_selected_conditional_scores = as.vector(sapply(as.numeric(reduced_selected), function(j){sqrt(abs(pdcor(X[,j],y, X[,setdiff(as.numeric(reduced_selected), c(j))])))}))
  results$final_selected_marginal_scores = sapply(as.numeric(reduced_selected), function(j){sqrt(dcor2d(X[,j],y))})
  results$final_pooled_score = dcor(X[, as.numeric(reduced_selected)],y)
  if(!is.null(final_selection_above_cutoff)){
    results$final_selection_above_cutoff_index = as.numeric(final_selection_above_cutoff)
    results$final_selection_above_cutoff_names = colnames(X)[as.numeric(final_selection_above_cutoff)]
    results$final_selection_above_cutoff_conditional_scores = as.vector(sapply(as.numeric(final_selection_above_cutoff), function(j){sqrt(abs(pdcor(X[,j],y, X[,setdiff(as.numeric(final_selection_above_cutoff), c(j))])))}))
    results$final_selection_above_cutoff_marginal_scores = sapply(as.numeric(final_selection_above_cutoff), function(j){sqrt(dcor2d(X[,j],y))})
    results$final_selection_above_cutoff_pooled_score = dcor(X[, as.numeric(final_selection_above_cutoff)],y)

    }

  if(perform_clustered_backward_selection){
    results$cluster_selected_index = as.numeric(cluster_selected)
    results$cluster_selected_names = colnames(X)[as.numeric(cluster_selected)]
    results$cluster_selected_marginal_scores = sapply(as.numeric(cluster_selected), function(j){sqrt(dcor2d(X[,j],y))})
  }


  results$initial_selection_index = as.numeric(first_selection)
  results$initial_selection_names = colnames(X)[as.numeric(first_selection)]
  results$initial_selected_marginal_scores = sapply(as.numeric(first_selection), function(j){sqrt(dcor2d(X[,j],y))})

  results$initial_retained_each_iter_index = as.numeric(retain_for_later)
  results$initial_retained_each_iter_names = colnames(X)[as.numeric(retain_for_later)]
  results$initial_retained_each_iter_marginal_scores = sapply(as.numeric(retain_for_later), function(j){sqrt(dcor2d(X[,j],y))})


  return(results)
}

#Compute dcor between each variable (column of X) specified by search_set and outcome y
#Computation time is reduced by computing the correlation for a randomly selected partitioning of the sample (size n) into groups of size by_n...
#The correlation is computed for each group and then averaged to obtain a final estimate...
#This reduces the complexity from n^2 to (by_n^2) * n/by_n = n*by_n.
compute_dcor = function(X, search_set, y, by_n){
  search_set = as.numeric(search_set)
  if(is.null(by_n)){
    return(apply(X[,search_set], 2, function(u){sqrt(dcor2d(u,y))}))
  }
  num_groups = max(round(nrow(X) / by_n),1)
  grp_assignment = sample(1:num_groups, size = nrow(X), replace = T)
  grp_scores = do.call(cbind, lapply(unique(grp_assignment), function(grp){
    row_id = which(grp_assignment==grp)
    compute_dcor(X[row_id,], search_set, y[row_id], NULL)
  }))
  return(rowMeans(grp_scores))

}

#Compute pdcor between each variable (column of X) specified by search_set and outcome y conditional on the variables specified by the set selected.
#Computation time is reduced by computing the correlation for a randomly selected partitioning of the sample (size n) into groups of size by_n...
#The correlation is computed for each group and then averaged to obtain a final estimate...
#This reduces the complexity from n^2 to (by_n^2) * n/by_n = n*by_n.
compute_pdcor = function(X, search_set, y, selected,  by_n, neg = F){
  search_set = as.numeric(search_set)
  selected = as.numeric(selected)

  if(is.null(by_n) & !neg){

    M = X[,search_set]
    if(!is.matrix(M)){
      M = as.matrix(M, ncol = 1)
    }
    return(apply(M, 2, function(u){sqrt(abs(pdcor(u,y, X[,selected])))}))
  }
  if(is.null(by_n) & neg){

    M = X[,search_set]
    if(!is.matrix(M)){
      M = as.matrix(M, ncol = 1)
    }
    return(apply(M, 2, function(u){pdcor(u,y, X[,selected])}))
  }

  num_groups = max(round(nrow(X) / by_n),1)
  grp_assignment = sample(1:num_groups, size = nrow(X), replace = T)
  grp_scores = do.call(cbind, lapply(unique(grp_assignment), function(grp){
    row_id = which(grp_assignment==grp)
    compute_pdcor(X[row_id,], search_set, y[row_id], selected,  by_n=NULL, neg = T)
  }))
  return(sqrt(abs(rowMeans(grp_scores))))

}

compute_pdcor_mat = function(X, y, z,  by_n, neg = F){
  if(!is.matrix(X)){
    X = as.matrix(X, ncol = 1)
  }

  if(is.null(by_n) & !neg){

    return(sqrt(abs(pdcor(X,y, z))))
  }
  if(is.null(by_n) & neg){

    return(pdcor(X,y, z))
  }
  num_groups = max(round(nrow(X) / by_n),1)

  grp_assignment = sample(1:num_groups, size = nrow(X), replace = T)
  grp_scores = unlist( lapply(unique(grp_assignment), function(grp){
    row_id = which(grp_assignment==grp)
    return(compute_pdcor_mat(X[row_id,], y[row_id], z[row_id,],  by_n=NULL, neg = T))
  }))

  return(sqrt(abs(mean(grp_scores))))

}

compute_dcor_2d = function(x, y, by_n){


  if(is.null(by_n)){

    return(sqrt(abs(dcor2d(x,y))))
  }
  num_groups = max(round(length(x) / by_n),1)
  grp_assignment = sample(1:num_groups, size = length(x), replace = T)
  grp_scores = unlist( lapply(unique(grp_assignment), function(grp){
    row_id = which(grp_assignment==grp)
    return(compute_dcor_2d(x[row_id], y[row_id],  by_n=NULL))
  }))

  return(mean(grp_scores))

}

#Clusters the variables (using hclust) in selected by their pearson or distance correlation.
#The number of clusters chosen is length(selected)/num_per_cluster.
clusterByCor = function(dat, num_per_cluster, selected, by_dcor = F, by_n, cluster_groups_only = F){

  if(num_per_cluster >= selected){
    if(cluster_groups_only){
      return(rep("1", length(selected)))
    }
    return(list(dat[,selected]))
  }

  if(is.null(num_per_cluster) | num_per_cluster==0){
    return(NULL)
  }

  if(length(selected)==1){
    if(cluster_groups_only){
      return("1")
    }
    return(list(as.matrix(dat[,selected], ncol=1)))
  }

  selected = as.numeric(selected)
  dat = as.data.frame(as.matrix(dat))
  dat = dat[,selected]

  if( is.vector(dat) | ncol(dat) ==1){

    dat = matrix(dat, ncol = 1)
  }



  colnames(dat) = c(selected)

  k = max(floor(length(selected)/num_per_cluster),1)
  if(by_dcor){
    cor_mat=as.matrix(outer(1:ncol(dat),1:ncol(dat),function(i,j){compute_dcor_2d(dat[,i],dat[,j], by_n)}))
  }
  else{

    cor_mat = cor(dat)
    if(any(is.na(cor_mat)) | any(is.infinite(cor_mat))){

    }
  }

  cor_dist = as.dist(1-abs(cor_mat))
  hclust_fit = hclust(cor_dist)
  keep_ind = cutree(hclust_fit, k)
  if(cluster_groups_only){
    return(keep_ind)
  }
  lstOfMats = lapply(split(as.data.frame(t(dat)), keep_ind), function(m){as.data.frame(t(m))})
  return(lstOfMats)
}



batch_forward_selection = function(X, y, num_to_select, search_set, filter_search_set_to,  compute_pdcor_by_n, selection_batch_size, dcor_retain_top_n, pdcor_retain_top_n, initial_max_selected, max_iterations, rejection_dcor_cutoff, rejection_pdcor_cutoff, increase_cutoff_by, min_proportion_reduce, verbose){
  retain_for_later = c()
  selected = c()
  for(i in 1:max_iterations){



    len = length(search_set)
    if(len <= num_to_select){
      min_proportion_reduce = 0
    }
    if(length(search_set)<=1){
      selected = c(selected, search_set)
      if(verbose){
        print("Search set has been reduced to 0. Stopping...")
      }
      break
    }
    if(length(selected)>= initial_max_selected){
      if(verbose){
        print("Limit initial_max_selected reached. Stopping...")
      }
      break
    }
    if(length(selected)==0){
      vals = compute_dcor(X, search_set,y,NULL) #apply(X[,search_set], 2, function(u){sqrt(dcor2d(u,y))})
      retain_for_later = c(retain_for_later,search_set[which(rank(-vals) <= dcor_retain_top_n)])
      if(verbose){
        expr1 <- quote(print("Top retained dcor values are: "))
        expr2<- quote(print(vals[which(rank(-vals) <= dcor_retain_top_n)]))
      }

    }
    else{
      vals = compute_pdcor(X, search_set,y,selected , compute_pdcor_by_n)    #apply(X[,search_set], 2, function(u){abs(pdcor(u,y, X[,selected]))})
      retain_for_later = c(retain_for_later,search_set[which(rank(-vals) <= pdcor_retain_top_n)])
      if(verbose){
        expr1 <- quote(print("Top retained pdcor values are: "))
        expr2<- quote(print(vals[which(rank(-vals) <= pdcor_retain_top_n)]))
      }
    }
    choose_ind = which(rank(-vals,ties.method = "random") <= selection_batch_size)
    choose_vals = vals[choose_ind]
    top_ind = search_set[choose_ind]


    if(!is.null(rejection_dcor_cutoff)){
      if(min(vals)> rejection_dcor_cutoff){
        keep = which(vals >= rejection_dcor_cutoff + increase_cutoff_by)

      }
      else{
        keep = which(vals >= rejection_dcor_cutoff)
      }
    }
    else{
      keep = 1:length(vals)
    }
    rejection_dcor_cutoff=rejection_pdcor_cutoff



    if(min_proportion_reduce>0 & !is.null(min_proportion_reduce) & (length(search_set)-length(keep))/(length(search_set)-selection_batch_size) < min_proportion_reduce)
    {

      keep = which(vals > quantile(vals,min_proportion_reduce ))

    }

    if(!is.null(filter_search_set_to) & length(search_set)-length(keep)> filter_search_set_to){
      keep = which(rank(vals) <= length(search_set) - filter_search_set_to)
    }
    reduce_to = search_set[keep]

    selected = c(selected, top_ind)
    search_set = setdiff(reduce_to, selected)

    if(verbose){
      print(paste0("Iter: ", i, ". Current number of batch-selected variables: ", length(selected), "."))
      print(paste0("Search set reduced to ", length(search_set), " variables."))
      eval(expr1)
      eval(expr2)
      qs = t(as.matrix(quantile(vals)))
      rownames(qs) = "(p)dcor quantiles"
      print(qs)
    }


  }

  first_selection = selected
  if(verbose){
    print(paste0("First selection process selected ", length(selected), " variables."))
  }
  output = list()
  output$selected = selected
  output$retain_for_later = retain_for_later
  return(output)
}



clustered_forward_selection = function(reduced_selected, X, y, num_to_select,  compute_pdcor_by_n , rejection_dcor_cutoff ,rejection_pdcor_cutoff, forward_selection_maximum_number_clusters, reduce_by_n_clusters, cluster_by_dcor , recompute_clusters  , verbose,fast_analysis ){


  check_names = reduced_selected
  reduce_by_n_clusters=1

  maximum_cluster_size = round(length(reduced_selected)/forward_selection_maximum_number_clusters)
  cluster_groups = as.vector(unlist(clusterByCor(X, maximum_cluster_size, reduced_selected, cluster_by_dcor,compute_pdcor_by_n, cluster_groups_only = T)))
  iter = max(20,ceiling(length(reduced_selected)/maximum_cluster_size) )
  search_set = reduced_selected
  selected = c()
  if(verbose){
    print("Performing cluster forward.")

  }
  for(i in 1:iter){

    len = length(selected)

    num_mat = length(unique(cluster_groups))


    if(len >= 2*num_to_select | length(search_set) <2){
      if( length(search_set) <2){
        selected = c(selected, search_set)
      }
      vals = compute_pdcor(X, search_set, y, selected, compute_pdcor_by_n)
      selected = c(selected, search_set[which(rank(-vals)<=max(10,num_to_select))])
      break
    }

    if(length(selected)==0){
      vals =  compute_dcor(X, search_set, y, compute_pdcor_by_n)
    }
    else{

      vals = compute_pdcor(X, search_set, y, selected, compute_pdcor_by_n)
    }
    vals = as.vector(vals)

    bestInGroup = function(grp, vals){
      choice = as.vector((vals == max(vals[cluster_groups ==grp]))  & (cluster_groups ==grp))

      choice = which(choice)


      if(length(choice)==0){


      }
      if(length(c(choice)) >1){
        choice = choice[1]
      }
      return(choice)
    }
    new_best_ind = unlist(sapply(unique(cluster_groups), bestInGroup, vals = vals))

    new_best = search_set[new_best_ind]
    selected = c(selected, new_best)
    if(verbose){
      print(quantile(vals[new_best_ind]))

    }
    search_set = setdiff(search_set,new_best)
    maximum_cluster_size = length(reduced_selected)/forward_selection_maximum_number_clusters
    cluster_groups = as.vector(unlist(clusterByCor(X, maximum_cluster_size, search_set, cluster_by_dcor,compute_pdcor_by_n, cluster_groups_only = T)))


  }
  output = list()
  output$selected = as.numeric(selected)

  if(!(all(reduced_selected %in% check_names))){
    stop("Something went wrong with clustering.")
  }
  return(output)



}





clustered_backward_selection = function(reduced_selected, X, y, num_to_select,  compute_pdcor_by_n , rejection_dcor_cutoff ,rejection_pdcor_cutoff, backward_selection_maximum_cluster_size, reduce_by_n_clusters, cluster_by_dcor , recompute_clusters  , verbose,fast_analysis, forward_selection_maximum_number_clusters ){
  check_names = reduced_selected
  maximum_cluster_size=backward_selection_maximum_cluster_size
  reduce_by_n_clusters=1
  mat_list = clusterByCor(X, maximum_cluster_size, reduced_selected, cluster_by_dcor,compute_pdcor_by_n)
  iter = ceiling(length(reduced_selected)/maximum_cluster_size) + 25
  entered = F
  final_selection_above_cutoff = NULL
  for(i in 1:iter){

    len = length(reduced_selected)
    num_mat = length(mat_list)

    if(len<=num_to_select){
      print("breakk")
      break
    }
    one_out_cor = function(j){
      return( compute_pdcor_mat(as.matrix(mat_list[[j]]), y, as.matrix(do.call(cbind, mat_list[-j])), compute_pdcor_by_n))           #abs(pdcor(as.matrix(mat_list[[j]]),y,as.matrix(do.call(cbind, mat_list[-j])))))
    }
    one_out_cor = Vectorize(one_out_cor)
    vals=one_out_cor(1:num_mat)

    min_ind = which(rank(vals, ties.method = "random") <= reduce_by_n_clusters)
    if(verbose){
      print(paste0("Minimum value is ", vals[min_ind]))
    }
    num_removed = sum(sapply(min_ind, function(i){ncol(mat_list[[i]])}))
    num=reduce_by_n_clusters
     n=0
    while(num_removed<backward_selection_maximum_cluster_size& n <=5){
     n = n+1
    min_ind = c(which(rank(vals, ties.method = "random") <= reduce_by_n_clusters + n))
     num_removed = sum(sapply(min_ind, function(i){ncol(mat_list[[i]])}))
    }

    if(len - num_removed < num_to_select){
      if(reduce_by_n_clusters==1){
        if(verbose){
          print("Next selection would have gone below num_to_select. Stopping.")

        }
        break
      }
      else{
        reduce_by_n_clusters = reduce_by_n_clusters -1
      }

    }

    if(max(vals[min_ind]) <= rejection_pdcor_cutoff | entered | (length(mat_list) >= num_to_select)){
      if(length(c(min_ind))>1){
        prev_names = unlist(lapply(mat_list[min_ind], colnames))
        removed_mat = as.matrix(do.call(cbind, mat_list[min_ind]))
        colnames(removed_mat) = as.character(prev_names)

      }
      else{
        removed_mat = mat_list[[min_ind]]
      }

      old_names = colnames(removed_mat)
      mat_list = mat_list[-min_ind]

      if(ncol(removed_mat) > maximum_cluster_size){

        inner_scores = apply(removed_mat, 2, function(x){compute_dcor_2d(x,y, 30)})
        keep_ind = which(rank(-inner_scores, ties.method = "random")<= (ncol(removed_mat) - maximum_cluster_size))
        removed_mat = removed_mat[, keep_ind]
        if(is.vector(removed_mat)){
          removed_mat = matrix(removed_mat, ncol=1)
        }
        colnames(removed_mat) = c(old_names[keep_ind])
        mat_list = c(mat_list, list(removed_mat))


      }
      reduced_selected = unlist(lapply(mat_list, colnames))
    }
    else if (!entered){
      if(verbose){
        print("Cluster with minimal pdcor is above rejection cutoff. Storing current selection and now performing cluster-based forward selection.")

      }
      final_selection_above_cutoff = reduced_selected
      entered=T


      break

      if(fast_analysis){
        compute_pdcor_by_n = max(75,compute_pdcor_by_n/2)
      }

      if(length(c(min_ind))>1){
        prev_names = unlist(lapply(mat_list[min_ind], colnames))
        removed_mat = as.matrix(do.call(cbind, mat_list[min_ind]))
        colnames(removed_mat) = as.character(prev_names)

      }
      else{
        removed_mat = mat_list[[min_ind]]
      }

      old_names = colnames(removed_mat)
      mat_list = mat_list[-min_ind]

      inner_scores = apply(removed_mat, 2, function(x){compute_dcor_2d(x,y, 30)})
      keep_ind = which(rank(-inner_scores, ties.method = "random")<= (ncol(removed_mat) - maximum_cluster_size))
      removed_mat = removed_mat[, keep_ind]
      if(is.vector(removed_mat)){
        removed_mat = matrix(removed_mat, ncol=1)
      }
      colnames(removed_mat) = c(old_names[keep_ind])
      mat_list = c(mat_list, list(removed_mat))
      reduced_selected = unlist(lapply(mat_list, colnames))

    }

    if(length(reduced_selected) < 2*num_to_select + 10){
      break
    }



    if(verbose){
      print(paste0("Iteration #", i, " of cluster-based reduction removed ", len - length(reduced_selected), " variables. ", "We retained ", length(reduced_selected), " variables."))

      qs = t(as.matrix(quantile(vals)))
      rownames(qs) = "(p)dcor quantiles"
      print(qs)
    }
    if(recompute_clusters){

      mat_list = clusterByCor(X, maximum_cluster_size, as.numeric(reduced_selected),cluster_by_dcor, compute_pdcor_by_n)
    }


  }


  if(entered){
    output_forward = clustered_forward_selection(reduced_selected, X, y, num_to_select,  compute_pdcor_by_n , rejection_dcor_cutoff ,rejection_pdcor_cutoff, forward_selection_maximum_number_clusters, reduce_by_n_clusters, cluster_by_dcor , recompute_clusters  , verbose,fast_analysis )
    reduced_selected = output_forward$selected
  }


  output = list()
  output$reduced_selected = reduced_selected
  output$final_selection_above_cutoff = final_selection_above_cutoff

  if(!(all(reduced_selected %in% check_names))){
    stop("Something went wrong with clustering.")
  }
  return(output)



}


backward_selection = function(entered, X,y, reduced_selected, num_to_select, compute_pdcor_by_n,rejection_dcor_cutoff, verbose, fast_analysis){
  final_selection_above_cutoff = NULL

  while(length(reduced_selected) > num_to_select){
    new_dat = X[,as.numeric(reduced_selected)]
    one_out_cor = function(j){
      return(compute_pdcor_mat(new_dat[,j], y, new_dat[,-j], compute_pdcor_by_n))               #abs(pdcor(new_dat[,j],y,new_dat[,-j])))
    }
    one_out_cor = Vectorize(one_out_cor)
    vals = one_out_cor(1:length(reduced_selected))
    reduced_selected = reduced_selected[-which.min(vals)]
    if(!entered & min(vals) >= rejection_dcor_cutoff){
      print("In one-left-out reduction process, no variables reached the rejection cutoff. Storing current variable selection and now continuing rejection process. ")
      final_selection_above_cutoff = reduced_selected
      entered = T
      if(fast_analysis){
        compute_pdcor_by_n = max(75,compute_pdcor_by_n/2)
      }
    }
    if(verbose){
      print(paste0("Iteration finished. Reduced selection set by one. We retained ", length(reduced_selected), " variables."))
      qs = t(as.matrix(quantile(vals)))
      rownames(qs) = "(p)dcor quantiles"
      print(qs)
    }

  }
  output = list()
  output$final_selection_above_cutoff = final_selection_above_cutoff
  output$reduced_selected = reduced_selected
  return(output)
}







