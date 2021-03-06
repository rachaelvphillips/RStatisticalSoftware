

#Selects num_to_select optimal variables for nonparametrically predicting y from the variables specified by search_set.

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
dcor_screen = function(X, y, num_to_select, search_set = 1:ncol(X), perform_filtering = T, filter_search_set_to = 10000,  compute_pdcor_by_n = 50, algorithm = c("leave_one_out", "forward_selection"), selection_batch_size = min(length(search_set)/20,12), dcor_retain_top_n = 5, pdcor_retain_top_n=1, initial_max_selected = ifelse(selection_batch_size==1, num_to_select, 15*num_to_select), max_iterations = min(length(search_set)/selection_batch_size,num_to_select), rejection_dcor_cutoff = 0.05, rejection_pdcor_cutoff = 0.1, increase_cutoff_by = 0.01, min_proportion_reduce = 0.05, reduce_clusters_of_size = 5, reduce_by_n_clusters = 1, cluster_by_dcor = F, recompute_clusters  = T, verbose = T){
  retain_for_later = c()
  algorithm = c(algorithm)[1]
  #Obtain initial selection set
  selected = c()
  forward_selection = algorithm == "forward_selection"
  if(forward_selection){
    selection_batch_size=1
    max_iterations = num_to_select
    reduce_clusters_of_size = 0
    min_proportion_reduce = -1
    rejection_dcor_cutoff=-1
  }
  if(verbose){
    print(paste0("Initial search set contains ", length(search_set), " variables"))
    print("Selecting initial set of variables and reducing search set...")
  }

  for(i in 1:max_iterations){
    print(" ")


    len = length(search_set)
    if(length(search_set)<=1){
      selected = c(selected, search_set)
      if(verbose){
        print("Search set has been reduced to 0. Stopping...")
      }
      break
    }
    if(length(selected)> initial_max_selected){
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

      top_ind = search_set[which(rank(-vals,ties.method = "random") <= selection_batch_size)]

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



      if(!forward_selection & !is.null(min_proportion_reduce) & (length(search_set)-length(keep))/(length(search_set)-selection_batch_size) < min_proportion_reduce)
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
    print(paste0("Reducing selection set by clusters of (average) size ", reduce_clusters_of_size, "..." ))
  }


  mat_list = clusterByCor(X, reduce_clusters_of_size, reduced_selected, cluster_by_dcor)
  iter = ceiling(length(selected)/reduce_clusters_of_size) + 10
  entered = F
  final_selection_above_cutoff = NULL
  for(i in 1:iter){

    len = length(reduced_selected)
    num_mat = length(mat_list)
    print(len)
    print(num_mat)
    if(len<=num_to_select){
      break
    }
    one_out_cor = function(j){
      return( compute_pdcor_mat(as.matrix(mat_list[[j]]), y, as.matrix(do.call(cbind, mat_list[-j])), compute_pdcor_by_n))           #abs(pdcor(as.matrix(mat_list[[j]]),y,as.matrix(do.call(cbind, mat_list[-j])))))
    }
    one_out_cor = Vectorize(one_out_cor)
    vals=one_out_cor(1:num_mat)

    min_ind = c(which(rank(-vals) <= reduce_by_n_clusters))

    num_removed = sum(sapply(min_ind, function(i){ncol(mat_list[[i]])}))
    num=reduce_by_n_clusters
    n=0
    while(num_removed<reduce_clusters_of_size & n <=5){
      n = n+1
      min_ind = c(which(rank(-vals) <= reduce_by_n_clusters + n))
      num_removed = sum(sapply(min_ind, function(i){ncol(mat_list[[i]])}))
    }

    if(len - num_removed < num_to_select){
      if(reduce_by_n_clusters==1){
        print("Next selection would have gone below num_to_select. Stopping.")
        break
      }
      else{
        reduce_by_n_clusters = reduce_by_n_clusters -1
      }

    }
    if(max(vals[min_ind]) <= rejection_dcor_cutoff | entered){
      mat_list = mat_list[-min_ind]
    }
    else if (!entered){
      print("Cluster with minimal pdcor is above rejection cutoff. Storing current selection and then continuing with removal until 2*num_to_select retained.")
      final_selection_above_cutoff = reduced_selected
      entered=T
      compute_pdcor_by_n=compute_pdcor_by_n/2
      print("Halving compute_pdcor_by_n for computational speed up.")
    }
    if(entered & length(reduced_selected) < 2*num_to_select){
      break
    }


    reduced_selected = unlist(lapply(mat_list, colnames))
    if(verbose){
      print(paste0("Iteration #", i, " of cluster-based reduction removed ", len - length(reduced_selected), " variables. ", "We retained ", length(reduced_selected), " variables."))

      qs = t(as.matrix(quantile(vals)))
      rownames(qs) = "(p)dcor quantiles"
      print(qs)
    }
    if(recompute_clusters){
      mat_list = clusterByCor(X, reduce_clusters_of_size, reduced_selected,cluster_by_dcor)
    }


  }
  if(verbose){
    print(paste0("Adding ", length(retain_for_later), " of the variables retained from each iteration of the first selection process."))
  }
  reduced_selected = union(retain_for_later, reduced_selected)
  if(verbose){
    print(paste0("Cluster-based reduction process retained ", length(reduced_selected), " variables."))
    print("Now reducing selection set by one variable at a time until desired number of variables reached.")
  }
  cluster_selected = reduced_selected

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
    }
    if(verbose){
      print(paste0("Iteration finished. Reduced selection set by one. We retained ", length(reduced_selected), " variables."))
      qs = t(as.matrix(quantile(vals)))
      rownames(qs) = "(p)dcor quantiles"
      print(qs)
    }

  }
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

  if(!is.null(final_selection_above_cutoff)){
    results$final_selection_above_cutoff_index = as.numeric(final_selection_above_cutoff)
    results$final_selection_above_cutoff_names = colnames(X)[as.numeric(final_selection_above_cutoff)]
    results$final_selection_above_cutoff_conditional_scores = as.vector(sapply(as.numeric(final_selection_above_cutoff), function(j){sqrt(abs(pdcor(X[,j],y, X[,setdiff(as.numeric(final_selection_above_cutoff), c(j))])))}))
    results$final_selection_above_cutoff_marginal_scores = sapply(as.numeric(final_selection_above_cutoff), function(j){sqrt(dcor2d(X[,j],y))})
  }


  results$cluster_selected_index = as.numeric(cluster_selected)
  results$cluster_selected_names = colnames(X)[as.numeric(cluster_selected)]
  results$cluster_selected_marginal_scores = sapply(as.numeric(cluster_selected), function(j){sqrt(dcor2d(X[,j],y))})

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
compute_pdcor = function(X, search_set, y, selected,  by_n){
  search_set = as.numeric(search_set)


  if(is.null(by_n)){
    M = X[,search_set]
    if(!is.matrix(M)){
      M = as.matrix(M, ncol = 1)
    }
    return(apply(M, 2, function(u){sqrt(abs(pdcor(u,y, X[,selected])))}))
  }
  num_groups = max(round(nrow(X) / by_n),1)
  grp_assignment = sample(1:num_groups, size = nrow(X), replace = T)
  grp_scores = do.call(cbind, lapply(unique(grp_assignment), function(grp){
    row_id = which(grp_assignment==grp)
    compute_pdcor(X[row_id,], search_set, y[row_id], selected,  by_n=NULL)
  }))
  return(rowMeans(grp_scores))

}

compute_pdcor_mat = function(X, y, z,  by_n){
  if(!is.matrix(X)){
    X = as.matrix(X, ncol = 1)
  }

  if(is.null(by_n)){

    return(sqrt(abs(pdcor(X,y, z))))
  }
  num_groups = max(round(nrow(X) / by_n),1)
  grp_assignment = sample(1:num_groups, size = nrow(X), replace = T)
  grp_scores = unlist( lapply(unique(grp_assignment), function(grp){
    row_id = which(grp_assignment==grp)
    return(compute_pdcor_mat(X[row_id,], y[row_id], z[row_id,],  by_n=NULL))
  }))

  return(mean(grp_scores))

}


compute_dcor_2d = function(x, y, by_n){


  if(is.null(by_n)){

    return(sqrt(abs(dcor2d(x,y))))
  }
  num_groups = max(round(nrow(X) / by_n),1)
  grp_assignment = sample(1:num_groups, size = nrow(X), replace = T)
  grp_scores = unlist( lapply(unique(grp_assignment), function(grp){
    row_id = which(grp_assignment==grp)
    return(compute_pdcor_mat(x[row_id], y[row_id],  by_n=NULL))
  }))

  return(mean(grp_scores))

}

#Clusters the variables (using hclust) in selected by their pearson or distance correlation.
#The number of clusters chosen is length(selected)/num_per_cluster.
clusterByCor = function(dat, num_per_cluster, selected, by_dcor = F, by_n){
  if(is.null(num_per_cluster) | num_per_cluster==0){
    return(NULL)
  }
  selected = as.numeric(selected)
  dat = dat[,selected]
  colnames(dat) = (selected)
  k = floor(length(selected)/num_per_cluster)
  if(by_dcor){
    cor_mat=as.matrix(outer(selected,selected,function(i,j){compute_dcor_2d(dat[,i],dat[,j], by_n)}))
  }
  else{
    cor_mat = cor(dat)
  }
  cor_dist = as.dist(1-abs(cor_mat))
  hclust_fit = hclust(cor_dist)
  keep_ind = cutree(hclust_fit, k)
  lstOfMats = lapply(split(as.data.frame(t(dat)), keep_ind), function(m){as.data.frame(t(m))})
  return(lstOfMats)
}









