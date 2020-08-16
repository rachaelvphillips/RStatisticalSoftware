generate_likelihood <- function(task, lrnr_binomial = Lrnr_glm$new(), lrnr_continuous = Lrnr_condensier$new() ){
  factor_list = c()
  for(name in names(task$npsem)){
    node <- task$npsem[[name]]

    if(name =="L0"){
      factor_list <- c(factor_list, name = LF_emp$new(name ))
    }
    else if(node$variable_type$type == "binomial"){
      factor_list <- c(factor_list, name = LF_fit_pooled$new(name,  lrnr_binomial))
    }
    else{
      factor_list <- c(factor_list, name = LF_fit_pooled$new(name, lrnr_continuous ))
    }


  }
  return(Likelihood_pooled$new(factor_list))
}

generate_likelihood_surv <- function(task, lrnr_binomial = Lrnr_glm$new(), lrnr_continuous = Lrnr_condensier$new() ){
  factor_list = c()
  for(name in names(task$npsem)){
    node <- task$npsem[[name]]

    if(name =="W"){
      factor_list <- c(factor_list, name = LF_emp$new(name ))
    }
    else if(node$variable_type$type == "binomial"){
      factor_list <- c(factor_list, name = LF_fit_pooled$new(name,  lrnr_binomial))
    }
    else{
      factor_list <- c(factor_list, name = LF_fit_pooled$new(name, lrnr_continuous ))
    }


  }
  return(Likelihood_pooled$new(factor_list))
}


generate_npsem <- function(data, baseline, time_dep_cov, trtment, outcome){

  num_time_points <- length(unique(data$t))
  node_names <-  get_node_names(num_time_points)
  node_list <- lapply(seq_along(node_names), get_node, node_names,baseline, trtment, time_dep_cov, outcome )
  return(node_list)
}


get_node = function(i, node_names, baseline, trtment, time_dep_cov, outcome){
  at_risk = NULL
  if(i==1 | i==2 |i == 3){
    #baseline = c(baseline, time_dep_cov)
  }
  if(i==1){
    parents = c()
    cov = baseline
  }
  else{
    parents = node_names[1:(i-1)]
    if(i%%2==0){
      cov = trtment
      at_risk = make_summary_measure_last_value("A_at_risk")
    } else {
      cov = time_dep_cov
    }
  }
  time =  as.numeric(stringr::str_match(node_names[i],"[^0-9].*([0-9]+)")[2])+1
  if(node_names[i] == "Y"){
    cov = outcome
    time = as.numeric(stringr::str_match(node_names[i-1],"[^0-9].*([0-9]+)")[2])+1
  }

  if(!is.null(baseline) & length(c(baseline)) >0){
    base_line_summary <- make_summary_measure_baseline(baseline)
  }
  else{
    base_line_summary <- c()
  }
  time_dep_summary <- make_summary_measure_last_value(c(trtment, time_dep_cov))
  if(i==1){
    summaries = make_summary_measure_NULL()
  }
  else if(i==2){
    summaries = base_line_summary
  }
  else if(i==3){
    summaries =  c( base_line_summary, make_summary_measure_last_value(trtment, strict_past = F))
  }
  else{
    summaries = c(base_line_summary, time_dep_summary)
  }

  return(define_lnode(node_names[i], cov, parents,time,summaries, at_risk_summary_function = at_risk))

}


get_node_names <- function(num_time_points){

  get_time_point = function(i){
    return(c(paste0("L",i), paste0("A",i)))
  }
  ltmle_nodes = c(unlist(lapply(seq_len(num_time_points)-1, get_time_point)), "Y")
  return(ltmle_nodes)
}









