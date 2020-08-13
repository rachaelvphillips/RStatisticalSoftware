generate_likelihood <- function(task, lrnr_binomial = Lrnr_glm$new(family = "binomial"), lrnr_continuous = Lrnr_condensier$new() ){
  factor_list = c()
  for(name in names(task$npsem)){
    node <- task$npsem[[name]]

    if(name =="L0"){
      factor_list <- c(factor_list, name = LF_emp$new(name ))
    }
    else if(node$variable_type$type == "binomial"){
      factor_list <- c(factor_list, name = LF_fit$new(name,  lrnr_binomial))
    }
    else{
      factor_list <- c(factor_list, name = LF_fit$new(name, lrnr_continuous ))
    }


  }
  return(Likelihood$new(factor_list))
}
