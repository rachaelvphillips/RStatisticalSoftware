
#' @importFrom resample colVars
#' @export
get_next_beta <- function(x,y, beta, last_risk, eps_max_scale = 1, last_eps = NULL, iter, weights){
  t = proc.time()

  link = (y - as.vector(x %*% beta))/ length(y)
  #x_summed = Matrix::colSums(x * link)
  if(iter==1){
    weights = 1

  }
  x_summed= Matrix::crossprod(x, link)
  grad = x_summed * beta
  can_grad = grad - (sum(grad * abs(beta))/ sum(beta^2))*abs(beta)

  if(is.null(last_eps) | (iter<50 & iter%%10 ==0) ){
    eps_max = 1/max(abs(can_grad))/eps_max_scale
    search_set = exp(seq(log(eps_max), log(.0000000001), length.out = 10))
  }
  else if (iter==200){
    eps_max = last_eps
    search_set = c( eps_max, eps_max/2, eps_max/4)
  }
  else{
    eps_max = last_eps
    #search_set = c(2*eps_max, eps_max, eps_max/2)
    search_set = c(eps_max)
    if(iter%%5!=1){
      return(list(Inf, Inf, (1 + eps_max*can_grad)*beta,  eps_max, weights))
    }
  }

  best_beta = NULL
  best_risk = last_risk
  best_eps = NULL
  t2 = proc.time()
  grid_results = lapply(search_set, function(eps){


    new_beta = (1 + eps*can_grad)*beta
    risk = mean((y - as.vector(x %*% new_beta))^2)
    if(risk < best_risk){

      best_eps <<- eps
      best_beta <<- new_beta
      best_risk <<- risk
    }
    return()
  })

  if(is.null(best_beta)){
    best_eps = min(search_set)/2
    if(iter>0){
      best_beta = get_next_beta(x,y, beta, last_risk, eps_max_scale = 1, last_eps = best_eps, iter=-1, weights)[[3]]

    }
    else if(iter <=-10 ){
      warning("Did you already converge? Current beta vector could not be improved.")

      best_beta = beta

    }
    else if(iter <0 ){
      print("here")
      best_beta = get_next_beta(x,y, beta, last_risk, eps_max_scale = 1, last_eps = best_eps, iter= iter -1, weights)[[3]]

    }

  }

  #EIC_norm = sqrt(sum((can_grad*weights)^2)/length(can_grad))
  EIC_norm = Inf

  return(list(EIC_norm, best_risk, best_beta,  best_eps, weights))
}


#' @export
run_descent_gaussian <- function(x,y, beta, max_iter = 250, verbose = F){
  thresh = min(1/sqrt(length(y))/log(length(y)),.1)
  iter = 1
  risks = c()
  norms = c()
  new_beta = beta
  risk = mean(( y- as.vector(x %*% beta))^2)
  eps_max_scale = 1
  best_eps = NULL
  weights = NULL
  while(iter < max_iter){

    new_beta_lst = get_next_beta(x,y, new_beta, risk, eps_max_scale, best_eps, iter, weights)
    if(iter ==1){
      weights = new_beta_lst[[5]]
    }
    if(iter%%5==1){
      last_risk = risk
      EIC_norm = new_beta_lst[[1]]
      norms = c(norms, EIC_norm)
      risk = new_beta_lst[[2]]
      best_eps = new_beta_lst[[4]]

      if(iter %% 50 ==0){
        if(verbose){
          print(c(risk, EIC_norm))
        }
      }
      risks = c(risks, risk)
      new_beta = new_beta_lst[[3]]


      if(EIC_norm < thresh){
        if(verbose){
          print(data.frame(cbind(risks, norms)))
          print("Converged!")
          print(iter)
        }
        return(new_beta)
      }
      if(abs(last_risk - risk)/last_risk < .0001 & (iter > 50)){

        if(verbose){
          print(data.frame(cbind(risks, norms)))
          print("Risk change too slow. Converged.")
          print(iter)
        }
        return(new_beta)
      }
    }
    iter = iter + 1
  }
  if(verbose){
    print(data.frame(cbind(risks)))
  }
  return(new_beta)
}
