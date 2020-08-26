
Lrnr_HALgrad <- R6Class(
  classname = "Lrnr_HALgrad",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(relRiskStop = 0.9, max_iter = 200, offset = NULL,...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c(
      "continuous"
    ),

    .train = function(task) {
      verbose <- getOption("sl3.verbose")
      params <- self$params
      MSE = private$MSE
      X = as.matrix(task$X)
      y = as.vector(task$Y)
      if(ncol(X) > 1) {
        beta <- coefficients(nnls::nnls(X, y))
        X = X %*% beta
      }

      x=X
      offset = params$offset

      relRiskStop = params$relRiskStop
      max_iter = params$max_iter

      #By default the algorithm begins with offsets f(x) = x and f(x) = bx where b is chosen to get the maximum value correctly.
      offset = c(offset, Vectorize(function(t){return(t)}), Vectorize(function(t){(quantile(y, 0.95))/quantile(x, 0.95)*(t) }) )


      wrapper = function(offset){
        private$GradientDescentUpdate(as.vector(X),as.vector(y),relRiskStop, max_iter,offset)
        }
      funcs = lapply(offset, wrapper)
      getMSE = function(func){
        return(MSE(func(X),y))
      }
      MSE_of_funcs = lapply(funcs, getMSE)
      best_func = funcs[[which.min(MSE_of_funcs)]]

      fit_object = list()
      fit_object$func = best_func


      return(fit_object)
    },

    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      X=task$X
      X = rowMeans(X)
      pred = self$fit_object$func
      predictions <- as.vector(sapply(X, pred))#rowMeans(apply(X, MARGIN =2, pred))
      return(predictions)
    },

makeFunction = function(f,grid){

  f =unlist(f)
  grid = unlist(grid)

  func = function(x){
    if(x < grid[1]){
      return(min(x,f[1]))
    }
    else if(x>grid[length(grid)]){
      return(max(x,f[length(grid)]))
    }
    return(f[which.min(abs(x-grid))])
  }
  func = Vectorize(func)
  return(func)
},
MSE = function(pred,y) {
  mean((pred - y) ^ 2)
},

GradientDescentUpdate = function(x,y, relRiskStop, iter, offset) {
  x = unlist(x)
  y=unlist(y)
  ngrid=100
  n = length(x)
  tau0 = min(x)
  tau1 = max(x)

  grid = sort(unique(round(
    unlist(seq(tau0, tau1, length.out = ngrid)), digits = 6)))
  ngrid = length(grid)
  #Store function as list
  f_offset = offset(grid)
  MSE = function(x){return(private$MSE(x,y))}


  evalF = function(v,f){
    return(f[which.min(abs(v-grid))])
  }
  evalF = Vectorize(evalF, vectorize.args = "v")

  updateOnce = function(f_offset, eps = 0.5 / tau1) {
    f_offset=as.vector(f_offset)
    map_to_obs = evalF(x, f=f_offset)  #observed vales
    #Construct clever matrix

    help = function(j_u){
      ind = as.numeric((x<grid[j_u]))
      result = unlist(abs(map_to_obs*ind +(1-ind)*f_offset[j_u]))
    }
    t = proc.time()
    P1= do.call(rbind,lapply(1:ngrid,help))
    P2=outer(unlist(f_offset), unlist(map_to_obs))/abs(f_offset[length(f_offset)])
    clevMat=P1-P2



    Residuals = y - map_to_obs

    perturb = clevMat %*% Residuals

    new_fs = list()
    risks = c()
    epLst = c(30*eps, 10*eps, 4*eps, -eps, eps*2, 0.5, seq(0, eps, length.out = 10))
    for (ep in epLst) {
      cur = f_offset + ep * perturb / n

      new_fs = c(new_fs, list(cur))

      newRisk = MSE(evalF(x, f = cur))

      risks = c(risks, newRisk)
    }


    keep = which.min(risks)

    return(list(epLst[keep],new_fs[[keep]]))
  }


  cur_MSE = Inf
  eps = 100/max(abs(x))
  for (i in 1:iter) {
    lst <- updateOnce(f_offset, eps)
    f_offset = lst[[2]]
    eps = lst[[1]]
    if (i > 0 & i %% 5 == 0) {

      mse = MSE(evalF(x, f_offset))
      if (mse < cur_MSE) {
        if (mse / cur_MSE > relRiskStop) {
          break
        }
        cur_MSE = mse
      }
      else{
        break
      }
    }

  }

  result_func = private$makeFunction(f_offset,grid)
  return(result_func)
}
  )
)
