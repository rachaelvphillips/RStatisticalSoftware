hal_basis <- function(max_degree = 2, bins = 350, smoothness_orders = 0, include_zero_order = F, include_lower_order = F, ...) {
  max_degree <- max_degree
  bins <- bins
  smoothness_orders <- smoothness_orders
  include_zero_order <- include_zero_order
  include_lower_order <- include_lower_order
  generator <- function(X) {
    X <-as.matrix(X)
    smoothness_orders <- rep(smoothness_orders[[1]], ncol(X))
    basis_list <- hal9001fast::enumerate_basis(X, max_degree, bins = rep(bins[1], ncol(X)), order_map = smoothness_orders, include_zero_order = include_zero_order, include_lower_order = include_lower_order )
    gen <- function(X) {
      X <-as.matrix(X)
      x_basis <- hal9001fast::make_design_matrix(X, basis_list)
      return(x_basis)
    }
    return(gen)
  }
  return(generator)
}
class(hal_basis) <- "basis_config"

fourier_basis <- function(nbasis = 50, max_degree = 2, ...) {
  nbasis <- nbasis
  max_degree <- max_degree
  generator <- function(X) {
    basis_list <- list()
    var_names <- c()
    print(colnames(X))
    for(name in colnames(X)) {
      sd_val <- sd(X[,name])
      if(sd_val < 1e-9) {
        next
      }
      var_names <- c(var_names, name)
      rangeval <- c(min(X[,name]) - sd_val/2, max(X[,name]) +sd_val/2)
      basis <- fda::create.fourier.basis(rangeval = rangeval, nbasis = nbasis)
      basis_list[[name]] <- basis

    }
    num_orig <- length(var_names)

    gen <- function(X) {
      k <- NULL
      data <- do.call(cbind, lapply(var_names, function(name) {
        out <- as.data.table(fda::eval.basis(X[,name], basis_list[[name]]))
        k <<- ncol(out)
        return(out)
      }))
      data <- matrix(as.vector(as.matrix(data)), ncol = num_orig)

      form <- formula(paste0("~.^", max_degree))
      data <- model.frame(form, data = as.data.frame(data))
      data <- model.matrix(form, data = data)

      data <- data[,-1]
      data <- do.call(cbind, unlist(apply(data, 2, function(v) {list(matrix(v, ncol = k))}), recursive = F))
      data <- as.matrix(data)
      return(data)
    }

    return(gen)
  }
  return(generator)

}

class(fourier_basis) <- "basis_config"

