fourier_basis <- function(nbasis = 50, max_degree = 2, unpenalized = NULL, ...) {
  nbasis <- nbasis
  max_degree <- max_degree
  generator <- function(X) {
    X <- as.matrix(X)
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
      colnames(data) <- var_names
      if(max_degree == 1) {
        form <- formula(paste0("~."))

      } else {
        form <- formula(paste0("~.^", max_degree))

      }
      data <- model.frame(form, data = as.data.frame(data))
      data <- model.matrix(form, data = data)

      data <- as.matrix(data[,-1, drop = F])
      var_names <- colnames(data)
      data <- do.call(cbind, unlist(apply(data, 2, function(v) {list(matrix(v, ncol = k))}), recursive = F))
      data <- as.matrix(data)
      new_names <- sapply(var_names, function(a) paste0(a, "_", 1:k))
      colnames(data) <- new_names
      return(data)
    }

    return(gen)
  }
  return(generator)

}
