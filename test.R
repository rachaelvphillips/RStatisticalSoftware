library(hal9001)
library(SGL)
n=1000
X = cbind(replicate(10, runif(n)))
Y = rowSums(X[,sample(1:ncol(X),4)])
fit = fit_hal(X,Y, max_degree = 1, return_x_basis = T)




x_basis = fit$x_basis
data = list(x = as.matrix(x_basis), y = as.vector(Y))
basis_list = fit$basis_list
groups = as.vector(sapply(basis_list, function(b){b$cols}))

sg = cvSGL(data=data)

