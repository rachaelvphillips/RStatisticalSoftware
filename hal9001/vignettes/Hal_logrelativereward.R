#
#
# #devtools::install_github('Larsvanderlaan/RStatisticalSoftware/hal9001', build_vignettes = FALSE)
# library(simcausal)
# library(hal9001)
# library(glmnet)
#
#
# D <- DAG.empty()
# D <- D +
#   #node("W1", distr = "rnorm") +
#   node("W2", distr = "rnorm", mean =0,  var = 2) +
#   node("pi", distr = "rconst",   const =  0.5) +
#   node("Abin", distr = "rbinom", size = 1, prob = pi) +
#   node("A", distr = "rconst", const = 2*Abin - 1 ) +
#   node("R", distr = "rbinom", size = 1,
#        prob = plogis(0.3*A + 0.4*W2  )) +
#   node("truth1", distr = "rconst",   const = plogis(0.3*1 + 0.4*W2  )) +
#   node("truth0", distr = "rconst",  const =  plogis(0.3*-1 + 0.4*W2  )) +
#   node("trueRatio", distr = "rconst", const = (truth1/truth0)) +
#   node("outcome", distr = "rconst", const = (A+1)/2) +
#   node("weights", distr = "rconst", const = R/(A*pi + (1-A)/2))
# setD <- set.DAG(D)
# dat <- sim(setD, n = 5000)
# # only grab ID, W's, A, T.tilde, Delta
# Wname <- grep("W", colnames(dat), value = TRUE)
# df <- dat[, c( Wname, "A", "R")]
# true = dat$trueRatio
# mean(dat$Abin)
# # The simulator will generate death at time 0.
# # our package only allow positive integer time, so I add one to all times
# head(df)
# head(true)
# sum(df$A==1)
# sum(df$R==1)
# sum(df$R==1 & df$A==1)
# dat
# pi=dat$pi
# compute_LRR <- function(data, pi,  max_degree= 2){
#   R  = data[,"R"]
#   A = data[,"A"]
#
#   X = data[,-which(colnames(data) %in% c("R", "A"))]
#
#   X = as.matrix(X)
#   A=A
#   pseudo_outcome = 1*(A+1)/2
#   weights = R/(A*pi + (1-A)/2)
#
#   print(dim(X))
#   fit = fit_hal(X=X, Y = pseudo_outcome, weights = weights, max_degree= max_degree,  family = "binomial", smoothness_orders = 0, lower.limits =-Inf,  num_bins = 1000,   return_x_basis = T)
#  # new_fit = glmnet(y=pseudo_outcome,x=fit$x_basis, lambda = fit$lambda_star, weights = weights )
#   predictions = predict(fit$glmnet, newx = fit$x_basis, response = "link", s = fit$lambda_star)
#   print(length(predictions))
#   ratio = exp(predictions)
#   print(length(ratio))
#   return(list(fit, ratio = ratio, avg = mean(ratio)))
# }
#
#
# compute_LRR_glm <- function(data, pi,  max_degree= 1){
#   R  = data[,"R"]
#   A = data[,"A"]
#
#   X = data[,-which(colnames(data) %in% c("R", "A")),drop=F]
#
#   X = as.matrix(cbind(X,X))
#
#   pseudo_outcome = (A+1)/2
#   w = R/(A*pi + (1-A)/2)
#
#
#   fit = cv.glmnet(x=X, y=pseudo_outcome, weights = w, family = "binomial", penalty = c(0,1))
#
#   predictions = predict(fit, newx = X, type = "link")
#   ratio = exp(predictions)
#
#   return(list(fit, ratio = ratio, avg = mean(ratio)))
# }
#
# compute_LRR_naive  <- function(data, pi, include_A = F, max_degree = 1){
#   Y = data[,"R"]
#   X = as.matrix(data[,-which(colnames(data)=="R")])
#   print(colnames(X))
#   fit = fit_hal(X=X, Y = Y, max_degree= max_degree, family = "binomial")
#   X_A0 = X
#   X_A0[,"A"] = -1
#   X_A1 = X
#   X_A1[,"A"] = 1
#   R0 = predict(fit, new_data = X_A0)
#   R1 = predict(fit, new_data = X_A1)
#   #basis = (fit$basis_list[which(fit$coefs!=0)])
#   #print(lapply(basis, function(b){b$cols}))
#
#   return(list(ratio = R1/R0,avg =  mean(R1/R0)))
# }
#
# compute_LRR_naive_glm  <- function(data, pi, include_A = F, max_degree = 1){
#
#   fit = glm("R~ W2 + A + W2*A",data = data, family = binomial())
#
#   data0 = data
#   data0[,"A"] = -1
#   data1 = data
#   data1[,"A"] = 1
#   R0 = predict(fit, data0, type = "response")
#   R1 = predict(fit, data1, type = "response")
#   #basis = (fit$basis_list[which(fit$coefs!=0)])
#   #print(lapply(basis, function(b){b$cols}))
#   return(list(ratio = R1/R0,avg =  mean(R1/R0)))
# }
#
# pred1 = compute_LRR_glm(df, pi)
# pred2 = compute_LRR(df, pi)
#
#
# f = pred2[[1]]
#
# pred3 =  compute_LRR_naive_glm(df,pi)
#
#
#
# mean((log(pred2$ratio) -log(true))^2)
# mean((log(pred3$ratio) - log(true))^2)
#
# mean(((pred2$ratio) -(true))^2)
# mean(((pred3$ratio) - (true))^2)
#
#
# plot(df$W2, )
#
#
# plot(df$W2,log(pred2$ratio))
# plot(df$W2,(pred2$ratio))
# plot(df$W2,(pred3$ratio))
# plot(df$W2,true)
# plot(df$W2,log(true))
#
# f = pred2[[1]]
# length(which(f$coefs!=0))
#
#
# median(pred2$ratio)
# median(pred3$ratio)
# median(true)
# pred2$avg
# pred3$avg
# mean(true)
#
# quantile(df$W2)
