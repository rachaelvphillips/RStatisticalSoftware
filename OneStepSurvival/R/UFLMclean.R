




#####
# @param T_tilde Vector of failure/censoring times
# @param Delta Indicator of whether event was failure
# @param A Vector of binary treatment assignments
# @param W Matrix of baseline covariates
# @param t_max Maximum time point to consider
# @param n Size of sample
# @param dNt Matrix (empirical hazard) where dNt[n, t] = 1(Delta_n ==1 & t = T_tilde)
#meaning that it encodes whether individual n had a failure at time t.
# @param aliveMat Matrix like above except 1 at time t if the individual is still alive (). Used to compute Risk quickly.
# @param targetPoints Vector of times which to target
# @param estimator Estimator object that stores current hazard/survival estimates
# @param originalEstimator Estimator that contains the intitial superlearner fits
# @param delta initial step size for onestep algorithm
# @param weights Weights for the norm of the EIC we wish to solve
# @param targetFunc Custom parameter of survival curve to target
# @param gradient Gradient of targetFunc to be used in delta method
#
#
#' @importFrom R6 R6Class
#' @importFrom ggplot2 ggplot
#' @importFrom glmnet glmnet
#' @importFrom stats coef
#' @importFrom stats plogis
#' @importFrom stats qlogis
#' @importFrom resample colVars
#####
#' @export
OneStep2 <- R6::R6Class("OneStep2",
                        private = list(
                          varsInit = list(),
                          debug = list(),
                          T_tilde = NULL,
                          Delta= NULL,
                          A= NULL,
                          W= NULL,
                          t_max= NULL,
                          n = NULL,
                          dNt = NULL,
                          aliveMat = NULL,
                          EIC_Indicator = NULL,
                          targetPoints= NULL,
                          estimator = NULL,
                          originalEstimator = NULL,
                          delta = 0.00005,
                          weights = c(),
                          weightsDiff = NULL,
                          weightsCustom=NULL,
                          EICWeights = NULL,
                          targetFunc = NULL,
                          gradient = NULL,
                          #Stores history of norms, risks, and psi
                          EIC_normLst = NULL,
                          EIC_meanLst = NULL,
                          boundBelow = 0.0005,
                          boundAbove = 0.99,
                          cleverMatLongLst = NULL,
                          EICMatlst = NULL,
                          EICMatObsLst = NULL,
                          risk = Inf,
                          historyNorm1 = c(),
                          historyRisk1 = c(),
                          historyNorm0 = c(),
                          historyRisk0 = c(),
                          historyRisk = c(),
                          historyNorm = c(),
                          historyPsi1 = c(),
                          historyPsi0 = c(),
                          historys= NULL,
                          converged = c(F,F),
                          lastRisk=Inf,
                          inner_prod = function(u,v, weights){
                            sum(weights*u*v)
                          },
                          norm = function(u, weights){
                            sqrt(sum(u^2*weights))
                          },
                          #Initializes dNt matrix
                          create_dNt  = function() {
                            max_t =  (private$t_max)
                            dNt <- matrix(0, nrow = length(private$A), ncol = private$t_max)

                            indices = 1:length(private$A)
                            change_indices = indices[(private$Delta ==1 & private$T_tilde <= max_t)]

                            changeStuff = function(n){
                              dNt[n, private$T_tilde[n]] <<- 1
                              return(NULL)
                            }
                            lapply(change_indices, changeStuff)

                            private$dNt = as.matrix(dNt)
                          },
                          #initializes aliveMat
                          create_EIC_indicator = function(){
                            Ttilde = private$T_tilde
                            helper = function(i,j){

                              ttilde = Ttilde[i]
                              if(ttilde <j){
                                return(0)
                              }
                              else{
                                return(1)
                              }}
                            helper = Vectorize(helper)
                            #tmp = matrix(1, nrow = nrow(cleverVecs), ncol = ncol(cleverVecs))
                            ind = outer(1:private$n, 1:private$t_max , FUN=helper)
                            private$EIC_Indicator = as.matrix(ind)
                          },
                          create_alive = function(){
                            x = private$dNt
                            max_t =  (private$t_max)
                            Ttilde = private$T_tilde
                            Delta = private$Delta

                            y = matrix(0, nrow = nrow(x), ncol = max_t)
                            setOnes = function(n){
                              if(Ttilde[n]> max_t){
                                y[n,1:(max_t)] <<- 1
                                return()
                              }
                              ttilde = Ttilde[n]
                              if(ttilde ==1 & Delta[n] ==1){
                                return()
                              }
                              else if(ttilde ==1 & Delta[n] ==0){
                                y[n,1] <<- 1
                                return()
                              }
                              else{
                                if(Delta[n] ==1){
                                  y[n,1:(ttilde-1)] <<- 1
                                  return()
                                }
                                else{
                                  #Assuming that censoring at time ttilde means we only know they survived up until ttilde -1
                                  y[n,1:(ttilde)] <<- 1
                                  return()
                                }
                              }

                            }

                            lapply(1: nrow(x), setOnes)
                            private$aliveMat = y
                          }

                        ),
                        public = list(
                          #T_tilde: vector of event times
                          #Delta: vector of indicators that =1 if failure event occured
                          #A: vector of treatment assignments
                          #W: matrix of covariates
                          #t_max: maximum time of the survival curve to consider
                          #targetPoints: time points of survival curve to target
                          initialize = function(T_tilde, Delta, A, W, t_max, targetPoints){
                            private$T_tilde = T_tilde
                            private$Delta = Delta
                            private$A = A
                            private$W = W
                            private$t_max = t_max
                            private$targetPoints = targetPoints
                            private$n = length(A)
                            private$create_dNt()
                            private$create_alive()
                            private$create_EIC_indicator()
                            private$cleverMatLongLst = list(matrix(0, nrow = private$n*t_max, ncol = length(targetPoints)),matrix(0, nrow = private$n*t_max, ncol = length(targetPoints)))
                          },
                          setDelta = function(delta){
                            private$delta = delta
                          },
                          getDebug = function(){
                            return(private$debug)
                          },
                          get_dNt  = function(long=F) {
                            if(long==F){
                              return(private$dNt)
                            }
                            else{
                              return(as.vector(private$dNt))
                            }
                          },
                          get_alive  = function(long=F) {
                            if(long==F){
                              return(private$aliveMat)
                            }
                            else{
                              return(as.vector(private$aliveMat))
                            }
                          },
                          bound = function(mat){
                            mat[mat < private$boundBelow] = private$boundBelow
                            mat[mat > private$boundAbove] = private$boundAbove
                            return(mat)
                          },
                          #Allows user to input their own fits of the censoring hazard, failure hazard, and conditional treatment probabilities.
                          #hazards must be a matrix with column indices being time and row indicies being the individual.
                          #g1W is the vector of P(A=1|W) for each individual.

                          addInitialEstimates = function(haz_failure_1, haz_failure_0, surv_censor_1, surv_censor_0, g1W){
                            haz_failure_1 = self$bound(haz_failure_1)
                            haz_failure_0 = self$bound(haz_failure_0)
                            surv_censor_1 = self$bound(surv_censor_1)
                            surv_censor_0 = self$bound(surv_censor_0)
                            g1W = self$bound(g1W)
                            print(data.frame(g1W))
                            EstFits  = list("hazard_failure_1" = haz_failure_1,
                                            "hazard_failure_0" = haz_failure_0,
                                            "surv_censor_1" = surv_censor_1,
                                            "surv_censor_0" = surv_censor_0,
                                            "g1W" = unlist(g1W))



                            private$estimator = Estimator$new(EstFits)

                            private$originalEstimator = Estimator$new(EstFits)$clone(deep=T)
                            self$getReady()
                          },
                          #Performs initial SuperLearner fit of the hazards and probabilities using a user-supplied library.
                          doInitialFit = function(sl_failure = c("SL.glm"),
                                                  sl_censoring = c("SL.glm"),
                                                  sl_treatment = c("SL.glm"),gtol = 1e-3){
                            source("sl_fit_survival.R")
                            sl_fits = initial_sl_fit(private$T_tilde,
                                                     private$Delta,
                                                     private$A,
                                                     private$W,
                                                     private$t_max,
                                                     sl_failure = sl_failure,
                                                     sl_censoring = sl_censoring,
                                                     sl_treatment =sl_treatment,
                                                     gtol = 1e-3)

                            private$estimator = Estimator$new(sl_fits)
                            private$originalEstimator = Estimator$new(sl_fits)$clone(deep=T)
                            self$getReady()


                          },
                          #targetFunc is a bivariate function of S1(t), S0(t) and
                          #gradient is the gradient of targetFunc (as a function of S1(t) and S0(t))
                          #Compute weights using initial fit.
                          getReady = function(){
                            self$recalculate()
                            vars1 = self$getColVars(1,"obs")
                            vars0=self$getColVars(0,"obs")
                            private$varsInit = list(vars0, vars1)
                            inv_vars1 = 1/vars1
                            inv_vars0 = 1/vars0


                            inv_vars0[inv_vars0 < 0.0005] = 0
                            print(paste0("Not targetting A =0: ",private$targetPoints[inv_vars0 < 0.005]))


                            inv_vars0[inv_vars1 < 0.0005] = 0
                            print(paste0("Not targetting A =1:", private$targetPoints[inv_vars1 < 0.005]))
                            private$EICWeights[[2]] = inv_vars1/sum(inv_vars1)
                            private$EICWeights[[1]] = inv_vars0/sum(inv_vars0)
                          },
                          setTargetMapping = function(targetFunc, gradient){

                            private$targetFunc = targetFunc
                            private$gradient = gradient
                          },
                          computeGradient = function(survValues){
                            private$gradient(survValues[1], survValues[2])
                          },
                          computeTarget = function(survValues){
                            #print(paste0("works", private$targetFunc))
                            private$targetFunc(survValues[1], survValues[2])
                          },
                          #Returns matrix where column t and row n is the influence curve
                          #for t-th component of parameter and person n evaluated at the given A
                          updateHistory = function(A_cf){

                            if(A_cf==1){
                              private$historyNorm1 = c(private$historyNorm1, self$getEICNorm(1))
                              private$historyRisk1 = c(private$historyRisk1, self$getRisk())
                              private$historyPsi1 = rbind(private$historyPsi1, self$computeEstimates(1))
                            }
                            else if(A_cf==0){
                              private$historyNorm0 = c(private$historyNorm0, self$getEICNorm(0))
                              private$historyRisk0 = c(private$historyRisk0, self$getRisk())
                              private$historyPsi0 = rbind(private$historyPsi0, self$computeEstimates(0))
                            }
                          },
                          setEICMeanNormFast = function(A_cf){
                            X0 = self$getColMeans(0, "obs")
                            X1 = self$getColMeans(1, "obs")
                            private$EIC_normLst = list(private$norm(X0,self$getEICWeights(0)), private$norm(X1, self$getEICWeights(1)))
                            private$EIC_meanLst = list(X0, X1)
                          },
                          #Sets the long clever matrix needed for computing EIC and update step
                          #Then computes the EICMean and Norm quickly. Then risk.
                          recalculate = function(){
                            self$setCleverMatLong()
                            #self$setEICMat()
                            self$setEICMeanNormFast()
                            self$computeRisk()

                          },
                          getEICMat = function(A_cf, full = F, A = "obs"){
                            return(self$computeEICMat(A_cf = A_cf, full = full, A = A)
                            if(full == T){
                              return(self$computeEICMat(A_cf = A_cf, full = full, A = A))
                            }
                            else if(A == "obs" & A_cf %in% c(0,1)){
                              return(private$EICMatObsLst[[A_cf+1]])
                            }
                            else if(A_cf %in% c(0,1) & A %in% c(0,1) &A!=A_cf){
                              return(matrix(0,nrow=private$n, ncol = length(private$targetPoints)))
                            }
                            else if(A_cf == "custom"){

                              S1 = self$computeEstimates(1)
                              S0 =  self$computeEstimates(0)
                              preds = matrix(c(S0,S1), ncol =2)
                              grad = apply(X=preds,MARGIN=1, FUN=self$computeGradient)
                              C0 = t(t(self$getEICMat(A_cf=0, F, A=A))*grad[1,])
                              C1 = t(t(self$getEICMat(A_cf=1,F, A=A))*grad[2,])
                              return(C1+C0)
                            }
                            else{
                              return(private$EICMatlst[[A+1]])
                            }
                          },

                          setEICMat = function(){
                            X0 = self$computeEICMat(0, full =F, A=0)
                            X1 = self$computeEICMat(1, full =F, A=1)
                            private$EICMatlst = list(X0, X1)
                            ind= private$A==1
                            observed1 = X1*ind
                            observed0 = X0*(1-ind)
                            private$EICMatObsLst  = list(observed0, observed1)

                          },
                          computeEICMat = function(A_cf, full = F,A="obs"){

                            return(do.call(cbind,(lapply(private$targetPoints,self$computeEIC, A_cf = A_cf, full = full, A=A))))
                          },

                          getColMeans = function(A_cf, A="obs"){

                            Qmat = as.matrix(self$getHazard(A))
                            dNt = as.matrix(self$get_dNt())
                            residual_mat = (dNt - Qmat)
                            ind = private$EIC_Indicator
                            cleverVecs =self$getCleverMatLong(A_cf=A_cf,A=A)

                            resVec = colSums(cleverVecs*as.vector((residual_mat)*(ind)/private$n))
                            return(resVec)



                          },
                          #Recomputes EIC from clever covariates.
                          #Seems to be faster than using the stored clever covariates.
                          #Only needs to be computed once at beggining for weights and at end for inference.
                          getColVars = function(A_cf, A = "obs"){
                            return(resample::colVars(self$computeEICMat(A_cf, full = F, A)))
                          },
                          computeEIC = function(t_tgt, A_cf, full = F, A="obs", mat =F ){

                            #cleverVecs1=matrix(self$getCleverMatLong(A_cf=A_cf,A=A)[,which(private$targetPoints ==t_tgt)], ncol = private$t_max, nrow = private$n)
                            cleverVecs = self$getCleverMatSpeed(A_cf=A_cf,t_tgt =t_tgt, A=A)

                            Qmat = as.matrix(self$getHazard(A))
                            dNt = as.matrix(self$get_dNt())
                            residual_mat = dNt - Qmat
                            ind = private$EIC_Indicator
                            #return(list(as.vector(cleverVecs1)*as.vector(residual_mat*ind), data.frame(cleverVecs*residual_mat*ind)))
                            part1 = rowSums(cleverVecs*residual_mat*ind)

                            # if(mat==T){
                            #   cleverVecs =as.vector(self$getCleverMatLong(A_cf=A_cf,A=A))
                            #   resVec = colSums(matrix(colMeans(matrix(cleverVecs*as.vector((residual_mat)*(ind)), nrow = private$n)), nrow = private$t_max))
                            #   return(resVec)
                            #  hi = apply(X=resVec, MARGIN=2, FUN = function(v){rowSums(matrix(v, ncol = private$t_max, nrow = private$n))})
                            #  return(hi)
                            #  }
                            #return((cleverVecs*residual_mat*ind))

                            if(!full){
                              return(part1)
                            }
                            Ttilde = private$T_tilde
                            t_max = private$t_max
                            est = private$estimator
                            if(A_cf %in% c(0,1)){
                              part2 = est$evalSurvival_failure(t_tgt, A_cf)- mean(unlist(est$evalSurvival_failure(t_tgt, A_cf)))
                            }
                            else if(A_cf == "diff"){
                              part2 = est$evalSurvival_failure(t_tgt, 1)- est$evalSurvival_failure(t_tgt, 0)-  (mean(unlist(est$evalSurvival_failure(t_tgt, 1))) - mean(unlist(est$evalSurvival_failure(t_tgt, 0))))
                            }
                            else if(A_cf == "custom"){

                              parts = as.matrix(cbind(est$evalSurvival_failure(t_tgt, 0)- mean(unlist(est$evalSurvival_failure(t_tgt, 0))), est$evalSurvival_failure(t_tgt, 1)- mean(unlist(est$evalSurvival_failure(t_tgt, 1)))))
                              S1 = self$getSurvival(1)[t_tgt]
                              S0 =  self$getSurvival(0)[t_tgt]

                              grad = self$computeGradient(c(S0, S1))
                              part2 = parts %*% grad
                            }
                            #Add support for diff

                            return(part1 + part2)

                          },
                          getCleverMatLong = function(A_cf, A){
                            if(A_cf %in% c(0,1)& A%in% c(0,1) &A!=A_cf){
                              X1 = self$getCleverMatLong(1,1)
                              return(matrix(0, nrow=nrow(X1), ncol=ncol(X1)))
                            }
                            if(A == "obs" & A_cf %in% c(0,1)){


                              ind = as.numeric(private$A==A_cf)


                              X = self$getCleverMatLong(A_cf,A_cf)
                              z=X*ind

                              return(as.matrix(z))

                            }
                            if(A_cf == "custom"){
                              S1 = self$computeEstimates(1)
                              S0 =  self$computeEstimates(0)
                              preds = matrix(c(S0,S1), ncol =2)
                              grad = apply(X=preds,MARGIN=1, FUN=self$computeGradient)
                              C0 = t(t(self$getCleverMatLong(A_cf=0, A=A))*grad[1,])
                              C1 = t(t(self$getCleverMatLong(A_cf=1, A=A))*grad[2,])
                              return(C0+C1)
                            }
                            if(A_cf == "diff"){

                              return(self$getCleverMatLong(1,A) - self$getCleverMatLong(0,A))
                            }

                            return( private$cleverMatLongLst[[A_cf+1]])


                          },
                          setCleverMatLong = function(){

                            private$cleverMatLongLst[[1]] =self$computeCleverMatLong(A_cf = 0, A=0)
                            private$cleverMatLongLst[[2]] =self$computeCleverMatLong(A_cf = 1, A=1)



                          },
                          computeCleverMatWide = function(A_cf, A="obs"){

                            return(as.matrix(do.call(cbind,lapply(1:private$t_max, self$getCleverMat, A_cf = A_cf, A = A))))
                          },

                          computeCleverMatLong = function(A_cf, A = "obs"){

                            result = as.matrix(do.call(cbind, lapply(private$targetPoints, FUN = function(t){as.vector(self$getCleverMatSpeed(A_cf=A_cf,t,A=A))})))

                            if(max(result)>1000){
                              print("Large value in clever cov")
                            }



                            return(as.matrix(result))
                          },

                          #Returns the matrix of clever covariates (for all target times) for time t
                          #SLOW, do not use.
                          getCleverMat = function(t, A_cf, A="obs",addIndic = F){
                            getForT= function(t) {
                              return(do.call(cbind,lapply(private$targetPoints, self$getCleverCovariate, A_cf = A_cf, t = t, A=A,addIndic = addIndic)))
                            }
                            return(getForT(t))
                          },


                          #Returns matrix of clever covariates at a specific target time
                          #Different than above. But using this is faster when constructing the long cleverMat matrix.
                          getCleverMatSpeed = function(A_cf, t_tgt,A=A){
                            if(A_cf == "custom"){
                              S1 = self$computeEstimates(1)
                              S0 =  self$computeEstimates(0)
                              preds = matrix(c(S0,S1), ncol =2)
                              grad = apply(X=preds,MARGIN=1, FUN=self$computeGradient)
                              C0 = self$getCleverMatSpeed(A_cf=0,t_tgt=t_tgt, A=A)*grad[1,]
                              C1 = self$getCleverMatSpeed(A_cf=1,t_tgt=t_tgt, A=A)*grad[2,]
                              return(C0+C1)
                            }
                            if(A_cf == "diff"){

                              return(self$getCleverMatSpeed(1,, t_tgt, A) - self$getCleverMatSpeed(0, t_Tgt, A))
                            }
                            est = private$estimator


                            if(A=="obs"){
                              A_ind = as.numeric(private$A==A_cf)
                            }
                            else {
                              A_ind = as.numeric(A==A_cf)
                            }

                            g_A = as.numeric(est$eval_g(A_cf))


                            G_censor =force(est$evalSurvival_censorMatLeft(A_cf))

                            S_t = force(est$evalSurvival_failureMat(A_cf))

                            S_tgt = est$evalSurvival_failure(t_tgt, A_cf)

                            clever_cov = -1 *(A_ind*S_tgt/g_A)/ (G_censor*S_t)
                            #clever_cov1 = -1 *A_ind*(S_tgt/S_t)/ (G_censor*g_A)
                            #print(max(abs(clever_cov1-clever_cov)))
                            if(t_tgt == private$t_max){
                              return(clever_cov)
                            }
                            clever_cov[,(t_tgt+1):ncol(clever_cov)]<-0


                            return(as.matrix(clever_cov))


                          },



                          #Returns the vector of clever covariates for targeting time t_tgt, time t, target parameter A_cf, evaluated at treatment A
                          #addIndic is used to compute the EIC. It takes into account the 1(Ttilde>=t) in the formula for the EIC.
                          #The indicator is absorbed into the clever covariate function for speed.
                          getCleverCovariate = function(t, A_cf, t_tgt, A="obs", addIndic = F){
                            if(t > t_tgt){return(rep(0,private$n))}
                            if(A_cf =="diff"){
                              return(self$getCleverCovariate(t, A_cf=1, t_tgt, A=A, addIndic = addIndic) -self$getCleverCovariate(t, A_cf=0, t_tgt,A=A, addIndic = addIndic) )
                            }
                            if(A_cf == "custom"){
                              S1 = mean(self$getSurvival(1)[,t_tgt])
                              S0 =  mean(self$getSurvival(0)[,t_tgt])
                              grad = self$computeGradient(c(S0, S1))
                              cleverCovariates = as.matrix(cbind(self$getCleverCovariate(t=t,A_cf = 0, t_tgt = t_tgt, A = A, addIndic = addIndic) , self$getCleverCovariate(t=t,A_cf = 1, t_tgt = t_tgt, A = A, addIndic = addIndic) ))

                              return(cleverCovariates %*% grad)
                            }


                            est = private$estimator
                            if(A=="obs"){
                              A_ind = as.numeric(private$A==A_cf)
                            }
                            else {
                              A_ind = as.numeric(A==A_cf)
                            }

                            g_A = as.numeric(est$eval_g(A_cf))


                            G_censor = as.numeric(est$evalSurvival_censorLeft(t, A_cf))
                            #print(data.frame(G_censor))
                            #speedup https://stackoverflow.com/questions/17308551/do-callrbind-list-for-uneven-number-of-column
                            S_t = est$evalSurvival_failure(t, A_cf)
                            S_tgt = est$evalSurvival_failure(t_tgt, A_cf)

                            clever_cov = -1 *(A_ind / (g_A*G_censor)) * (S_tgt/S_t)
                            if(addIndic){
                              clever_cov = clever_cov * as.numeric((private$T_tilde>=t ))
                            }

                            return(clever_cov)
                          },

                          #Performs one iteration of OneStep targeting the parameter S_1 - S_0 directly.
                          updateDiff = function(){
                            delta = private$delta


                            ptm <- proc.time()


                            EIC_mat = self$getEICMat(A_cf=1,A="obs")-self$getEICMat(A_cf=0,A="obs")


                            #mean_var = Rfast2::colmeansvars(as.matrix(EIC_mat), parallel = F)
                            EIC_mean = colMeans(EIC_mat)
                            EIC_var = resample::colVars(EIC_mat)
                            EIC_var[EIC_var<1e-4] = 1e-4

                            inv_EIC_var = 1/EIC_var
                            weights = private$weightsDiff
                            if(is.null(weights)){
                              private$weightsDiff = inv_EIC_var/sum(inv_EIC_var)
                            }


                            #weights = rep(1/length(EIC_mean), length(EIC_mean))


                            curRisk = self$computeRisk()




                            EIC_norm = private$norm(EIC_mean, weights)
                            print(paste("hi",weights,EIC_norm,curRisk))



                            if(EIC_norm < 1/(sqrt(private$n) * log(private$n))){
                              print("Already Converged... returning...")
                              #return()


                            }


                            cleverMats1 = lapply(1:private$t_max, self$getCleverMat, A_cf = "diff",A=1)
                            cleverMats0 = lapply(1:private$t_max, self$getCleverMat, A_cf = "diff",A=0)

                            offset1 = data.frame(stats::qlogis(self$getHazard(A_cf=1, long = F)))
                            offset0 = data.frame(stats::qlogis(self$getHazard(A_cf=0, long = F)))


                            getDirection = function(){
                              dir = EIC_mean /weights
                              norm_dir = private$norm(dir, weights )
                              return(delta*dir/norm_dir)
                            }



                            updateAtT_1 = function(t, dir){
                              clever_t1 = cleverMats1[[t]]
                              offset_t1 = as.vector(offset1[,t])
                              update_t1 = stats::plogis(as.vector(offset_t1 + (clever_t1) %*%  (dir))) #Rfast::mat.mult(as.matrix(clever_t), as.matrix(dir))))
                              return(update_t1)
                            }
                            updateAtT_0 = function(t, dir){
                              clever_t0 = cleverMats0[[t]]
                              offset_t0 = as.vector(offset0[,t])
                              update_t0 = stats::plogis(as.vector(offset_t0 + (clever_t0) %*%  (dir))) #Rfast::mat.mult(as.matrix(clever_t), as.matrix(dir))))
                              return(update_t0)
                            }

                            dir = getDirection()
                            oldHaz1 = self$getHazard(1)
                            oldHaz0 = self$getHazard(0)
                            new_haz_1 = as.matrix(do.call(cbind, lapply(1:private$t_max, updateAtT_1, dir = dir)))
                            new_haz_0 = as.matrix(do.call(cbind, lapply(1:private$t_max, updateAtT_0, dir = dir)))

                            newRisk = self$computeUpdatedRisk(new_haz_0, new_haz_1)

                            lastRisk =curRisk
                            if(newRisk > lastRisk){
                              print(newRisk)
                              print(lastRisk)
                              print(paste0("Warning: Risk increased... Step rejected... Halving stepsize...new delta=", private$delta/2))
                              private$delta = delta/2

                              if(private$delta <1e-8){
                                print(paste0("Convergence failed... Delta too small...delta=", private$delta))
                                #private$converged[A_cf+1]=T
                                return()
                              }
                              #self$updateOnce(A_cf)
                              return("redo")

                            }
                            else{

                              #private$delta = deltanew
                            }

                            private$estimator$update(haz = as.matrix(new_haz_0), 0)

                            private$estimator$update(haz = as.matrix(new_haz_1), 1)

                            newRisk=self$computeRisk()

                            if(newRisk > curRisk){
                              print("Warning: Risk increased!")




                            }
                            else{
                              #print(paste0("risk diff", newRisk- curRisk))
                            }




                          },
                          ##############
                          #Computes the criterion change in Psi squared divided by change in Risk to evaluate sizes of delta steps.
                          getQuotient = function(haz, A_cf){


                            surv = data.frame(t(apply(haz[,1:max(private$targetPoints)], 1, function(x){cumprod(1-x)})))
                            surv = surv[,private$targetPoints]

                            newPsi = colMeans(surv)
                            curPsi = self$computeEstimates(A_cf)

                            diff = abs(newPsi-curPsi)
                            max_diff = max(diff)
                            PsiChange = mean((diff)^2)
                            if(A_cf ==1){
                              Risk = self$computeUpdatedRisk(newHazard1 = haz)
                            }
                            else{
                              Risk = self$computeUpdatedRisk(newHazard0 = haz)
                            }

                            curRisk = self$getRisk()
                            diff = -(Risk - curRisk)
                            return(c(max_diff, Risk - curRisk, Risk, PsiChange/diff))
                          },

                          ####################
                          #Targets the parameter S_0 or S_1 depending on whether A_cf =0 or 1.
                          #Targetting both S_0 and S_1 requires running the below function twice for each value.
                          updateOnce = function(A_cf, search_step = T, byRisk = F,recalc=T, stop = T, debug = F,  maxDelta = 0.005){
                            t=proc.time()
                            delta = private$delta
                            if(private$converged[A_cf+1] == T){
                              print("Already converged... returning...")
                              #return()
                            }
                            if(is.null(private$EICWeights)){
                              self$recalculate()
                            }
                            weights = self$getEICWeights(A_cf)
                            EIC_mean = self$getEICMean(A_cf)
                            EIC_norm = self$getEICNorm(A_cf)
                            curRisk = self$getRisk()

                            if(EIC_norm < 1/(sqrt(private$n) * log(private$n))){
                              print(EIC_norm)
                              print("We Converged!")
                              return() }

                            cleverMats1 = self$getCleverMatLong(A_cf=A_cf,A=A_cf)
                            haz = self$getHazard(A_cf=A_cf, long = F)
                            haz = self$bound(haz)

                            offset = stats::qlogis(haz)





                            getDirection = function(delt){
                              n=private$n
                              len=length(EIC_mean)
                              convergedIndex = (abs(EIC_mean)/ sqrt(private$varsInit[[A_cf+1]]) ) < 1/n
                              #dir = EIC_mean /abs(weights)
                              dir = EIC_mean/ sqrt(private$varsInit[[A_cf+1]])

                              norm_dir = private$norm(dir,1/len)
                              dir[convergedIndex] = 0

                              #dir = EIC_mean/weights
                              #norm_dir = private$norm(dir, weights)
                              if(debug){
                                #print(data.frame(dir/norm_dir))

                              }

                              return(as.vector((dir/norm_dir)*delt))
                            }


                            updateAllT = function(dir){
                              dir = as.vector(dir)
                              cleverAll = cleverMats1
                              offset_All = as.vector(offset)

                              update_All = stats::plogis(offset_All + cleverAll %*%  (dir))
                              return(matrix(update_All, nrow = private$n, ncol = private$t_max))
                            }

                            if(debug){
                              drawDirection = function(n){
                                v=sample(c(-1,1), length(private$targetPoints), replace=T)
                                return(delta*v/private$norm(v,1))
                              }
                              directions = lapply(1:30, drawDirection)
                              hazards = lapply(directions, updateAllT)
                              if(A_cf==1){
                                risks = unlist(lapply(hazards, self$computeUpdatedRisk, newHazard0=NULL))
                              }
                              else{
                                risks=unlist(lapply(hazards, self$computeUpdatedRisk, newHazard1=NULL))
                              }
                              private$debug = list(data.frame(do.call(cbind,directions)), data.frame(risks - self$getRisk()))


                            }
                            #Extremely inefficient???
                            getOptHazard = function(){
                              deltas = c(delta/10, delta/4, delta/2, delta, delta*2, 4*delta, 10*delta)
                              deltas = deltas[deltas<maxDelta]
                              directions = lapply(deltas, getDirection)
                              hazards = lapply(directions, updateAllT)
                              if(byRisk){
                                if(A_cf==1){
                                  risks =unlist(lapply(hazards, self$computeUpdatedRisk, newHazard0 =NULL))
                                  ind = which.min(risks)
                                }
                                else{
                                  risks =unlist(lapply(hazards, self$computeUpdatedRisk, newHazard1 =NULL))
                                  ind = which.min(risks)
                                }
                                return(list(risks[ind],hazards[[ind]]))

                              }

                              quotientStats = as.matrix(do.call(rbind, lapply(hazards, self$getQuotient, A_cf = A_cf)))
                              rownames(quotientStats) = deltas
                              colnames(quotientStats) = c("ChangePsi_max", "risk change", "-Likelihood", "Quotient")
                              if(debug){
                                print(data.frame(quotientStats))
                              }

                              ind = which.max(quotientStats[,4])

                              new_haz = hazards[[ind]]
                              newRisk = quotientStats[ind,3]
                              newPsiChange = quotientStats[ind,1]
                              deltanew = deltas[ind]
                              private$delta = deltanew

                              return(list(new_haz, newRisk, newPsiChange))
                            }

                            if(T){
                              if(search_step & !byRisk){
                                lst= getOptHazard()
                                new_haz = lst[[1]]
                                newRisk = lst[[2]]
                                newPsiChange = lst[[3]]
                              }
                              else if(search_step & byRisk){
                                lst= getOptHazard()
                                new_haz = lst[[2]]
                                newRisk = lst[[1]]

                                newPsiChange = 1
                              }
                              else{
                                newPsiChange=1
                                new_haz = (updateAllT(getDirection(delta)))


                                if(A_cf ==1){
                                  newRisk  =self$computeUpdatedRisk(newHazard1 = new_haz)
                                }
                                else{
                                  newRisk  = self$computeUpdatedRisk(newHazard0 = new_haz)
                                }
                              }
                            }
                            #new_haz = updateAllT(getDirection(delta))
                            #print(paste("LGOGLGLGLO",self$computeUpdatedRisk(newHazard1 = new_haz) ))
                            #newPsiChange = 1
                            new_haz[new_haz<1e-5] = 1e-5
                            new_haz[new_haz>1-1e-5] = 1-1e-5


                            if(newRisk > curRisk){
                              print(paste0("Warning: Risk increased... Step rejected... Halving stepsize...new delta=", private$delta/2))
                              #private$delta = delta/2
                              if(private$delta <1e-16){
                                print(paste0("Convergence failed... Delta too small...delta=", private$delta))

                                return()
                              }
                              if(stop){
                                private$delta = delta/2
                                return("redo")
                              }
                              #return("redo")

                            }
                            t =proc.time()
                            private$estimator$update(haz = as.matrix(new_haz), A_cf)

                            t = proc.time()
                            if(recalc){
                              #If recalculate is T (e.g. updateOnce isn't called through update function)
                              #then will recalculate new state of estimator using new fit and update history
                              #Between consecutive updateOnce calls with same arg A_cf, recalculate() must be called.
                              self$recalculate()
                              self$updateHistory(A_cf)
                              newNorm = self$getEICNorm(A_cf)
                              if(EIC_norm < 1/(sqrt(private$n) * log(private$n))){
                                print("Converged!")
                              }
                            }

                            newNorm = self$getEICNorm(A_cf)

                            if(newPsiChange < 1e-4){
                              #private$delta = delta*2
                              #print(paste0("Convergence too slow... Doubling step size... new delta=", private$delta))
                              print(paste0("Convergence too slow...  either you selected your step size too small or you have converged"))

                            }

                            return()

                          },
                          #Update function that performs a maximum amount of iterations given by iter with initial delta value delta.
                          #A_cf specifies whether to target S_1, S_0, their difference, or the custom parameter.
                          update = function(iter = 10, delta = NULL, A_cf = NULL){
                            if(!is.null(delta)){
                              private$delta = delta
                            }
                            tmp = private$delta

                            helper = function(A_cf,k){
                              if(k==0){
                                print("Too many recursions... Iteration failed...")
                                return()
                              }
                              if(A_cf =="diff"){
                                self$updateDiff()
                              }
                              else if(A_cf == "custom"){
                                self$updateCustom()
                              }
                              else{
                                x = self$updateOnce(A_cf=A_cf)
                              }

                              #if(!is.matrix(x)  & is.character(x) & x=="redo"){
                              #helper(A_cf, k-1)
                              #}
                            }
                            if(is.null(A_cf)){
                              replicate(iter,helper(A_cf = 1, 2))
                              private$delta = tmp
                              replicate(iter,helper(A_cf = 0, 2))

                            }
                            else if(A_cf ==1 | A_cf ==0){
                              replicate(iter,helper(A_cf = A_cf, 2))

                            }
                            else if(A_cf == "custom"){
                              replicate(iter,helper(A_cf = "custom", 2))
                            }



                            private$delta =delta
                            return(self$getHistory())
                          },
                          #Computes risk from hazard function. This is only called in the update functions for speed.
                          #Allows us to compute Risk without needing to construct an estimator object.
                          #We do not want to construct an estimator function if the risk increases, this allows us to do this.
                          computeUpdatedRisk = function(newHazard0=NULL, newHazard1=NULL){
                            if(is.null(newHazard0)){
                              newHazard0 = as.matrix(self$getHazard(0))
                            }
                            else if(is.null(newHazard1)){
                              newHazard1 = as.matrix(self$getHazard(1))
                            }


                            A = private$A
                            ind = A==1
                            haz =newHazard1*(ind) + newHazard0*(1-ind)
                            dNt = as.matrix(private$dNt)
                            alive = as.matrix(private$aliveMat)

                            haz = self$bound(haz)
                            loss = rowSums(as.matrix(-dNt* log(haz) - (alive)* log(1 - haz)))
                            mean(loss)
                          },
                          #Computes risk of current state of estimator using the estimator object.
                          getRisk = function(){
                            return(private$risk)
                          },
                          computeRisk = function(){

                            haz = as.matrix(self$getHazard("obs"))
                            dNt = as.matrix(private$dNt)
                            alive = as.matrix(private$aliveMat)
                            haz = self$bound(haz)
                            loss = rowSums(as.matrix(-dNt* log(haz) - (alive)* log(1 - haz)))
                            risk = mean(loss)
                            private$risk = risk

                            return(risk)
                          },
                          getEICWeights = function(A_cf){
                            #return(1/length(private$targetPoints))
                            return(private$EICWeights[[A_cf+1]])
                          },
                          setEICWeights = function(){

                            doOnce = function(A_cf){
                              EIC_mat = as.matrix(self$getEICMat(A_cf,A="obs"))
                              #mean_var = Rfast2::colmeansvars(as.matrix(EIC_mat), parallel = F)
                              EIC_mean = colMeans(EIC_mat)
                              EIC_var = resample::colVars(EIC_mat)
                              EIC_var[EIC_var<1e-4] = 1e-4
                              inv_EIC_var = 1/EIC_var
                              print(inv_EIC_var)
                              inv_EIC_var[inv_EIC_var < 0.0005] = 0

                              print(paste0("Not targetting ",A_cf, private$targetPoints[inv_EIC_var < 0.005]))
                              weights = inv_EIC_var/sum(inv_EIC_var)
                              #weights = rep(1/length(EIC_mean), length(EIC_mean))
                              return(weights)
                            }

                            private$EICWeights = list(doOnce(0), doOnce(1))
                          },
                          getEICNorm = function(A_cf){
                            return(private$EIC_normLst[[A_cf+1]])
                          },
                          getEICMean = function(A_cf){
                            private$EIC_meanLst[[A_cf+1]]
                          },
                          # setEICNorm = function(A_cf){
                          #   X0 = self$computeEICNorm(0)
                          #   X1 = self$computeEICNorm(1)
                          #   private$EIC_normLst = list(X0[[2]], X1[[2]])
                          #   private$EIC_meanLst = list(X0[[1]], X1[[1]])
                          # },
                          #
                          # #Computes current (weighted) norm of EIC for targetParameter A_cf
                          # computeEICNorm = function(A_cf){
                          #   EIC_mat = as.matrix((self$getEICMat(A_cf=A_cf, A="obs", full=F)))
                          #
                          #   if(is.null(private$EICWeights)){
                          #     self$setEICWeights()
                          #   }
                          #   weights = self$getEICWeights(A_cf)
                          #   EIC_mean = colMeans(EIC_mat)
                          #
                          #   EIC_norm = private$norm(EIC_mean, weights)
                          #
                          #   return(list(EIC_mean,EIC_norm))
                          # },
                          #Computes the current estimate for target Parameter A_cf
                          computeEstimates = function(A_cf){
                            if(A_cf == "diff"){
                              return(self$computeEstimates(1) - self$computeEstimates(0))
                            }
                            if(A_cf == "custom"){

                              survEstimates = as.matrix(cbind(self$computeEstimates(0), self$computeEstimates(1)))
                              return(apply(X = survEstimates, FUN=self$computeTarget, MARGIN = 1))
                            }
                            results= lapply(private$targetPoints, private$estimator$evalSurvival_failure, A_cf = A_cf)
                            df = data.frame(do.call(cbind, results))
                            colnames(df) = paste0(private$targetPoints, "_", A_cf)
                            return(colMeans(df))
                          },
                          #Computes conditional survival curve estimate given A = A_cf
                          getSurvival = function(A_cf){
                            mat = private$estimator$evalSurvival_failureMat()

                            return(mat)
                          },

                          #Computes conditional hazard evaluated at A=1, A=0, or the observed values of A for each individual.
                          getHazard = function(A_cf, long = F){
                            if(A_cf == "obs"){
                              haz1 = private$estimator$evalHazard_failureMat(1)
                              haz0 = private$estimator$evalHazard_failureMat(0)
                              ind = as.numeric(private$A==1)


                              haz = ind*haz1 + (1-ind)*haz0

                            }
                            else{
                              haz = private$estimator$evalHazard_failureMat(A_cf)
                            }
                            if(long){
                              return(as.vector(haz))
                            }
                            return(haz)
                          },
                          #Returns originalEstimator
                          getOriginalEstimator = function(){
                            return(private$originalEstimator)
                          },
                          #Returns current estimator

                          #Resets the estimator back to the initial estimator (e.g. superlearner fit)
                          reset = function(){
                            tmp = private$originalEstimator$clone(deep=T)
                            private$estimator = private$originalEstimator
                            private$originalEstimator = tmp

                          },

                          #Returns history for S_1 and S_0
                          getHistory = function(){
                            dg = data.frame()
                            df = data.frame()
                            if(length(private$historyNorm0)>0){
                              dg = data.frame(cbind((private$historyNorm0), (private$historyRisk0)))

                              colnames(dg) = c("EIC_norms_0",  "Risks_0")
                              rownames(dg) = 1:length(dg[,1])-1
                            }
                            if(length(private$historyNorm1)>0){
                              df = data.frame(cbind((private$historyNorm1), (private$historyRisk1)))

                              colnames(df) = c("EIC_norms_1",  "Risks_1")
                              rownames(df) = 1:length(df[,1])-1
                            }




                            return(list(dg,df, data.frame(private$historyPsi1 ), data.frame(private$historyPsi0)))
                          },
                          plot = function(A_cf, x_lab = "Time after beggining monotherapy (months)", y_lab = "Probability of survival", main = "Counterfactual survival curves for time until death for PD1 vs PDL1 treatments.", fill_lab = "Treatment"){
                            if(A_cf == "diff"){
                              results_CI = self$getCI(A_cf)
                              results_CI$lower = results_CI[, paste0("lower_bound","_", A_cf)]
                              results_CI$upper = results_CI[, paste0("upper_bound","_", A_cf)]

                              g1 = ggplot2::ggplot(data =  results_CI, aes_string("t",  paste0("prediction","_", A_cf))) +
                                geom_point(aes_string(x = "t", y =  paste0("prediction","_", A_cf)), legend=  F,  xlab="X", ylab="Y", colour=alpha('red')) +
                                geom_point(aes_string(x = "t", y =  paste0("lower_bound","_", A_cf)), legend=  F,  xlab="X", ylab="Y", colour=alpha('red', 0.2)) +
                                geom_point(aes_string(x = "t", y =  paste0("upper_bound","_", A_cf)), legend=  F,  xlab="X", ylab="Y", colour=alpha('red', 0.2)) +
                                geom_smooth(aes(colour = "red"), se = F)+
                                geom_ribbon(aes(ymin = lower, ymax = upper, fill = A_cf), alpha= 0.2, color = NA) +
                                scale_fill_manual(values=c( "blue"))  + scale_x_continuous(breaks = private$targetPoints) +
                                geom_hline(yintercept=0, linetype="dashed", color = "red", size=1)

                              return(g1)
                            }
                            else if(A_cf == "custom"){
                              results_CI = self$getCI(A_cf)
                              results_CI$lower = results_CI[, paste0("lower_bound","_", A_cf)]
                              results_CI$upper = results_CI[, paste0("upper_bound","_", A_cf)]


                              g1 = ggplot2::ggplot(data =  results_CI, aes_string("t",  paste0("prediction","_", A_cf))) +
                                geom_point(aes_string(x = "t", y =  paste0("prediction","_", A_cf)), legend=  F,  xlab="X", ylab="Y", colour=alpha('red')) +
                                geom_point(aes_string(x = "t", y =  paste0("lower_bound","_", A_cf)), legend=  F,  xlab="X", ylab="Y", colour=alpha('red', 0.2)) +
                                geom_point(aes_string(x = "t", y =  paste0("upper_bound","_", A_cf)), legend=  F,  xlab="X", ylab="Y", colour=alpha('red', 0.2)) +
                                geom_smooth(aes(colour = "red"), se = F)+
                                geom_ribbon(aes(ymin = lower, ymax = upper, fill = A_cf), alpha= 0.2, color = NA) +
                                scale_fill_manual(values=c( "blue")) + scale_x_continuous(breaks = private$targetPoints)
                              return(g1)
                            }

                            results_CI_a = self$getCI(A_cf)
                            results_CI_a$lower = results_CI_a[, paste0("lower_bound","_", A_cf)]
                            results_CI_a$upper = results_CI_a[, paste0("upper_bound","_", A_cf)]
                            results_CI_a$pred = results_CI_a[, paste0("prediction","_", A_cf)]

                            g1 = ggplot2::ggplot(data =  results_CI_a, aes_string("t",  paste0("prediction","_", A_cf)))+ scale_x_continuous(breaks = private$targetPoints) +
                              geom_point(aes_string(x = "t", y =  paste0("prediction","_", A_cf)), legend=  F,   colour=alpha('red')) +
                              geom_point(aes_string(x = "t", y =  paste0("lower_bound","_", A_cf)), legend=  F,  colour=alpha('red', 0.2)) +
                              geom_point(aes_string(x = "t", y =  paste0("upper_bound","_", A_cf)), legend=  F,  colour=alpha('red', 0.2)) +
                              geom_smooth(aes(colour = "red"),se = F)+
                              geom_ribbon(aes(ymin = lower, ymax = upper, fill = A_cf), alpha= 0.2, color = NA) +
                              scale_x_continuous(breaks = private$targetPoints)
                            #

                            results_CI_b = self$getCI(1-A_cf)
                            results_CI_b$lower = results_CI_b[, paste0("lower_bound","_", 1-A_cf)]
                            results_CI_b$upper = results_CI_b[, paste0("upper_bound","_", 1-A_cf)]
                            results_CI_b$pred = results_CI_b[, paste0("prediction","_", 1-A_cf)]
                            g2 = ggplot2::ggplot(data = results_CI_b, aes_string("t",  paste0("prediction","_", 1-A_cf))) +
                              geom_point(aes_string(x = "t", y =  paste0("prediction","_", 1-A_cf)), legend=  F,   colour=alpha('blue')) +
                              geom_point(aes_string(x = "t", y =  paste0("lower_bound","_", 1-A_cf)), legend=  F, colour=alpha('blue', 0.2)) +
                              geom_point(aes_string(x = "t", y =  paste0("upper_bound","_", 1-A_cf)), legend=  F,  ,colour=alpha('blue', 0.2)) +
                              geom_smooth(aes(colour = "red"), se = F)+
                              geom_ribbon(aes(ymin = lower, ymax = upper, fill = 1-A_cf), alpha= 0.2, color = NA) +
                              scale_x_continuous(breaks = private$targetPoints)



                            x1 = results_CI_a

                            x0 = results_CI_b
                            x1$trt = "1"
                            x0$trt = "0"
                            names = c("t", "lower", "pred", "upper", "trt")
                            x1 = x1[,names]
                            x0 = x0[,names]



                            dat =data.frame(rbind(x1, x0))


                            g3 = ggplot2::ggplot(data =  dat, aes(x = t, y = pred, fill = trt)) + labs(fill = fill_lab ) + geom_point()+     geom_smooth( se = F)+ geom_ribbon(aes(ymin = lower, ymax = upper), alpha= 0.2, color = NA) +
                              scale_x_continuous(breaks = private$targetPoints) + xlab(x_lab) + ylab(y_lab) + ggtitle(main)

                            #g3 = ggplot2::ggplot(data =  dat, aes(x = t, y = pred, color = trt, fill = trt))  + geom_point()+ geom_smooth()+ geom_ribbon(aes(ymin = lower, ymax = upper, fill = trt,color, alpha= 0.05))
                            return(list(g1,g2,g3))
                          },
                          getCI = function(A_cf){
                            #est = private$estimator
                            surv_1 = self$computeEstimates(A_cf) #as.vector(colMeans(est$evalSurvival_failureMat(A_cf)))

                            EICMat = data.frame(self$computeEICMat(A_cf, full = T))
                            colnames(EICMat) = NULL
                            std_err = compute_simultaneous_ci(as.matrix(EICMat))

                            upper_bound1 <- surv_1 + 1.96 * std_err
                            lower_bound1 <- surv_1 - 1.96 * std_err


                            interv = private$targetPoints
                            results =data.frame(cbind(interv, lower_bound1, surv_1,upper_bound1))
                            colnames(results) = c("t", paste0("lower_bound","_", A_cf), paste0("prediction","_", A_cf), paste0("upper_bound","_", A_cf))
                            return(results)

                          }

                        )
)



#' @export
Estimator = R6::R6Class("Estimator",
                        private = list(
                          hazards_failure = NULL,
                          surv_failure = NULL,
                          #Survival censor corresponds to left evaluation
                          surv_censor = NULL,
                          g1W = NULL,
                          getSurv = function(ind){
                            haz = private$hazards_failure[[ind+1]]
                            surv = as.matrix(t(apply(haz, MARGIN=1, function(v){cumprod(1-v)})))
                            surv[surv<0.01] = 0.01
                            return(surv)
                          }


                        ),


                        public = list(
                          initialize = function(lstOfHazardFits){


                            private$hazards_failure = lstOfHazardFits[c("hazard_failure_0", "hazard_failure_1")]
                            private$surv_censor = lstOfHazardFits[c("surv_censor_0", "surv_censor_1")]
                            ncols=ncol(private$surv_censor[[1]])
                            nrows=nrow(private$surv_censor[[1]])
                            private$surv_censor[[1]] = cbind(rep(1,nrows), private$surv_censor[[1]][,1:(ncols-1)])
                            private$surv_censor[[2]] = cbind(rep(1,nrows), private$surv_censor[[2]][,1:(ncols-1)])
                            private$surv_failure = list(private$getSurv(0),private$getSurv(1))
                            private$g1W = lstOfHazardFits[["g1W"]]


                          },

                          evalHazard_failureMat = function(A_cf=NULL){

                            return(private$hazards_failure[[A_cf+1]])
                          },

                          evalSurvival_failureMat = function(A_cf=NULL){

                            return(private$surv_failure[[A_cf+1]])
                          },
                          evalSurvival_censorMatLeft = function(A_cf=NULL){

                            return(private$surv_censor[[A_cf+1]])
                          },
                          evalHazard_failure = function(t, A_cf=NULL){

                            return(private$hazards_failure[[A_cf+1]][,t])
                          },

                          evalSurvival_failure = function(t, A_cf=NULL){

                            return(private$surv_failure[[A_cf+1]][,t])
                          },
                          evalSurvival_censorLeft = function(t, A_cf=NULL){

                            return(private$surv_censor[[A_cf+1]][,t])
                          },
                          eval_g = function(A_cf = 1){
                            if(A_cf ==1){
                              return(private$g1W)
                            }
                            else{
                              return(1 - private$g1W)
                            }
                          },
                          update = function(haz, A_cf){
                            private$hazards_failure[[A_cf+1]] = haz
                            private$surv_failure[[A_cf+1]] = private$getSurv(A_cf)
                          }

                        )
)

#' @export
expit = function(x){
  return(exp(x)/(1+exp(x)))
}

#' @export
logit = function(x){
  return(log(x/(1-x)))
}

#' @export
compute_simultaneous_ci <- function(eic_fit) {
  rownames(eic_fit) = NULL
  colnames(eic_fit) = NULL
  # compute the value to +- around the Psi_n
  n <- nrow(eic_fit)
  sigma_squared <- stats::cov(eic_fit)
  sigma <- stats::cor(eic_fit)
  # impute when the variance are zero
  sigma_squared[is.na(sigma_squared)] <- 1e-10
  sigma[is.na(sigma)] <- 1e-10

  variance_marginal <- diag(sigma_squared)
  q <- compute_q(corr = sigma, B = 1e3, alpha = 0.05)
  return(sqrt(variance_marginal) / sqrt(n) * q)
}

#' @export
compute_q <- function(corr, B = 1e3, alpha = 0.05) {
  dim <- nrow(corr)
  z <- apply(
    abs(MASS::mvrnorm(B, mu = rep(0, dim), Sigma = corr)), 1, max
  )
  return(as.numeric(stats::quantile(z, 1 - alpha)))
}


