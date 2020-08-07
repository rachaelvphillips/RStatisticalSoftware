n=1000

#Always do zero order so that the l1 penalization is Vnorm (even if continuous).
#Do binning for computation.


screen_vars_hal <- function(X, Y, fit = NULL, bins = 750){

  if(is.null(fit)){
    Xq = quantizer(X, bins = bins)
    fit = fit_hal(X=Xq, Y=Y, max_degree = 1)
    #fit = fit_hal(X=X, Y=Y, max_degree = 1, basis_list = fit$basis_list, cv_select = F, lambda = fit$lambda_star)
  }



  basis = fit$basis_list
  coefs = fit$coefs
  keep = which(fit$coefs!=0)
  coefs = coefs[keep]
  basis = basis[keep]

  cols = sapply(basis, function(term){term$cols
  })

  colsSelected = unique(cols)

  keep_basis = which(fit$coefs !=0)
  getMatchCol = function(col, basis){
    which(sapply(basis, function(b){b$cols==col}))
  }

  filter_by_dcor = function(col, thresh = 0.05, other = NULL){
    if(!is.null(other)){
      test =pdcor.test(X[,col],Y, X[,setdiff(other,col)])
    }
    test = dcorT.test(X[,col],Y)
    return(test$p.value <= thresh)
  }

  keep_cols = sapply(colsSelected, filter_by_dcor)
  print(sum(1-keep_cols))
  colsSelected = colsSelected[keep_cols]

  counts = table(cols)

  getVnorm = function(col){
    index = which(colsSelected==col)
    sum(abs(coefs[index]))
  }

  Vnorms = sapply(colsSelected, getVnorm)
  names(Vnorms) = colsSelected





  return(list(counts, Vnorms,  colsSelected))

}




filter_by_dcor = function(col, thresh = 0.1, other = NULL){
  if(!is.null(other)){
    test =pdcor.test(X[,col],Y, X[,setdiff(other,col)], R=50)
  }
  else{
    test = dcorT.test(X[,col],Y)
  }
  return(test$p.value <= thresh )
}




