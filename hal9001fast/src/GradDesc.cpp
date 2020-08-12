// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include "hal9001fast_types.h"
using namespace Rcpp;

// [[Rcpp::export]]
SpMat calc_canon_grad(SpMat& X, const SpMat& Y, const SpMat& beta) {
  SpMat link;


  link  =  Y - X * beta;



  return(link);
}
