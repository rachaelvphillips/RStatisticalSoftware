// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppNumerical.h>
using namespace Numer;
using namespace Rcpp;
// [[Rcpp::export]]
int eval_basis(NumericVector x, const List& basis) {
  IntegerVector cols = as<IntegerVector>(basis["cols"]);
  NumericVector cutoffs = as<NumericVector>(basis["cutoffs"]);
  int p = cols.length();
  for (int i = 0; i<p; i++) {
    double obs = x[cols[i] - 1]; // using 1-indexing for basis columns
    if(!(obs > cutoffs[i])) {
      return(0);
    }
  }
  return(1);
}


/*** R
eval_basis(c(2,2,3), list("cols" = c(1,2), "cutoffs" = c(2,1)))
*/

class Centered_Basis: public Func
{
private:
  List basis;
  NumericVector val;
public:
  Centered_Basis(NumericVector val, List& basis) : basis(basis), val(val){}

  double operator()(const double& x) const
  {
    NumericVector z = val;
    z.push_front(x);
    return eval_basis(z, basis);
  }
};
// [[Rcpp::export]]
double integrate(NumericVector val, List& basis, double lower, double upper){

  Centered_Basis f(val, basis);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);
  return res;
};


