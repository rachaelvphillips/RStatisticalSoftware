// [[Rcpp::depends(RcppEigen)]]
#ifndef HAL9001_TYPES_H
#define HAL9001_TYPES_H

#include <RcppEigen.h>
using namespace Rcpp;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseMatrix<int> IntSpMat;
typedef Eigen::Map<SpMat> MSpMat;
typedef MSpMat::InnerIterator MInIterMat;

typedef SpMat::InnerIterator InIterMat;
//typedef SpMat::InnerVectorReturnType InVec;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;


class OperatorInterface {
public:
  virtual double operator() (SpVec& v, int index = 0) {}
  virtual double operator() (double val, int index = 0) {}
  virtual SpVec operator() (SpVec& v_in) {}
  virtual NumericVector operator() (NumericVector& v_in){}
  virtual void update_step_size(double step_size_) = 0;
  virtual void update_cache_vars(double old_x_i, double new_x_i, int index) = 0;
  virtual void update_cache_vars(SpVec& x, int rank, int num_threads) = 0;
};

class prox_l1 : public OperatorInterface {
public:
  double weight;
  double step_size;
  //Applies soft threshold to vector at given index
  double operator() (NumericVector& v, int index = 0) {
    double val;
    double threshold = step_size * weight;
    val = v[index];
    if (val > threshold) {
      return val - threshold;
    } else if (val < -threshold) {
      return val + threshold;
    } else {
      return 0.;
    }
  }
  //Applies soft threshold to value
  double operator() (double val, int index = 0) {
    double threshold = step_size * weight;
    if (val > threshold) {
      return val - threshold;
    } else if (val < -threshold) {
      return val + threshold;
    } else {
      return 0.;
    }
  }
  // full update operator
  // Performs full soft thresholding on vector v_in
  SpVec operator() (SpVec& v_in) {
    
    
    int len = v_in.size();
    SpVec v_out(len);
    double threshold = step_size * weight;
    double val;
    int i;
    for (InIterVec it(v_in); it; ++it)
    {
      val = it.value(); // == vec[ it.index() ]
      i = it.index();
      if (val > threshold) {
        v_out.insert(i) = val - threshold;
      } else if (val < -threshold) {
        v_out.insert(i) = val + threshold;
      }
      
    }
    return(v_out);
  }
  
  NumericVector operator() (NumericVector& v_in) {
    
    
    
    double threshold = step_size * weight;
    double val;
    int len = v_in.size();
    
    for (int i=0; i < len; ++i)
    {
      val = v_in[i]; // == vec[ it.index() ]
      
      if (val > threshold) {
        v_in[i] = val - threshold;
      } else if (val < -threshold) {
        v_in[i] = val + threshold;
      }
      else{
        v_in[i] =0;
      }
      
    }
    return(v_in);
  }
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
  }
  void update_cache_vars(SpVec* x, int rank, int num_threads) {
  }
  
  void update_step_size(double step_size_) {
    step_size = step_size_;
  }
  
  prox_l1 (double step_size_, double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
  }
  prox_l1 () {
    step_size = 0.;
    weight = 1.;
  }
};



template <typename Mat>
class forward_grad_for_square_loss : public OperatorInterface {
public:
  const SpMat A;
  NumericVector b, Atx;  
  const SpMat At; // this is for the sync-parallel stuff
  double weight;
  double step_size;
  
  double operator() (Vector* x, int index) {
   return(0)
  }
  
  void operator() (Vector* v_in, Vector* v_out) {
    int m = A->rows(), n = A->cols();
    Vector Atv_in(m, 0.);
    trans_multiply(*A, *v_in, Atv_in);
    add(Atv_in, *b, -1.);
    Vector temp(n, 0.);
    multiply(*A, Atv_in, temp);
    scale(temp, -step_size * weight);
    add(temp, *v_in);
    for (int i = 0; i < m; i++) {
      (*v_out)[i] = temp[i];
    }
  }
  
  double operator() (double val, int index = 1) {
    return DOUBLE_MARKER;
  }
  
  void update_step_size(double step_size_) {
    step_size = step_size_;
  }  
  
  void update_cache_vars(double old_x_i, double new_x_i, int index) {
    add(Atx, A, index, -old_x_i + new_x_i);
  }
  
  void update_cache_vars(Vector* x, int rank, int num_threads){
    if(At == nullptr){
      std::cout << "No matrix transpose available, aborting" << std::endl ;
      assert(At != nullptr);
    }
    int m = At->rows(); //y=A'*x
    int block_size = m/num_threads;
    int start_idx = rank*(block_size);
    int end_idx = (rank == num_threads-1)? m : start_idx+block_size;
    for(int iter=start_idx; iter != end_idx; ++iter){
      (*Atx)[iter]=dot(At, x, iter);
    }
  }
  
  forward_grad_for_square_loss () {
    step_size = 0.;
    weight = 1.;
  }
  forward_grad_for_square_loss (double step_size_, double weight_ = 1.) {
    step_size = step_size_;
    weight = weight_;
  }
  
  forward_grad_for_square_loss (SpMat& A_, NumericVector& b_, NumericVector& Atx_,
                                double step_size_ = 1., double weight_ = 1.,
                                SpMat& At_) {
    step_size = step_size_;
    weight = weight_;
    A = A_;
    b = b_;
   
    At = At_;
  }
};

#endif //HAL9001_TYPES_H
