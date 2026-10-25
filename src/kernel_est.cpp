#define ARMA_NO_DEBUG
#define ARMA_USE_BLAS

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

// define the K(u) cumulative kernel function
static inline double Kernel(int Ktype, double x){
  if(Ktype == 0) return R::pnorm(x, 0.0, 1.0, 1, 0); // Gauss kernel
  else{ // if(Ktype == 1)
    if(fabs(x) <= 1) return 0.5 + 3.0*(x - x*x*x/3.0)/4.0; // Epanechnikov kernel
    else if(x > 1) return 1.0;
    else return 0.0;
  }
  // else return 0;
}

// define the k(u) kernel function
static inline double kernel(int Ktype, double x){
  if(Ktype == 0) return R::dnorm(x, 0.0, 1.0, 0); // Gauss kernel
  else{ // if(Ktype == 1)
    if(fabs(x) <= 1) return 0.75*(1 - x*x); // Epanechnikov kernel
    else return 0.0;
  }
  // else return 0;
}

// define the F_kernel
// [[Rcpp::export]]
NumericVector cdf_kernel_C(NumericVector x, NumericVector X, int Ktype,
                           double bwd){
  // x: points to compute CDF
  // X: data
  // Ktype: type of kernel, 0 - Gaussian, 1 - Epanechnikov
  // bwd: bandwidth
  int n = X.size(), m = x.size();
  double uij = 0.0;
  NumericVector num(m);
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      uij = (x[i] - X[j])/bwd;
      num[i] += Kernel(Ktype, uij);
    }
  }
  return num/n;
}

// define the f_kernel
// [[Rcpp::export]]
NumericVector pdf_kernel_C(NumericVector x, NumericVector X, int Ktype,
                           double bwd){
  // x: points to compute CDF
  // X: data
  // Ktype: type of kernel, 0 - Gaussian, 1 - Epanechnikov
  // bwd: bandwidth
  int n = X.size(), m = x.size();
  double uij = 0.0;
  NumericVector num(m);
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      uij = (x[i] - X[j])/bwd;
      num[i] += kernel(Ktype, uij);
    }
  }
  return num/n/bwd;
}

// define Simpson's rule integral
// static inline
// [[Rcpp::export]]
double simpson(arma::vec fx, int nx, double hx){
  // nx must be an even number - the number of intervals
  double res = 0.0;
  int pl, pr, pmid;
  for(int i = 0; i < nx/2; i++){
    pr = 2*(i + 1);
    pl = pr - 2;
    pmid = pr - 1;
    res += fx[pl] + 4.0*fx[pmid] + fx[pr];
  }
  return res*hx/3.0;
}

// define the leave-one-out F_kernel
// [[Rcpp::export]]
double cdf_loocv_kernel_C(NumericVector x, NumericVector X, int Ktype,
                          double bwd, double hx){
  int n = X.size(), m = x.size();
  double uij = 0.0;
  arma::mat w(m, n, fill::zeros);
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      uij = (x[i] - X[j])/bwd;
      w(i, j) = Kernel(Ktype, uij);
    }
  }
  arma::vec nF_kernel(m, fill::zeros);
  nF_kernel = sum(w, 1);
  arma::mat F_loo(m, n, fill::zeros);
  for(int i = 0; i < n; i++){
    F_loo.col(i) = nF_kernel - w.col(i);
  }
  F_loo = F_loo/(n - 1);
  arma::vec res_cv(n, fill::zeros);
  arma::vec ff(m, fill::zeros);
  for (int i = 0; i < n; i++){
    // arma::vec ff(m, fill::zeros);
    for(int j = 0; j < m; j++){
      if(x[j] >= X[i]) ff[j] = pow(1.0 - F_loo(j, i), 2);
      else ff[j] = pow(F_loo(j, i), 2);
    }
    res_cv[i] = simpson(ff, m - 1, hx);
  }
  return sum(res_cv)/n;
}

// define the LOOCV bandwidth
// [[Rcpp::export]]
List cv_bwd_C(NumericVector x, NumericVector X, int Ktype, double hx,
              NumericVector bwd_seq){
  arma::vec bwds = as<vec>(bwd_seq);
  double bwd_pts = bwds.size();
  arma::vec cv_out(bwd_pts, fill::zeros);
  for(int i = 0; i < bwd_pts; i++){
    cv_out[i] = cdf_loocv_kernel_C(x, X, Ktype, bwds[i], hx);
  }
  double id_h_cv = cv_out.index_min();
  List out;
  out["bwd_seq"] = wrap(bwds);
  out["cv_out"] = wrap(cv_out);
  out["h_opt"] = bwds[id_h_cv];
  return out;
}

