#include <Rcpp.h>
#include <cmath>

// Tensor trapezoidal rule on x <= y, with half of each diagonal weight.
// [[Rcpp::export]]
double arocs_integral_grid_C(Rcpp::NumericVector f1,
                             Rcpp::NumericVector f2,
                             Rcpp::NumericVector f3, double dx) {
  const int n = f1.size();
  if (n < 2 || f2.size() != n || f3.size() != n || !R_finite(dx) || dx <= 0)
    Rcpp::stop("Invalid integration grid.");
  long double total = 0;
  for (int i = 0; i < n; ++i) {
    if (i % 128 == 0) Rcpp::checkUserInterrupt();
    double wi = (i == 0 || i == n - 1) ? 0.5 : 1.0;
    for (int j = i; j < n; ++j) {
      double wj = (j == 0 || j == n - 1) ? 0.5 : 1.0;
      double value = std::hypot(std::hypot(f1[i] * f3[j], f2[i] * f3[j]),
                                f1[i] * f2[j]);
      total += wi * wj * (i == j ? 0.5 : 1.0) * value;
    }
  }
  return static_cast<double>(total * dx * dx);
}
