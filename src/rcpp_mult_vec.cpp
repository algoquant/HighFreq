#include <Rcpp.h>
using namespace Rcpp;

// The function rcpp_mult_vec() multiplies two vectors
//' @export
// [[Rcpp::export]]
NumericVector rcpp_mult_vec(NumericVector x, NumericVector y) {
  return x * y;
}
