#include <Rcpp.h>
using namespace Rcpp;

//' The function rcpp_mult() multiplies two numbers.
//' @export
// [[Rcpp::export]]
double rcpp_mult(double x, double y) {
  return x * y;
}
