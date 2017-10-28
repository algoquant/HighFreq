#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. 
// You can source this function into an R session using the 
// function Rcpp::sourceCpp() 
// (or via the Source button on the editor toolbar).
// Learn more about Rcpp at:
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/

// The function rcpp_ou_proc() simulates an Ornstein-Uhlenbeck process
//' @export
// [[Rcpp::export]]
NumericVector rcpp_ou_proc(int len_gth, double eq_price, double vol_at, double the_ta, NumericVector r_norm) {
  NumericVector price_s(len_gth);
  NumericVector re_turns(len_gth);
  price_s[0] = eq_price;
  for (int i = 1; i < len_gth; ++i) {
    re_turns[i] = the_ta*(eq_price - price_s[i-1]) + vol_at*r_norm[i-1];
    price_s[i] = price_s[i-1] * exp(re_turns[i]);
  }
  return price_s;
}
