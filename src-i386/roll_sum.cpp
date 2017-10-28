#include <Rcpp.h>
using namespace Rcpp;

// The function roll_sum() calculates the rolling sum over a vector
//' @export
// [[Rcpp::export]]
NumericVector roll_sum(NumericVector re_turns, int look_back) {
  int len_gth = re_turns.size();
  NumericVector rolling_sum(len_gth);
  
  rolling_sum[look_back-1] = sum(re_turns[Range(0, (look_back-1))]);
  for (int i = look_back; i < len_gth; ++i) {
//    re_turns[i] = the_ta*(eq_price - rolling_sum[i-1]) + vol_at*r_norm[i-1];
    rolling_sum[i] = rolling_sum[i-1] + re_turns[i] - re_turns[i-look_back];
  }
  
  return rolling_sum;
}
