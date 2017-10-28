#include <Rcpp.h>
using namespace Rcpp;

// The function roll_wsum() calculates the weighted rolling sum over a vector
//' @export
// [[Rcpp::export]]
NumericVector roll_wsum(NumericVector re_turns, NumericVector wei_ghts) {
  int len_gth = re_turns.size();
  int look_back = wei_ghts.size();
  NumericVector rolling_sum(len_gth);
  
  rolling_sum[look_back-1] = sum(wei_ghts * re_turns[Range(0, (look_back-1))]);
  for (int i = look_back; i < len_gth; ++i) {
//    re_turns[i] = the_ta*(eq_price - rolling_sum[i-1]) + vol_at*r_norm[i-1];
    rolling_sum[i] = sum(wei_ghts * re_turns[Range(i-look_back+1, i)]);
  }
  
  return rolling_sum;
}
