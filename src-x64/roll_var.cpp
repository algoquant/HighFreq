#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]

NumericVector roll_var(NumericVector re_turns, int look_back) {
  int len_gth = re_turns.size();
  NumericVector var_iance(len_gth);
  NumericVector mean_s(len_gth);
  
  var_iance[look_back-1] = sum(re_turns[Range(0, (look_back-1))]);
  for (int i = look_back; i < len_gth; ++i) {
//    var_iance[i] = var_iance[i-1] * exp(re_turns[i]);
    var_iance[i] = var_iance[i] + re_turns[i] - re_turns[i-look_back];
  }
  
  return var_iance;
}
