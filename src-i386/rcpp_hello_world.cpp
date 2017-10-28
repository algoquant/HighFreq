
#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
CharacterVector rcpp_hello_world() {
  
  CharacterVector x = CharacterVector::create( "hello", "world" )  ;
  //    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
  //    List z            = List::create( x, y ) ;
  
  return x ;
}
