////////////////////////////
// Test functions to be included in HighFreq.cpp
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/HighFreq/sandbox/test_fun.cpp")


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
arma::mat calc_invmatr(arma::mat& matrixv) {
  
  return arma::inv(matrixv);
  
}  // end calc_invmatr

