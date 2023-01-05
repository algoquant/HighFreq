////////////////////////////
// Test functions to be included in HighFreq.cpp
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/HighFreq/sandbox/test_fun.cpp")


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // include RcppArmadillo header file
using namespace arma; // use Armadillo C++ namespace


// [[Rcpp::export]]
arma::mat calc_invmat(arma::mat& matrixv) {
  
  return arma::inv(matrixv);
  
}  // end calc_invmat



// Calculate the alpha and beta regression coefficients
// without using the Moore-Penrose inverse arma::pinv().
// [[Rcpp::export]]
arma::mat calc_betas(const arma::mat& response, 
                     const arma::mat& predictor) {
  
  // arma::uword nrows = predictor.n_rows;
  // arma::uword ncols = predictor.n_cols;
  // Response and Predictor means
  arma::mat respmeans = arma::mean(response, 0);
  arma::mat predmeans = arma::mean(predictor, 0);
  arma::mat predz = predictor.each_row() - predmeans;
  
  // Initialize the variables
  arma::mat covarespred = arma::trans(response.each_row()-respmeans)*predz;
  arma::mat covarpred = arma::trans(predz)*predz;
  arma::mat betas = covarespred*arma::inv(covarpred);
  arma::mat alphas = respmeans - arma::dot(betas, predmeans);
  
  return arma::join_rows(alphas, betas);
  
}  // end calc_betas

