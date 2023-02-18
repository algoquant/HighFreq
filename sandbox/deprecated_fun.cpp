////////////////////////////
// Deprecated functions or old versions of functions excluded from HighFreq.cpp
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/HighFreq/sandbox/deprecated_fun.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
using namespace arma;
// Use STL
using namespace std;


////////////////////////////////////////////////////////////
//' Calculate the rolling sums over a single-column \emph{time series} or a
//' single-column \emph{matrix} using \emph{Rcpp}.
//' 
//' @param \code{tseries} A single-column \emph{time series} or a single-column
//'   \emph{matrix}.
//' 
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of elements of data used for calculating the sum.
//'
//' @param \code{weightv} A single-column \emph{matrix} of weights.
//'
//' @return A single-column \emph{matrix} of the same length as the argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{roll_vec()} calculates a single-column \emph{matrix} of
//'   rolling sums, over a single-column \emph{matrix} of data, using fast
//'   \emph{Rcpp} \code{C++} code.  The function \code{roll_vec()} is several
//'   times faster than \code{rutils::roll_sum()} which uses vectorized \code{R}
//'   code.
//'
//' @examples
//' \dontrun{
//' # Define a single-column matrix of returns
//' retsp <- zoo::coredata(na.omit(rutils::etfenv$returns$VTI))
//' # Calculate rolling sums over 11-period look-back intervals
//' sum_rolling <- HighFreq::roll_vec(retsp, look_back=11)
//' # Compare HighFreq::roll_vec() with rutils::roll_sum()
//' all.equal(HighFreq::roll_vec(retsp, look_back=11), 
//'          rutils::roll_sum(retsp, look_back=11), 
//'          check.attributes=FALSE)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_vec(retsp, look_back=11),
//'   Rcode=rutils::roll_sum(retsp, look_back=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_vec(const arma::mat& tseries, 
                  arma::uword look_back, 
                  const arma::colvec& weightv) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat sumr(nrows, 1);
  
  // Warmup period
  sumr[0] = tseries[0];
  for (arma::uword it = 1; it < look_back; it++) {
    sumr(it) = sumr(it-1) + tseries(it);
  }  // end for
  
  // Remaining period
  for (arma::uword it = look_back; it < nrows; it++) {
    sumr(it) = sumr(it-1) + tseries(it) - tseries(it-look_back);
  }  // end for
  
  return sumr;
  
}  // end roll_vec




////////////////////////////////////////////////////////////
//' Calculate the rolling weighted sums over a single-column \emph{time series}
//' or a single-column \emph{matrix} using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A single-column \emph{time series} or a single-column
//'   \emph{matrix}.
//' 
//' @param \code{weightv} A single-column \emph{matrix} of weights.
//'
//' @return A single-column \emph{matrix} of the same length as the argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{roll_vecw()} calculates the rolling weighted sums of a
//'   single-column \emph{matrix} over its past values (a convolution with the
//'   single-column \emph{matrix} of weights), using \code{RcppArmadillo}. It
//'   performs a similar calculation as the standard \code{R} function
//'   \cr\code{stats::filter(x=series, filter=weightv, method="convolution",
//'   sides=1)}, but it's over \code{6} times faster, and it doesn't produce any
//'   \code{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Define a single-column matrix of returns
//' retsp <- zoo::coredata(na.omit(rutils::etfenv$returns$VTI))
//' # Create simple weights
//' weightv <- c(1, rep(0, 10))
//' # Calculate rolling weighted sums
//' weighted <- HighFreq::roll_vecw(tseries=retsp, weightv=weightv)
//' # Compare with original
//' all.equal(zoo::coredata(retsp), weighted, check.attributes=FALSE)
//' # Second example
//' # Create exponentially decaying weights
//' weightv <- exp(-0.2*1:11)
//' weightv <- weightv/sum(weightv)
//' # Calculate rolling weighted sums
//' weighted <- HighFreq::roll_vecw(tseries=retsp, weightv=weightv)
//' # Calculate rolling weighted sums using filter()
//' filtered <- stats::filter(x=retsp, filter=weightv, method="convolution", sides=1)
//' # Compare both methods
//' all.equal(filtered[-(1:11)], weighted[-(1:11)], check.attributes=FALSE)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_vecw(tseries=retsp, weightv=weightv),
//'   Rcode=stats::filter(x=retsp, filter=weightv, method="convolution", sides=1),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_vecw(const arma::mat& tseries, 
                   const arma::colvec& weightv) {
  
  arma::uword nrows = tseries.n_rows;
  arma::uword look_back = weightv.n_rows;
  arma::mat sumr(nrows, 1);
  arma::mat weightr = arma::reverse(weightv);
  // arma::mat weightr = weightv;
  
  // Warmup period
  sumr.rows(0, look_back-2) = tseries.rows(0, look_back-2);
  
  // Remaining periods
  for (arma::uword it = look_back-1; it < nrows; it++) {
    sumr(it) = arma::dot(weightr, tseries.rows(it-look_back+1, it));
  }  // end for
  
  return sumr;
  
}  // end roll_vecw

