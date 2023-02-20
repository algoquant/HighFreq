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

// run_zscores() is deprecated because run_reg() is the new version
////////////////////////////////////////////////////////////
//' Calculate the z-scores of trailing regressions of streaming \emph{time
//' series} of returns.
//' 
//' @param \code{respv} A single-column \emph{time series} or a single-column
//'   \emph{matrix} of response data.
//' 
//' @param \code{predv} A \emph{time series} or a \emph{matrix} of predictor
//'   data.
//' 
//' @param \code{lambda} A decay factor which multiplies past
//'   estimates.
//'   
//' @param \code{demean} A \emph{Boolean} specifying whether the \emph{z-scores}
//'   should be de-meaned (the default is \code{demean = TRUE}).
//'
//' @return A \emph{matrix} with the z-scores, betas, and the variances of the
//'   predictor data.
//'
//' @details
//'   The function \code{run_zscores()} calculates the vectors of \emph{betas}
//'   \eqn{\beta_t} and the residuals \eqn{\epsilon_t} of trailing regressions
//'   by recursively weighting the current estimates with past estimates, using
//'   the decay factor \eqn{\lambda}:
//'   \deqn{
//'     \sigma^2_t = (1-\lambda) p^2_t + \lambda \sigma^2_{t-1}
//'   }
//'   \deqn{
//'     cov_t = (1-\lambda) r_t p_t + \lambda cov_{t-1}
//'   }
//'   \deqn{
//'     \beta_t = (1-\lambda) \frac{cov_t}{\sigma^2_t} + \lambda \beta_{t-1}
//'   }
//'   \deqn{
//'     \epsilon_t = (1-\lambda) (r_t - \beta_t p_t) + \lambda \epsilon_{t-1}
//'   }
//'   Where \eqn{cov_t} is the vector of covariances between the
//'   response and predictor returns, at time \eqn{t};
//'   \eqn{\sigma^2_t} is the predictor variance,
//'   and \eqn{r_t} and \eqn{p_t} are the streaming returns of the response
//'   and predictor data.
//' 
//'   The above formulas for \eqn{\sigma^2} and \eqn{cov} are
//'   approximate because they don't subtract the means before squaring the
//'   returns.  But they're very good approximations for daily returns.
//' 
//'   The matrices \eqn{\sigma^2}, \eqn{cov}, and \eqn{\beta} have the
//'   same number of rows as the input argument \code{predv}.
//'
//'   If the argument \code{demean = TRUE} (the default) then the
//'   \emph{z-scores} \eqn{z_t} are calculated as equal to the residuals
//'   \eqn{\epsilon_t} minus their means \eqn{\bar{\epsilon}}, divided by their
//'   volatilities \eqn{\varsigma_t}:
//'   \deqn{
//'     z_t = \frac{\epsilon_t - \bar{\epsilon}}{\varsigma_t}
//'   }
//'   If the argument \code{demean = FALSE} then the \emph{z-scores} are
//'   only divided by their volatilities without subtracting their means:
//'   \deqn{
//'     z_t = \frac{\epsilon_t}{\varsigma_t}
//'   }
//' 
//'   The above recursive formulas are convenient for processing live streaming
//'   data because they don't require maintaining a buffer of past data.
//'   The formulas are equivalent to a convolution with exponentially decaying
//'   weights, but they're much faster to calculate.
//'   Using exponentially decaying weights is more natural than using a sliding
//'   look-back interval, because it gradually "forgets" about the past data.
//' 
//'   The value of the decay factor \eqn{\lambda} must be in the range between
//'   \code{0} and \code{1}.
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the trailing \emph{z-score} values have
//'   a stronger dependence on past data.  This is equivalent to a long
//'   look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the trailing \emph{z-score} values
//'   have a weaker dependence on past data.  This is equivalent to a short
//'   look-back interval.
//' 
//'   The function \code{run_zscores()} returns multiple columns of data. If the
//'   matrix \code{predv} has \code{n} columns then \code{run_zscores()}
//'   returns a matrix with \code{2n+1} columns.  The first column contains the
//'   \emph{z-scores}, and the remaining columns contain the \emph{betas} and
//'   the \emph{variances} of the predictor data.
//' 
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' retsp <- na.omit(rutils::etfenv$returns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' respv <- retsp[, 1]
//' # Predictor matrix equals VTI and IEF returns
//' predv <- retsp[, -1]
//' # Calculate the trailing z-scores
//' lambda <- 0.9
//' zscores <- HighFreq::run_zscores(respv=respv, predv=predv, lambda=lambda)
//' # Plot the trailing z-scores
//' datav <- cbind(cumsum(respv), zscores[, 1])
//' colnames(datav) <- c("XLF", "zscores")
//' colnamev <- colnames(datav)
//' dygraphs::dygraph(datav, main="Z-Scores of XLF Versus VTI and IEF") %>%
//'   dyAxis("y", label=colnamev[1], independentTicks=TRUE) %>%
//'   dyAxis("y2", label=colnamev[2], independentTicks=TRUE) %>%
//'   dySeries(name=colnamev[1], axis="y", label=colnamev[1], strokeWidth=2, col="blue") %>%
//'   dySeries(name=colnamev[2], axis="y2", label=colnamev[2], strokeWidth=2, col="red")
//' }
//' 
 arma::mat run_zscores(const arma::mat& respv, 
                       const arma::mat& predv,
                       double lambda, // Decay factor which multiplies the past values 
                       bool demean = true) {
   
   arma::uword nrows = predv.n_rows;
   arma::uword ncols = predv.n_cols;
   // arma::mat var1 = arma::square(tseries.col(0));
   arma::mat varv = arma::square(predv);
   arma::mat betas = arma::ones(nrows, ncols);
   arma::mat zscores = arma::ones(nrows, 1);
   arma::mat vars = arma::ones(nrows, 1);
   arma::mat meanm = arma::zeros(nrows, 1);
   double lambda1 = 1-lambda;
   
   // Multiply each column of predictor by the response
   arma::mat covars = predv;
   covars.each_col() %= respv;
   
   // Perform loop over the rows
   for (arma::uword it = 1; it < nrows; it++) {
     // Calculate the z-score as the weighted sum of products of returns.
     // cout << "Calculating vars: " << it << endl;
     varv.row(it) = lambda1*varv.row(it) + lambda*varv.row(it-1);  // this is wrong
     // cout << "Calculating covars: " << it << endl;
     covars.row(it) = lambda1*covars.row(it) + lambda*covars.row(it-1);
     // cout << "Calculating betas: " << it << endl;
     betas.row(it) = lambda1*covars.row(it)/vars.row(it) + lambda*betas.row(it-1);
     // cout << "Calculating zscores: " << it << endl;
     zscores.row(it) = lambda1*(respv.row(it) - arma::dot(betas.row(it), predv.row(it))) + lambda*zscores.row(it-1);
     // Calculate the mean and variance of the z-scores.
     meanm.row(it) = lambda1*zscores.row(it) + lambda*meanm.row(it-1);
     vars.row(it) = lambda1*arma::square(zscores.row(it) - zscores.row(it-1)) + lambda*vars.row(it-1);
   }  // end for
   
   if (demean)
     return arma::join_rows((zscores - meanm)/arma::sqrt(vars), betas, vars);
   else
     return arma::join_rows(zscores/arma::sqrt(vars), betas, vars);
   
 }  // end run_zscores
 
 
 
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

