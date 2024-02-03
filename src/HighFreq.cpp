// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
#include <vector>
// Create hooks for RcppArmadillo
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
// Use STL
using namespace std;
// For eigen solver SymEigsSolver
using namespace arma::newarp;

////////////////////////////////////////////////////////////
// Rcpp and RcppArmadillo functions for package HighFreq
////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////
// Functions miscellaneous
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
//' Create a named list of model parameters that can be passed into regression
//' and machine learning functions.
//' 
//' @param \code{method} A \emph{character string} specifying the type of
//'   regression model (the default is \code{method = "least_squares"}).
//'   
//' @param \code{intercept} A \emph{Boolean} specifying whether an intercept
//'   term should be added to the predictor (the default is \code{intercept =
//'   TRUE}).
//'
//' @param \code{singmin} A \emph{numeric} threshold level for discarding
//'   small \emph{singular values} in order to regularize the inverse of the
//'   predictor matrix (the default is \code{1e-5}).
//'   
//' @param \code{dimax} An \emph{integer} equal to the number of \emph{singular
//'   values} used for calculating the \emph{reduced inverse} of the
//'   predictor matrix (the default is \code{dimax = 0} - standard matrix
//'   inverse using all the \emph{singular values}).
//'   
//' @param \code{confl} The confidence level for calculating the quantiles of
//'   returns (the default is \code{confl = 0.75}).
//'
//' @param \code{alpha} The shrinkage intensity of \code{returns} (with values
//'   between \code{0} and \code{1} - the default is \code{0}).
//'   
//' @return A named list of model parameters that can be passed into regression
//' and machine learning functions.
//'   
//' @details
//'   The function \code{param_reg()} creates a named list of model parameters
//'   that can be passed into regression and machine learning functions.  For
//'   example into the functions \code{calc_reg()} and \code{roll_reg()}.
//'   
//'   The function \code{param_reg()} simplifies the creation of regression
//'   parameter lists.  The users can create a parameter list with the default
//'   values, or they can specify custom parameter values.
//'
//' @examples
//' \dontrun{
//' # Create a default list of regression parameters
//' controlv <- HighFreq::param_reg()
//' unlist(controlv)
//' # Create a custom list of regression parameters
//' controlv <- HighFreq::param_reg(intercept=FALSE, method="regular", dimax=4)
//' unlist(controlv)
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List param_reg(std::string regmod = "least_squares",  // Type of regression model
                     bool intercept = true,  // Add unit intercept column to the predictor matrix?
                     double singmin = 1e-5, // Threshold level for discarding small singular values
                     arma::uword dimax = 0, // Number of eigen vectors for dimension reduction
                     std::string residscale = "none",  // Method for scaling the residuals
                     double confl = 0.1, // Confidence level for calculating the quantiles of returns
                     double alpha = 0.0) {  // Shrinkage intensity of returns
  
  Rcpp::List controlv = Rcpp::List::create(Rcpp::Named("regmod") = regmod,
                                           Rcpp::Named("intercept") = intercept,
                                           Rcpp::Named("singmin") = singmin,
                                           Rcpp::Named("dimax") = dimax,
                                           Rcpp::Named("residscale") = residscale,
                                           Rcpp::Named("confl") = confl,
                                           Rcpp::Named("alpha") = alpha);
  
  return controlv;
  
}  // end param_reg


////////////////////////////////////////////////////////////
//' Create a named list of model parameters that can be passed into portfolio
//' optimization functions.
//' 
//' @param \code{method} A \emph{character string} specifying the method for
//'   calculating the portfolio weights (the default is \code{method =
//'   "sharpem"}).
//'   
//' @param \code{singmin} A \emph{numeric} threshold level for discarding
//'   small \emph{singular values} in order to regularize the inverse of the
//'   \code{covariance matrix} of \code{returns} (the default is \code{1e-5}).
//'   
//' @param \code{dimax} An \emph{integer} equal to the number of \emph{singular
//'   values} used for calculating the \emph{reduced inverse} of the
//'   \code{covariance matrix} of \code{returns} matrix (the default is
//'   \code{dimax = 0} - standard matrix inverse using all the \emph{singular
//'   values}).
//'   
//' @param \code{confl} The confidence level for calculating the quantiles of
//'   returns (the default is \code{confl = 0.75}).
//'
//' @param \code{alpha} The shrinkage intensity of \code{returns} (with values
//'   between \code{0} and \code{1} - the default is \code{0}).
//' 
//' @param \code{rankw} A \emph{Boolean} specifying whether the weights should
//'   be ranked (the default is \code{rankw = FALSE}).
//'
//' @param \code{centerw} A \emph{Boolean} specifying whether the weights should
//'   be centered (the default is \code{centerw = FALSE}).
//'
//' @param \code{scalew} A \emph{character string} specifying the method for
//'   scaling the weights (the default is \code{scalew = "voltarget"}).
//'
//' @param \code{voltarget} A \emph{numeric} volatility target for scaling the
//'   weights (the default is \code{0.01})
//' 
//' 
//' @return A named list of model parameters that can be passed into portfolio
//'   optimization functions.
//'
//' @details
//'   The function \code{param_portf()} creates a named list of model parameters
//'   that can be passed into portfolio optimization functions.  For example
//'   into the functions \code{calc_weights()} and \code{back_test()}.
//'   See the function \code{calc_weights()} for more details.
//'   
//'   The function \code{param_portf()} simplifies the creation of portfolio
//'   optimization parameter lists.  The users can create a parameter list with
//'   the default values, or they can specify custom parameter values.
//'
//' @examples
//' \dontrun{
//' # Create a default list of portfolio optimization parameters
//' controlv <- HighFreq::param_portf()
//' unlist(controlv)
//' # Create a custom list of portfolio optimization parameters
//' controlv <- HighFreq::param_portf(method="regular", dimax=4)
//' unlist(controlv)
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List param_portf(std::string method = "sharpem",  // Type of portfolio optimization model
                       double singmin = 1e-5, // Threshold level for discarding small singular values
                       arma::uword dimax = 0, // Number of eigen vectors for dimension reduction
                       double confl = 0.1, // Confidence level for calculating the quantiles of returns
                       double alpha = 0.0, // Shrinkage intensity of returns
                       bool rankw = false, // Should the weights be ranked?
                       bool centerw = false, // Should the weights be centered?
                       std::string scalew = "voltarget", // Method for scaling the weights
                       double voltarget = 0.01) { // Volatility target for scaling the weights

  Rcpp::List controlv = Rcpp::List::create(Rcpp::Named("method") = method,
                                           Rcpp::Named("singmin") = singmin,
                                           Rcpp::Named("dimax") = dimax,
                                           Rcpp::Named("confl") = confl,
                                           Rcpp::Named("alpha") = alpha,
                                           Rcpp::Named("rankw") = rankw,
                                           Rcpp::Named("centerw") = centerw,
                                           Rcpp::Named("scalew") = scalew,
                                           Rcpp::Named("voltarget") = voltarget);
  
  return controlv;
  
}  // end param_portf



////////////////////////////////////////////////////////////
//' Apply a lag to a single-column \emph{time series} or a \emph{vector} 
//' using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A single-column \emph{time series} or a
//'   \emph{vector}.
//'
//' @param \code{lagg} An \emph{integer} equal to the number of periods to lag.
//'   (The default is \code{lagg = 1}.)
//'
//' @param \code{pad_zeros} \emph{Boolean} argument: Should the output be padded
//'   with zeros? (The default is \code{pad_zeros = TRUE}.)
//'
//' @return A column \emph{vector} with the same number of elements as the input
//'   time series.
//'
//' @details
//'   The function \code{lag_vec()} applies a lag to the input \emph{time
//'   series} \code{tseries} by shifting its elements by the number equal to the
//'   argument \code{lagg}.  For positive \code{lagg} values, the elements are
//'   shifted forward in time (down), and for negative \code{lagg} values they
//'   are shifted backward (up).
//'   
//'   The output \emph{vector} is padded with either zeros (the default), or
//'   with data from \code{tseries}, so that it has the same number of element
//'   as \code{tseries}.
//'   If the \code{lagg} is positive, then the first element is copied and added
//'   upfront.
//'   If the \code{lagg} is negative, then the last element is copied and added
//'   to the end.
//'   
//'   As a rule, if \code{tseries} contains returns data, then the output
//'   \emph{matrix} should be padded with zeros, to avoid data snooping.
//'   If \code{tseries} contains prices, then the output \emph{matrix} should
//'   be padded with the prices.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random returns
//' retp <- rnorm(1e6)
//' # Compare lag_vec() with rutils::lagit()
//' all.equal(drop(HighFreq::lag_vec(retp)), 
//'   rutils::lagit(retp))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::lag_vec(retp),
//'   Rcode=rutils::lagit(retp),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec lag_vec(const arma::vec& tseries, 
                  arma::sword lagg = 1, 
                  bool pad_zeros = true) {
  
  arma::uword nrows = (tseries.n_elem - 1);
  
  if (lagg > 0) {
    if (pad_zeros) {
      // Pad front with zeros
      return arma::join_cols(arma::zeros<vec>(lagg), 
                             tseries.subvec(0, nrows-lagg));
    } else {
      // Pad front with first element of tseries
      return arma::join_cols(arma::repelem(tseries.subvec(0, 0), lagg, 1), 
                             tseries.subvec(0, nrows-lagg));
    }  // end if
  } else {
    if (pad_zeros) {
      // Pad back with zeros
      return arma::join_cols(tseries.subvec(-lagg, nrows), 
                             arma::zeros<vec>(-lagg));
    } else {
      // Pad back with last element of tseries
      return arma::join_cols(tseries.subvec(-lagg, nrows), 
                             arma::repelem(tseries.subvec(nrows, nrows), -lagg, 1));
    }  // end if
  }  // end if
  
}  // end lag_vec




////////////////////////////////////////////////////////////
//' Apply a lag to the rows of a \emph{time series} or a \emph{matrix} using
//' \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lagg} An \emph{integer} equal to the number of periods to lag
//'   (the default is \code{lagg = 1}).
//'
//' @param \code{pad_zeros} \emph{Boolean} argument: Should the output be padded
//'   with zeros? (The default is \code{pad_zeros = TRUE}.)
//'   
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{lagit()} applies a lag to the input \emph{matrix} by
//'   shifting its rows by the number equal to the argument \code{lagg}. For
//'   positive \code{lagg} values, the rows are shifted \emph{forward} (down),
//'   and for negative \code{lagg} values they are shifted \emph{backward} (up).
//'   
//'   The output \emph{matrix} is padded with either zeros (the default), or
//'   with rows of data from \code{tseries}, so that it has the same dimensions
//'   as \code{tseries}.
//'   If the \code{lagg} is positive, then the first row is copied and added
//'   upfront.
//'   If the \code{lagg} is negative, then the last row is copied and added
//'   to the end.
//'   
//'   As a rule, if \code{tseries} contains returns data, then the output
//'   \emph{matrix} should be padded with zeros, to avoid data snooping.
//'   If \code{tseries} contains prices, then the output \emph{matrix} should
//'   be padded with the prices.
//'
//' @examples
//' \dontrun{
//' # Create a matrix of random returns
//' retp <- matrix(rnorm(5e6), nc=5)
//' # Compare lagit() with rutils::lagit()
//' all.equal(HighFreq::lagit(retp), rutils::lagit(retp))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::lagit(retp),
//'   Rcode=rutils::lagit(retp),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat lagit(const arma::mat& tseries, 
                arma::sword lagg = 1, 
                bool pad_zeros = true) {
  
  arma::uword nrows = (tseries.n_rows - 1);
  arma::uword ncols = tseries.n_cols;
  
  if (lagg > 0) {
    // Positive lag
    if (pad_zeros) {
      // Pad front with zeros
      return arma::join_cols(arma::zeros<mat>(lagg, ncols), 
                             tseries.rows(0, nrows-lagg));
    } else {
      // Pad front with first element of tseries
      return arma::join_cols(arma::repmat(tseries.rows(0, 0), lagg, 1), 
                             tseries.rows(0, nrows-lagg));
    }  // end if
  } else {
    // Negative lag
    if (pad_zeros) {
      // Pad back with zeros
      return arma::join_cols(tseries.rows(-lagg, nrows), 
                             arma::zeros<mat>(-lagg, ncols));
    } else {
      // Pad back with last element of tseries
      return arma::join_cols(tseries.rows(-lagg, nrows), 
                             arma::repmat(tseries.rows(nrows, nrows), -lagg, 1));
    }  // end if
  }  // end if
  
  // Old code below
  // if (lagg > 0)
  //   // Positive lag
  //   return arma::join_cols(arma::repelem(tseries.row(0), lagg, 1), 
  //                          tseries.rows(0, nrows-lagg));
  // else
  //   // Negative lag
  //   return arma::join_cols(tseries.rows(-lagg, nrows), 
  //                          arma::repelem(tseries.row(nrows), -lagg, 1));
  
}  // end lagit




////////////////////////////////////////////////////////////
//' Calculate the differences between the neighboring elements of a
//' single-column \emph{time series} or a \emph{vector}.
//' 
//' @param \code{tseries} A single-column \emph{time series} or a \emph{vector}.
//' 
//' @param \code{lagg} An \emph{integer} equal to the number of time periods to
//'   lag when calculating the differences (the default is \code{lagg = 1}).
//'   
//' @param \code{pad_zeros} \emph{Boolean} argument: Should the output
//'   \emph{vector} be padded (extended) with zeros, in order to return a
//'   \emph{vector} of the same length as the input? (the default is
//'   \code{pad_zeros = TRUE})
//'
//' @return A column \emph{vector} containing the differences between the
//'   elements of the input vector.
//'
//' @details
//'   The function \code{diff_vec()} calculates the differences between the
//'   input \emph{time series} or \emph{vector} and its lagged version.
//'   
//'   The argument \code{lagg} specifies the number of lags.  For example, if
//'   \code{lagg=3} then the differences will be taken between each element
//'   minus the element three time periods before it (in the past).  The default
//'   is \code{lagg = 1}.
//' 
//'   The argument \code{pad_zeros} specifies whether the output \emph{vector}
//'   should be padded (extended) with zeros at the front, in order to
//'   return a \emph{vector} of the same length as the input.  The default is
//'   \code{pad_zeros = TRUE}. The padding operation can be time-consuming,
//'   because it requires the copying the data in memory.
//'   
//'   The function \code{diff_vec()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, which makes it several times faster than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random returns
//' retp <- rnorm(1e6)
//' # Compare diff_vec() with rutils::diffit()
//' all.equal(drop(HighFreq::diff_vec(retp, lagg=3, pad=TRUE)),
//'   rutils::diffit(retp, lagg=3))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::diff_vec(retp, lagg=3, pad=TRUE),
//'   Rcode=rutils::diffit(retp, lagg=3),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec diff_vec(const arma::vec& tseries, arma::uword lagg = 1, bool pad_zeros = true) {
  
  arma::uword length = (tseries.n_elem - 1);
  
  if (pad_zeros)
    // Pad the output with zeros at the front
    return (tseries - arma::join_cols(tseries.subvec(0, lagg - 1), 
                                      tseries.subvec(0, length - lagg)));
  else
    // Don't pad the output
    return (tseries.subvec(lagg, length) - tseries.subvec(0, length - lagg));
  
}  // end diff_vec




////////////////////////////////////////////////////////////
//' Calculate the row differences of a \emph{time series} or a \emph{matrix}
//' using \emph{RcppArmadillo}.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lagg} An \emph{integer} equal to the number of rows (time
//'   periods) to lag when calculating the differences (the default is
//'   \code{lagg = 1}).
//'   
//' @param \code{pad_zeros} \emph{Boolean} argument: Should the output
//'   \emph{matrix} be padded (extended) with zero values, in order to return a
//'   \emph{matrix} with the same number of rows as the input? (the default is
//'   \code{pad_zeros = TRUE})
//'
//' @return A \emph{matrix} containing the differences between the rows of the
//'   input \emph{matrix} \code{tseries}.
//'
//' @details
//'   The function \code{diffit()} calculates the differences between the rows
//'   of the input \emph{matrix} \code{tseries} and its lagged version.
//'   
//'   The argument \code{lagg} specifies the number of lags applied to the rows
//'   of the lagged version of \code{tseries}. 
//'   For positive \code{lagg} values, the lagged version of \code{tseries} has
//'   its rows shifted \emph{forward} (down) by the number equal to \code{lagg}
//'   rows. For negative \code{lagg} values, the lagged version of
//'   \code{tseries} has its rows shifted \emph{backward} (up) by the number
//'   equal to \code{-lagg} rows.
//'   For example, if \code{lagg=3} then the lagged version will have its rows
//'   shifted down by \code{3} rows, and the differences will be taken between
//'   each row minus the row three time periods before it (in the past). The
//'   default is \code{lagg = 1}.
//' 
//'   The argument \code{pad_zeros} specifies whether the output \emph{matrix}
//'   should be padded (extended) with zero values in order to return a
//'   \emph{matrix} with the same number of rows as the input \code{tseries}.
//'   The default is \code{pad_zeros = TRUE}. If \code{pad_zeros = FALSE} then
//'   the return \emph{matrix} has a smaller number of rows than the input
//'   \code{tseries}. The padding operation can be time-consuming, because it
//'   requires the copying the data in memory.
//'   
//'   The function \code{diffit()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, which makes it several times faster than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Create a matrix of random data
//' datav <- matrix(sample(15), nc=3)
//' # Calculate differences with lagged rows
//' HighFreq::diffit(datav, lagg=2)
//' # Calculate differences with advanced rows
//' HighFreq::diffit(datav, lagg=-2)
//' # Compare HighFreq::diffit() with rutils::diffit()
//' all.equal(HighFreq::diffit(datav, lagg=2), 
//'   rutils::diffit(datav, lagg=2), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::diffit(datav, lagg=2),
//'   Rcode=rutils::diffit(datav, lagg=2),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat diffit(const arma::mat& tseries, arma::sword lagg = 1, bool pad_zeros = true) {
  
  arma::uword nrows = (tseries.n_rows-1);
  arma::uword ncols = tseries.n_cols;
  
  if (lagg > 0) {
    // Positive lag
    // Matrix difference without padding
    arma::mat diffmat = (tseries.rows(lagg, nrows) - tseries.rows(0, nrows - lagg));
    if (pad_zeros) {
      // Pad diffmat with zeros at the front
      return arma::join_cols(arma::zeros<mat>(lagg, ncols), diffmat);
    } else {
      // Don't pad the output
      return diffmat;
    }  // end if pad_zeros
  } else {
    // Negative lag
    // Matrix difference without padding
    arma::mat diffmat = (tseries.rows(0, nrows + lagg) - tseries.rows(-lagg, nrows));
    if (pad_zeros) {
      // Pad diffmat with zeros at the back
      return arma::join_cols(diffmat, arma::zeros<mat>(-lagg, ncols));
    } else {
      // Don't pad the output
      return diffmat;
    }  // end if pad_zeros
  }  // end if lagg
  
}  // end diffit




////////////////////////////////////////////////////////////
//' Calculate a vector of end points that divides an integer time sequence of
//' time periods into equal time intervals.
//'
//' @param \code{length} An \emph{integer} equal to the length of the time
//'   sequence to be divided into equal intervals.
//'   
//' @param \code{step} The number of time periods in each interval between
//'   neighboring end points (the default is \code{step = 1}).
//' 
//' @param \code{stub} An \emph{integer} equal to the first non-zero end point
//'   (the default is \code{stub = 0}).
//'
//' @param \code{stubs} A \emph{Boolean} specifying whether to include stub
//'   intervals (the default is \code{stubs = TRUE}).
//'
//' @return A vector of equally spaced \emph{integers} representing the end
//'   points.
//'
//' @details
//'   The end points are a vector of integers which divide the sequence of time
//'   periods of length equal to \code{length} into equally spaced time
//'   intervals.
//'   The number of time periods between neighboring end points is equal to the
//'   argument \code{step}.
//'   If a whole number of intervals doesn't fit over the whole sequence, then
//'   \code{calc_endpoints()} adds a stub interval at the end.
//'   A stub interval is one where the number of periods between neighboring end
//'   points is less than the argument \code{step}.
//'   
//'   If \code{stubs = TRUE} (the default) then the first end point is 
//'   equal to \code{0} (since indexing in \code{C++} code starts at \code{0}).
//'   The first non-zero end point is equal to \code{step} or \code{stub} (if
//'   it's not zero).
//'   If \code{stub = 0} (the default) then the first end point is equal to
//'   \code{0} (even if \code{stubs = FALSE}).
//'   If \code{stubs = TRUE} (the default) then the last end point is always
//'   equal to \code{length-1}.
//'   The argument \code{stub} should be less than the \code{step}: \code{stub <
//'   step}.
//'   
//'   If \code{step = 1} and \code{stub = 0} (the default), then the vector of
//'   end points is simply equal to:
//'   \deqn{
//'     \{ 0, 1, 2, ..., length - 1 \}
//'   }
//'
//'   If \code{stub = 0} (the default) and \code{stubs = TRUE} (the default)
//'   then the vector of end points is equal to:
//'   \deqn{
//'     \{ 0, step, 2*step, ..., length - 1 \}
//'   }
//'   
//'   If \code{stub = 0} (the default) and \code{stubs = FALSE} then the vector
//'   of end points is equal to:
//'   \deqn{
//'     \{ 0, step, 2*step, ..., n*step \}
//'   }
//'   
//'   If \code{stub > 0} and \code{stubs = TRUE} (the default), then the vector
//'   of end points is equal to:
//'   \deqn{
//'     \{ 0, stub, stub + step, ..., length - 1 \}
//'   }
//'   
//'   For example, the end points for \code{length = 20}, divided into intervals
//'   of \code{step = 5} are equal to: \code{0, 5, 10, 15, 19}.
//'   
//'   If \code{stub = 1} then the first non-zero end point is equal to \code{1}
//'   and the end points are equal to: \code{0, 1, 6, 11, 16, 19}.
//'   The stub interval at the beginning is equal to \code{2} (including
//'   \code{0} and \code{1}).
//'   The stub interval at the end is equal to \code{3 = 19 - 16}.
//'   
//'   The end points for \code{length = 21} divided into intervals of length
//'   \code{step = 5}, with \code{stub = 0}, are equal to: \code{0, 5, 10, 15,
//'   20}.
//'   The beginning interval is equal to \code{5}.
//'   The end interval is equal to \code{5 = 20 - 15}.
//'   
//'   If \code{stub = 1} then the first non-zero end point is equal to \code{1}
//'   and the end points are equal to: \code{0, 1, 6, 11, 16, 20}.
//'   The beginning stub interval is equal to \code{2}.
//'   The end stub interval is equal to \code{4 = 20 - 16}.
//'   
//'   The function \code{calc_endpoints()} is similar to the function
//'   \code{rutils::calc_endpoints()} from package
//'   \href{https://github.com/algoquant/rutils}{rutils}.
//'   
//'   But the end points are shifted by \code{-1} compared to \code{R} code
//'   because indexing starts at \code{0} in \code{C++} code, while it starts at
//'   \code{1} in \code{R} code. So if \code{calc_endpoints()} is used in
//'   \code{R} code then \code{1} should be added to it.
//'   
//' @examples
//' # Calculate the end points without a stub interval
//' HighFreq::calc_endpoints(length=20, step=5)
//' # Calculate the end points with a final stub interval
//' HighFreq::calc_endpoints(length=23, step=5)
//' # Calculate the end points with initial and final stub intervals
//' HighFreq::calc_endpoints(length=20, step=5, stub=2)
//'
//' @export
// [[Rcpp::export]]
arma::uvec calc_endpoints(arma::uword length,  // The length of the sequence
                          arma::uword step = 1,  // The number of periods between neighboring end points
                          arma::uword stub = 0,  // The first non-zero end point
                          bool stubs = true) {  // Include stub intervals?

  // Number of initial end points
  arma::uword numpts = length / step + 3;
  // Define the end points
  arma::uvec endd;
  endd.zeros(numpts);
  // Define the last end point
  arma::uword lastp = length - 1;

  // Calculate the initial end points - including extra end points at the end
  if (stub == 0) {
    for (arma::uword it = 0; it < numpts; ++it) {
      endd(it) = it*step;
    }  // end for
  } else if ((stub > 0) & (stubs)) {
    for (arma::uword it = 1; it < numpts; ++it) {
      endd(it) = stub + (it-1)*step;
    }  // end for
  } else {
    for (arma::uword it = 0; it < numpts; ++it) {
      endd(it) = stub + it*step;
    }  // end for
  }  // end if
  // std::cout << "endd = " << arma::conv_to<rowvec>::from(endd) << std::endl;
  
  // arma::uvec endd = arma::regspace<uvec>(stub, step, lastp + step);
  // Find the index of the largest element of endd which is less than lastp
  arma::uword endpp = 0;
  for (arma::uword it = 0; endd(it) < lastp; ++it) {
    endpp++;
  }  // end for
  // std::cout << "endpp = " << endpp << std::endl;
  
  // Trim the initial end points at the end - remove extra end points at the end
  // Subset endd to include the smallest element of endd which is equal to or greater than lastp
  endd = endd.subvec(0, endpp);

  // Set the stub intervals at the end
  if (stubs) {
    // Include stub intervals
    // Set the last end point to lastp - last element of endd
    endd[endpp] = lastp;
  } else {
    // Do not include the last end point - no stub interval at the end
    // Exclude the last element greater than lastp
    if (endd[endpp] > lastp) {
      endd = endd.subvec(0, endpp-1);
    }  // end if
  }  // end if
  
  return endd;
  
  // Old code below
  // if ((stub == 0) & (remainp == 0)) {
  //   std::cout << "(stub == 0) & (remainp == 0)" << std::endl;
  //   // No stub interval
  //   endd = arma::regspace<uvec>(0, step, length);
  //   endd.back() = lastp;
  //   // endd.transform([](arma::uword val) {return (val - 1);});
  //   // endd.front() = 0;
  // } else if ((stub == 0) & (remainp > 0)) {
  //   std::cout << "(stub == 0) & (remainp > 0)" << std::endl;
  //   // Add stub interval at end
  //   endd = arma::regspace<uvec>(0, step, length + step);
  //   // endd.transform([](arma::uword val) {return (val - 1);});
  //   // endd.front() = 0;
  //   endd.back() = lastp;
  // } else if ((stub > 0) & (remainp == 0)) {
  //   std::cout << "(stub > 0) & (remainp == 0)" << std::endl;
  //   // Add initial stub interval equal to stub
  //   arma::uvec endm = arma::regspace<uvec>(stub, step, length);
  //   endm.back() = lastp;
  //   endd.zeros(numpts+1);
  //   endd.subvec(1, numpts) = endm;
  // } else if ((stub > 0) & (remainp > 0)) {
  //   std::cout << "(stub > 0) & (remainp > 0)" << std::endl;
  //   // Add initial stub interval equal to stub
  //   arma::uvec endm = arma::regspace<uvec>(stub, step, lastp + step);
  //   endm.back() = lastp;
  //   endd.zeros(numpts+1);
  //   endd.subvec(1, numpts) = endm;
  // }  // end if
  
}  // end calc_endpoints




////////////////////////////////////////////////////////////
//' Calculate a vector of start points by lagging (shifting) a vector of end
//' points.
//'
//' @param \code{endd} An \emph{integer} vector of end points.
//'   
//' @param \code{lookb} The length of the look-back interval, equal to the
//'   lag (shift) applied to the end points.
//'   
//' @return An \emph{integer} vector with the same number of elements as the
//'   vector \code{endd}.
//'
//' @details
//'   The start points are equal to the values of the vector \code{endd} lagged
//'   (shifted) by an amount equal to \code{lookb}.  In addition, an extra
//'   value of \code{1} is added to them, to avoid data overlaps.  The lag
//'   operation requires appending a beginning warmup interval containing zeros,
//'   so that the vector of start points has the same length as the \code{endd}.
//'   
//'   For example, consider the end points for a vector of length \code{25}
//'   divided into equal intervals of length \code{5}: \code{4, 9, 14, 19, 24}.
//'   (In \code{C++} the vector indexing starts at \code{0} not \code{1}, so
//'   it's shifted by \code{-1}.)
//'   Then the start points for \code{lookb = 2} are equal to: \code{0, 0, 
//'   5, 10, 15}.  The differences between the end points minus the
//'   corresponding start points are equal to \code{9}, except for the warmup
//'   interval.
//'   
//' @examples
//' # Calculate end points
//' endd <- HighFreq::calc_endpoints(length=55, step=5)
//' # Calculate start points corresponding to the end points
//' startp <- HighFreq::calc_startpoints(endd, lookb=5)
//'
//' @export
// [[Rcpp::export]]
arma::uvec calc_startpoints(arma::uvec endd, arma::uword lookb) {
  
  arma::uword numpts = endd.n_elem;
  arma::uvec startp = arma::join_cols(arma::zeros<uvec>(lookb), 
                                      endd.subvec(0, numpts - lookb - 1) + 1);
  
  return startp;
  
}  // end calc_startpoints



////////////////////////////////////////////////////////////
//' Count the number of consecutive \code{TRUE} elements in a Boolean vector,
//' and reset the count to zero after every \code{FALSE} element.
//' 
//' @param \code{tseries} A \emph{Boolean vector} of data.
//'
//' @return An \emph{integer vector} of the same length as the argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{roll_count()} calculates the number of consecutive
//'   \code{TRUE} elements in a Boolean vector, and it resets the count to zero
//'   after every \code{FALSE} element.
//'   
//'   For example, the Boolean vector {\code{FALSE}, \code{TRUE}, \code{TRUE},
//'   \code{FALSE}, \code{FALSE}, \code{TRUE}, \code{TRUE}, \code{TRUE},
//'   \code{TRUE}, \code{TRUE}, \code{FALSE}}, is translated into {\code{0},
//'   \code{1}, \code{2}, \code{0}, \code{0}, \code{1}, \code{2}, \code{3},
//'   \code{4}, \code{5}, \code{0}}.
//'   
//' @examples
//' \dontrun{
//' # Calculate the number of consecutive TRUE elements
//' drop(HighFreq::roll_count(c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)))
//' }
//' @export
// [[Rcpp::export]]
arma::uvec roll_count(const arma::uvec& tseries) {
  
  arma::uword length = tseries.n_elem;
  arma::uvec count_true(length);
  
  // Initialize count
  count_true(0) = tseries(0);
  // Loop over tseries
  for (arma::uword it = 1; it < length; it++) {
    if (tseries(it))
      // Add count number
      count_true(it) = count_true(it-1) + 1;
    else
      // Reset count to zero
      count_true(it) = tseries(it);
  }  // end for
  
  return count_true;
  
}  // end roll_count



////////////////////////////////////////////////////////////
//' Calculate the run length encoding of a single-column \emph{time series},
//' \emph{matrix}, or a \emph{vector}.
//' 
//' @param \code{tseries} A single-column \emph{time series}, \emph{matrix}, or
//'   a \emph{vector}.
//'
//' @return A \emph{list} with two \emph{vectors}: a \emph{vector} of encoded
//'   data and an \emph{integer vector} of data counts (repeats).
//'
//' @details
//'   The function \code{encode_it()} calculates the run length encoding of a
//'   single-column \emph{time series}, \emph{matrix}, or a \emph{vector}.
//'   
//'   The run length encoding of a \emph{vector} consists of two \emph{vectors}:
//'   a \emph{vector} of encoded data (consecutive data values) and of an
//'   \emph{integer vector} of the data counts (the number of times the same
//'   value repeats in succession).
//'   
//'   Run length encoding (RLE) is a data compression algorithm which encodes
//'   the data in two \emph{vectors}: the consecutive data values and their
//'   counts.  If a data value occurs several times in succession then it is
//'   recorded only once and its corresponding count is equal to the number of
//'   times it occurs. Run-length encoding is different from a contingency
//'   table.
//'   
//' @examples
//' \dontrun{
//' # Create a vector of data
//' datav <- sample(5, 31, replace=TRUE)
//' # Calculate the run length encoding of datav
//' HighFreq::encode_it(datav)
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List encode_it(arma::vec tseries) {
  
  // Define vector of encoded data
  std::vector<double> codev;
  codev.reserve(tseries.size());
  // Define vector of data counts (repeats)
  std::vector<int> countv;
  countv.reserve(tseries.size());
  
  // Define iterators
  std::vector<int>::reverse_iterator revit = countv.rbegin();
  
  // Initialize the data
  // Copy the first element of tseries to previous value
  double preval = tseries(0);
  // std::cout << "Initialize the data, preval = " << preval << std::endl;
  // Copy the previous value to the encoded data vector
  codev.push_back(preval);
  // Set the first counter vector value to 0
  countv.push_back(0);
  revit = countv.rbegin();
  // std::cout << "Initialize the data, revit = " << *revit << std::endl;
  
  // Perform loop over tseries
  for (auto inpit: tseries) {
    // std::cout << "for loop, inpit = " << *inpit << std::endl;
    if (preval == inpit) {
      // Data was repeated - increment the current counter vector value by 1
      (*revit)++;
      // std::cout << "Data was repeated, inpit = " << inpit << std::endl;
      // std::cout << "Data was repeated, preval = " << preval << std::endl;
    } else {
      // Data was not repeated
      // std::cout << "Data was not repeated, inpit = " << inpit << std::endl;
      // std::cout << "Data was not repeated, preval = " << preval << std::endl;
      // Copy the tseries value to the previous value
      preval = inpit;
      // Copy the previous value to the encoded data vector
      codev.push_back(preval);
      // Set the next counter vector value to 1
      countv.push_back(1);
      // Reset the counter iterator
      revit = countv.rbegin();
    }  // end if
    
  }  // end for
  
  // Return a list with the encoded data and the counter vector
  return Rcpp::List::create(
    Rcpp::_["data"] = codev,
    Rcpp::_["counts"] = countv
  );
  
}  // end encode_it



////////////////////////////////////////////////////////////
//' Calculate the \emph{vector} of data from its run length encoding.
//' 
//' @param \code{encodel} A \emph{list} with two \emph{vectors}: a \emph{numeric
//'   vector} of encoded data and an \emph{integer vector} of data counts
//'   (repeats).
//'
//' @return A \code{numeric vector}.
//' 
//' @details
//'   The function \code{decode_it()} the \emph{vector} of data from its run
//'   length encoding.
//'   
//'   The run length encoding of a \emph{vector} consists of two \emph{vectors}:
//'   a \emph{numeric vector} of encoded data (consecutive data values) and of
//'   an \emph{integer vector} of the data counts (the number of times the same
//'   value repeats in succession).
//'   
//'   Run length encoding (RLE) is a data compression algorithm which encodes
//'   the data in two \emph{vectors}: the consecutive data values and their
//'   counts.  If a data value occurs several times in succession then it is
//'   recorded only once and its corresponding count is equal to the number of
//'   times it occurs. Run-length encoding is different from a contingency
//'   table.
//'   
//' @examples
//' \dontrun{
//' # Create a vector of data
//' datav <- sample(5, 31, replace=TRUE)
//' # Calculate the run length encoding of datav
//' rle <- HighFreq::encode_it(datav)
//' # Decode the data from its run length encoding
//' decodev <- HighFreq::decode_it(rle)
//' all.equal(datav, decodev)
//' }
//' 
//' @export
// [[Rcpp::export]]
std::vector<double> decode_it(Rcpp::List encodel) {
  
  // Extract vector of encoded data
  std::vector<double> codev = encodel["data"];
  // Extract vector of data counts (repeats)
  std::vector<int> countv = encodel["counts"];
  // Define the output decoded vector
  std::vector<double> decodev;
  decodev.reserve(std::accumulate(countv.begin(), countv.end(), 0));
  
  // Perform loop over the codev and countv vectors
  for (std::size_t it = 0; it < codev.size(); it++) {
    for (int j = 0; j < countv[it]; j++) {
      decodev.push_back(codev[it]);
    }  // end for
  }  // end for
  
  return decodev;
  
}  // end decode_it



////////////////////////////////////////////////////////////
//' Calculate the ranks of the elements of a single-column \emph{time series},
//' \emph{matrix}, or a \emph{vector} using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A single-column \emph{time series}, \emph{matrix}, or
//'   a \emph{vector}.
//'
//' @return An \emph{integer vector} with the ranks of the elements of the
//'   \code{tseries}.
//'
//' @details
//'   The function \code{calc_ranks()} calculates the ranks of the elements of a
//'   single-column \emph{time series}, \emph{matrix}, or a \emph{vector}. 
//'   
//'   The permutation index is an integer vector which sorts a given vector into
//'   ascending order. 
//'   The permutation index of the permutation index is the \emph{reverse}
//'   permutation index, because it sorts the vector from ascending order back
//'   into its original unsorted order.
//'   The ranks of the elements are equal to the \emph{reverse} permutation
//'   index. The function \code{calc_ranks()} calculates the \emph{reverse}
//'   permutation index.
//'   
//'   The ranks produced by \code{calc_ranks()} start at zero, following the 
//'   \code{C++} convention.
//'   
//'   The \code{Armadillo} function \code{arma::sort_index()} calculates the
//'   permutation index which sorts a given vector into an ascending order.
//'   Applying the function \code{arma::sort_index()} twice:\cr
//'   \code{arma::sort_index(arma::sort_index())},\cr
//'   calculates the \emph{reverse} permutation index to sort the vector from
//'   ascending order back into its original unsorted order.
//'   
//'   The function \code{calc_ranks()} calls the \code{Armadillo} function
//'   \code{arma::sort_index()} twice to calculate the \emph{reverse}
//'   permutation index, to sort the vector from ascending order back into its
//'   original unsorted order.
//'   
//' @examples
//' \dontrun{
//' # Create a vector of data
//' datav <- rnorm(1e3)
//' # Calculate the ranks of the elements using R code and RcppArmadillo
//' all.equal(rank(datav), drop(HighFreq::calc_ranks(datav))+1)
//' # Compare the speed of R code with RcppArmadillo
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcode=rank(datav),
//'   Rcpp=calc_ranks(datav),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::uvec calc_ranks(arma::vec tseries) {
  
  return (arma::sort_index(arma::sort_index(tseries)));
  
}  // end calc_ranks



////////////////////////////////////////////////////////////
//' Calculate the ranks of the elements of a single-column \emph{time series},
//' \emph{matrix}, or a \emph{vector} using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A single-column \emph{time series}, \emph{matrix}, or
//'   a \emph{vector}.
//'
//' @return An \emph{integer vector} with the ranks of the elements of
//'   \code{tseries}.
//'
//' @details
//'   The function \code{calc_ranks_stl()} calculates the ranks of the elements
//'   of a single-column \emph{time series}, \emph{matrix}, or a \emph{vector}.
//'   The function \code{calc_ranks_stl()} is slightly faster than the function
//'   \code{calc_ranks()}.
//'   
//'   The permutation index is an integer vector which sorts a given vector into
//'   ascending order. 
//'   The permutation index of the permutation index is the \emph{reverse}
//'   permutation index, because it sorts the vector from ascending order back
//'   into its original unsorted order.
//'   The ranks of the elements are equal to the \emph{reverse} permutation
//'   index. The function \code{calc_ranks()} calculates the \emph{reverse}
//'   permutation index.
//'   
//'   The ranks produced by \code{calc_ranks_stl()} start at zero, following the
//'   \code{C++} convention.
//'
//'   The \code{STL} \code{C++} function \code{std::sort()} sorts a vector into
//'   ascending order. It can also be used to calculate the permutation index
//'   which sorts the vector into an ascending order.
//'   
//'   The function \code{calc_ranks_stl()} calls the function \code{std::sort()}
//'   twice:
//'   First, it calculates the permutation index which sorts the vector
//'   \code{tseries} into ascending order.
//'   Second, it calculates the permutation index of the permutation index,
//'   which are the ranks (the \emph{reverse} permutation index) of the vector
//'   \code{tseries}.
//' 
//' @examples
//' \dontrun{
//' # Create a vector of data
//' datav <- rnorm(1e3)
//' # Calculate the ranks of the elements using R code and RcppArmadillo
//' all.equal(rank(datav), drop(HighFreq::calc_ranks_stl(datav))+1)
//' # Compare the speed of R code with RcppArmadillo
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcode=rank(datav),
//'   Rcpp=HighFreq::calc_ranks_stl(datav),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::uvec calc_ranks_stl(arma::vec tseries) {
  
  size_t ndata = tseries.size();
  // size_t ndata = sizeof(tseries);
  
  // Define index of integers along tseries
  arma::uvec indeks(ndata);
  // Define the ranks of the vector elements
  arma::uvec ranks(ndata);
  // Fill the vectors with a sequence of consecutive integers
  // The function std::iota() is simiar to the R function seq_along()
  std::iota(indeks.begin(), indeks.end(), 0);
  std::iota(ranks.begin(), ranks.end(), 0);
  
  // Calculate the permutation index by sorting the sequence of consecutive integers 
  // according to the order of tseries.
  std::sort(indeks.begin(), indeks.end(), 
            // Lambda comparison function defines sorting order.
            // The brackets [] are used to pass in variables from the outer scope of the lambda function.
            // The "&" passes the outer scope variables by reference.
            [&tseries](int i1, int i2) {return tseries[i1] < tseries[i2];});
  
  // Calculate the ranks (inverse permutation index) by sorting the sequence of consecutive integers 
  // according to the order of the permutation index indeks.
  std::sort(ranks.begin(), ranks.end(), 
            // Lambda comparison function defines sorting order.
            // The brackets [] are used to pass in variables from the outer scope of the lambda function.
            // The "&" passes the outer scope variables by reference.
            [&indeks](int i1, int i2) {return indeks[i1] < indeks[i2];});
  
  return ranks;
  
}  // end calc_ranks_stl



// The function remove_dup() removes consecutive duplicate elements 
// from the input vector of strings.
// It doesn't remove all duplicate elements.  
// It doesn't remove duplicate elements which don't neighbor each other.
// It uses the STL algorithm std::unique().
// [[Rcpp::export]]
std::vector<std::string> remove_dup(std::vector<std::string> stringv) {
  
  std::size_t ndata = stringv.size();
  // Define vector iterator
  std::vector<std::string>::iterator stringit;
  // Define vector of output strings
  std::vector<std::string> outv = stringv;
  outv.reserve(ndata);
  
  // Remove consecutive duplicate elements
  stringit = std::unique(outv.begin(), outv.end());
  // Resize the output vector
  outv.resize(std::distance(outv.begin(), stringit));
  
  return outv;
  
}  // end remove_dup




////////////////////////////////////////////////////////////
// Functions for matrix algebra
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
//' Multiply element-wise the rows or columns of a \emph{matrix} times a
//' \emph{vector}.
//' 
//' @param \code{vector} A \emph{numeric} \emph{vector}.
//' 
//' @param \code{matrix} A \emph{numeric} \emph{matrix}.
//' 
//' @param \code{byrow} A \emph{Boolean} argument: if \code{TRUE} then multiply
//'   the rows of \code{matrix} by \code{vector}, otherwise multiply the columns
//'   (the default is \code{byrow = TRUE}.)
//' 
//' @return A \emph{matrix} equal to the product of the elements of
//'   \code{matrix} times \code{vector}, with the same dimensions as the
//'   argument \code{matrix}.
//' 
//' @details
//'   The function \code{mult_mat()} multiplies element-wise the rows or columns
//'   of a \emph{matrix} times a \emph{vector}.
//'
//'   If \code{byrow = TRUE} (the default), then function \code{mult_mat()}
//'   multiplies the rows of the argument \code{matrix} times the argument
//'   \code{vector}.
//'   Otherwise it multiplies the columns of \code{matrix}.
//' 
//'   In \code{R}, \emph{matrix} multiplication is performed by columns.
//'   Performing multiplication by rows is often required, for example when
//'   multiplying asset returns by portfolio weights.
//'   But performing multiplication by rows requires explicit loops in \code{R},
//'   or it requires \emph{matrix} transpose.  And both are slow.
//'
//'   The function \code{mult_mat()} uses \code{RcppArmadillo} \code{C++}
//'   code, so when multiplying large \emph{matrix} columns it's several times
//'   faster than vectorized \code{R} code, and it's even much faster compared
//'   to \code{R} when multiplying the \emph{matrix} rows.
//' 
//'   The function \code{mult_mat()} performs loops over the \emph{matrix} rows
//'   and columns using the \code{Armadillo} operators \code{each_row()} and
//'   \code{each_col()}, instead of performing explicit \code{for()} loops (both
//'   methods are equally fast).
//'   
//' @examples
//' \dontrun{
//' # Create vector and matrix data
//' matrixv <- matrix(round(runif(25e4), 2), nc=5e2)
//' vectorv <- round(runif(5e2), 2)
//' 
//' # Multiply the matrix rows using R
//' matrixr <- t(vectorv*t(matrixv))
//' # Multiply the matrix rows using C++
//' matrixp <- HighFreq::mult_mat(vectorv, matrixv, byrow=TRUE)
//' all.equal(matrixr, matrixp)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_mat(vectorv, matrixv, byrow=TRUE),
//'     Rcode=t(vectorv*t(matrixv)),
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//'     
//' # Multiply the matrix columns using R
//' matrixr <- vectorv*matrixv
//' # Multiply the matrix columns using C++
//' matrixp <- HighFreq::mult_mat(vectorv, matrixv, byrow=FALSE)
//' all.equal(matrixr, matrixp)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_mat(vectorv, matrixv, byrow=FALSE),
//'     Rcode=vectorv*matrixv,
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat mult_mat(arma::vec vectorv,
                   arma::mat matrixv,
                   bool byrow = true) {
  
  arma::uword ndata = vectorv.n_elem;
  arma::uword nrows = matrixv.n_rows;
  arma::uword ncols = matrixv.n_cols;
  
  if (byrow && (ndata == ncols)) {
    // Multiply every row of matrixv by vectorv
    matrixv = matrixv.each_row() % vectorv.t();
  } else if (!byrow && (ndata == nrows)) {
    // Multiply every column of matrixv by vectorv
    matrixv = matrixv.each_col() % vectorv;
  } else {
    // Do nothing
    cout << "Nothing done: Vector length is neither equal to the number of columns nor to the rows of the matrix!" << std::endl;
  }  // end if
  
  return matrixv;
  
}  // end mult_mat



////////////////////////////////////////////////////////////
//' Multiply the rows or columns of a \emph{matrix} times a \emph{vector},
//' element-wise and in place, without copying the data in memory.
//' 
//' @param \code{vectorv} A \emph{numeric} \emph{vector}.
//' 
//' @param \code{matrixv} A \emph{numeric} \emph{matrix}.
//' 
//' @param \code{byrow} A \emph{Boolean} argument: if \code{TRUE} then multiply
//'   the rows of \code{matrixv} by \code{vectorv}, otherwise multiply the columns
//'   (the default is \code{byrow = TRUE}.)
//' 
//' @return Void (no return value - modifies the data in place).
//' 
//' @details
//'   The function \code{mult_mat_ref()} multiplies the rows or columns of a
//'   \emph{matrix} times a \emph{vector}, element-wise and in place, without
//'   copying the data in memory.
//'
//'   It accepts a \emph{pointer} to the argument \code{matrixv}, and it
//'   overwrites the old \code{matrix} values with the new values. It performs
//'   the calculation in place, without copying the \emph{matrix} in memory,
//'   which can significantly increase the computation speed for large matrices.
//'
//'   If \code{byrow = TRUE} (the default), then function \code{mult_mat_ref()}
//'   multiplies the rows of the argument \code{matrixv} times the argument
//'   \code{vectorv}.
//'   Otherwise it multiplies the columns of \code{matrixv}.
//' 
//'   In \code{R}, \emph{matrix} multiplication is performed by columns.
//'   Performing multiplication by rows is often required, for example when
//'   multiplying asset returns by portfolio weights.
//'   But performing multiplication by rows requires explicit loops in \code{R},
//'   or it requires \emph{matrix} transpose.  And both are slow.
//'
//'   The function \code{mult_mat_ref()} uses \code{RcppArmadillo} \code{C++}
//'   code, so when multiplying large \emph{matrix} columns it's several times
//'   faster than vectorized \code{R} code, and it's even much faster compared
//'   to \code{R} when multiplying the \emph{matrix} rows.
//' 
//'   The function \code{mult_mat_ref()} performs loops over the \emph{matrix}
//'   rows and columns using the \code{Armadillo} operators \code{each_row()}
//'   and \code{each_col()}, instead of performing explicit \code{for()} loops
//'   (both methods are equally fast).
//'   
//' @examples
//' \dontrun{
//' # Create vector and matrix data
//' matrixv <- matrix(round(runif(25e4), 2), nc=5e2)
//' vectorv <- round(runif(5e2), 2)
//' 
//' # Multiply the matrix rows using R
//' matrixr <- t(vectorv*t(matrixv))
//' # Multiply the matrix rows using C++
//' HighFreq::mult_mat_ref(vectorv, matrixv, byrow=TRUE)
//' all.equal(matrixr, matrixv)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_mat_ref(vectorv, matrixv, byrow=TRUE),
//'     Rcode=t(vectorv*t(matrixv)),
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//'     
//' # Multiply the matrix columns using R
//' matrixr <- vectorv*matrixv
//' # Multiply the matrix columns using C++
//' HighFreq::mult_mat_ref(vectorv, matrixv, byrow=FALSE)
//' all.equal(matrixr, matrixv)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_mat_ref(vectorv, matrixv, byrow=FALSE),
//'     Rcode=vectorv*matrixv,
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
void mult_mat_ref(arma::vec vectorv,
                  arma::mat matrixv,
                  bool byrow = true) {
  
  arma::uword ndata = vectorv.n_elem;
  arma::uword nrows = matrixv.n_rows;
  arma::uword ncols = matrixv.n_cols;
  
  if (byrow && (ndata == ncols)) {
    // Multiply every row of matrixv by vectorv
    matrixv.each_row() %= vectorv.t();
  } else if (!byrow && (ndata == nrows)) {
    // Multiply every column of matrixv by vectorv
    matrixv.each_col() %= vectorv;
  } else {
    // Do nothing
    cout << "Nothing done: Vector length is neither equal to the number of columns nor to the rows of the matrix!" << std::endl;
  }  // end if
  
  // return matrixv;
  
}  // end mult_mat_ref




////////////////////////////////////////////////////////////
//' Calculate the eigen decomposition of a square, symmetric matrix using
//' \code{RcppArmadillo}.
//' 
//' @param \code{matrixv} A square, symmetric matrix.
//'
//' @param \code{eigenval} A \emph{vector} of eigen values.
//' 
//' @param \code{eigenvec} A \emph{matrix} of eigen vectors.
//' 
//' @return Void (no return value - passes the eigen values and eigen vectors by
//'   reference).
//'
//' @details
//'   The function \code{calc_eigen()} calculates the eigen decomposition of a
//'   square, symmetric matrix using \code{RcppArmadillo}.  It calls the
//'   \code{Armadillo} function \code{arma::eig_sym()} to calculate the eigen
//'   decomposition. 
//'   
//'   For small matrices, the function \code{calc_eigen()} is several times
//'   faster than the \code{R} function \code{eigen()}, since
//'   \code{calc_eigen()} has no overhead in \code{R} code. But for large
//'   matrices, they are about the same, since both call \code{C++} code.
//'
//' @examples
//' \dontrun{
//' # Create random positive semi-definite matrix
//' matrixv <- matrix(runif(25), nc=5)
//' matrixv <- t(matrixv) %*% matrixv
//' # Calculate the eigen decomposition using RcppArmadillo
//' eigenval <- numeric(5) # Allocate eigen values
//' eigenvec <- matrix(numeric(25), nc=5) # Allocate eigen vectors
//' HighFreq::calc_eigen(matrixv, eigenval, eigenvec)
//' # Calculate the eigen decomposition using R
//' eigenr <- eigen(matrixv)
//' # Compare the eigen decompositions
//' all.equal(eigenr$values, drop(eigenval))
//' all.equal(abs(eigenr$vectors), abs(eigenvec))
//' # Compare the speed of Rcpp with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_eigen(matrixv, eigenval, eigenvec),
//'   Rcode=eigen(matrixv),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
void calc_eigen(const arma::mat& matrixv,
                arma::vec& eigenval, // Eigen values
                arma::mat& eigenvec) { // Eigen vectors
  
  // arma::mat eigenvec;
  // arma::vec eigenval;
  arma::eig_sym(eigenval, eigenvec, matrixv);
  // Reverse the order of elements from largest eigenvalue to smallest, similar to R
  eigenval = arma::flipud(eigenval);
  eigenvec = arma::fliplr(eigenvec);
  
  // Reverse the order of elements from largest eigenvalue to smallest, similar to R
  // return Rcpp::List::create(Rcpp::Named("values") = arma::flipud(eigenval),
  //                           Rcpp::Named("vectors") = arma::fliplr(eigenvec));
  
}  // end calc_eigen



////////////////////////////////////////////////////////////
//' Calculate the partial eigen decomposition of a dense symmetric matrix using
//' \code{RcppArmadillo}.
//' 
//' @param \code{matrixv} A square matrix.
//'
//' @param \code{neigen} An \emph{integer} equal to the number of eigenvalues
//'   to be calculated.
//'
//' @return A list with two elements: a \emph{vector} of eigenvalues (named
//'   "values"), and a \emph{matrix} of eigenvectors (named "vectors").
//'
//' @details
//'   The function \code{calc_eigenp()} calculates the partial eigen
//'   decomposition (the lowest order principal components, with the largest
//'   eigenvalues) of a dense matrix using RcppArmadillo.  It calls the internal
//'   \code{Armadillo} eigen solver \code{SymEigsSolver} in the namespace
//'   \code{arma::newarp} to calculate the partial eigen decomposition.
//'   
//'   The eigen solver \code{SymEigsSolver} uses the Implicitly Restarted
//'   Lanczos Method (IRLM) which was adapted from the
//'   \href{https://en.wikipedia.org/wiki/ARPACK}{ARPACK} library. The eigen
//'   solver \code{SymEigsSolver} was implemented by
//'   \href{https://github.com/yixuan/arpack-arma}{Yixuan Qiu}.
//'   
//'   The function \code{arma::eigs_sym()} also calculates the partial eigen
//'   decomposition using the eigen solver \code{SymEigsSolver}, but it only
//'   works for sparse matrices which are not standard R matrices.
//'   
//'   For matrices smaller than \code{100} rows, the function
//'   \code{calc_eigenp()} is slower than the function \code{calc_eigen()} which
//'   calculates the full eigen decomposition.  But it's faster for very large
//'   matrices.
//'
//' @examples
//' \dontrun{
//' # Create random positive semi-definite matrix
//' matrixv <- matrix(runif(100), nc=10)
//' matrixv <- t(matrixv) %*% matrixv
//' # Calculate the partial eigen decomposition
//' neigen <- 5
//' eigenp <- HighFreq::calc_eigenp(matrixv, neigen)
//' # Calculate the eigen decomposition using RcppArmadillo
//' eigenval <- numeric(10) # Allocate eigen values
//' eigenvec <- matrix(numeric(100), nc=10) # Allocate eigen vectors
//' HighFreq::calc_eigen(matrixv, eigenval, eigenvec)
//' # Compare the eigen decompositions
//' all.equal(eigenp$values[1:neigen], eigenval[1:neigen])
//' all.equal(abs(eigenp$vectors), abs(eigenvec[, 1:neigen]))
//' # Compare the speed of partial versus full decomposition
//' summary(microbenchmark(
//'   partial=HighFreq::calc_eigenp(matrixv, neigen),
//'   full=HighFreq::calc_eigen(matrixv, eigenval, eigenvec),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List calc_eigenp(arma::mat& matrixv, const arma::uword& neigen) {
  
  arma::uword nrows = matrixv.n_rows;
  arma::mat eigenvec;
  arma::vec eigenval;
  
  // Construct matrix operation object using the wrapper class DenseGenMatProd
  DenseGenMatProd<double> matop(matrixv);
  
  // Construct eigen solver object, requesting the largest three eigenvalues
  SymEigsSolver <double, EigsSelect::LARGEST_ALGE, DenseGenMatProd<double>> eigs(matop, neigen, nrows);
  
  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute();
  
  // Retrieve results
  if (nconv > 0) {
    eigenval = eigs.eigenvalues();
    eigenvec = eigs.eigenvectors();
  }  // end if
  
  // Reverse the order of elements from largest eigenvalue to smallest, similar to R
  return Rcpp::List::create(Rcpp::Named("values") = arma::flipud(eigenval),
                            Rcpp::Named("vectors") = arma::fliplr(eigenvec));
  
}  // end calc_eigenp




////////////////////////////////////////////////////////////
//' Calculate the \emph{reduced inverse} of a symmetric \emph{matrix} of
//' data using eigen decomposition.
//' 
//' @param \code{matrixv} A symmetric \emph{matrix} of data.
//'   
//' @param \code{dimax} An \emph{integer} equal to the number of \emph{eigen
//'   values} used for calculating the \emph{reduced inverse} of the matrix
//'   \code{matrixv} (the default is \code{dimax = 0} - standard matrix inverse
//'   using all the \emph{eigen values}).
//' 
//' @param \code{singmin} A \emph{numeric} threshold level for discarding
//'   small \emph{eigen values} in order to regularize the inverse of the matrix
//'   \code{matrixv} (the default is \code{0.0}).
//'
//' @return A \emph{matrix} equal to the \emph{reduced inverse} of the
//'   matrix \code{matrixv}.
//'
//' @details
//'   The function \code{calc_inv()} calculates the \emph{reduced inverse}
//'   of the matrix \code{matrixv} using eigen decomposition.
//'   
//'   The function \code{calc_inv()} first performs eigen decomposition of the
//'   matrix \code{matrixv}.
//'   The eigen decomposition of a matrix \eqn{\strong{C}} is defined as the
//'   factorization:
//'   \deqn{
//'     \strong{C} = \strong{O} \, \Sigma \, \strong{O}^T
//'   } Where \eqn{\strong{O}} is the matrix of \emph{eigen vectors} and
//'   \eqn{\Sigma} is a diagonal matrix of \emph{eigen values}.
//'   
//'   The inverse \eqn{\strong{C}^{-1}} of the matrix \eqn{\strong{C}} can be
//'   calculated from the eigen decomposition as:
//'   \deqn{
//'     \strong{C}^{-1} = \strong{O} \, \Sigma^{-1} \, \strong{O}^T
//'   }
//'   
//'   The \emph{reduced inverse} of the matrix \eqn{\strong{C}} is obtained
//'   by removing \emph{eigen vectors} with very small \emph{eigen values}:
//'   \deqn{
//'     \strong{C}^{-1} = \strong{O}_{dimax} \, \Sigma^{-1}_{dimax} \, \strong{O}^T_{dimax}
//'   }
//'   Where \eqn{\strong{O}_{dimax}} is the matrix of \emph{eigen vectors} 
//'   that correspond to the largest \emph{eigen values} \eqn{\Sigma_{dimax}}. 
//'   
//'   The function \code{calc_inv()} applies regularization to the matrix
//'   inverse using the arguments \code{dimax} and \code{singmin}.
//'   
//'   The function \code{calc_inv()} applies regularization by discarding the
//'   smallest \emph{eigen values} \eqn{\Sigma_i} that are less than the
//'   threshold level \code{singmin} times the sum of all the \emph{eigen
//'   values}: \deqn{\Sigma_i < eigen\_thresh \cdot (\sum{\Sigma_i})}
//'   
//'   It also discards additional \emph{eigen vectors} so that only the highest
//'   order \emph{eigen vectors} remain, up to order \code{dimax}.
//'   It calculates the \emph{reduced inverse} from the eigen decomposition
//'   using only the largest \emph{eigen values} up to \code{dimax}.  For
//'   example, if
//'   \code{dimax = 3} then it only uses the \code{3} highest order \emph{eigen
//'   vectors}, with the largest \emph{eigen values}. This has the effect of
//'   dimension reduction.
//'   
//'   If the matrix \code{matrixv} has a large number of small \emph{eigen
//'   values}, then the number of remaining \emph{eigen values} may be less than
//'   \code{dimax}.
//'   
//' @examples
//' \dontrun{
//' # Calculate ETF returns
//' retp <- na.omit(rutils::etfenv$returns[, c("VTI", "TLT", "DBC")])
//' # Calculate covariance matrix
//' covmat <- cov(retp)
//' # Calculate matrix inverse using RcppArmadillo
//' invmat <- HighFreq::calc_inv(covmat)
//' # Calculate matrix inverse in R
//' invr <- solve(covmat)
//' all.equal(invmat, invr, check.attributes=FALSE)
//' # Calculate reduced inverse using RcppArmadillo
//' invmat <- HighFreq::calc_inv(covmat, dimax=3)
//' # Calculate reduced inverse using eigen decomposition in R
//' eigend <- eigen(covmat)
//' dimax <- 1:3
//' invr <- eigend$vectors[, dimax] %*% (t(eigend$vectors[, dimax])/eigend$values[dimax])
//' # Compare RcppArmadillo with R
//' all.equal(invmat, invr)
//' }
//' 
//' @examples
// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& matrixv, 
                   arma::uword dimax = 0, // Number of eigen vectors for dimension reduction
                   double singmin = 0.0) { // Threshold for discarding small eigen values
  
  // Allocate Eigen variables
  arma::uword ncols = matrixv.n_cols;
  arma::vec eigenval(ncols, fill::zeros); // Eigen values
  arma::mat eigenvec(ncols, ncols, fill::zeros); // Eigen vectors
  // Calculate the eigen decomposition - the eigenvalues are in ascending order
  arma::eig_sym(eigenval, eigenvec, matrixv);
  // Calculate the number of non-small singular values
  arma::uword neigen = arma::sum(eigenval > singmin*arma::sum(eigenval));
  
  // If no regularization then set dimax to neigen
  if (dimax == 0) {
    // Set dimax
    dimax = neigen;
  } else {
    // Adjust dimax
    dimax = std::min(dimax, neigen);
  }  // end if
  
  // Remove all small singular values
  eigenval = eigenval.subvec(ncols-dimax, ncols-1);
  eigenvec = eigenvec.cols(ncols-dimax, ncols-1);
  
  // Calculate the reduced inverse from the eigen decomposition
  return eigenvec*arma::diagmat(1/eigenval)*eigenvec.t();
  
}  // end calc_inv



////////////////////////////////////////////////////////////
//' Calculate the \emph{reduced inverse} of a \emph{matrix} of data using
//' Singular Value Decomposition (\emph{SVD}).
//' 
//' @param \code{matrixv} A \emph{matrix} of data.
//'   
//' @param \code{dimax} An \emph{integer} equal to the number of \emph{singular
//'   values} used for calculating the \emph{reduced inverse} of the matrix
//'   \code{matrixv} (the default is \code{dimax = 0} - standard matrix inverse
//'   using all the \emph{singular values}).
//' 
//' @param \code{singmin} A \emph{numeric} threshold level for discarding
//'   small \emph{singular values} in order to regularize the inverse of the
//'   matrix \code{matrixv} (the default is \code{0.0}).
//'
//' @return A \emph{matrix} equal to the \emph{reduced inverse} of the
//'   matrix \code{matrixv}.
//'
//' @details
//'   The function \code{calc_invsvd()} calculates the \emph{reduced
//'   inverse} of the matrix \code{matrixv} using Singular Value Decomposition
//'   (\emph{SVD}).
//'   
//'   The function \code{calc_invsvd()} first performs Singular Value
//'   Decomposition (\emph{SVD}) of the matrix \code{matrixv}.
//'   The \emph{SVD} of a matrix \eqn{\strong{C}} is defined as the
//'   factorization:
//'   \deqn{
//'     \strong{C} = \strong{U} \, \Sigma \, \strong{V}^T
//'   }
//'   Where \eqn{\strong{U}} and \eqn{\strong{V}} are the left and right
//'   \emph{singular matrices}, and \eqn{\Sigma} is a diagonal matrix of
//'   \emph{singular values}.
//'   
//'   The inverse \eqn{\strong{C}^{-1}} of the matrix \eqn{\strong{C}} can be
//'   calculated from the \emph{SVD} matrices as:
//'   \deqn{
//'     \strong{C}^{-1} = \strong{V} \, \Sigma^{-1} \, \strong{U}^T
//'   }
//'   
//'   The \emph{reduced inverse} of the matrix \eqn{\strong{C}} is obtained
//'   by removing \emph{singular vectors} with very small \emph{singular values}:
//'   \deqn{
//'     \strong{C}^{-1} = \strong{V}_n \, \Sigma_n^{-1} \, \strong{U}_n^T
//'   }
//'   Where \eqn{\strong{U}_n}, \eqn{\strong{V}_n} and \eqn{\Sigma_n} are the
//'   \emph{SVD} matrices with the rows and columns corresponding to very small
//'   \emph{singular values} removed.
//'   
//'   The function \code{calc_invsvd()} applies regularization to the matrix
//'   inverse using the arguments \code{dimax} and \code{singmin}.
//'   
//'   The function \code{calc_invsvd()} applies regularization by discarding the
//'   smallest \emph{singular values} \eqn{\sigma_i} that are less than the
//'   threshold level \code{singmin} times the sum of all the
//'   \emph{singular values}: \deqn{\sigma_i < eigen\_thresh \cdot
//'   (\sum{\sigma_i})}
//'   
//'   It also discards additional \emph{singular vectors} so that only the
//'   highest order \emph{singular vectors} remain, up to order \code{dimax}. It
//'   calculates the \emph{reduced inverse} from the \emph{SVD} matrices
//'   using only the largest \emph{singular values} up to order \code{dimax}.
//'   For example, if \code{dimax = 3} then it only uses the \code{3} highest
//'   order \emph{singular vectors}, with the largest \emph{singular values}.
//'   This has the effect of dimension reduction.
//'   
//'   If the matrix \code{matrixv} has a large number of small \emph{singular
//'   values}, then the number of remaining \emph{singular values} may be less
//'   than \code{dimax}.
//'   
//' @examples
//' \dontrun{
//' # Calculate ETF returns
//' retp <- na.omit(rutils::etfenv$returns[, c("VTI", "TLT", "DBC")])
//' # Calculate covariance matrix
//' covmat <- cov(retp)
//' # Calculate matrix inverse using RcppArmadillo
//' invmat <- HighFreq::calc_invsvd(covmat)
//' # Calculate matrix inverse in R
//' invr <- solve(covmat)
//' all.equal(invmat, invr, check.attributes=FALSE)
//' # Calculate reduced inverse using RcppArmadillo
//' invmat <- HighFreq::calc_invsvd(covmat, dimax=3)
//' # Calculate reduced inverse from SVD in R
//' svdec <- svd(covmat)
//' dimax <- 1:3
//' invr <- svdec$v[, dimax] %*% (t(svdec$u[, dimax])/svdec$d[dimax])
//' # Compare RcppArmadillo with R
//' all.equal(invmat, invr)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_invsvd(const arma::mat& matrixv, 
                      arma::uword dimax = 0, // The number of singular vectors for dimension reduction
                      double singmin = 0.0) { // Threshold for discarding small singular values
  
  // Allocate SVD variables
  arma::vec svdval;  // Singular values
  arma::mat svdu, svdv;  // Singular matrices
  // Calculate the SVD - the singular values are in ascending order
  arma::svd(svdu, svdval, svdv, matrixv);
  // Calculate the number of non-small singular values
  arma::uword svdnum = arma::sum(svdval > singmin*arma::sum(svdval));
  
  // If no regularization then set dimax to (svdnum - 1)
  if (dimax == 0) {
    // Set dimax
    dimax = svdnum - 1;
  } else {
    // Adjust dimax
    dimax = std::min(dimax - 1, svdnum - 1);
  }  // end if
  
  // Remove all small singular values
  svdval = svdval.subvec(0, dimax);
  svdu = svdu.cols(0, dimax);
  svdv = svdv.cols(0, dimax);
  
  // Calculate the reduced inverse from the SVD decomposition
  return svdv*arma::diagmat(1/svdval)*svdu.t();
  
}  // end calc_invsvd



////////////////////////////////////////////////////////////
//' Calculate the approximate inverse of a square \emph{matrix} recursively
//' using the Schulz formula (without copying the data in memory).
//' 
//' @param \code{matrixv} A \emph{matrix} of data to be inverted.
//' 
//' @param \code{invmat} A \emph{matrix} of data equal to the starting point for
//'   the recursion.
//' 
//' @param \code{niter} An \emph{integer} equal to the number of recursion
//'   iterations.
//' 
//' @return No return value.
//'
//' @details
//'   The function \code{calc_invrec()} calculates the approximate inverse
//'   \eqn{\strong{A}^{-1}} of a square \emph{matrix} \eqn{\strong{A}}
//'   recursively using the Schulz formula:
//'   \deqn{
//'     \strong{A}_{i+1}^{-1} \leftarrow 2 \strong{A}_i^{-1} - \strong{A}_i^{-1} \strong{A} \strong{A}_i^{-1}
//'   }
//'   The Schulz formula is repeated \code{niter} times.
//'   The Schulz formula is useful for updating the inverse when the matrix
//'   \eqn{\strong{A}} changes only slightly.  For example, for updating the
//'   inverse of the covariance matrix as it changes slowly over time.
//'
//'   The function \code{calc_invrec()} accepts a \emph{pointer} to the argument
//'   \code{invmat} (which is the initial value of the inverse matrix for the
//'   recursion), and it overwrites the old inverse matrix values with the
//'   updated inverse values.
//'
//'   The function \code{calc_invrec()} performs the calculation in place,
//'   without copying the matrix in memory, which can significantly increase the
//'   computation speed for large matrices.
//'   
//'   The function \code{calc_invrec()} doesn't return a value.
//'   The function \code{calc_invrec()} performs the calculations using
//'   \code{C++} \code{Armadillo} code.
//'   
//' @examples
//' \dontrun{
//' # Calculate a random matrix
//' matrixv <- matrix(rnorm(100), nc=10)
//' # Define the initial value of the inverse matrix
//' invmat <- solve(matrixv) + matrix(rnorm(100, sd=0.1), nc=10)
//' # Calculate the inverse in place using RcppArmadillo
//' HighFreq::calc_invrec(matrixv, invmat, 3)
//' # Multiply the matrix times its inverse
//' multmat <- matrixv %*% invmat
//' round(multmat, 4)
//' # Calculate the sum of the off-diagonal elements
//' sum(multmat[upper.tri(multmat)])
//' # Compare RcppArmadillo with R
//' all.equal(invmat, solve(matrixv))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'    rcode=solve(matrixv),
//'    cppcode=HighFreq::calc_invrec(matrixv, invmat, 3),
//'    times=10))[, c(1, 4, 5)]
//' }
//' 
//' @export
// [[Rcpp::export]]
void calc_invrec(const arma::mat& matrixv, 
                arma::mat& invmat, 
                arma::uword niter=1) {
  
  // Calculate the inverse of matrixv recursively
  for (arma::uword it = 0; it < niter; it++) {
    invmat = 2*invmat - invmat*matrixv*invmat;
  }  // end for
  
}  // end calc_invrec



////////////////////////////////////////////////////////////
//' Calculate the inverse of a square \emph{matrix} in place, without copying
//' the data in memory.
//' 
//' @param \code{matrixv} A \emph{matrix} of data to be inverted.  (The argument
//'   is interpreted as a \emph{pointer} to a \emph{matrix}, and it is
//'   overwritten with the inverse matrix.)
//' 
//' @return No return value.
//'
//' @details
//'   The function \code{calc_invref()} calculates the inverse of a square
//'   \emph{matrix} in place, without copying the data in memory. It accepts a
//'   \emph{pointer} to the argument \code{matrixv} (which is the \code{matrix}
//'   to be inverted), and it overwrites the old matrix values with the inverse
//'   matrix values. It performs the calculation in place, without copying the
//'   data in memory, which can significantly increase the computation speed for
//'   large matrices.
//'
//'   The function \code{calc_invref()} doesn't return a value.
//'   The function \code{calc_invref()} calls the \code{C++} \code{Armadillo}
//'   function \code{arma::inv()} to calculate the matrix inverse.
//'   
//' @examples
//' \dontrun{
//' # Calculate a random matrix
//' matrixv <- matrix(rnorm(100), nc=10)
//' # Copy matrixv to a matrix in a different memory location
//' invmat <- matrixv + 0
//' # Calculate the inverse in place using RcppArmadillo
//' HighFreq::calc_invref(invmat)
//' # Multiply the matrix times its inverse
//' multmat <- matrixv %*% invmat
//' round(multmat, 4)
//' # Calculate the sum of the off-diagonal elements
//' sum(multmat[upper.tri(multmat)])
//' # Compare RcppArmadillo with R
//' all.equal(invmat, solve(matrixv))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'    rcode=solve(matrixv),
//'    cppcode=calc_invref(matrixv),
//'    times=10))[, c(1, 4, 5)]
//' }
//' 
//' @export
// [[Rcpp::export]]
void calc_invref(arma::mat& matrixv) {
  
  matrixv = arma::inv(matrixv);
  
}  // end calc_invref



////////////////////////////////////////////////////////////
//' Standardize (center and scale) the columns of a \emph{time series} of data
//' in place, without copying the data in memory, using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A \emph{time series} or \emph{matrix} of data.
//' 
//' @param \code{center} A \emph{Boolean} argument: if \code{TRUE} then center
//'   the columns so that they have zero mean or median (the default is
//'   \code{TRUE}).
//' 
//' @param \code{scale} A \emph{Boolean} argument: if \code{TRUE} then scale the
//'   columns so that they have unit standard deviation or MAD (the default is
//'   \code{TRUE}).
//' 
//' @param \code{use_median} A \emph{Boolean} argument: if \code{TRUE} then the 
//'   centrality (central tendency) is calculated as the \emph{median} and the 
//'   dispersion is calculated as the \emph{median absolute deviation}
//'   (\emph{MAD}) (the default is \code{FALSE}).
//'   If \code{use_median = FALSE} then the centrality is calculated as the
//'   \emph{mean} and the dispersion is calculated as the \emph{standard
//'   deviation}.
//'
//' @return Void (no return value - modifies the data in place).
//'
//' @details
//'   The function \code{calc_scale()} standardizes (centers and scales) the
//'   columns of a \emph{time series} of data in place, without copying the data
//'   in memory, using \code{RcppArmadillo}.
//'
//'   If the arguments \code{center} and \code{scale} are both \code{TRUE} and
//'   \code{use_median} is \code{FALSE} (the defaults), then \code{calc_scale()}
//'   performs the same calculation as the standard \code{R} function
//'   \code{scale()}, and it calculates the centrality (central tendency) as the
//'   \emph{mean} and the dispersion as the \emph{standard deviation}.
//'
//'   If the arguments \code{center} and \code{scale} are both \code{TRUE} (the
//'   defaults), then \code{calc_scale()} standardizes the data.
//'   If the argument \code{center} is \code{FALSE} then \code{calc_scale()}
//'   only scales the data (divides it by the standard deviations).
//'   If the argument \code{scale} is \code{FALSE} then \code{calc_scale()}
//'   only demeans the data (subtracts the means).
//'   
//'   If the argument \code{use_median} is \code{TRUE}, then it calculates the
//'   centrality as the \emph{median} and the dispersion as the \emph{median
//'   absolute deviation} (\emph{MAD}).
//'
//'   If the number of rows of \code{tseries} is less than \code{3} then it
//'   does nothing and \code{tseries} is not scaled.
//'   
//'   The function \code{calc_scale()} accepts a \emph{pointer} to the argument
//'   \code{tseries}, and it overwrites the old data with the standardized data.
//'   It performs the calculation in place, without copying the data in memory,
//'   which can significantly increase the computation speed for large time
//'   series.
//'
//'   The function \code{calc_scale()} uses \code{RcppArmadillo} \code{C++}
//'   code, so on a typical time series it can be over \emph{10} times faster
//'   than the function \code{scale()}.
//'   
//' @examples
//' \dontrun{
//' # Calculate a time series of returns
//' retp <- zoo::coredata(na.omit(rutils::etfenv$returns[, c("IEF", "VTI")]))
//' # Demean the returns
//' demeaned <- apply(retp, 2, function(x) (x-mean(x)))
//' HighFreq::calc_scale(retp, scale=FALSE)
//' all.equal(demeaned, retp, check.attributes=FALSE)
//' # Calculate a time series of returns
//' retp <- zoo::coredata(na.omit(rutils::etfenv$returns[, c("IEF", "VTI")]))
//' # Standardize the returns
//' retss <- scale(retp)
//' HighFreq::calc_scale(retp)
//' all.equal(retss, retp, check.attributes=FALSE)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcode=scale(retp),
//'   Rcpp=HighFreq::calc_scale(retp),
//'   times=100))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
void calc_scale(arma::mat& tseries, 
                bool center = true, 
                bool scale = true, 
                bool use_median = false) {
  
  arma::uword nrows = tseries.n_rows;
  arma::uword ncols = tseries.n_cols;
  
  // Perform a loop over the columns
  if (nrows > 2) {
    if (scale and center) {
      if (use_median) {
        for (arma::uword it = 0; it < ncols; it++) {
          tseries.col(it) = (tseries.col(it) - center*arma::median(tseries.col(it)));
          tseries.col(it) = tseries.col(it)/arma::median(arma::abs(tseries.col(it)));
        }  // end for
      } else {
        for (arma::uword it = 0; it < ncols; it++) {
          tseries.col(it) = (tseries.col(it) - center*arma::mean(tseries.col(it)));
          tseries.col(it) = tseries.col(it)/arma::stddev(tseries.col(it));
        }  // end for
      }  // end if
    } else if (scale and (not center)) {
      if (use_median) {
        for (arma::uword it = 0; it < ncols; it++) {
          tseries.col(it) = tseries.col(it)/arma::median(arma::abs(tseries.col(it)));
        }  // end for
      } else {
        for (arma::uword it = 0; it < ncols; it++) {
          tseries.col(it) = tseries.col(it)/arma::stddev(tseries.col(it));
        }  // end for
      }  // end if
    } else if ((not scale) and center) {
      if (use_median) {
        for (arma::uword it = 0; it < ncols; it++) {
          tseries.col(it) = (tseries.col(it) - center*arma::median(tseries.col(it)));
        }  // end for
      } else {
        for (arma::uword it = 0; it < ncols; it++) {
          tseries.col(it) = (tseries.col(it) - center*arma::mean(tseries.col(it)));
        }  // end for
      }  // end if
    }  // end if
  }  // end if
  
}  // end calc_scale




////////////////////////////////////////////////////////////
// Functions for rolling aggregations
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
//' Aggregate a time series of data into a single bar of \emph{OHLC} data.
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} with multiple
//'   columns of data.
//'
//' @return A \emph{matrix} containing a single row, with the \emph{open},
//'   \emph{high}, \emph{low}, and \emph{close} values, and also the total
//'   \emph{volume} (if provided as either the second or fifth column of
//'   \code{tseries}).
//'
//' @details
//'   The function \code{agg_ohlc()} aggregates a time series of data into a
//'   single bar of \emph{OHLC} data. It can accept either a single column of
//'   data or four columns of \emph{OHLC} data.
//'   It can also accept an additional column containing the trading volume.
//'   
//' The function \code{agg_ohlc()} calculates the \emph{open} value as equal to
//' the \emph{open} value of the first row of \code{tseries}.
//'   The \emph{high} value as the maximum of the \emph{high} column of
//'   \code{tseries}.
//'   The \emph{low} value as the minimum of the \emph{low} column of
//'   \code{tseries}.
//'   The \emph{close} value as the \emph{close} of the last row of
//'   \code{tseries}.
//'   The \emph{volume} value as the sum of the \emph{volume} column of
//'   \code{tseries}.
//'
//'   For a single column of data, the \emph{open}, \emph{high}, \emph{low}, and
//'   \emph{close} values are all the same.
//'
//' @examples
//' \dontrun{
//' # Define matrix of OHLC data
//' ohlc <- coredata(rutils::etfenv$VTI[, 1:5])
//' # Aggregate to single row matrix
//' ohlcagg <- HighFreq::agg_ohlc(ohlc)
//' # Compare with calculation in R
//' all.equal(drop(ohlcagg),
//'   c(ohlc[1, 1], max(ohlc[, 2]), min(ohlc[, 3]), ohlc[NROW(ohlc), 4], sum(ohlc[, 5])), 
//'   check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat agg_ohlc(const arma::mat& tseries) {
  
  arma::uword nrows = tseries.n_rows;
  arma::uword ncols = tseries.n_cols;
  
  // Number of output columns
  arma::uword numohlc = ncols;
  if (ncols < 4)
    // Add volume column for non-OHLC data
    numohlc = 4 + ncols - 1;
  // Allocate output matrix
  arma::mat ohlc(1, numohlc);
  
  if (ncols < 4) {
    // Aggregate time series into a single bar of OHLC data.
    ohlc(0, 0) = tseries(0, 0);
    ohlc(0, 1) = arma::max(tseries.col(0));
    ohlc(0, 2) = arma::min(tseries.col(0));
    ohlc(0, 3) = tseries(nrows-1, 0);
    if (ncols == 2) {
      // Aggregate volume data.
      ohlc(0, 4) = arma::sum(tseries.col(1));
    }  // end if
  } else {
    // Aggregate OHLC time series into a single bar of OHLC data.
    ohlc(0, 0) = tseries(0, 0);
    ohlc(0, 1) = arma::max(tseries.col(1));
    ohlc(0, 2) = arma::min(tseries.col(2));
    ohlc(0, 3) = tseries(nrows-1, 3);
    if (ncols == 5) {
      // Aggregate volume data.
      ohlc(0, 4) = arma::sum(tseries.col(4));
    }  // end if
  }  // end if
  
  return ohlc;
  
}  // end agg_ohlc




////////////////////////////////////////////////////////////
//' Aggregate a time series to an \emph{OHLC} time series with lower
//' periodicity.
//'
//' Given a time series of prices at a higher periodicity (say seconds), it
//' calculates the \emph{OHLC} prices at a lower periodicity (say minutes).
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} with multiple
//'   columns of data.
//'   
//' @param \emph{endd} An \emph{integer vector} of end points.
//'
//' @return A \emph{matrix} with \emph{OHLC} data, with the number of rows equal
//'   to the number of \emph{endd} minus one.
//'   
//' @details
//'   The function \code{roll_ohlc()} performs a loop over the end points
//'   \emph{endd}, along the rows of the data \code{tseries}. At each end point,
//'   it selects the past rows of the data \code{tseries}, starting at the first
//'   bar after the previous end point, and then calls the function
//'   \code{agg_ohlc()} on the selected data \code{tseries} to calculate the
//'   aggregations.
//'   
//'   The function \code{roll_ohlc()} can accept either a single column of data
//'   or four columns of \emph{OHLC} data.
//'   It can also accept an additional column containing the trading volume.
//'
//'   The function \code{roll_ohlc()} performs a similar aggregation as the
//'   function \code{to.period()} from package
//'   \href{https://cran.r-project.org/web/packages/xts/index.html}{xts}.
//'
//' @examples
//' \dontrun{
//' # Define matrix of OHLC data
//' ohlc <- rutils::etfenv$VTI[, 1:5]
//' # Define end points at 25 day intervals
//' endd <- HighFreq::calc_endpoints(NROW(ohlc), step=25)
//' # Aggregate over endd:
//' ohlcagg <- HighFreq::roll_ohlc(tseries=ohlc, endd=endd)
//' # Compare with xts::to.period()
//' ohlcagg_xts <- .Call("toPeriod", ohlc, as.integer(endd+1), TRUE, NCOL(ohlc), FALSE, FALSE, colnames(ohlc), PACKAGE="xts")
//' all.equal(ohlcagg, coredata(ohlcagg_xts), check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_ohlc(const arma::mat& tseries, arma::uvec endd) {
  
  // arma::uword nrows = tseries.n_rows;
  arma::uword ncols = tseries.n_cols;
  // Number of output rows
  arma::uword numpts = endd.size();
  // Number of output columns
  arma::uword numohlc = 5; // With volume column
  if ((ncols == 1) || (ncols == 4))
    numohlc = 4; // No volume column
  // Allocate output matrix
  arma::mat ohlcagg(numpts-1, numohlc, fill::zeros);
  
  // Perform loop over the end points
  for (arma::uword it = 1; it < numpts; it++) {
    // cout << "it: " << it << endl;
    // Aggregate the OHLC
    ohlcagg.row(it-1) = agg_ohlc(tseries.rows(endd(it-1)+1, endd(it)));
  }  // end for
  
  // Return the aggregations
  return ohlcagg;
  
}  // end roll_ohlc



////////////////////////////////////////////////////////////
//' Calculate the rolling convolutions (weighted sums) of a \emph{time series}
//' with a single-column \emph{matrix} of weights.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//' 
//' @param \code{weightv} A single-column \emph{matrix} of weights.
//'
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{roll_conv()} calculates the convolutions of the
//'   \emph{matrix} columns with a single-column \emph{matrix} of weights.  It
//'   performs a loop over the \emph{matrix} rows and multiplies the past
//'   (higher) values by the weights.  It calculates the rolling weighted sums
//'   of the past data.
//'   
//'   The function \code{roll_conv()} uses the \code{Armadillo} function
//'   \code{arma::conv2()}. It performs a similar calculation to the standard
//'   \code{R} function \cr\code{filter(x=tseries, filter=weightv,
//'   method="convolution", sides=1)}, but it's over \code{6} times faster, and
//'   it doesn't produce any leading \code{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Calculate a time series of returns
//' retp <- na.omit(rutils::etfenv$returns[, c("IEF", "VTI")])
//' # Create simple weights equal to a 1 value plus zeros
//' weightv <- c(1, rep(0, 10))
//' # Calculate rolling weighted sums
//' retf <- HighFreq::roll_conv(retp, weightv)
//' # Compare with original
//' all.equal(coredata(retp), retf, check.attributes=FALSE)
//' # Second example
//' # Calculate exponentially decaying weights
//' weightv <- exp(-0.2*(1:11))
//' weightv <- weightv/sum(weightv)
//' # Calculate rolling weighted sums
//' retf <- HighFreq::roll_conv(retp, weightv)
//' # Calculate rolling weighted sums using filter()
//' retc <- filter(x=retp, filter=weightv, method="convolution", sides=1)
//' # Compare both methods
//' all.equal(retc[-(1:11), ], retf[-(1:11), ], check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_conv(const arma::mat& tseries, const arma::colvec& weightv) {
  
  arma::uword lookb = weightv.n_rows-2;
  arma::uword nrows = tseries.n_rows-1;
  
  // Calculate the convolutions
  arma::mat convmat = arma::conv2(tseries, weightv, "full");
  
  // Copy the warmup period
  convmat.rows(0, lookb) = tseries.rows(0, lookb);
  
  return convmat.rows(0, nrows);
  
}  // end roll_conv




////////////////////////////////////////////////////////////
//' Calculate the rolling sums over a \emph{time series} or a \emph{matrix}
//' using \emph{Rcpp}.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lookb} The length of the look-back interval, equal to the
//'   number of data points included in calculating the rolling sum (the default
//'   is \code{lookb = 1}).
//'   
//' @param \code{weightv} A single-column \emph{matrix} of weights (the default
//'   is \code{weightv = 0}).
//'
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{roll_sum()} calculates the rolling \emph{weighted} sums
//'   over the columns of the data \code{tseries}.
//'   
//'   If the argument \code{weightv} is equal to zero (the default), then the
//'   function \code{roll_sum()} calculates the simple rolling sums of the 
//'   \emph{time series} data \eqn{p_t} over the look-back interval \eqn{\Delta}:
//'   \deqn{
//'     \bar{p}_t = \sum_{j=(t-\Delta+1)}^{t} p_j
//'   }
//'   
//'   If the \code{weightv} argument has the same number of rows as the argument
//'   \code{tseries}, then the function \code{roll_sum()} calculates rolling 
//'   \emph{weighted} sums of the \emph{time series} data \eqn{p_t} in two steps.
//'   
//'   It first calculates the rolling sums of the products of the weights
//'   \eqn{w_t} times the \emph{time series} data \eqn{p_t} over the look-back
//'   interval \eqn{\Delta}:
//'   \deqn{
//'     \bar{w}_t = \sum_{j=(t-\Delta+1)}^{t} w_j
//'   }
//'   \deqn{
//'     \bar{p}^w_t = \sum_{j=(t-\Delta+1)}^{t} w_j p_j
//'   }
//'   
//'   It then calculates the rolling \emph{weighted} sums \eqn{\bar{p}_t} as the
//'   ratio of the sum products of the weights and the data, divided by the sums
//'   of the weights:
//'   \deqn{
//'     \bar{p}_t = \frac{\bar{p}^w_t}{\bar{w}_t}
//'   }
//'   
//'   The function \code{roll_sum()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{tseries}.
//' 
//'   The function \code{roll_sum()} is written in \code{C++} \code{Armadillo}
//'   code, so it's much faster than equivalent \code{R} code.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' retp <- na.omit(rutils::etfenv$returns[, c("VTI", "IEF")])
//' # Define parameters
//' lookb <- 11
//' # Calculate rolling sums and compare with rutils::roll_sum()
//' sumc <- HighFreq::roll_sum(retp, lookb)
//' sumr <- rutils::roll_sum(retp, lookb)
//' all.equal(sumc, coredata(sumr), check.attributes=FALSE)
//' # Calculate rolling sums using R code
//' sumr <- apply(zoo::coredata(retp), 2, cumsum)
//' sumlag <- rbind(matrix(numeric(2*lookb), nc=2), sumr[1:(NROW(sumr) - lookb), ])
//' sumr <- (sumr - sumlag)
//' all.equal(sumc, sumr, check.attributes=FALSE)
//' # Calculate weights equal to the trading volumes
//' weightv <- quantmod::Vo(rutils::etfenv$VTI)
//' weightv <- weightv[zoo::index(retp)]
//' # Calculate rolling weighted sums
//' sumc <- HighFreq::roll_sum(retp, lookb, 1/weightv)
//' # Plot dygraph of the weighted sums
//' datav <- cbind(retp$VTI, sumc[, 1])
//' colnames(datav) <- c("VTI", "Weighted")
//' endd <- rutils::calc_endpoints(datav, interval="weeks")
//' dygraphs::dygraph(cumsum(datav)[endd], main=colnames(datav)) %>% 
//'   dyOptions(colors=c("blue", "red"), strokeWidth=2) %>% 
//'   dyLegend(width=300)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_sum(const arma::mat& tseries, 
                   arma::uword lookb = 1,
                   const arma::colvec& weightv = 0) {
  
  arma::uword startp;
  arma::uword nrows = tseries.n_rows;
  arma::uword nweights = weightv.n_elem;
  arma::mat sumr(nrows, tseries.n_cols);
  
  if (nweights == 1) { // Calculate sums without weights
    // Warmup period
    sumr.row(0) = tseries.row(0);
    for (arma::uword it = 1; it < lookb; it++) {
      sumr.row(it) = sumr.row(it-1) + tseries.row(it);
    }  // end for
    // Perform loop over the remaining rows
    for (arma::uword it = lookb; it < nrows; it++) {
      startp = it - lookb;
      sumr.row(it) = sumr.row(it-1) + tseries.row(it) - tseries.row(startp);
    }  // end for
    // Old, clever code - doesn't produce correct warmup values
    // Calculate the cumulative sum
    // arma::mat cumsumv = arma::cumsum(tseries, 0);
    // Return the differences of the cumulative sum
    // sumr = diffit(cumsumv, lookb, true);
  } else if (nweights == nrows) { // Calculate sums with weights
    // Calculate means with weights
    double weightr = weightv(0);
    double weightrr;
    sumr.row(0) = weightr*tseries.row(0);
    // Warmup period
    for (arma::uword it = 1; it < lookb; it++) {
      weightrr = weightr;
      // Sum up the weighted means
      weightr = weightr + weightv(it);
      sumr.row(it) = sumr.row(it-1) + weightv(it)*tseries.row(it);
      // Divide by the mean weight
      sumr.row(it-1) = sumr.row(it-1)/weightrr;
    }  // end for
    // Perform loop over the remaining rows
    for (arma::uword it = lookb; it < nrows; it++) {
      startp = it - lookb;
      weightrr = weightr;
      // Update the weighted means - add the new element and subtract the old one
      weightr = weightr + weightv(it) - weightv(startp);
      sumr.row(it) = sumr.row(it-1) + weightv(it)*tseries.row(it) - weightv(startp)*tseries.row(startp);
      // Divide by the mean weight
      sumr.row(it-1) = sumr.row(it-1)/weightrr;
    }  // end for
    sumr.row(nrows-1) = sumr.row(nrows-1)/weightr;
  }  // end if
  
  return sumr;
  
}  // end roll_sum




////////////////////////////////////////////////////////////
//' Calculate the rolling sums at the end points of a \emph{time series} or a
//' \emph{matrix}.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endd} An \emph{integer} vector of end points (the default is
//'   \code{endd = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{lookb} The number of end points in the look-back interval
//'   (the default is \code{lookb = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points.
//'
//' @return A \emph{matrix} with the same number of columns as the input time
//'   series \code{tseries}, and the number of rows equal to the number of end
//'   points.
//'   
//' @details
//'   The function \code{roll_sumep()} calculates the rolling sums at the end
//'   points of the \emph{time series} \code{tseries}.
//'   
//'   The function \code{roll_sumep()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, which makes it several times faster than \code{R} code.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' retp <- na.omit(rutils::etfenv$returns[, c("VTI", "IEF")])
//' # Define end points at 25 day intervals
//' endd <- HighFreq::calc_endpoints(NROW(retp), step=25)
//' # Define start points as 75 day lag of end points
//' startp <- HighFreq::calc_startpoints(endd, lookb=3)
//' # Calculate rolling sums using Rcpp
//' sumc <- HighFreq::roll_sumep(retp, startp=startp, endd=endd)
//' # Calculate rolling sums using R code
//' sumr <- sapply(1:NROW(endd), function(ep) {
//' colSums(retp[(startp[ep]+1):(endd[ep]+1), ])
//'   })  # end sapply
//' sumr <- t(sumr)
//' all.equal(sumc, sumr, check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_sumep(const arma::mat& tseries, 
                     arma::uvec startp = 0, 
                     arma::uvec endd = 0, 
                     arma::uword step = 1, 
                     arma::uword lookb = 1, 
                     arma::uword stub = 0) {
  
  // Allocate end points
  arma::uword nrows = tseries.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endd) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endd;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by lookb
    startpts = calc_startpoints(endpts, lookb);
    // Old code for startpts
    // Start points equal to end points lagged by lookb - without adding +1
    // arma::uword numpts = endpts.n_elem;
    // arma::uvec startpts = arma::join_cols(arma::zeros<uvec>(lookb), 
    //                                        endpts.subvec(0, numpts - lookb - 1));
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  // Allocate sums matrix
  arma::uword numpts = endpts.n_elem;
  arma::mat sums = arma::zeros(numpts, tseries.n_cols);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate sums
    if (endpts(ep) > startpts(ep)) {
      sums.row(ep) = arma::sum(tseries.rows(startpts(ep), endpts(ep)));
    }  // end if
  }  // end for
  
  return sums;
  
  // Old code using arma::cumsum() - sums exclude startpts
  // Calculate cumulative sums at end points.
  // arma::mat cumsumv = arma::cumsum(tseries, 0);
  // arma::mat cum_start = cumsumv.rows(startpts);
  // arma::mat cum_end = cumsumv.rows(endpts);
  
  // Return the differences of the cumulative returns
  // return diffit(cumsumv, 1, true);
  // return (cum_end - cum_start);
  
}  // end roll_sumep




////////////////////////////////////////////////////////////
//' Calculate the rolling weighted sums over a \emph{time series} or a
//' \emph{matrix} using \emph{Rcpp}.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{endd} An \emph{integer} vector of end points (the default is
//'   \code{endd = NULL}).
//'   
//' @param \code{lookb} The length of the look-back interval, equal to the
//'   number of data points included in calculating the rolling sum (the default
//'   is \code{lookb = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = NULL}).
//' 
//' @param \code{weightv} A single-column \emph{matrix} of weights (the default
//'   is \code{weightv = NULL}).
//'
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{roll_sumw()} calculates the rolling weighted sums over
//'   the columns of the data \code{tseries}.
//' 
//'   The function \code{roll_sumw()} calculates the rolling weighted sums as
//'   convolutions of the columns of \code{tseries} with the \emph{column
//'   vector} of weights using the \code{Armadillo} function
//'   \code{arma::conv2()}.  It performs a similar calculation to the standard
//'   \code{R} function \cr\code{stats::filter(x=retp, filter=weightv,
//'   method="convolution", sides=1)}, but it can be many times faster, and it
//'   doesn't produce any leading \code{NA} values.
//'   
//'   The function \code{roll_sumw()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{tseries}.
//' 
//'   The arguments \code{weightv}, \code{endd}, and \code{stub} are
//'   optional.
//'   
//'   If the argument \code{weightv} is not supplied, then simple sums are
//'   calculated, not weighted sums.
//'   
//'   If either the \code{stub} or \code{endd} arguments are supplied,
//'   then the rolling sums are calculated at the end points. 
//'   
//'   If only the argument \code{stub} is supplied, then the end points are
//'   calculated from the \code{stub} and \code{lookb} arguments. The first
//'   end point is equal to \code{stub} and the end points are spaced
//'   \code{lookb} periods apart.
//'   
//'   If the arguments \code{weightv}, \code{endd}, and \code{stub} are
//'   not supplied, then the sums are calculated over a number of data points
//'   equal to \code{lookb}.
//'   
//'   The function \code{roll_sumw()} is also several times faster than
//'   \code{rutils::roll_sum()} which uses vectorized \code{R} code.
//'   
//'   Technical note:
//'   The function \code{roll_sumw()} has arguments with default values equal to
//'   \code{NULL}, which are implemented in \code{Rcpp} code.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Calculate historical returns
//' retp <- na.omit(rutils::etfenv$returns[, c("VTI", "IEF")])
//' # Define parameters
//' lookb <- 11
//' # Calculate rolling sums and compare with rutils::roll_sum()
//' sumc <- HighFreq::roll_sum(retp, lookb)
//' sumr <- rutils::roll_sum(retp, lookb)
//' all.equal(sumc, coredata(sumr), check.attributes=FALSE)
//' # Calculate rolling sums using R code
//' sumr <- apply(zoo::coredata(retp), 2, cumsum)
//' sumlag <- rbind(matrix(numeric(2*lookb), nc=2), sumr[1:(NROW(sumr) - lookb), ])
//' sumr <- (sumr - sumlag)
//' all.equal(sumc, sumr, check.attributes=FALSE)
//' 
//' # Calculate rolling sums at end points
//' stubv <- 21
//' sumc <- HighFreq::roll_sumw(retp, lookb, stub=stubv)
//' endd <- (stubv + lookb*(0:(NROW(retp) %/% lookb)))
//' endd <- endd[endd < NROW(retp)]
//' sumr <- apply(zoo::coredata(retp), 2, cumsum)
//' sumr <- sumr[endd+1, ]
//' sumlag <- rbind(numeric(2), sumr[1:(NROW(sumr) - 1), ])
//' sumr <- (sumr - sumlag)
//' all.equal(sumc, sumr, check.attributes=FALSE)
//' 
//' # Calculate rolling sums at end points - pass in endd
//' sumc <- HighFreq::roll_sumw(retp, endd=endd)
//' all.equal(sumc, sumr, check.attributes=FALSE)
//' 
//' # Create exponentially decaying weights
//' weightv <- exp(-0.2*(1:11))
//' weightv <- matrix(weightv/sum(weightv), nc=1)
//' # Calculate rolling weighted sum
//' sumc <- HighFreq::roll_sumw(retp, weightv=weightv)
//' # Calculate rolling weighted sum using filter()
//' retc <- filter(x=retp, filter=weightv, method="convolution", sides=1)
//' all.equal(sumc[-(1:11), ], retc[-(1:11), ], check.attributes=FALSE)
//' 
//' # Calculate rolling weighted sums at end points
//' sumc <- HighFreq::roll_sumw(retp, endd=endd, weightv=weightv)
//' all.equal(sumc, retc[endd+1, ], check.attributes=FALSE)
//' 
//' # Create simple weights equal to a 1 value plus zeros
//' weightv <- matrix(c(1, rep(0, 10)), nc=1)
//' # Calculate rolling weighted sum
//' weighted <- HighFreq::roll_sumw(retp, weightv=weightv)
//' # Compare with original
//' all.equal(coredata(retp), weighted, check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_sumw(const arma::mat& tseries,
                    Rcpp::Nullable<Rcpp::IntegerVector> endd = R_NilValue, 
                    arma::uword lookb = 1,
                    Rcpp::Nullable<int> stub = R_NilValue, 
                    Rcpp::Nullable<Rcpp::NumericVector> weightv = R_NilValue) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat cumsumv;
  
  // Calculate the rolling sums
  if (weightv.isNotNull()) {
    // Coerce weights from Rcpp to Armadillo vector
    arma::vec weighta = Rcpp::as<vec>(weightv);
    arma::uword nweights = weighta.n_elem;
    // Calculate the weighted averages as convolutions
    cumsumv = arma::conv2(tseries, weighta, "full");
    // Copy the warmup period
    // cout << "nweights = " << nweights << endl;
    cumsumv.rows(0, nweights-2) = tseries.rows(0, nweights-2);
    cumsumv = cumsumv.rows(0, nrows-1);
    // cout << "cumsumv.n_rows = " << cumsumv.n_rows << endl;
  } else {
    // Calculate the cumulative sum
    cumsumv = arma::cumsum(tseries, 0);
  }  // end if
  
  
  // Declare empty end points
  arma::uvec endpts;
  // Update end points
  if (endd.isNotNull()) {
    // Copy endd
    endpts = Rcpp::as<uvec>(endd);
  } else if (stub.isNotNull()) {
    // Calculate end points with stub
    endpts = arma::regspace<uvec>(Rcpp::as<uword>(stub), lookb, nrows + lookb);
    endpts = endpts.elem(find(endpts < nrows));
  }  // end if
  
  
  // Subset the rolling sums according the end points
  if (endpts.is_empty() && weightv.isNotNull()) {
    // Do nothing
    // Return the weighted averages (convolutions) at each point
    // return cumsumv;
  } else if (endpts.is_empty() && !weightv.isNotNull()) {
    // Return unweighted rolling sums at each point
    cumsumv = diffit(cumsumv, lookb, true);
  } else if (!endpts.is_empty() && weightv.isNotNull()) {
    // Return the weighted averages (convolutions) at end points
    cumsumv = cumsumv.rows(endpts);
  } else if (!endpts.is_empty() && !weightv.isNotNull()) {
    // Return the unweighted rolling sums at end points
    cumsumv = cumsumv.rows(endpts);
    cumsumv = diffit(cumsumv, 1, true);
  }  // end if
  
  return cumsumv;
  
}  // end roll_sumw




////////////////////////////////////////////////////////////
// Functions for online aggregations of streaming data
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
//' Calculate the exponential moving average (EMA) of streaming \emph{time
//' series} data using an online recursive formula.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lambda} A decay factor which multiplies past estimates.
//'   
//' @param \code{weightv} A single-column \emph{matrix} of weights.
//'
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{run_mean()} calculates the exponential moving average
//'   (EMA) of the streaming \emph{time series} data \eqn{p_t} by recursively
//'   weighting present and past values using the decay factor \eqn{\lambda}. If
//'   the \code{weightv} argument is equal to zero, then the function
//'   \code{run_mean()} simply calculates the exponentially weighted moving
//'   average value of the streaming \emph{time series} data \eqn{p_t}:
//'   \deqn{
//'     \bar{p}_t = \lambda \bar{p}_{t-1} + (1-\lambda) p_t = (1-\lambda) \sum_{j=0}^{n} \lambda^j p_{t-j}
//'   }
//'   
//'   Some applications require applying additional weight factors, like for
//'   example the volume-weighted average price indicator (VWAP).  Then the
//'   streaming prices can be multiplied by the streaming trading volumes.
//'   
//'   If the argument \code{weightv} has the same number of rows as the argument
//'   \code{tseries}, then the function \code{run_mean()} calculates the
//'   exponential moving average (EMA) in two steps.
//'   
//'   First it calculates the trailing mean weights \eqn{\bar{w}_t}:
//'   \deqn{
//'     \bar{w}_t = \lambda \bar{w}_{t-1} + (1-\lambda) w_t
//'   }
//'   
//'   Second it calculates the trailing mean products \eqn{\bar{w p}_t} of the
//'   weights \eqn{w_t} and the data \eqn{p_t}:
//'   \deqn{
//'     \bar{w p}_t = \lambda \bar{w p}_{t-1} + (1-\lambda) w_t p_t
//'   }
//'   Where \eqn{p_t} is the streaming data, \eqn{w_t} are the streaming
//'   weights, \eqn{\bar{w}_t} are the trailing mean weights, and \eqn{\bar{w p}_t}
//'   are the trailing mean products of the data and the weights.
//'   
//'   The trailing mean weighted value \eqn{\bar{p}_t} is equal to the ratio of the
//'   data and weights products, divided by the mean weights:
//'   \deqn{
//'     \bar{p}_t = \frac{\bar{w p}_t}{\bar{w}_t}
//'   }
//' 
//'   The above online recursive formulas are convenient for processing live
//'   streaming data because they don't require maintaining a buffer of past
//'   data.
//'   The formulas are equivalent to a convolution with exponentially decaying
//'   weights, but they're much faster to calculate.
//'   Using exponentially decaying weights is more natural than using a sliding
//'   look-back interval, because it gradually "forgets" about the past data.
//'   
//'   The value of the decay factor \eqn{\lambda} must be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the trailing mean values have a stronger
//'   dependence on past data.  This is equivalent to a long look-back
//'   interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the trailing mean values have a
//'   weaker dependence on past data.  This is equivalent to a short look-back
//'   interval.
//' 
//'   The function \code{run_mean()} performs the same calculation as the
//'   standard \code{R} function\cr\code{stats::filter(x=series, filter=lambda,
//'   method="recursive")}, but it's several times faster.
//' 
//'   The function \code{run_mean()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{tseries}.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical prices
//' ohlc <- rutils::etfenv$VTI
//' closep <- quantmod::Cl(ohlc)
//' # Calculate the trailing means
//' lambdaf <- 0.9 # Decay factor
//' meanv <- HighFreq::run_mean(closep, lambda=lambdaf)
//' # Calculate the trailing means using R code
//' pricef <- (1-lambdaf)*filter(closep, 
//'   filter=lambdaf, init=as.numeric(closep[1, 1])/(1-lambdaf), 
//'   method="recursive")
//' all.equal(drop(meanv), unclass(pricef), check.attributes=FALSE)
//' 
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::run_mean(closep, lambda=lambdaf),
//'   Rcode=filter(closep, filter=lambdaf, init=as.numeric(closep[1, 1])/(1-lambdaf), method="recursive"),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//'   
//' # Calculate weights equal to the trading volumes
//' weightv <- quantmod::Vo(ohlc)
//' # Calculate the exponential moving average (EMA)
//' meanw <- HighFreq::run_mean(closep, lambda=lambdaf, weightv=weightv)
//' # Plot dygraph of the EMA
//' datav <- xts(cbind(meanv, meanw), zoo::index(ohlc))
//' colnames(datav) <- c("means trailing", "means weighted")
//' dygraphs::dygraph(datav, main="Trailing Means") %>%
//'   dyOptions(colors=c("blue", "red"), strokeWidth=2) %>%
//'   dyLegend(show="always", width=300)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_mean(const arma::mat& tseries, 
                   double lambda, // Decay factor which multiplies the past values 
                   const arma::colvec& weightv = 0) {
  
  arma::uword nrows = tseries.n_rows;
  arma::uword nweights = weightv.n_elem;
  arma::mat meanm(nrows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  if (nweights == 1) {
    meanm.row(0) = tseries.row(0);
    // Calculate means without weights
    for (arma::uword it = 1; it < nrows; it++) {
      // Calculate the means using the decay factor
      meanm.row(it) = lambda*meanm.row(it-1) + lambda1*tseries.row(it);
    }  // end for
  } else if (nweights == nrows) {
    // Calculate means with weights
    double meanw = weightv(0);
    double meanww;
    meanm.row(0) = meanw*tseries.row(0);
    for (arma::uword it = 1; it < nrows; it++) {
      // Calculate the means using the decay factor
      meanww = meanw;
      meanw = lambda*meanw + lambda1*weightv(it);
      meanm.row(it) = lambda*meanm.row(it-1) + lambda1*weightv(it)*tseries.row(it);
      // Divide by the mean weight
      meanm.row(it-1) = meanm.row(it-1)/meanww;
    }  // end for
    meanm.row(nrows-1) = meanm.row(nrows-1)/meanw;
  }  // end if
  
  return meanm;
  
}  // end run_mean



////////////////////////////////////////////////////////////
//' Calculate the trailing maximum values of streaming \emph{time series} data
//' using an online recursive formula.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lambda} A decay factor which multiplies past estimates.
//'   
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{run_max()} calculates the trailing maximum values of
//'   streaming \emph{time series} data by recursively weighting present and past
//'   values using the decay factor \eqn{\lambda}.
//'
//'   It calculates the trailing maximum values \eqn{p^{max}_t} of the streaming
//'   data \eqn{p_t} as follows:
//'   \deqn{
//'     p^{max}_t = max(p_t, \lambda p^{max}_{t-1} + (1-\lambda) p_t)
//'   }
//'   The first term in the sum is the maximum value multiplied by the decay
//'   factor \eqn{\lambda}, so that the past maximum value is gradually
//'   "forgotten". The second term pulls the maximum value to the current value
//'   \eqn{p_t}.
//'   
//'   The value of the decay factor \eqn{\lambda} must be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the past maximum values persist
//'   for longer.  This is equivalent to a long look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the past maximum values
//'   decay quickly, and the trailing maximum depends on the more recent
//'   streaming data.  This is equivalent to a short look-back interval.
//' 
//'   The above formula can also be expressed as:
//'   \deqn{
//'     p^{max}_t = \lambda max(p_t, p^{max}_{t-1}) + (1-\lambda) p_t
//'   }
//'   The first term is the maximum value multiplied by the decay factor
//'   \eqn{\lambda}, so that the past maximum value is gradually "forgotten".
//'   The second term pulls the maximum value to the current value \eqn{p_t}.
//'   
//'   The above recursive formula is convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//' 
//'   The function \code{run_max()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{tseries}.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical prices
//' closep <- zoo::coredata(quantmod::Cl(rutils::etfenv$VTI))
//' # Calculate the trailing maximums
//' lambdaf <- 0.9 # Decay factor
//' pricmax <- HighFreq::run_max(closep, lambda=lambdaf)
//' # Plot dygraph of VTI prices and trailing maximums
//' datav <- cbind(quantmod::Cl(rutils::etfenv$VTI), pricmax)
//' colnames(datav) <- c("prices", "max")
//' colnamev <- colnames(datav)
//' dygraphs::dygraph(datav, main="VTI Prices and Trailing Maximums") %>%
//'   dySeries(label=colnamev[1], strokeWidth=2, col="blue") %>%
//'   dySeries(label=colnamev[2], strokeWidth=2, col="red")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_max(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat maxv = arma::zeros(nrows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  // Perform loop over the rows
  maxv.row(0) = tseries.row(0);
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the max using the decay factor
    maxv.row(it) = lambda*arma::max(tseries.row(it), maxv.row(it-1)) + lambda1*tseries.row(it);
    // Alternative formula for the same:
    // maxv.row(it) = arma::max(tseries.row(it), lambda*maxv.row(it-1) + lambda1*tseries.row(it));
  }  // end for
  
  return maxv;
  
}  // end run_max




////////////////////////////////////////////////////////////
//' Calculate the trailing minimum values of streaming \emph{time series} data
//' using an online recursive formula.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lambda} A decay factor which multiplies past estimates.
//'   
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{run_min()} calculates the trailing minimum values of
//'   streaming \emph{time series} data by recursively weighting present and past
//'   values using the decay factor \eqn{\lambda}.
//'
//'   It calculates the trailing minimum values \eqn{p^{min}_t} of the streaming
//'   data \eqn{p_t} as follows:
//'   \deqn{
//'     p^{min}_t = min(p_t, \lambda p^{min}_{t-1} + (1-\lambda) p_t)
//'   }
//'   The first term in the sum is the minimum value multiplied by the decay
//'   factor \eqn{\lambda}, so that the past minimum value is gradually
//'   "forgotten". The second term pulls the minimum value to the current value
//'   \eqn{p_t}.
//'   
//'   The value of the decay factor \eqn{\lambda} must be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the past minimum values persist
//'   for longer.  This is equivalent to a long look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the past minimum values
//'   decay quickly, and the trailing minimum depends on the more recent
//'   streaming data.  This is equivalent to a short look-back interval.
//' 
//'   The above formula can also be expressed as:
//'   \deqn{
//'     p^{min}_t = \lambda min(p_t, p^{min}_{t-1}) + (1-\lambda) p_t
//'   }
//'   The first term is the minimum value multiplied by the decay factor
//'   \eqn{\lambda}, so that the past minimum value is gradually "forgotten".
//'   The second term pulls the minimum value to the current value \eqn{p_t}.
//'   
//'   The above recursive formula is convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//' 
//'   The function \code{run_min()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{tseries}.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical prices
//' closep <- zoo::coredata(quantmod::Cl(rutils::etfenv$VTI))
//' # Calculate the trailing minimums
//' lambdaf <- 0.9 # Decay factor
//' pricmin <- HighFreq::run_min(closep, lambda=lambdaf)
//' # Plot dygraph of VTI prices and trailing minimums
//' datav <- cbind(quantmod::Cl(rutils::etfenv$VTI), pricmin)
//' colnames(datav) <- c("prices", "min")
//' colnamev <- colnames(datav)
//' dygraphs::dygraph(datav, main="VTI Prices and Trailing Minimums") %>%
//'   dySeries(label=colnamev[1], strokeWidth=1, col="blue") %>%
//'   dySeries(label=colnamev[2], strokeWidth=1, col="red")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_min(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat minv = arma::zeros(nrows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  // Perform loop over the rows
  minv.row(0) = tseries.row(0);
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the min using the decay factor
    minv.row(it) = lambda*arma::min(tseries.row(it), minv.row(it-1)) + lambda1*tseries.row(it);
    // Alternative formula for the same:
    // minv.row(it) = arma::min(tseries.row(it), lambda*minv.row(it-1) + lambda1*tseries.row(it));
  }  // end for
  
  return minv;
  
}  // end run_min



////////////////////////////////////////////////////////////
//' Calculate the trailing variance of streaming \emph{time series} of data
//' using an online recursive formula.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//' 
//' @param \code{lambda} A decay factor which multiplies past estimates.
//'   
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{run_var()} calculates the trailing variance of
//'   streaming \emph{time series} of data \eqn{r_t}, by recursively weighting
//'   the past variance estimates \eqn{\sigma^2_{t-1}}, with the squared
//'   differences of the data minus its trailing means \eqn{(r_t -
//'   \bar{r}_t)^2}, using the decay factor \eqn{\lambda}:
//'   \deqn{
//'     \bar{r}_t = \lambda \bar{r}_{t-1} + (1-\lambda) r_t
//'   }
//'   \deqn{
//'     \sigma^2_t = \lambda \sigma^2_{t-1} + (1-\lambda) (r_t - \bar{r}_t)^2
//'   }
//'   Where \eqn{r_t} are the streaming data, \eqn{\bar{r}_t} are the trailing
//'   means, and \eqn{\sigma^2_t} are the trailing variance estimates.
//' 
//'   The above online recursive formulas are convenient for processing live
//'   streaming data because they don't require maintaining a buffer of past
//'   data.
//'   The formulas are equivalent to a convolution with exponentially decaying
//'   weights, but they're much faster to calculate.
//'   Using exponentially decaying weights is more natural than using a sliding
//'   look-back interval, because it gradually "forgets" about the past data.
//' 
//'   The value of the decay factor \eqn{\lambda} must be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the trailing variance values have a
//'   stronger dependence on past data.  This is equivalent to a long
//'   look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the trailing variance values have a
//'   weaker dependence on past data.  This is equivalent to a short look-back
//'   interval.
//' 
//'   The function \code{run_var()} performs the same calculation as the
//'   standard \code{R} function\cr\code{stats::filter(x=series,
//'   filter=weightv, method="recursive")}, but it's several times faster.
//' 
//'   The function \code{run_var()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{tseries}.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' retp <- zoo::coredata(na.omit(rutils::etfenv$returns$VTI))
//' # Calculate the trailing variance
//' lambdaf <- 0.9 # Decay factor
//' vars <- HighFreq::run_var(retp, lambda=lambdaf)
//' # Calculate centered returns
//' retc <- (retp - HighFreq::run_mean(retp, lambda=lambdaf))
//' # Calculate the trailing variance using R code
//' retc2 <- (1-lambdaf)*filter(retc^2, filter=lambdaf, 
//'   init=as.numeric(retc[1, 1])^2/(1-lambdaf), 
//'   method="recursive")
//' all.equal(vars, unclass(retc2), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::run_var(retp, lambda=lambdaf),
//'   Rcode=filter(retc^2, filter=lambdaf, init=as.numeric(retc[1, 1])^2/(1-lambdaf), method="recursive"),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_var(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::uword ncols = tseries.n_cols;
  arma::mat meanm = arma::zeros(nrows, ncols);
  arma::mat vars = arma::zeros(nrows, ncols);
  double lambda1 = 1-lambda;
  
  // Perform loop over the rows
  meanm.row(0) = tseries.row(0);
  vars.row(0) = arma::square(tseries.row(1) - tseries.row(0));
  // vars.row(0) = 0;
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the means using the decay factor
    meanm.row(it) = lambda*meanm.row(it-1) + lambda1*tseries.row(it);
    // Variance is the weighted sum of the past variance and the square of the data minus its mean
    vars.row(it) = lambda*vars.row(it-1) + lambda1*arma::square(tseries.row(it) - meanm.row(it));
  }  // end for
  
  return vars;
  
}  // end run_var



////////////////////////////////////////////////////////////
//' Calculate the trailing means, volatilities, and z-scores of a streaming
//' \emph{time series} of data using an online recursive formula.
//' 
//' @param \code{tseries} A single \emph{time series} or a single column
//'   \emph{matrix} of data.
//' 
//' @param \code{lambda} A decay factor which multiplies past estimates.
//'   
//' @return A \emph{matrix} with three columns (means, volatilities, and
//'   z-scores) and the same number of rows as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{run_zscores()} calculates the trailing means,
//'   volatilities, and z-scores of a single streaming \emph{time series} of
//'   data \eqn{r_t}, by recursively weighting the past variance estimates
//'   \eqn{\sigma^2_{t-1}}, with the squared differences of the data minus its
//'   trailing means \eqn{(r_t - \bar{r}_t)^2}, using the decay factor
//'   \eqn{\lambda}:
//'   \deqn{
//'     \bar{r}_t = \lambda \bar{r}_{t-1} + (1-\lambda) r_t
//'   }
//'   \deqn{
//'     \sigma^2_t = \lambda \sigma^2_{t-1} + (1-\lambda) (r_t - \bar{r}_t)^2
//'   }
//'   \deqn{
//'     z_t = \frac{r_t - \bar{r}_t}{\sigma_t}
//'   }
//'   Where \eqn{r_t} are the streaming data, \eqn{\bar{r}_t} are the trailing
//'   means, \eqn{\sigma^2_t} are the trailing variance estimates, and \eqn{z_t}
//'   are the z-scores.
//'   
//'   The above online recursive formulas are convenient for processing live
//'   streaming data because they don't require maintaining a buffer of past
//'   data.
//'   The formulas are equivalent to a convolution with exponentially decaying
//'   weights, but they're much faster to calculate.
//'   Using exponentially decaying weights is more natural than using a sliding
//'   look-back interval, because it gradually "forgets" about the past data.
//' 
//'   The value of the decay factor \eqn{\lambda} must be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the trailing variance values have a
//'   stronger dependence on past data.  This is equivalent to a long
//'   look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the trailing variance values have a
//'   weaker dependence on past data.  This is equivalent to a short look-back
//'   interval.
//' 
//'   The function \code{run_zscores()} returns a \emph{matrix} with three
//'   columns (means, volatilities, and z-scores) and the same number of rows as
//'   the input argument \code{tseries}.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical VTI log prices
//' pricev <- log(na.omit(rutils::etfenv$prices$VTI))
//' # Calculate the trailing variance and z-scores of prices
//' lambdaf <- 0.9 # Decay factor
//' zscores <- HighFreq::run_zscores(pricev, lambda=lambdaf)
//' datav <- cbind(pricev, zscores[, 3])
//' colnames(datav) <- c("VTI", "Z-Scores")
//' colnamev <- colnames(datav)
//' dygraphs::dygraph(datav, main="VTI Prices and Z-scores") %>%
//'    dyAxis("y", label=colnamev[1], independentTicks=TRUE) %>%
//'    dyAxis("y2", label=colnamev[2], independentTicks=TRUE) %>%
//'    dySeries(axis="y", label=colnamev[1], strokeWidth=2, col="blue") %>%
//'    dySeries(axis="y2", label=colnamev[2], strokeWidth=2, col="red") %>%
//'    dyLegend(show="always", width=300)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_zscores(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  // arma::uword ncols = tseries.n_cols;
  arma::mat meanm = arma::zeros(nrows, 1);
  arma::mat vars = arma::zeros(nrows, 1);
  double lambda1 = 1-lambda;
  
  // Perform loop over the rows
  meanm(0) = tseries(0);
  vars(0) = pow(tseries(1) - tseries(0), 2);
  // vars(0) = 0;
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the means using the decay factor
    meanm(it) = lambda*meanm(it-1) + lambda1*tseries(it);
    // Variance is the weighted sum of the past variance and the square of the data minus its mean
    vars(it) = lambda*vars(it-1) + lambda1*pow(tseries(it) - meanm(it), 2);
  }  // end for
  
  vars = sqrt(vars);
  arma::mat zscores = (tseries - meanm)/vars;
  
  return arma::join_rows(meanm, vars, zscores);
  
}  // end run_zscores



////////////////////////////////////////////////////////////
//' Calculate the trailing variance of streaming \emph{OHLC} price data using an
//' online recursive formula.
//' 
//' @param \code{ohlc} A \emph{time series} or a \emph{matrix} with \emph{OHLC}
//'   price data.
//'   
//' @param \code{lambda} A decay factor which multiplies past estimates.
//'
//' @return A single-column \emph{matrix} of variance estimates, with the same
//'   number of rows as the input \code{ohlc} price data.
//'
//' @details
//'   The function \code{run_var_ohlc()} calculates a single-column
//'   \emph{matrix} of variance estimates of streaming \emph{OHLC} price data.
//'   
//'   The function \code{run_var_ohlc()} calculates the variance from the
//'   differences between the \emph{Open}, \emph{High}, \emph{Low}, and
//'   \emph{Close} prices, using the \emph{Yang-Zhang} range volatility
//'   estimator:
//'   \deqn{
//'     \sigma^2_t = (1-\lambda) ((O_t - C_{t-1})^2 + 0.134 (C_t - O_t)^2 + 
//'     0.866 ((H_i - O_i) (H_i - C_i) + (L_i - O_i) (L_i - C_i))) + 
//'     \lambda \sigma^2_{t-1}
//'   }
//'   It recursively weighs the current variance estimate with the past
//'   estimates \eqn{\sigma^2_{t-1}}, using the decay factor \eqn{\lambda}.
//'
//'   The above recursive formula is convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//'   The formula is equivalent to a convolution with exponentially decaying
//'   weights, but it's faster.
//'   Using exponentially decaying weights is more natural than using a sliding
//'   look-back interval, because it gradually "forgets" about the past data.
//'   
//'   The function \code{run_var_ohlc()} does not calculate the logarithm of
//'   the prices.
//'   So if the argument \code{ohlc} contains dollar prices then
//'   \code{run_var_ohlc()} calculates the dollar variance.
//'   If the argument \code{ohlc} contains the log prices then
//'   \code{run_var_ohlc()} calculates the percentage variance.
//'   
//'   The function \code{run_var_ohlc()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, which makes it several times faster than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Extract the log OHLC prices of VTI
//' ohlc <- log(rutils::etfenv$VTI)
//' # Calculate the trailing variance
//' vart <- HighFreq::run_var_ohlc(ohlc, lambda=0.8)
//' # Calculate the rolling variance
//' varoll <- HighFreq::roll_var_ohlc(ohlc, lookb=5, method="yang_zhang", scale=FALSE)
//' datav <- cbind(vart, varoll)
//' colnames(datav) <- c("trailing", "rolling")
//' colnamev <- colnames(datav)
//' datav <- xts::xts(datav, index(ohlc))
//' # dygraph plot of VTI trailing versus rolling volatility
//' dygraphs::dygraph(sqrt(datav[-(1:111), ]), main="Trailing and Rolling Volatility of VTI") %>%
//'   dyOptions(colors=c("red", "blue"), strokeWidth=2) %>%
//'   dyLegend(show="always", width=300)
//' # Compare the speed of trailing versus rolling volatility
//' library(microbenchmark)
//' summary(microbenchmark(
//'   trailing=HighFreq::run_var_ohlc(ohlc, lambda=0.8),
//'   rolling=HighFreq::roll_var_ohlc(ohlc, lookb=5, method="yang_zhang", scale=FALSE),
//'   times=10))[, c(1, 4, 5)]
//' }
//' @export
// [[Rcpp::export]]
arma::mat run_var_ohlc(const arma::mat& ohlc, 
                       double lambda) {
  
  // Allocate variance matrix
  arma::uword nrows = ohlc.n_rows;
  arma::mat vars = arma::zeros(nrows, 1);
  double lambda1 = 1-lambda;
  double coeff = 0.134;

  // Calculate all the different intra-day and day-over-day returns 
  // (differences of OHLC prices)
  arma::mat openp = ohlc.col(0);
  arma::mat highp = ohlc.col(1);
  arma::mat lowp = ohlc.col(2);
  arma::mat closep = ohlc.col(3);
  arma::mat opcl(closep.n_rows, 1);
  opcl = (openp - lagit(closep, 1, false));
  arma::mat clop = (closep - openp);
  arma::mat clhi = (closep - highp);
  arma::mat cllow = (closep - lowp);
  arma::mat hilow = (highp - lowp);
  arma::mat hiop = (highp - openp);
  arma::mat lowop = (lowp - openp);
  
  // Perform loop over the rows
  vars.row(0) = arma::square(opcl.row(0)) + coeff*arma::square(clop.row(0)) +
    (coeff-1)*(clhi.row(0)*hiop.row(0) + cllow.row(0)*lowop.row(0));
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the variance as the weighted sum of squared returns minus the squared means
    vars.row(it) = lambda1*(arma::square(opcl.row(it)) + coeff*arma::square(clop.row(it)) +
      (coeff-1)*(clhi.row(it)*hiop.row(it) + cllow.row(it)*lowop.row(it))) + lambda*vars.row(it-1);
  }  // end for
  
  return vars;
  
}  // end run_var_ohlc



////////////////////////////////////////////////////////////
//' Calculate the correlation matrix from the covariance matrix.
//' 
//' @param \code{covmat} A \emph{matrix} of covariances.
//' 
//' @return Void (no return value - modifies the covariance matrix in place).
//' 
//' @details
//'   The function \code{push_cov2cor()} calculates the correlation matrix from
//'   the covariance matrix, in place, without copying the data in memory.
//'   
//'   The function \code{push_cov2cor()} accepts a \emph{pointer} to the
//'   covariance matrix, and it overwrites it with the correlation matrix.
//'
//'   The function \code{push_cov2cor()} is written in \code{RcppArmadillo}
//'   \code{C++} so it's much faster than \code{R} code.
//'   
//' @examples
//' \dontrun{
//' # Calculate a time series of returns
//' retp <- na.omit(rutils::etfenv$returns[, c("IEF", "VTI", "DBC")])
//' # Calculate the covariance matrix of returns
//' covmat <- cov(retp)
//' # Calculate the correlation matrix of returns
//' push_cov2cor(covmat)
//' all.equal(covmat, cor(retp))
//' }
//' 
//' @export
// [[Rcpp::export]]
void push_cov2cor(arma::mat& covmat) {
  
  arma::vec volv = arma::sqrt(covmat.diag());
  covmat = covmat/(volv*arma::trans(volv));
  
}  // end push_cov2cor



////////////////////////////////////////////////////////////
//' Update the trailing covariance matrix of streaming asset returns,
//' with a row of new returns using an online recursive formula.
//' 
//' @param \code{retsn} A \emph{vector} of new asset returns.
//' 
//' @param \code{covmat} A trailing covariance \emph{matrix} of asset returns.
//' 
//' @param \code{meanv} A \emph{vector} of trailing means of asset returns.
//' 
//' @param \code{lambdacov} A decay factor which multiplies the past covariance.
//' 
//' @return Void (no return value - modifies the trailing covariance matrix
//'   and the return means in place).
//' 
//' @details
//'   The function \code{push_covar()} updates the trailing covariance matrix of
//'   streaming asset returns, with a row of new returns.  It updates the
//'   covariance matrix in place, without copying the data in memory.
//'   
//'   The streaming asset returns \eqn{r_t} contain multiple columns and the
//'   parameter \code{retsn} represents a single row of \eqn{r_t} - the asset
//'   returns at time \eqn{t}.  The elements of the vectors \code{retsn} and 
//'   \code{meanv} represent single rows of data with multiple columns.
//'   
//'   The function \code{push_covar()} accepts \emph{pointers} to the arguments
//'   \code{covmat} and \code{meanv}, 
//'   and it overwrites the old values with the new values. It performs the
//'   calculation in place, without copying the data in memory, which can
//'   significantly increase the computation speed for large matrices.
//'
//'   First, the function \code{push_covar()} updates the trailing means
//'   \eqn{\bar{r}_t} of the streaming asset returns \eqn{r_t} by recursively
//'   weighting present and past values using the decay factor \eqn{\lambda}:
//'   \deqn{
//'     \bar{r}_t = \lambda \bar{r}_{t-1} + (1-\lambda) r_t
//'   }
//'   This recursive formula is equivalent to the exponentially weighted moving
//'   average of the streaming asset returns \eqn{r_t}.
//'
//'   It then calculates the demeaned returns:
//'   \deqn{
//'     \hat{r}_t = r_t - \bar{r}_t
//'   }
//'   
//'   Finally, it updates the trailing covariance matrix of the returns:
//'   \deqn{
//'     {cov}_t = \lambda {cov}_{t-1} + (1-\lambda) \hat{r}^T_t \hat{r}_t
//'   }
//'   
//'   The decay factor \eqn{\lambda} determines the strength of the updates,
//'   with smaller \eqn{\lambda} values giving more weight to the new data. If
//'   the asset returns are not stationary, then applying more weight to the new
//'   returns reduces the bias of the trailing covariance matrix, but it also
//'   increases its variance. Simulation can be used to find the value of the
//'   \eqn{\lambda} parameter to achieve the best bias-variance tradeoff.
//'   
//'   The function \code{push_covar()} is written in \code{RcppArmadillo}
//'   \code{C++} so it's much faster than \code{R} code.
//'   
//' @examples
//' \dontrun{
//' # Calculate a time series of returns
//' retp <- na.omit(rutils::etfenv$returns[, c("IEF", "VTI", "DBC")])
//' # Calculate the returns without last row
//' nrows <- NROW(retp)
//' retss <- retp[-nrows]
//' # Calculate the covariance of returns
//' meanv <- colMeans(retss)
//' covmat <- cov(retss)
//' # Update the covariance of returns
//' HighFreq::push_covar(retsn=retp[nrows], covmat=covmat, meanv=meanv, lambdacov=0.9)
//' }
//' 
//' @export
// [[Rcpp::export]]
void push_covar(const arma::rowvec& retsn, // Row of new asset returns
                arma::mat& covmat,  // Covariance matrix
                arma::rowvec& meanv, // Trailing means of the returns
                const double& lambdacov) { // Covariance decay factor
  
  double lambda1 = 1-lambdacov;
  
  // Update the means of the returns
  meanv = lambdacov*meanv + lambda1*retsn;
  // Calculate the centered returns
  arma::rowvec datav = (retsn - meanv);
  
  // Update the covariance of the returns
  covmat = lambdacov*covmat + lambda1*arma::trans(datav)*datav;
  
}  // end push_covar



////////////////////////////////////////////////////////////
//' Update the trailing eigen values and eigen vectors of streaming asset return
//' data, with a row of new returns.
//' 
//' @param \code{retsn} A \emph{vector} of new asset returns.
//' 
//' @param \code{covmat} A trailing covariance \emph{matrix} of asset returns.
//' 
//' @param \code{eigenval} A \emph{vector} of eigen values.
//' 
//' @param \code{eigenvec} A \emph{matrix} of eigen vectors.
//' 
//' @param \code{eigenret} A \emph{vector} of eigen portfolio returns.
//' 
//' @param \code{meanv} A \emph{vector} of trailing means of asset returns.
//' 
//' @param \code{lambdacov} A decay factor which multiplies the past covariance.
//' 
//' @return Void (no return value - modifies the trailing eigen values, eigen
//'   vectors, the eigen portfolio returns, and the return means in place).
//' 
//' @details
//'   The function \code{push_eigen()} updates the trailing eigen values, eigen
//'   vectors, and the eigen portfolio returns of streaming asset returns, with
//'   a row of new data.  It updates the eigenelements in place, without copying
//'   the data in memory.
//'   
//'   The streaming asset returns \eqn{r_t} contain multiple columns and the
//'   parameter \code{retsn} represents a single row of \eqn{r_t} - the asset
//'   returns at time \eqn{t}.  The elements of the vectors \code{retsn},
//'   \code{eigenret}, and \code{meanv} represent single rows of data with
//'   multiple columns.
//'   
//'   The function \code{push_eigen()} accepts \emph{pointers} to the arguments
//'   \code{eigenval}, \code{eigenval}, \code{eigenvec}, \code{meanv}, and
//'   \code{eigenret}, and it overwrites the old values with the new values. It
//'   performs the calculation in place, without copying the data in memory,
//'   which can significantly increase the computation speed for large matrices.
//'
//'   First, the function \code{push_eigen()} calls the function
//'   \code{HighFreq::push_covar()} to update the trailing covariance matrix of
//'   streaming asset returns, with a row of new returns.  It updates the
//'   covariance matrix in place, without copying the data in memory.
//'   
//'   It then calls the \code{Armadillo} function \code{arma::eig_sym} to
//'   calculate the eigen decomposition of the trailing covariance matrix.
//'   
//'   The function \code{push_eigen()} calculates the eigen portfolio returns by
//'   multiplying the scaled asset returns times the eigen vectors
//'   \eqn{\strong{v}_{t-1}}:
//'   \deqn{
//'     r^{eigen}_t = \strong{v}_{t-1} \frac{r_t}{\sigma_{t-1}}
//'   }
//'   Where \eqn{\strong{v}_{t-1}} is the matrix of previous eigen vectors that
//'   are passed by reference through the parameter \code{eigenvec}. The eigen
//'   returns \eqn{r^{eigen}_t} are the returns of the eigen portfolios, with
//'   weights equal to the eigen vectors \eqn{\strong{v}_{t-1}}. The eigen
//'   weights are applied to the asset returns scaled by their volatilities.
//'   The eigen returns \eqn{r^{eigen}_t} are passed by reference through the
//'   parameter \code{eigenret}. 
//'   
//'   The decay factor \eqn{\lambda} determines the strength of the updates,
//'   with smaller \eqn{\lambda} values giving more weight to the new data. If
//'   the asset returns are not stationary, then applying more weight to the new
//'   returns reduces the bias of the trailing covariance matrix, but it also
//'   increases its variance. Simulation can be used to find the value of the
//'   \eqn{\lambda} parameter to achieve the best bias-variance tradeoff.
//'   
//'   The function \code{push_eigen()} is written in \code{RcppArmadillo}
//'   \code{C++} so it's much faster than \code{R} code.
//'   
//' @examples
//' \dontrun{
//' # Calculate a time series of returns
//' retp <- na.omit(rutils::etfenv$returns[, c("IEF", "VTI", "DBC")])
//' # Calculate the returns without the last row
//' nrows <- NROW(retp)
//' retss <- retp[-nrows]
//' # Calculate the previous covariance of returns
//' meanv <- colMeans(retss)
//' covmat <- cov(retss)
//' # Update the covariance of returns
//' eigenret <- numeric(NCOL(retp))
//' HighFreq::push_eigen(retsn=retp[nrows], covmat=covmat, 
//'   eigenval=eigenval, eigenvec=eigenvec, 
//'   eigenret=eigenret, meanv=meanv, lambdacov=0.9)
//' }
//' 
//' @export
// [[Rcpp::export]]
void push_eigen(const arma::rowvec& retsn, // Row of new asset returns
                arma::mat& covmat,  // Covariance matrix
                arma::vec& eigenval, // Eigen values
                arma::mat& eigenvec, // Eigen vectors
                arma::rowvec& eigenret, // Row of eigen portfolio returns
                arma::rowvec& meanv, // Trailing means of the returns
                const double& lambdacov) { // Covariance decay factor
  
  // Scale the returns by their volatility
  arma::rowvec varv = arma::trans(covmat.diag());
  varv.replace(0, 1);
  arma::rowvec retsc = retsn/arma::sqrt(varv);
  // Calculate the eigen portfolio returns - the products of the previous eigen vectors times the scaled returns
  eigenret = retsc*eigenvec;
  // Update the covariance matrix
  push_covar(retsc, covmat, meanv, lambdacov);
  // Calculate the eigen decomposition
  arma::eig_sym(eigenval, eigenvec, covmat);
  
}  // end push_eigen



////////////////////////////////////////////////////////////
//' Update the trailing eigen values and eigen vectors of streaming asset return
//' data, with a row of new returns, using the \emph{SGA} algorithm.
//' 
//' @param \code{retsn} A \emph{vector} of new asset returns.
//' 
//' @param \code{eigenval} A \emph{vector} of eigen values.
//' 
//' @param \code{eigenvec} A \emph{matrix} of eigen vectors.
//' 
//' @param \code{eigenret} A \emph{vector} of eigen portfolio returns.
//' 
//' @param \code{meanv} A \emph{vector} of trailing means of asset returns.
//' 
//' @param \code{varv} A \emph{vector} of the trailing asset variances.
//' 
//' @param \code{lambda} A decay factor which multiplies the past mean and
//'   variance.
//' 
//' @param \code{gamma} A \emph{numeric} gain factor which multiplies the past
//'   eigenelements.
//' 
//' @return Void (no return value - modifies the trailing eigen values, eigen
//'   vectors, the return means, and the return variances in place).
//' 
//' @details
//'   The function \code{push_sga()} updates the trailing eigen values, eigen
//'   vectors, and the eigen portfolio returns of streaming asset returns, with
//'   a row of new data, using the \emph{SGA} algorithm. It updates the
//'   eigenelements in place, without copying the data in memory.
//'   
//'   The streaming asset returns \eqn{r_t} contain multiple columns and the
//'   parameter \code{retsn} represents a single row of \eqn{r_t} - the asset
//'   returns at time \eqn{t}.  The elements of the vectors \code{retsn},
//'   \code{meanv}, and \code{varv} represent single rows of data with multiple
//'   columns.
//'   
//'   The function \code{push_sga()} accepts \emph{pointers} to the arguments
//'   \code{eigenval}, \code{eigenvec}, \code{meanv}, and \code{varv},
//'   and it overwrites the old values with the new values. It performs the
//'   calculation in place, without copying the data in memory, which can
//'   significantly increase the computation speed for large matrices.
//'
//'   First, the function \code{push_sga()} updates the trailing means
//'   \eqn{\bar{r}_t} and variances \eqn{\sigma^2_t} of the streaming asset
//'   returns \eqn{r_t} by recursively weighting present and past values
//'   using the decay factor \eqn{\lambda}:
//'   \deqn{
//'     \bar{r}_t = \lambda \bar{r}_{t-1} + (1-\lambda) r_t
//'   }
//'   \deqn{
//'     \sigma^2_t = \lambda \sigma^2_{t-1} + (1-\lambda) (r_t - \bar{r}_t)^2
//'   }
//'   The past values \eqn{\bar{r}_{t-1}} and \eqn{\sigma^2_{t-1}} are passed in
//'   by reference through the variables \code{meanv} and \code{varv}. The
//'   updated values are then passed out by reference.
//'
//'   These recursive formulas are equivalent to the exponentially weighted
//'   moving averages of the streaming asset returns \eqn{r_t}.
//'
//'   It then calculates a vector of the eigen portfolio returns:
//'   \deqn{
//'     r^{eigen}_t = \strong{v}_{t-1} \frac{r_t}{\sigma_{t-1}}
//'   }
//'   Where \eqn{\strong{v}_{t-1}} is the matrix of previous eigen vectors that
//'   are passed by reference through the parameter \code{eigenvec}. The eigen
//'   returns \eqn{r^{eigen}_t} are the returns of the eigen portfolios, with
//'   weights equal to the eigen vectors \eqn{\strong{v}_{t-1}}. The eigen
//'   weights are applied to the asset returns scaled by their volatilities.
//'   The eigen returns \eqn{r^{eigen}_t} are passed by reference through the
//'   parameter \code{eigenret}. 
//'   
//'   The function \code{push_sga()} then standardizes the columns of the new
//'   returns:
//'   \deqn{
//'     \hat{r}_t = \frac{r_t - \bar{r}_t}{\sigma_t}
//'   }
//'   
//'   Finally, the vector of eigen values \eqn{\Lambda_{j, t}} and the matrix of
//'   eigen vectors \eqn{\strong{v}_{j, t}} (\eqn{j} is the column index) are
//'   then updated using the \emph{SGA} algorithm:
//'   \deqn{
//'     \Lambda_{j, t} = (1-\gamma) \Lambda_{j, t-1} + \gamma \phi_{j, t-1}
//'   }
//'   \deqn{
//'     \strong{v}_{j, t} = \strong{v}_{j, t-1} + \gamma \phi_{j, t-1} (\hat{r}_{t} - \phi_{j, t-1} \strong{v}_{j, t-1} - 2 \sum_{i=1}^{j-1} \phi_{i, t-1} \strong{v}_{i, t-1})
//'   }
//'   Where \eqn{\phi_{j, t-1} = \hat{r}_{t} \strong{v}_{j, t-1}} are the matrix
//'   products of the new data times the previous eigen vectors. 
//'   
//'   The gain factor \eqn{\gamma} determines the strength of the updates, with
//'   larger \eqn{\gamma} values giving more weight to the new data. If the
//'   asset returns are not stationary, then applying more weight to the new
//'   returns reduces the bias of the trailing eigen vectors, but it also
//'   increases their variance. Simulation can be used to find the value of the
//'   \eqn{\gamma} parameter to achieve the best bias-variance tradeoff.
//'   
//'   A description of the \emph{SGA} algorithm can be found in the package
//'   \href{https://cran.r-project.org/web/packages/onlinePCA/index.html}{onlinePCA} and in the 
//'   \href{https://paperswithcode.com/paper/online-principal-component-analysis-in-high}{Online PCA paper}.
//'   
//'   The function \code{push_sga()} is written in \code{RcppArmadillo}
//'   \code{C++} code and it calls the \code{Armadillo} function
//'   \code{arma::qr_econ()} to perform the QR decomposition, to calculate the
//'   eigen vectors.
//'   
//' @examples
//' \dontrun{
//' # Calculate a time series of returns
//' retp <- na.omit(rutils::etfenv$returns[, c("IEF", "VTI", "DBC")])
//' # Calculate the covariance of returns without the last row
//' nrows <- NROW(retp)
//' retss <- retp[-nrows]
//' HighFreq::calc_scale(retss)
//' meanv <- colMeans(retss)
//' varv <- sapply(retss, var)
//' covmat <- cov(retss)
//' ncols <- NCOL(retss)
//' # Calculate the eigen decomposition using RcppArmadillo
//' eigenval <- numeric(ncols) # Allocate eigen values
//' eigenvec <- matrix(numeric(ncols^2), nc=ncols) # Allocate eigen vectors
//' HighFreq::calc_eigen(covmat, eigenval, eigenvec)
//' # Update the eigen decomposition using SGA
//' eigenret <- numeric(NCOL(retp))
//' HighFreq::push_sga(retsn=retp[nrows], 
//'   eigenval=eigenval, eigenvec=eigenvec, 
//'   eigenret=eigenret, meanv=meanv, varv=varv, lambda=0.9, gamma=0.1)
//' }
//' 
//' @export
// [[Rcpp::export]]
void push_sga(const arma::rowvec& retsn, // Row of new asset returns
              arma::rowvec& eigenval, // Eigen values
              arma::mat& eigenvec, // Eigen vectors
              arma::rowvec& eigenret, // Row of eigen portfolio returns
              arma::rowvec& meanv, // Trailing means of the returns
              arma::rowvec& varv, // Trailing variances of the returns
              const double& lambda, // Decay factor which multiplies the past mean and variance
              const double& gamma) { // Gain factor which multiplies the past eigenelements
  
  double lambda1 = 1-lambda;
  
  // Calculate the eigen portfolio returns - the products of the previous eigen vectors times the scaled returns
  arma::rowvec volv = arma::sqrt(varv);
  eigenret = (retsn/volv)*eigenvec;
  
  // Update the mean and variance of the returns
  meanv = lambda*meanv + lambda1*retsn;
  varv = lambda*varv + lambda1*arma::square(retsn-meanv);
  // Calculate the standardized returns
  volv = arma::sqrt(varv);
  arma::rowvec datav = (retsn-meanv)/volv;
  // arma::rowvec datav = (retsn-meanv);
  // std::cout << "datav: " << std::endl << datav << std::endl;
  // std::cout << "arma::trans(datav): " << std::endl << arma::trans(datav) << std::endl;
  
  // Calculate the phis - the products of the eigen vectors times the standardized returns
  // arma::rowvec vecd = arma::trans(eigenvec)*arma::trans(datav);
  arma::rowvec vecd = datav*eigenvec;
  // std::cout << "vecd: " << std::endl << arma::trans(vecd) << std::endl;
  // std::cout << "vecd: " << std::endl << vecd << std::endl;
  // Calculate the updated eigen values
  arma::rowvec neigenval = (1-gamma)*eigenval + gamma*arma::square(vecd);
  // std::cout << "eigenval: " << std::endl << eigenval << std::endl;
  // std::cout << "neigenval: " << std::endl << neigenval << std::endl;
  
  // Calculate diagonal matrix of gain factors
  arma::uword ncols = eigenvec.n_rows;
  arma::mat dgamma = arma::mat(ncols, ncols, arma::fill::zeros);
  dgamma.diag().fill(gamma);
  
  // Perform QR decomposition of Q to get the eigen vectors W
  arma::mat W, R, Q = eigenvec;
  Q += arma::trans(datav)*vecd*dgamma;
  // std::cout << "Q: " << std::endl << Q << std::endl;
  arma::qr_econ(W, R, Q);
  
  // Sort the eigen values and eigen vectors
  arma::uvec sorti = arma::sort_index(neigenval, "descend");
  // std::cout << "sorti: " << std::endl << sorti << std::endl;
  eigenval = neigenval.cols(sorti);
  eigenvec = W.cols(sorti);
  
}  // end push_sga



////////////////////////////////////////////////////////////
//' Calculate the trailing covariances of two streaming \emph{time series} of
//' returns using an online recursive formula.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} with two
//'   columns of returns data.
//' 
//' @param \code{lambda} A decay factor which multiplies past estimates.
//'   
//' @return A \emph{matrix} with five columns of data: the trailing covariances,
//'   the variances, and the mean values of the two columns of the argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{run_covar()} calculates the trailing covariances of two
//'   streaming \emph{time series} of returns, by recursively weighting the past
//'   covariance estimates \eqn{{cov}_{t-1}}, with the products of their
//'   returns minus their means, using the decay factor \eqn{\lambda}:
//'   \deqn{
//'     \bar{x}_t = \lambda \bar{x}_{t-1} + (1-\lambda) x_t
//'   }
//'   \deqn{
//'     \bar{y}_t = \lambda \bar{y}_{t-1} + (1-\lambda) y_t
//'   }
//'   \deqn{
//'     \sigma^2_{x t} = \lambda \sigma^2_{x t-1} + (1-\lambda) (x_t - \bar{x}_t)^2
//'   }
//'   \deqn{
//'     \sigma^2_{y t} = \lambda \sigma^2_{y t-1} + (1-\lambda) (y_t - \bar{y}_t)^2
//'   }
//'   \deqn{
//'     {cov}_t = \lambda {cov}_{t-1} + (1-\lambda) (x_t - \bar{x}_t) (y_t - \bar{y}_t)
//'   }
//'   Where \eqn{{cov}_t} is the trailing covariance estimate at time \eqn{t},
//'   \eqn{\sigma^2_{x t}}, \eqn{\sigma^2_{y t}}, \eqn{\bar{x}_t} and
//'   \eqn{\bar{x}_t} are the trailing variances and means of the returns, and
//'   \eqn{x_t} and \eqn{y_t} are the two streaming returns data.
//' 
//'   The above online recursive formulas are convenient for processing live
//'   streaming data because they don't require maintaining a buffer of past
//'   data. The formulas are equivalent to a convolution with exponentially
//'   decaying weights, but they're much faster to calculate.
//'   Using exponentially decaying weights is more natural than using a sliding
//'   look-back interval, because it gradually "forgets" about the past data.
//' 
//'   The value of the decay factor \eqn{\lambda} must be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the trailing covariance values have a
//'   stronger dependence on past data.  This is equivalent to a long
//'   look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the trailing covariance values have
//'   a weaker dependence on past data.  This is equivalent to a short
//'   look-back interval.
//' 
//'   The function \code{run_covar()} returns five columns of data: the trailing 
//'   covariances, the variances, and the mean values of the two columns of the
//'   argument \code{tseries}.  This allows calculating the trailing
//'   correlations, betas, and alphas.
//' 
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' retp <- zoo::coredata(na.omit(rutils::etfenv$returns[, c("IEF", "VTI")]))
//' # Calculate the trailing covariance
//' lambdaf <- 0.9 # Decay factor
//' covars <- HighFreq::run_covar(retp, lambda=lambdaf)
//' # Calculate the trailing correlation
//' correl <- covars[, 1]/sqrt(covars[, 2]*covars[, 3])
//' # Calculate the trailing covariance using R code
//' nrows <- NROW(retp)
//' retm <- matrix(numeric(2*nrows), nc=2)
//' retm[1, ] <- retp[1, ]
//' retd <- matrix(numeric(2*nrows), nc=2)
//' covarr <- numeric(nrows)
//' covarr[1] <- retp[1, 1]*retp[1, 2]
//' for (it in 2:nrows) {
//'   retm[it, ] <- lambdaf*retm[it-1, ] + (1-lambdaf)*(retp[it, ])
//'   retd[it, ] <- retp[it, ] - retm[it, ]
//'   covarr[it] <- lambdaf*covarr[it-1] + (1-lambdaf)*retd[it, 1]*retd[it, 2]
//' } # end for
//' all.equal(covars[, 1], covarr, check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_covar(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::uword ncols = tseries.n_cols;
  arma::mat meanm(nrows, ncols);
  arma::mat meand(1, ncols); // Centered row
  arma::mat vars(nrows, ncols);
  arma::mat covar(nrows, 1);
  double lambda1 = 1-lambda;

  // Perform loop over the rows
  meanm.row(0) = tseries.row(0);
  vars.row(0) = arma::square(tseries.row(0));
  covar.row(0) = tseries(0, 0)*tseries(0, 1);
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the mean as the weighted sum
    meanm.row(it) = lambda*meanm.row(it-1) + lambda1*tseries.row(it);
    meand = tseries.row(it) - meanm.row(it);
    // Calculate the covariance as the weighted sum of products of returns
    vars.row(it) = lambda*vars.row(it-1) + lambda1*arma::square(meand);
    covar.row(it) = lambda*covar.row(it-1) + lambda1*(meand(0)*meand(1));
  }  // end for

  return arma::join_rows(covar, vars, meanm);

  // Slower code below - because push_covar() calls arma::trans()
  // 
  // arma::uword nrows = tseries.n_rows;
  // arma::uword ncols = tseries.n_cols;
  // arma::mat vars(nrows, ncols);
  // arma::mat covar(nrows, 1);
  // arma::rowvec meanv = tseries.row(0);
  // arma::mat covmat = arma::trans(tseries.row(0))*tseries.row(0);
  // 
  // // Copy the covariance data
  // vars.row(0) = arma::trans(covmat.diag());
  // covar.row(0) = covmat(0, 1);
  // 
  // // Perform loop over the rows
  // for (arma::uword it = 1; it < nrows; it++) {
  //   // Update the covariance matrix
  //   push_covar(tseries.row(it), covmat, meanv, lambda);
  //   // Copy the covariance data
  //   vars.row(it) = arma::trans(covmat.diag());
  //   covar.row(it) = covmat(0, 1);
  // }  // end for
  // 
  // return arma::join_rows(covar, vars);
  
}  // end run_covar



////////////////////////////////////////////////////////////
//' Calculate the trailing autocovariance of a \emph{time series} of returns
//' using an online recursive formula.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} with a single
//'   column of returns data.
//' 
//' @param \code{lambda} A decay factor which multiplies past estimates.
//'   
//' @param \code{lagg} An \emph{integer} equal to the number of periods to lag.
//'   (The default is \code{lagg = 1}.)
//'
//' @return A \emph{matrix} with three columns of data: the trailing
//'   autocovariances, the variances, and the mean values of the argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{run_autocovar()} calculates the trailing
//'   autocovariance of a streaming \emph{time series} of returns, by
//'   recursively weighting the past autocovariance estimates \eqn{{cov}_{t-1}},
//'   with the products of their returns minus their means, using the decay
//'   factor \eqn{\lambda}:
//'   \deqn{
//'     \bar{x}_t = \lambda \bar{x}_{t-1} + (1-\lambda) x_t
//'   }
//'   \deqn{
//'     \sigma^2_{t} = \lambda \sigma^2_{t-1} + (1-\lambda) (x_t - \bar{x}_t)^2
//'   }
//'   \deqn{
//'     {cov}_t = \lambda {cov}_{t-1} + (1-\lambda) (x_t - \bar{x}_t) (x_{t-l} - \bar{x}_{t-l})
//'   }
//'   Where \eqn{{cov}_t} is the trailing autocovariance estimate at time
//'   \eqn{t}, with \code{lagg=l}.
//'   And \eqn{\sigma^2_{t}} and \eqn{\bar{x}_t} are the trailing variances and
//'   means of the streaming data.
//' 
//'   The above online recursive formulas are convenient for processing live
//'   streaming data because they don't require maintaining a buffer of past
//'   data. The formulas are equivalent to a convolution with exponentially
//'   decaying weights, but they're much faster to calculate.
//'   Using exponentially decaying weights is more natural than using a sliding
//'   look-back interval, because it gradually "forgets" about the past data.
//' 
//'   The value of the decay factor \eqn{\lambda} must be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the trailing covariance values have a
//'   stronger dependence on past data.  This is equivalent to a long
//'   look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the trailing covariance values have
//'   a weaker dependence on past data.  This is equivalent to a short
//'   look-back interval.
//' 
//'   The function \code{run_autocovar()} returns three columns of data: the
//'   trailing autocovariances, the variances, and the mean values of the
//'   argument \code{tseries}.  This allows calculating the trailing
//'   autocorrelations.
//' 
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' retp <- zoo::coredata(na.omit(rutils::etfenv$returns$VTI))
//' # Calculate the trailing autocovariance
//' lambdaf <- 0.9 # Decay factor
//' lagg <- 3
//' covars <- HighFreq::run_autocovar(retp, lambda=lambdaf, lagg=lagg)
//' # Calculate the trailing autocorrelation
//' correl <- covars[, 1]/covars[, 2]
//' # Calculate the trailing autocovariance using R code
//' nrows <- NROW(retp)
//' retm <- numeric(nrows)
//' retm[1] <- retp[1, ]
//' retd <- numeric(nrows)
//' covarr <- numeric(nrows)
//' covarr[1] <- retp[1, ]^2
//' for (it in 2:nrows) {
//'   retm[it] <- lambdaf*retm[it-1] + (1-lambdaf)*(retp[it])
//'   retd[it] <- retp[it] - retm[it]
//'   covarr[it] <- lambdaf*covarr[it-1] + (1-lambdaf)*retd[it]*retd[max(it-lagg, 1)]
//' } # end for
//' all.equal(covarr, covars[, 1])
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_autocovar(const arma::mat& tseries, 
                        double lambda, 
                        arma::uword lagg = 1) {
 
 arma::uword nrows = tseries.n_rows;
 arma::uword ncols = tseries.n_cols;
 arma::mat meanm(nrows, ncols);
 arma::mat meand(nrows, ncols); // Centered series
 arma::mat vars(nrows, ncols);
 arma::mat covar(nrows, 1);
 double lambda1 = 1-lambda;
 
 // Initialize the first rows
 meanm.row(0) = tseries.row(0);
 vars.row(0) = arma::square(tseries.row(0));
 covar.row(0) = arma::square(tseries.row(0));
 
 // Warmup period
 if (lagg > 1) {
   for (arma::uword it = 1; it < lagg; it++) {
     // Calculate the mean as the weighted sum
     meanm.row(it) = lambda*meanm.row(it-1) + lambda1*tseries.row(it);
     meand.row(it) = tseries.row(it) - meanm.row(it);
     // Calculate the covariance as the weighted sum of products of returns
     vars.row(it) = lambda*vars.row(it-1) + lambda1*arma::square(meand.row(it));
     covar.row(it) = lambda*covar.row(it-1) + lambda1*(meand.row(it)*meand.row(0));
   }  // end for
 }  // end if
 
 // Perform loop over the remaining rows
 for (arma::uword it = lagg; it < nrows; it++) {
   // Calculate the mean as the weighted sum
   meanm.row(it) = lambda*meanm.row(it-1) + lambda1*tseries.row(it);
   meand.row(it) = tseries.row(it) - meanm.row(it);
   // Calculate the covariance as the weighted sum of products of returns
   vars.row(it) = lambda*vars.row(it-1) + lambda1*arma::square(meand.row(it));
   covar.row(it) = lambda*covar.row(it-1) + lambda1*(meand.row(it)*meand.row(it - lagg));
 }  // end for
 
 return arma::join_rows(covar, vars, meanm);
 
}  // end run_autocovar



////////////////////////////////////////////////////////////
//' Perform regressions on the streaming \emph{time series} of response and
//' predictor data, and calculate the regression coefficients, the residuals,
//' and the forecasts, using online recursive formulas.
//' 
//' @param \code{respv} A single-column \emph{time series} or a single-column
//'   \emph{matrix} of response data.
//' 
//' @param \code{predm} A \emph{time series} or a \emph{matrix} of predictor
//'   data.
//' 
//' @param \code{lambda} A decay factor which multiplies past estimates.
//'   
//' @param \code{controlv} A \emph{list} of model parameters (see Details).
//'
//'   
//' @return A \emph{matrix} with the regression coefficients, the residuals, and
//'   the forecasts (in that order - see details), with the same number of rows
//'   as the predictor argument \code{predm}.
//'
//' @details
//'   The function \code{run_reg()} performs regressions on the streaming \emph{time
//'   series} of response \eqn{r_t} and predictor \eqn{p_t} data:
//'   \deqn{
//'     r_t = \beta_t p_t + \epsilon_t
//'   }
//'   Where \eqn{\beta_t} are the trailing regression coefficients and
//'   \eqn{\epsilon_t} are the residuals.
//'   
//'   It recursively updates the covariance matrix \eqn{{cov}_t} between the
//'   response and the predictor data, and the covariance matrix
//'   \eqn{{cov}_{pt}} between the predictors, using the decay factor
//'   \eqn{\lambda}:
//'   \deqn{
//'     {cov}_t = \lambda {cov}_{t-1} + (1-\lambda) r^T_t p_t
//'   }
//'   \deqn{
//'     {cov}_{p t} = \lambda {cov}_{p (t-1)} + (1-\lambda) p^T_t p_t
//'   }
//'   
//'   It calculates the regression coefficients \eqn{\beta_t} as equal to the
//'   covariance matrix between the response and the predictor data
//'   \eqn{{cov}_t}, divided by the covariance matrix between the predictors
//'   \eqn{{cov}_{pt}}:
//'   \deqn{
//'     \beta_t = {cov}_t \, {cov}^{-1}_{p t}
//'   }
//'   
//'   It calculates the residuals \eqn{\epsilon_t} as the difference between the
//'   response \eqn{r_t} minus the fitted values \eqn{\beta_t p_t}:
//'   \deqn{
//'     \epsilon_t = r_t - \beta_t p_t
//'   }
//'   
//'   And the residual variance \eqn{\sigma^2_t} as:
//'   \deqn{
//'     \bar{\epsilon}_t = \lambda \bar{\epsilon}_{t-1} + (1-\lambda) \epsilon_t
//'   }
//'   \deqn{
//'     \sigma^2_t = \lambda \sigma^2_{t-1} + (1-\lambda) (\epsilon_t - \bar{\epsilon}_t)^2
//'   }
//' 
//'   It then calculates the regression forecasts \eqn{f_t}, as equal to the
//'   past regression coefficients \eqn{\beta_{t-1}} times the current predictor
//'   data \eqn{p_t}:
//'   \deqn{
//'     f_t = \beta_{t-1} p_t
//'   }
//' 
//'   It finally calculates the forecast errors as the difference between the
//'   response minus the regression forecasts: \eqn{r_t - f_t}.
//' 
//'   The coefficient matrix \eqn{\beta} and the residuals \eqn{\epsilon} have
//'   the same number of rows as the predictor argument \code{predm}.
//'
//'   The function \code{run_reg()} accepts a list of regression model
//'   parameters through the argument \code{controlv}.
//'   The argument \code{controlv} contains the parameters \code{regmod} and
//'   \code{residscale}.
//'   Below is a description of how these parameters work.
//'   The list of model parameters can be created using the function
//'   \code{param_reg()}.  
//'
//'   The number of regression coefficients is equal to the number of columns of
//'   the predictor matrix \code{n}.
//'   If the predictor matrix contains a unit intercept column then the first
//'   regression coefficient is equal to the alpha value \eqn{\alpha}.
//'
//'   If \code{regmod = "least_squares"} (the default) then it performs the
//'   standard least squares regression.  This is currently the only option.
//' 
//'   The \emph{residuals} and the the \emph{forecast errors} may be scaled by
//'   their volatilities to obtain the \emph{z-scores}. 
//'   The default is \code{residscale = "none"} - no scaling.
//'   If the argument \code{residscale = "scale"} then the \emph{residuals}
//'   \eqn{\epsilon_t} are divided by their volatilities \eqn{\sigma_t}
//'   without subtracting their means:
//'   \deqn{
//'     \epsilon_t = \frac{\epsilon_t}{\sigma_t}
//'   }
//'   If the argument \code{residscale = "standardize"} then the residual means
//'   \eqn{\bar{\epsilon}} are subtracted from the \emph{residuals}, and then
//'   they are divided by their volatilities \eqn{\sigma_t}:
//'   \deqn{
//'     \epsilon_t = \frac{\epsilon_t - \bar{\epsilon}}{\sigma_t}
//'   }
//'   Which are equal to the \emph{z-scores}.
//'   
//'   The \emph{forecast errors} are also scaled in the same way as the
//'   \emph{residuals}, according to the argument\code{residscale}.
//' 
//'   The above online recursive formulas are convenient for processing live
//'   streaming data because they don't require maintaining a buffer of past
//'   data.
//'   The above recursive formulas are equivalent to a convolution with
//'   exponentially decaying weights, but they're much faster to calculate.
//'   Using exponentially decaying weights is more natural than using a sliding
//'   look-back interval, because it gradually "forgets" about the past data.
//'
//'   The value of the decay factor \eqn{\lambda} must be in the range between
//'   \code{0} and \code{1}.
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, so the trailing values have a greater
//'   dependence on past data.  This is equivalent to a long look-back
//'   interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, so the trailing values have a weaker
//'   dependence on past data.  This is equivalent to a short look-back
//'   interval.
//' 
//'   The function \code{run_reg()} returns multiple columns of data, with the
//'   same number of rows as the predictor matrix \code{predm}. If the predictor
//'   matrix \code{predm} has \code{n} columns then \code{run_reg()} returns a
//'   matrix with \code{n+2} columns.
//'   The first \code{n} columns contain the regression coefficients (with the
//'   first column equal to the alpha value \eqn{\alpha}).
//'   The last \code{2} columns are the regression residuals and the forecast
//'   errors.
//' 
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' retp <- na.omit(rutils::etfenv$returns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' respv <- retp[, 1]
//' # Predictor matrix equals VTI and IEF returns
//' predm <- retp[, -1]
//' # Add unit intercept column to the predictor matrix
//' predm <- cbind(rep(1, NROW(predm)), predm)
//' # Calculate the trailing regressions
//' lambdaf <- 0.9 # Decay factor
//' # Create a list of regression parameters
//' controlv <- HighFreq::param_reg(residscale="standardize")
//' regs <- HighFreq::run_reg(respv=respv, predm=predm, lambda=lambda, controlv=controlv)
//' # Plot the trailing residuals
//' datav <- cbind(cumsum(respv), regs[, NCOL(regs)])
//' colnames(datav) <- c("XLF", "residuals")
//' colnamev <- colnames(datav)
//' dygraphs::dygraph(datav["2008/2009"], main="Residuals of XLF Versus VTI and IEF") %>%
//'   dyAxis("y", label=colnamev[1], independentTicks=TRUE) %>%
//'   dyAxis("y2", label=colnamev[2], independentTicks=TRUE) %>%
//'   dySeries(axis="y", strokeWidth=2, col="blue") %>%
//'   dySeries(axis="y2", strokeWidth=2, col="red") %>%
//'   dyLegend(show="always", width=300)
//' 
//' # Calculate the trailing regressions using R code
//' lambda1 <- (1-lambdaf)
//' respv <- zoo::coredata(respv)
//' predm <- zoo::coredata(predm)
//' nrows <- NROW(predm)
//' ncols <- NCOL(predm)
//' covrespred <- respv[1, ]*predm[1, ]
//' covpred <- outer(predm[1, ], predm[1, ])
//' betas <- matrix(numeric(nrows*ncols), nc=ncols)
//' betas[1, ] <- covrespred %*% MASS::ginv(covpred)
//' resids <- numeric(nrows)
//' residm <- 0
//' residv <- 0
//' for (it in 2:nrows) {
//'  covrespred <- lambdaf*covrespred + lambda1*respv[it, ]*predm[it, ]
//'  covpred <- lambdaf*covpred + lambda1*outer(predm[it, ], predm[it, ])
//'  betas[it, ] <- covrespred %*% MASS::ginv(covpred)
//'  resids[it] <- respv[it, ] - (betas[it, ] %*% predm[it, ])
//'  residm <- lambdaf*residm + lambda1*resids[it]
//'  residv <- lambdaf*residv + lambda1*(resids[it] - residm)^2
//'  resids[it] <- (resids[it] - residm)/sqrt(residv)
//' } # end for
//' # Compare values, excluding warmup period
//' all.equal(regs[-(1:1e3), ], cbind(betas, resids)[-(1:1e3), ], check.attributes=FALSE)
//'  
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_reg(const arma::mat& respv, // Response vector
                  const arma::mat& predm, // Predictor matrix
                  double lambda, // Decay factor which multiplies the past values
                  Rcpp::List controlv) { // List of regression model parameters

  // Add unit intercept column to the predictor matrix
  // bool intercept = Rcpp::as<int>(controlv["intercept"]);
  arma::uword nrows = predm.n_rows;
  // arma::mat predm = predm; // Predictor matrix with intercept column
  // if (intercept)
  //   predm = arma::join_rows(ones(nrows), predm);
  
  arma::uword ncols = predm.n_cols;
  arma::mat covrespred = arma::zeros(1, ncols); // Covariance between the response and predictor
  arma::mat covpred = arma::zeros(ncols, ncols); // Covariance between the predictors
  arma::mat betas = arma::ones(nrows, ncols); // Betas
  arma::mat fcasts = arma::zeros(nrows, 1); // Forecasts
  arma::mat fcastm = arma::zeros(nrows, 1); // Forecast error means
  arma::mat fcastv = arma::ones(nrows, 1); // Forecast error variance
  arma::mat resids = arma::zeros(nrows, 1); // Residuals
  arma::mat residm = arma::zeros(nrows, 1); // Residual means
  arma::mat residv = arma::ones(nrows, 1); // Residual variance
  double lambda1 = 1-lambda;
  
  // Initialize the variables
  // cout << "Initializing the variables" << endl;
  covrespred = respv.row(0)*predm.row(0);
  covpred = arma::trans(predm.row(0))*predm.row(0);
  betas.row(0) = covrespred*arma::inv(covpred);
  resids.row(0) = respv.row(0) - arma::dot(betas.row(0), predm.row(0));
  residm.row(0) = resids.row(0);
  residv.row(0) = arma::square(resids.row(0));
  fcasts.row(0) = resids.row(0);
  // fcastm.row(0) = fcasts.row(0);
  // fcastv.row(0) = arma::square(fcasts.row(0));
  // Perform loop over the rows
  for (arma::uword it = 1; it < nrows; it++) {
    // Update the covariance between the response and predictor
    // cout << "Updating the covariances: " << it << endl;
    covrespred = lambda*covrespred + lambda1*respv.row(it)*predm.row(it);
    covpred = lambda*covpred + lambda1*arma::trans(predm.row(it))*predm.row(it);
    // cout << "Calculating betas: " << it << endl;
    // Update the betas and alphas
    betas.row(it) = covrespred*arma::inv(covpred);
    // Calculate the betas and alphas approximately assuming predictors are orthogonal - old method - only slightly faster
    // residv.row(it) = lambda*residv.row(it-1) + lambda1*arma::square(predd);
    // betas.row(it) = lambda*betas.row(it-1) + lambda1*covrespred/residv.row(it);
    // Calculate the residuals
    // cout << "Calculating residuals: " << it << endl;
    resids.row(it) = respv.row(it) - arma::dot(betas.row(it), predm.row(it));
    // resids.row(it) = -arma::dot(betas.row(it), predm.row(it));
    // Calculate the mean and variance of the residuals
    residm.row(it) = lambda*residm.row(it-1) + lambda1*resids.row(it);
    residv.row(it) = lambda*residv.row(it-1) + lambda1*arma::square(resids.row(it) - residm.row(it));
    // cout << "Calculating forecasts: " << it << endl;
    // Calculate the forecasts
    // fcasts.row(it) = arma::dot(betas.row(it-1), predm.row(it));
    // Calculate the forecast errors
    fcasts.row(it) = respv.row(it) - arma::dot(betas.row(it-1), predm.row(it));
    // Calculate the mean and variance of the forecast errors
    fcastm.row(it) = lambda*fcastm.row(it-1) + lambda1*fcasts.row(it);
    fcastv.row(it) = lambda*fcastv.row(it-1) + lambda1*arma::square(fcasts.row(it) - fcastm.row(it));
  }  // end for
  
  // Type of residual scaling - the default is "none"
  std::string residscale = Rcpp::as<std::string>(controlv["residscale"]);

  // Scale the forecast errors and the residuals
  if (residscale == "scale") {
    // Divide the residuals by their volatility
    resids = resids/arma::sqrt(residv);
    fcasts = fcasts/arma::sqrt(fcastv);
  } else if (residscale == "standardize") {
    // Center the residuals and divide them by their volatility
    resids = (resids - residm)/arma::sqrt(residv);
    fcasts = (fcasts - fcastm)/arma::sqrt(fcastv);
  }  // end if
  
  return arma::join_rows(betas, resids, fcasts);
  
}  // end run_reg




////////////////////////////////////////////////////////////
// Functions for statistics
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
// Define C++ enum type for different methods of regularization,
// methodsfor calculating variance, skewness, kurtosis, covariance, regression, 
// and matrix inverse.
enum methodenum {moment, least_squares, quantile, nonparametric, regular, sharpem, 
              maxsharpe, maxsharpemed, minvarlin, minvarquad, kellym, robustm, 
              sumsq, sumone, voltarget, voleqw};

// Map string to C++ enum type for switch statement.
// This is needed because Rcpp can't map C++ enum type to R variable SEXP.
methodenum calc_method(std::string method) {
  if (method == "moment" || method == "m") 
    return methodenum::moment;
  else if (method == "least_squares" || method == "l")
    return methodenum::least_squares;
  else if (method == "quantile" || method == "q")
    return methodenum::quantile;
  else if (method == "nonparametric" || method == "n")
    return methodenum::nonparametric;
  else if (method == "regular")
    return methodenum::regular;
  else if (method == "sharpem")
    return methodenum::sharpem;
  else if (method == "maxsharpe")
    return methodenum::maxsharpe;
  else if (method == "maxsharpemed")
    return methodenum::maxsharpemed;
  else if (method == "minvarlin")
    return methodenum::minvarlin;
  else if (method == "minvarquad")
    return methodenum::minvarquad;
  else if (method == "kellym")
    return methodenum::kellym;
  else if (method == "robustm")
    return methodenum::robustm;
  else if (method == "sumsq")
    return methodenum::sumsq;
  else if (method == "voltarget")
    return methodenum::voltarget;
  else 
    return methodenum::moment;
}  // end calc_method



////////////////////////////////////////////////////////////
//' Calculate the mean (location) of the columns of a \emph{time series} or a
//' \emph{matrix} using \code{RcppArmadillo}.
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'
//' @param \code{method} A \emph{character string} specifying the type of the
//'   mean (location) model (the default is \code{method = "moment"} - see
//'   Details).
//'
//' @param \code{confl} The confidence level for calculating the quantiles of
//'   returns (the default is \code{confl = 0.75}).
//'
//' @return A single-row matrix with the mean (location) of the columns of
//'   \code{tseries}.
//'
//' @details
//'   The function \code{calc_mean()} calculates the mean (location) values of
//'   the columns of the \emph{time series} \code{tseries} using \code{C++}
//'   \code{RcppArmadillo} code.
//'
//'   If \code{method = "moment"} (the default) then \code{calc_mean()}
//'   calculates the location as the mean - the first moment of the data.
//'
//'   If \code{method = "quantile"} then it calculates the location
//'   \eqn{\bar{r}} as the average of the quantiles as follows:
//'   \deqn{
//'     \bar{r} = \frac{q_{\alpha} + q_{1-\alpha}}{2}
//'   }
//'   Where \eqn{\alpha} is the confidence level for calculating the quantiles
//'   (argument \code{confl}).
//'
//'   If \code{method = "nonparametric"} then it calculates the location as the
//'   median.
//'   
//'   The code examples below compare the function \code{calc_mean()} with the
//'   mean (location) calculated using \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' retp <- na.omit(rutils::etfenv$returns[, c("XLP", "VTI")])
//' # Calculate the column means in RcppArmadillo
//' HighFreq::calc_mean(retp)
//' # Calculate the column means in R
//' sapply(retp, mean)
//' # Compare the values
//' all.equal(drop(HighFreq::calc_mean(retp)), 
//'   sapply(retp, mean), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_mean(retp),
//'   Rcode=sapply(retp, mean),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the quantile mean (location)
//' HighFreq::calc_mean(retp, method="quantile", confl=0.9)
//' # Calculate the quantile mean (location) in R
//' colSums(sapply(retp, quantile, c(0.9, 0.1), type=5))
//' # Compare the values
//' all.equal(drop(HighFreq::calc_mean(retp, method="quantile", confl=0.9)), 
//'   colSums(sapply(retp, quantile, c(0.9, 0.1), type=5)), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_mean(retp, method="quantile", confl=0.9),
//'   Rcode=colSums(sapply(retp, quantile, c(0.9, 0.1), type=5)),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the column medians in RcppArmadillo
//' HighFreq::calc_mean(retp, method="nonparametric")
//' # Calculate the column medians in R
//' sapply(retp, median)
//' # Compare the values
//' all.equal(drop(HighFreq::calc_mean(retp, method="nonparametric")), 
//'   sapply(retp, median), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_mean(retp, method="nonparametric"),
//'   Rcode=sapply(retp, median),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_mean(const arma::mat& tseries,
                    std::string method = "moment", 
                    double confl = 0.75) {
  
  // Apply different calculation methods for location
  switch(calc_method(method)) {
  case methodenum::moment: {  // moment
    return arma::mean(tseries);
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-confl, confl};
    arma::mat quantiles = arma::quantile(tseries, levels);
    return (quantiles.row(0) + quantiles.row(1))/2;
  }  // end quantile
  case methodenum::nonparametric: {  // nonparametric
    return arma::median(tseries);
  }  // end nonparametric
  default : {
    cout << "Warning: Invalid method parameter" << endl;
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end default
  }  // end switch
  
}  // end calc_mean




////////////////////////////////////////////////////////////
//' Calculate the variance of a single-column \emph{time series} or a
//' \emph{vector} using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A single-column \emph{time series} or a \emph{vector}.
//'
//' @return A \emph{numeric} value equal to the variance of the \emph{vector}.
//'
//' @details
//'   The function \code{calc_varvec()} calculates the variance of a
//'   \emph{vector} using \code{RcppArmadillo} \code{C++} code, so it's
//'   significantly faster than the \code{R} function \code{var()}.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random returns
//' retp <- rnorm(1e6)
//' # Compare calc_varvec() with standard var()
//' all.equal(HighFreq::calc_varvec(retp), var(retp))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_varvec(retp),
//'   Rcode=var(retp),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
double calc_varvec(const arma::vec& tseries) {
  
  return arma::var(tseries);
  
}  // end calc_varvec




////////////////////////////////////////////////////////////
//' Calculate the dispersion (variance) of the columns of a \emph{time series}
//' or a \emph{matrix} using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'   
//' @param \code{method} A \emph{character string} specifying the type of the
//'   dispersion model (the default is \code{method = "moment"} - see Details).
//'    
//' @param \code{confl} The confidence level for calculating the quantiles of
//'   returns (the default is \code{confl = 0.75}).
//'
//' @return A row vector equal to the dispersion of the columns of the matrix
//'   \code{tseries}.
//'
//' @details
//'   The function \code{calc_var()} calculates the dispersion of the
//'   columns of a \emph{time series} or a \emph{matrix} of data using
//'   \code{RcppArmadillo} \code{C++} code.
//'   
//'   The dispersion is a measure of the variability of the data.  Examples of
//'   dispersion are the variance and the Median Absolute Deviation (\emph{MAD}).
//'
//'   If \code{method = "moment"} (the default) then \code{calc_var()}
//'   calculates the dispersion as the second moment of the data (the variance).
//'   Then \code{calc_var()} performs the same calculation as the function
//'   \code{colVars()} from package
//'   \href{https://cran.r-project.org/web/packages/matrixStats/index.html}{matrixStats},
//'   but it's much faster because it uses \code{RcppArmadillo} \code{C++} code.
//'
//'   If \code{method = "quantile"} then it calculates the dispersion as the
//'   difference between the quantiles as follows:
//'   \deqn{
//'     \sigma = q_{\alpha} - q_{1-\alpha}
//'   }
//'   Where \eqn{\alpha} is the confidence level for calculating the quantiles.
//'   
//'   If \code{method = "nonparametric"} then it calculates the dispersion as the
//'   Median Absolute Deviation (\emph{MAD}):
//'   \deqn{
//'     MAD = median(abs(x - median(x)))
//'   }
//'   It also multiplies the \emph{MAD} by a factor of \code{1.4826}, to make it
//'   comparable to the standard deviation.
//'
//'   If \code{method = "nonparametric"} then \code{calc_var()} performs the
//'   same calculation as the function \code{stats::mad()}, but it's much faster
//'   because it uses \code{RcppArmadillo} \code{C++} code.
//'
//'   If the number of rows of \code{tseries} is less than \code{3} then it
//'   returns zeros.
//'   
//' @examples
//' \dontrun{
//' # Calculate VTI and XLF returns
//' retp <- na.omit(rutils::etfenv$returns[, c("VTI", "XLF")])
//' # Compare HighFreq::calc_var() with standard var()
//' all.equal(drop(HighFreq::calc_var(retp)), 
//'   apply(retp, 2, var), check.attributes=FALSE)
//' # Compare HighFreq::calc_var() with matrixStats
//' all.equal(drop(HighFreq::calc_var(retp)), 
//'   matrixStats::colVars(retp), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with matrixStats and with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_var(retp),
//'   matrixStats=matrixStats::colVars(retp),
//'   Rcode=apply(retp, 2, var),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Compare HighFreq::calc_var() with stats::mad()
//' all.equal(drop(HighFreq::calc_var(retp, method="nonparametric")), 
//'   sapply(retp, mad), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with stats::mad()
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_var(retp, method="nonparametric"),
//'   Rcode=sapply(retp, mad),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_var(const arma::mat& tseries,
                   std::string method = "moment", 
                   double confl = 0.75) {
  
  arma::uword ncols = tseries.n_cols;
  // Return zeros if not enough data
  if (tseries.n_rows < 3) {
    return arma::zeros<rowvec>(ncols);
  }  // end if
  
  // Apply different calculation methods for dispersion
  switch(calc_method(method)) {
  case methodenum::moment: {  // moment
    return arma::var(tseries);
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-confl, confl};
    arma::mat quantiles = arma::quantile(tseries, levels);
    return (quantiles.row(1) - quantiles.row(0));
  }  // end quantile
  case methodenum::nonparametric: {  // MAD
    arma::mat medians = arma::median(tseries);
    arma::mat mads(1, ncols);
    // Loop over columns of tseries
    for (arma::uword it = 0; it < ncols; it++) {
      mads.col(it) = arma::median(arma::abs(tseries.col(it) - arma::as_scalar(medians.col(it))));
    }  // end for
    // tseries.each_row() -= arma::median(tseries, 0);
    // return 1.4826*arma::median(arma::abs(tseries), 0);
    return 1.4826*mads;
  }  // end nonparametric
  default : {
    cout << "Warning: Invalid method parameter" << endl;
    return arma::zeros<rowvec>(ncols);
  }  // end default
  }  // end switch
  
}  // end calc_var



////////////////////////////////////////////////////////////
//' Calculate the covariance matrix of the columns of a \emph{time series}
//' using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'   
//' @param \code{method} A \emph{character string} specifying the type of the
//'   covariance model (the default is \code{method = "moment"} - see Details).
//'    
//' @param \code{confl} The confidence level for calculating the quantiles of
//'   returns (the default is \code{confl = 0.75}).
//'
//' @return A square matrix with the covariance coefficients of the columns of
//'   the \emph{time series} \code{tseries}.
//'
//' @details
//'   The function \code{calc_covar()} calculates the covariance matrix of the
//'   columns of a \emph{time series} or a \emph{matrix} of data using
//'   \code{RcppArmadillo} \code{C++} code.
//'   The covariance is a measure of the codependency of the data.
//'
//'   If \code{method = "moment"} (the default) then \code{calc_covar()}
//'   calculates the covariance as the second co-moment:
//'   \deqn{
//'     \sigma_{xy} = \frac{1}{n-1} \sum_{i=1}^n (x_i - \bar{x}) (y_i - \bar{y})
//'   }
//'   Then \code{calc_covar()} performs the same calculation as the \code{R}
//'   function \code{stats::cov()}.
//'
//'   If \code{method = "quantile"} then it calculates the covariance as the
//'   difference between the quantiles as follows:
//'   \deqn{
//'     \mu = q_{\alpha} - q_{1-\alpha}
//'   }
//'   Where \eqn{\alpha} is the confidence level for calculating the quantiles.
//'   
//'   If \code{method = "nonparametric"} then it calculates the covariance as the
//'   Median Absolute Deviation (\emph{MAD}):
//'   \deqn{
//'     MAD = median(abs(x - median(x)))
//'   }
//'   It also multiplies the \emph{MAD} by a factor of \code{1.4826}, to make it
//'   comparable to the standard deviation.
//'
//'   If \code{method = "nonparametric"} then \code{calc_covar()} performs the
//'   same calculation as the function \code{stats::mad()}, but it's much faster
//'   because it uses \code{RcppArmadillo} \code{C++} code.
//'
//'   If the number of rows of \code{tseries} is less than \code{3} then it
//'   returns zeros.
//'   
//' @examples
//' \dontrun{
//' # Calculate VTI and XLF returns
//' retp <- na.omit(rutils::etfenv$returns[, c("VTI", "XLF")])
//' # Compare HighFreq::calc_covar() with standard var()
//' all.equal(drop(HighFreq::calc_covar(retp)), 
//'   cov(retp), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with matrixStats and with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_covar(retp),
//'   Rcode=cov(retp),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Compare HighFreq::calc_covar() with stats::mad()
//' all.equal(drop(HighFreq::calc_covar(retp, method="nonparametric")), 
//'   sapply(retp, mad), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with stats::mad()
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_covar(retp, method="nonparametric"),
//'   Rcode=sapply(retp, mad),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_covar(const arma::mat& tseries,
                     std::string method = "moment", 
                     double confl = 0.75) {
  
  arma::uword ncols = tseries.n_cols;
  // Return zeros if not enough data
  // if (tseries.n_rows < 3) {
  //   return arma::zeros<rowvec>(ncols);
  // }  // end if
  
  // Apply different calculation methods for covariance
  switch(calc_method(method)) {
  case methodenum::moment: {  // moment
    return arma::cov(tseries);
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-confl, confl};
    arma::mat quantiles = arma::quantile(tseries, levels);
    return (quantiles.row(1) - quantiles.row(0));
  }  // end quantile
  case methodenum::nonparametric: {  // MAD
    arma::mat medians = arma::median(tseries);
    arma::mat mads(1, ncols);
    // Loop over columns of tseries
    for (arma::uword it = 0; it < ncols; it++) {
      mads.col(it) = arma::median(arma::abs(tseries.col(it) - arma::as_scalar(medians.col(it))));
    }  // end for
    // tseries.each_row() -= arma::median(tseries, 0);
    // return 1.4826*arma::median(arma::abs(tseries), 0);
    return 1.4826*mads;
  }  // end nonparametric
  default : {
    cout << "Warning: Invalid method parameter" << endl;
    return arma::zeros<rowvec>(ncols);
  }  // end default
  }  // end switch
  
}  // end calc_covar




////////////////////////////////////////////////////////////
//' Calculate the variance of returns aggregated over the end points. 
//'
//' @param \code{pricev} A \emph{time series} or a \emph{matrix} of prices.
//'
//' @param \code{step} The number of time periods in each interval between
//'   neighboring end points (the default is \code{step = 1}).
//' 
//' @return The variance of aggregated returns.
//'
//' @details
//'   The function \code{calc_var_ag()} calculates the variance of returns
//'   aggregated over the end points.
//'
//'   It first calculates the end points spaced apart by the number of periods
//'   equal to the argument \code{step}.  Then it calculates the aggregated
//'   returns by differencing the prices \code{pricev} calculated at the end
//'   points. Finally it calculates the variance of the returns.
//'
//'   The choice of the first end point is arbitrary, so \code{calc_var_ag()}
//'   calculates the different end points for all the possible starting points.
//'   It then calculates the variance values for all the different end points
//'   and averages them.
//'
//'   The aggregated volatility \eqn{\sigma_t} increases with the length of the
//'   aggregation interval \eqn{\Delta t}.
//'   The aggregated volatility increases as the length of the aggregation
//'   interval \eqn{\Delta t} raised to the power of the \emph{Hurst exponent}
//'   \eqn{H}:
//'     \deqn{
//'       \sigma_t = \sigma {\Delta t}^H
//'     }
//'   Where \eqn{\sigma} is the daily return volatility.
//' 
//'   The function \code{calc_var_ag()} can therefore be used to calculate the
//'   \emph{Hurst exponent} from the variance ratio.
//'
//' @examples
//' \dontrun{
//' # Calculate the prices
//' closep <- na.omit(rutils::etfenv$prices[, c("XLP", "VTI")])
//' closep <- log(closep)
//' # Calculate the variance of daily returns
//' calc_var_ag(prices, step=1)
//' # Calculate the variance using R
//' sapply(rutils::diffit(closep), var)
//' # Calculate the variance of returns aggregated over 21 days
//' calc_var_ag(prices, step=21)
//' # The variance over 21 days is approximately 21 times the daily variance
//' 21*calc_var_ag(prices, step=1)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_var_ag(const arma::mat& pricev, 
                      arma::uword step = 1) {
  
  if (step == 1)
    // Calculate the variance without aggregations.
    return arma::var(diffit(pricev, 1, false));
  else {
    // Allocate aggregations, end points, and variance.
    arma::uword nrows = pricev.n_rows;
    arma::mat aggs;
    arma::uvec endd;
    // The number of rows is equal to step so that loop works for stub=0
    arma::mat vars(step, pricev.n_cols);
    // Perform loop over the stubs
    for (arma::uword stub = 0; stub < step; stub++) {
      endd = calc_endpoints(nrows, step, stub, false);
      // endd = arma::regspace<uvec>(stub, step, nrows + step);
      // endd = endd.elem(find(endd < nrows));
      aggs = pricev.rows(endd);
      vars.row(stub) = arma::var(diffit(aggs, 1, false));
    }  // end for
    return mean(vars);
  }  // end if
  
}  // end calc_var_ag




////////////////////////////////////////////////////////////
//' Calculate the variance of returns from \emph{OHLC} prices using different
//' price range estimators.
//'
//' @param \code{ohlc} A \emph{time series} or a \emph{matrix} of \emph{OHLC}
//'   prices.
//'   
//' @param \code{method} A \emph{character string} representing the price range
//'   estimator for calculating the variance.  The estimators include:
//'   \itemize{
//'     \item "close" close-to-close estimator,
//'     \item "rogers_satchell" Rogers-Satchell estimator,
//'     \item "garman_klass" Garman-Klass estimator,
//'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
//'     \item "yang_zhang" Yang-Zhang estimator,
//'    }
//'    (The default is the \code{method = "yang_zhang"}.)
//'    
//' @param \code{closel} A \emph{vector} with the lagged \emph{close} prices
//'   of the \emph{OHLC time series}.  This is an optional argument. (The
//'   default is \code{closel = 0}).
//'   
//' @param \code{scale} \emph{Boolean} argument: Should the returns be divided
//'   by the time index, the number of seconds in each period? (The default is
//'   \code{scale = TRUE}).
//'
//' @param \code{index} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument (the default is \code{index = 0}).
//'   
//' @return A single \emph{numeric} value equal to the variance of the
//'   \emph{OHLC time series}.
//'
//' @details
//'   The function \code{calc_var_ohlc()} calculates the variance from all the
//'   different intra-day and day-over-day returns (defined as the differences
//'   of \emph{OHLC} prices), using several different variance estimation
//'   methods.
//'
//'   The function \code{calc_var_ohlc()} does not calculate the logarithm of
//'   the prices.
//'   So if the argument \code{ohlc} contains dollar prices then
//'   \code{calc_var_ohlc()} calculates the dollar variance.
//'   If the argument \code{ohlc} contains the log prices then
//'   \code{calc_var_ohlc()} calculates the percentage variance.
//'
//'   The default \code{method} is \emph{"yang_zhang"}, which theoretically
//'   has the lowest standard error among unbiased estimators.
//'   The methods \emph{"close"}, \emph{"garman_klass_yz"}, and
//'   \emph{"yang_zhang"} do account for \emph{close-to-open} price jumps, while
//'   the methods \emph{"garman_klass"} and \emph{"rogers_satchell"} do not
//'   account for \emph{close-to-open} price jumps.
//'
//'   If \code{scale} is \code{TRUE} (the default), then the returns are
//'   divided by the differences of the time index (which scales the variance to
//'   the units of variance per second squared). This is useful when calculating
//'   the variance from minutes bar data, because dividing returns by the
//'   number of seconds decreases the effect of overnight price jumps. If the
//'   time index is in days, then the variance is equal to the variance per day
//'   squared.
//'   
//'   If the number of rows of \code{ohlc} is less than \code{3} then it
//'   returns zero.
//'   
//'   The optional argument \code{index} is the time index of the \emph{time
//'   series} \code{ohlc}. If the time index is in seconds, then the
//'   differences of the index are equal to the number of seconds in each time
//'   period.  If the time index is in days, then the differences are equal to
//'   the number of days in each time period.
//'   
//'   The optional argument \code{closel} are the lagged \emph{close} prices
//'   of the \emph{OHLC time series}.  Passing in the lagged \emph{close} prices
//'   speeds up the calculation, so it's useful for rolling calculations.
//'   
//'   The function \code{calc_var_ohlc()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, and it's over \code{10} times faster than
//'   \code{calc_var_ohlc_r()}, which is implemented in \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Extract the log OHLC prices of SPY
//' ohlc <- log(HighFreq::SPY)
//' # Extract the time index of SPY prices
//' indeks <- c(1, diff(xts::.index(ohlc)))
//' # Calculate the variance of SPY returns, with scaling of the returns
//' HighFreq::calc_var_ohlc(ohlc, 
//'  method="yang_zhang", scale=TRUE, index=indeks)
//' # Calculate variance without accounting for overnight jumps
//' HighFreq::calc_var_ohlc(ohlc, 
//'  method="rogers_satchell", scale=TRUE, index=indeks)
//' # Calculate the variance without scaling the returns
//' HighFreq::calc_var_ohlc(ohlc, scale=FALSE)
//' # Calculate the variance by passing in the lagged close prices
//' closel <- HighFreq::lagit(ohlc[, 4])
//' all.equal(HighFreq::calc_var_ohlc(ohlc), 
//'   HighFreq::calc_var_ohlc(ohlc, closel=closel))
//' # Compare with HighFreq::calc_var_ohlc_r()
//' all.equal(HighFreq::calc_var_ohlc(ohlc, index=indeks), 
//'   HighFreq::calc_var_ohlc_r(ohlc))
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_var_ohlc(ohlc),
//'   Rcode=HighFreq::calc_var_ohlc_r(ohlc),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
double calc_var_ohlc(const arma::mat& ohlc, 
                     std::string method = "yang_zhang", 
                     arma::colvec closel = 0, 
                     bool scale = true, // Divide the returns by time index
                     arma::colvec index = 0) {
  
  arma::uword nrows = ohlc.n_rows;
  double coeff = 0.34/(1.34 + (nrows+1)/(nrows-1));
  
  if (nrows < 3) {
    // Return zero if not enough data
    return 0;
  }  // end if
  
  if (!scale || (index.n_rows == 1)) {
    index = arma::ones(nrows);
    // cout << "ohlc.n_rows = " << nrows << endl;
    // cout << "index.n_rows = " << index.n_rows << endl;
  }  // end if
  
  // Calculate all the different intra-day and day-over-day returns 
  // (differences of OHLC prices)
  arma::mat openp = ohlc.col(0);
  arma::mat highp = ohlc.col(1);
  arma::mat lowp = ohlc.col(2);
  arma::mat closep = ohlc.col(3);
  arma::mat opcl = arma::zeros(nrows, 1);
  if (closel.n_rows == 1) {
    closel = arma::join_cols(closep.row(0), closep.rows(0, nrows-2));
    opcl = (openp - closel)/index;
  } else {
    opcl = (openp - closel)/index;
  }  // end if
  arma::mat clop = (closep - openp)/index;
  arma::mat clhi = (closep - highp)/index;
  arma::mat clow = (closep - lowp)/index;
  arma::mat hilow = (highp - lowp)/index;
  arma::mat hiop = (highp - openp)/index;
  arma::mat lowop = (lowp - openp)/index;
  
  if (method == "close") {
    // cout << "Calc method is Close" << endl;
    return arma::as_scalar(arma::var(arma::diff(closep)));
  } else if (method == "rogers_satchell") {
    // cout << "Calc method is Rogers-Satchell" << endl;
    return -(arma::dot(clhi, hiop) +
             arma::dot(clow, lowop))/nrows;
  } else if (method == "garman_klass") {
    // cout << "Calc method is Garman-Klass" << endl;
    return (0.5*arma::dot(hilow, hilow) -
            (2*log(2)-1)*arma::dot(clop, clop))/nrows;
  } else if (method == "garman_klass_yz") {
    // cout << "Calc method is Garman-Klass-YZ" << endl;
    return arma::as_scalar((0.5*arma::dot(hilow, hilow) -
                           (2*log(2)-1)*arma::dot(clop, clop))/nrows + 
                           arma::var(opcl));
  } else if (method == "yang_zhang") {
    // cout << "Calc method is Yang-Zhang" << endl;
    return arma::as_scalar(arma::var(opcl) + coeff*arma::var(clop) +
                           (coeff-1)*(arma::dot(clhi, hiop) + 
                           arma::dot(clow, lowop))/nrows);
  } else {
    cout << "Wrong calc method!" << endl;
    return 1;
  }  // end if
  
  // cout << "Calc method is " << method << endl;
  
}  // end calc_var_ohlc




////////////////////////////////////////////////////////////
//' Calculate the variance of aggregated \emph{OHLC} prices using different
//' price range estimators.
//'
//' @param \code{ohlc} A \emph{time series} or a \emph{matrix} of \emph{OHLC}
//'   prices.
//'
//' @param \code{step} The number of time periods in each interval between
//'   neighboring end points.
//' 
//' @param \code{method} A \emph{character string} representing the price range
//'   estimator for calculating the variance.  The estimators include:
//'   \itemize{
//'     \item "close" close-to-close estimator,
//'     \item "rogers_satchell" Rogers-Satchell estimator,
//'     \item "garman_klass" Garman-Klass estimator,
//'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
//'     \item "yang_zhang" Yang-Zhang estimator,
//'    }
//'    (The default is the \code{method = "yang_zhang"}.)
//'    
//' @param \code{closel} A \emph{vector} with the lagged \emph{close} prices
//'   of the \emph{OHLC time series}.  This is an optional argument. (The
//'   default is \code{closel = 0}).
//'   
//' @param \code{scale} \emph{Boolean} argument: Should the returns be divided
//'   by the time index, the number of seconds in each period? (The default is
//'   \code{scale = TRUE}).
//'
//' @param \code{index} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument (the default is \code{index = 0}).
//'   
//' @return The variance of aggregated \emph{OHLC} prices.
//'
//' @details
//'   The function \code{calc_var_ohlc_ag()} calculates the variance of
//'   \emph{OHLC} prices aggregated over the end points.
//'
//'   It first calculates the end points spaced apart by the number of periods
//'   equal to the argument \code{step}.  Then it aggregates the \emph{OHLC}
//'   prices to the end points. Finally it calculates the variance of the
//'   aggregated \emph{OHLC} prices.
//'
//'   The choice of the first end point is arbitrary, so \code{calc_var_ohlc_ag()}
//'   calculates the different end points for all the possible starting points.
//'   It then calculates the variance values for all the different end points
//'   and averages them.
//'
//'   The aggregated volatility \eqn{\sigma_t} increases with the length of the
//'   aggregation interval \eqn{\Delta t}.
//'   The aggregated volatility increases as the length of the aggregation
//'   interval \eqn{\Delta t} raised to the power of the \emph{Hurst exponent}
//'   \eqn{H}:
//'     \deqn{
//'       \sigma_t = \sigma {\Delta t}^H
//'     }
//'   Where \eqn{\sigma} is the daily return volatility.
//' 
//'   The function \code{calc_var_ohlc_ag()} can therefore be used to calculate
//'   the \emph{Hurst exponent} from the variance ratio.
//'
//' @examples
//' \dontrun{
//' # Calculate the log ohlc prices
//' ohlc <- log(rutils::etfenv$VTI)
//' # Calculate the daily variance of percentage returns
//' calc_var_ohlc_ag(ohlc, step=1)
//' # Calculate the variance of returns aggregated over 21 days
//' calc_var_ohlc_ag(ohlc, step=21)
//' # The variance over 21 days is approximately 21 times the daily variance
//' 21*calc_var_ohlc_ag(ohlc, step=1)
//' }
//' 
//' @export
// [[Rcpp::export]]
double calc_var_ohlc_ag(const arma::mat& ohlc,
                        arma::uword step, 
                        std::string method = "yang_zhang", 
                        arma::colvec closel = 0, 
                        bool scale = true, // Divide the returns by time index
                        arma::colvec index = 0) {
  
  if (step == 1)
    // Calculate the variance without aggregations.
    return calc_var_ohlc(ohlc, method, closel, scale, index);
  else {
    // Allocate aggregations, end points, and variance.
    arma::uword nrows = ohlc.n_rows;
    arma::mat aggs;
    arma::uvec endd;
    arma::mat vars(step, 1);
    // Perform loop over the stubs
    for (arma::uword stub = 0; stub < step; stub++) {
      endd = calc_endpoints(nrows, step, stub, false);
      aggs = roll_ohlc(ohlc, endd);
      vars.row(stub) = calc_var_ohlc(aggs, method, closel, scale, index);
    }  // end for
    return arma::as_scalar(mean(vars));
  }  // end if
  
}  // end calc_var_ohlc_ag



////////////////////////////////////////////////////////////
//' Calculate the skewness of the columns of a \emph{time series} or a
//' \emph{matrix} using \code{RcppArmadillo}.
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'
//' @param \code{method} A \emph{character string} specifying the type of the
//'   skewness model (the default is \code{method = "moment"} - see Details).
//'
//' @param \code{confl} The confidence level for calculating the quantiles of
//'   returns (the default is \code{confl = 0.75}).
//'
//' @return A single-row matrix with the skewness of the columns of
//'   \code{tseries}.
//'
//' @details
//'   The function \code{calc_skew()} calculates the skewness of the columns of
//'   a \emph{time series} or a \emph{matrix} of data using \code{C++}
//'   \code{RcppArmadillo} code.
//'
//'   If \code{method = "moment"} (the default) then \code{calc_skew()}
//'   calculates the skewness as the third moment of the data.
//'
//'   If \code{method = "quantile"} then it calculates the skewness
//'   \eqn{\varsigma} from the differences between the quantiles of the data as
//'   follows:
//'   \deqn{
//'     \varsigma = \frac{q_{\alpha} + q_{1-\alpha} - 2 q_{0.5}}{q_{\alpha} - q_{1-\alpha}}
//'   }
//'   Where \eqn{\alpha} is the confidence level for calculating the quantiles.
//'
//'   If \code{method = "nonparametric"} then it calculates the skewness as the
//'   difference between the mean of the data minus its median, divided by the
//'   standard deviation.
//'   
//'   If the number of rows of \code{tseries} is less than \code{3} then it
//'   returns zeros.
//'   
//'   The code examples below compare the function \code{calc_skew()} with the
//'   skewness calculated using \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Define a single-column time series of returns
//' retp <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the moment skewness
//' HighFreq::calc_skew(retp)
//' # Calculate the moment skewness in R
//' calc_skewr <- function(x) {
//'   x <- (x-mean(x))
//'   sum(x^3)/var(x)^1.5/NROW(x)
//' }  # end calc_skewr
//' all.equal(HighFreq::calc_skew(retp), 
//'   calc_skewr(retp), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(retp),
//'   Rcode=calc_skewr(retp),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the quantile skewness
//' HighFreq::calc_skew(retp, method="quantile", confl=0.9)
//' # Calculate the quantile skewness in R
//' calc_skewq <- function(x, a = 0.75) {
//'   	quantiles <- quantile(x, c(1-a, 0.5, a), type=5)
//'   	(quantiles[3] + quantiles[1] - 2*quantiles[2])/(quantiles[3] - quantiles[1])
//' }  # end calc_skewq
//' all.equal(drop(HighFreq::calc_skew(retp, method="quantile", confl=0.9)), 
//'   calc_skewq(retp, a=0.9), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(retp, method="quantile"),
//'   Rcode=calc_skewq(retp),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the nonparametric skewness
//' HighFreq::calc_skew(retp, method="nonparametric")
//' # Compare HighFreq::calc_skew() with R nonparametric skewness
//' all.equal(drop(HighFreq::calc_skew(retp, method="nonparametric")), 
//'   (mean(retp)-median(retp))/sd(retp), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(retp, method="nonparametric"),
//'   Rcode=(mean(retp)-median(retp))/sd(retp),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_skew(const arma::mat& tseries,
                    std::string method = "moment", 
                    double confl = 0.75) {
  
  // Return zeros if not enough data
  if (tseries.n_rows < 3) {
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end if
  
  // Apply different calculation methods for skew
  switch(calc_method(method)) {
  case methodenum::moment: {  // moment
    arma::uword nrows = tseries.n_rows;
    arma::uword ncols = tseries.n_cols;
    arma::mat meanm = arma::mean(tseries);
    arma::mat vars = arma::var(tseries);
    arma::mat skewness(1, ncols);
    // Center the columns of tseries
    // tseries.each_row() -= meanm;
    // return arma::sum(arma::pow(tseries, 3))/arma::pow(vars, 1.5)/nrows;
    for (arma::uword it = 0; it < ncols; it++) {
      skewness.col(it) = arma::sum(arma::pow(tseries.col(it) - arma::as_scalar(meanm.col(it)), 3))/arma::pow(vars.col(it), 1.5)/nrows;
    }  // end for
    return skewness;
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-confl, 0.5, confl};
    arma::mat quantiles = arma::quantile(tseries, levels);
    return (quantiles.row(2) + quantiles.row(0) - 2*quantiles.row(1))/(quantiles.row(2) - quantiles.row(0));
  }  // end quantile
  case methodenum::nonparametric: {  // nonparametric
    return (arma::mean(tseries) - arma::median(tseries))/arma::stddev(tseries);
  }  // end nonparametric
  default : {
    cout << "Warning: Invalid method parameter" << endl;
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end default
  }  // end switch
  
}  // end calc_skew



////////////////////////////////////////////////////////////
//' Calculate the kurtosis of the columns of a \emph{time series} or a
//' \emph{matrix} using \code{RcppArmadillo}.
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'
//' @param \code{method} A \emph{character string} specifying the type of the
//'   kurtosis model (the default is \code{method = "moment"} - see Details).
//'
//' @param \code{confl} The confidence level for calculating the quantiles of
//'   returns (the default is \code{confl = 0.75}).
//'
//' @return A single-row matrix with the kurtosis of the columns of
//'   \code{tseries}.
//'
//' @details
//'   The function \code{calc_kurtosis()} calculates the kurtosis of the columns
//'   of the \emph{matrix} \code{tseries} using \code{RcppArmadillo} \code{C++}
//'   code.
//'
//'   If \code{method = "moment"} (the default) then \code{calc_kurtosis()}
//'   calculates the fourth moment of the data.
//'   But it doesn't center the columns of \code{tseries} because that requires
//'   copying the matrix \code{tseries} in memory, so it's time-consuming.
//'
//'   If \code{method = "quantile"} then it calculates the skewness
//'   \eqn{\kappa} from the differences between the quantiles of the data as
//'   follows:
//'   \deqn{
//'     \kappa = \frac{q_{\alpha} - q_{1-\alpha}}{q_{0.75} - q_{0.25}}
//'   }
//'   Where \eqn{\alpha} is the confidence level for calculating the quantiles.
//'
//'   If \code{method = "nonparametric"} then it calculates the kurtosis as the
//'   difference between the mean of the data minus its median, divided by the
//'   standard deviation.
//'   
//'   If the number of rows of \code{tseries} is less than \code{3} then it
//'   returns zeros.
//'   
//'   The code examples below compare the function \code{calc_kurtosis()} with the
//'   kurtosis calculated using \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Define a single-column time series of returns
//' retp <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the moment kurtosis
//' HighFreq::calc_kurtosis(retp)
//' # Calculate the moment kurtosis in R
//' calc_kurtr <- function(x) {
//'   x <- (x-mean(x))
//'   sum(x^4)/var(x)^2/NROW(x)
//' }  # end calc_kurtr
//' all.equal(HighFreq::calc_kurtosis(retp), 
//'   calc_kurtr(retp), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_kurtosis(retp),
//'   Rcode=calc_kurtr(retp),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the quantile kurtosis
//' HighFreq::calc_kurtosis(retp, method="quantile", confl=0.9)
//' # Calculate the quantile kurtosis in R
//' calc_kurtq <- function(x, a=0.9) {
//'   	quantiles <- quantile(x, c(1-a, 0.25, 0.75, a), type=5)
//'   	(quantiles[4] - quantiles[1])/(quantiles[3] - quantiles[2])
//' }  # end calc_kurtq
//' all.equal(drop(HighFreq::calc_kurtosis(retp, method="quantile", confl=0.9)), 
//'   calc_kurtq(retp, a=0.9), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_kurtosis(retp, method="quantile"),
//'   Rcode=calc_kurtq(retp),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the nonparametric kurtosis
//' HighFreq::calc_kurtosis(retp, method="nonparametric")
//' # Compare HighFreq::calc_kurtosis() with R nonparametric kurtosis
//' all.equal(drop(HighFreq::calc_kurtosis(retp, method="nonparametric")), 
//'   (mean(retp)-median(retp))/sd(retp), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_kurtosis(retp, method="nonparametric"),
//'   Rcode=(mean(retp)-median(retp))/sd(retp),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_kurtosis(const arma::mat& tseries,
                        std::string method = "moment", 
                        double confl = 0.75) {
  
  // Return zeros if not enough data
  if (tseries.n_rows < 3) {
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end if
  
  // Apply different calculation methods for kurtosis
  switch(calc_method(method)) {
  case methodenum::moment: {  // Fourth moment
    arma::uword nrows = tseries.n_rows;
    arma::uword ncols = tseries.n_cols;
    arma::mat meanm = arma::mean(tseries);
    arma::mat vars = arma::var(tseries);
    arma::mat kurtosis(1, ncols);
    // Don't center the columns of tseries because that requires copying the matrix of data, so it's time-consuming
    // Loop over columns of tseries
    for (arma::uword it = 0; it < ncols; it++) {
      kurtosis.col(it) = arma::sum(arma::pow(tseries.col(it) - arma::as_scalar(meanm.col(it)), 4))/arma::pow(vars.col(it), 2)/nrows;
    }  // end for
    // tseries.each_row() -= meanm;
    // return arma::sum(arma::pow(tseries, 4))/arma::pow(vars, 2)/nrows;
    return kurtosis;
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-confl, 0.25, 0.75, confl};
    arma::mat quantiles = arma::quantile(tseries, levels);
    return (quantiles.row(3) - quantiles.row(0))/(quantiles.row(2) - quantiles.row(1));
  }  // end quantile
  case methodenum::nonparametric: {  // nonparametric
    return (arma::mean(tseries) - arma::median(tseries))/arma::stddev(tseries);
  }  // end nonparametric
  default : {
    cout << "Warning: Invalid method parameter" << endl;
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end default
  }  // end switch
  
}  // end calc_kurtosis




////////////////////////////////////////////////////////////
//' Calculate the Hurst exponent from the volatility ratio of aggregated returns.
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of log prices.
//'
//' @param \code{aggv} A \emph{vector} of aggregation intervals.
//' 
//' @return The Hurst exponent calculated from the volatility ratio of
//'   aggregated returns.  If \code{tseries} contains multiple columns, then the
//'   function \code{calc_hurst()} returns a single-row matrix of Hurst
//'   exponents.
//'
//' @details
//'   The function \code{calc_hurst()} calculates the Hurst exponent from the
//'   ratios of the volatilities of aggregated returns.
//'   
//'   An aggregation interval is equal to the number of time periods between the
//'   neighboring aggregation end points.
//'
//'   The aggregated volatility \eqn{\sigma_t} increases with the length of the
//'   aggregation interval \eqn{\Delta t}.
//'   The aggregated volatility increases as the length of the aggregation
//'   interval \eqn{\Delta t} raised to the power of the \emph{Hurst exponent}
//'   \eqn{H}:
//'     \deqn{
//'       \sigma_t = \sigma {\Delta t}^H
//'     }
//'   Where \eqn{\sigma} is the daily return volatility.
//' 
//'   For a single aggregation interval \eqn{\Delta t}, the \emph{Hurst
//'   exponent} \eqn{H} is equal to the logarithm of the ratio of the
//'   volatilities divided by the logarithm of the aggregation interval
//'   \eqn{\Delta t}:
//'     \deqn{
//'       H = \frac{\log{\sigma_t} - \log{\sigma}}{\log{\Delta t}}
//'     }
//' 
//'   For a \emph{vector} of aggregation intervals \eqn{\Delta t_i}, the
//'   \emph{Hurst exponent} \eqn{H} is equal to the regression slope between the
//'   logarithms of the aggregated volatilities \eqn{\sigma_i} versus the
//'   logarithms of the aggregation intervals \eqn{\Delta t_i}:
//'     \deqn{
//'       H = \frac{\code{cov}(\log{\sigma_i}, \log{\Delta t_i})}{\code{var}(\log{\Delta t_i})}
//'     }
//' 
//'   The function \code{calc_hurst()} calls the function \code{calc_var_ag()}
//'   to calculate the variance of aggregated returns \eqn{\sigma^2_t}.
//' 
//' @examples
//' \dontrun{
//' # Calculate the log prices
//' closep <- na.omit(rutils::etfenv$prices[, c("XLP", "VTI")])
//' closep <- log(closep)
//' # Calculate the Hurst exponents for a 21 day aggregation interval
//' HighFreq::calc_hurst(prices, aggv=21)
//' # Calculate the Hurst exponents for a vector of aggregation intervals
//' aggv <- seq.int(from=3, to=35, length.out=9)^2
//' HighFreq::calc_hurst(prices, aggv=aggv)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_hurst(const arma::mat& tseries, 
                     const arma::vec& aggv) {
  
  // If only single agg value then calculate the Hurst exponent from a single data point
  if (aggv.n_rows == 1) {
    return 0.5*arma::log(calc_var_ag(tseries, aggv(0))/calc_var_ag(tseries, 1))/log(aggv(0));
  }  // end if
  
  // Allocate the objects
  arma::uword nrows = aggv.n_rows;
  arma::mat volv(nrows, tseries.n_cols, fill::zeros);
  
  // Calculate the log volatility at the agg points
  for (arma::uword it = 0; it < nrows; it++) {
    volv.row(it) = 0.5*arma::log(calc_var_ag(tseries, aggv(it)));
  }  // end for
  
  // Calculate the log of the agg points
  arma::mat agglog = arma::log(aggv);
  
  // Calculate the Hurst exponent from the regression slopes
  arma::mat varagg = arma::var(agglog);
  return (arma::cov(volv, agglog).t())/varagg(0);
  
}  // end calc_hurst



////////////////////////////////////////////////////////////
//' Calculate the Hurst exponent from the volatility ratio of aggregated
//' \emph{OHLC} prices.
//'
//' @param \code{ohlc} A \emph{time series} or a \emph{matrix} of \emph{OHLC}
//'   prices.
//'
//' @param \code{step} The number of time periods in each interval between
//'   neighboring end points.
//' 
//' @param \code{method} A \emph{character string} representing the price range
//'   estimator for calculating the variance.  The estimators include:
//'   \itemize{
//'     \item "close" close-to-close estimator,
//'     \item "rogers_satchell" Rogers-Satchell estimator,
//'     \item "garman_klass" Garman-Klass estimator,
//'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
//'     \item "yang_zhang" Yang-Zhang estimator,
//'    }
//'    (The default is the \code{method = "yang_zhang"}.)
//'    
//' @param \code{closel} A \emph{vector} with the lagged \emph{close} prices
//'   of the \emph{OHLC time series}.  This is an optional argument. (The
//'   default is \code{closel = 0}).
//'   
//' @param \code{scale} \emph{Boolean} argument: Should the returns be divided
//'   by the time index, the number of seconds in each period? (The default is
//'   \code{scale = TRUE}).
//'
//' @param \code{index} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument (the default is \code{index = 0}).
//'   
//' @return The Hurst exponent calculated from the volatility ratio of
//'   aggregated \emph{OHLC} prices.
//'
//' @details
//' The function \code{calc_hurst_ohlc()} calculates the Hurst exponent from the
//' ratios of the variances of aggregated \emph{OHLC} prices.
//'
//'   The aggregated volatility \eqn{\sigma_t} increases with the length of the
//'   aggregation interval \eqn{\Delta t}.
//'   The aggregated volatility increases as the length of the aggregation
//'   interval \eqn{\Delta t} raised to the power of the \emph{Hurst exponent}
//'   \eqn{H}:
//'     \deqn{
//'       \sigma_t = \sigma {\Delta t}^H
//'     }
//'   Where \eqn{\sigma} is the daily return volatility.
//' 
//'   The \emph{Hurst exponent} \eqn{H} is equal to the logarithm of the ratio
//'   of the volatilities divided by the logarithm of the time interval
//'   \eqn{\Delta t}:
//'     \deqn{
//'       H = \frac{\log{\sigma_t} - \log{\sigma}}{\log{\Delta t}}
//'     }
//' 
//'   The function \code{calc_hurst_ohlc()} calls the function
//'   \code{calc_var_ohlc_ag()} to calculate the aggregated variance
//'   \eqn{\sigma^2_t}.
//' 
//' @examples
//' \dontrun{
//' # Calculate the log ohlc prices
//' ohlc <- log(rutils::etfenv$VTI)
//' # Calculate the Hurst exponent from 21 day aggregations
//' calc_hurst_ohlc(ohlc, step=21)
//' }
//' 
//' @export
// [[Rcpp::export]]
double calc_hurst_ohlc(const arma::mat& ohlc,
                       arma::uword step, 
                       std::string method = "yang_zhang", 
                       arma::colvec closel = 0, 
                       bool scale = true, // Divide the returns by time index
                       arma::colvec index = 0) {
  
  return 0.5*log(calc_var_ohlc_ag(ohlc, step, method, closel, scale, index)/
                 calc_var_ohlc_ag(ohlc, 1, method, closel, scale, index))/log(step);
  
}  // end calc_hurst_ohlc




////////////////////////////////////////////////////////////
//' Perform multivariate linear regression using least squares and return a
//' named list of regression coefficients, their t-values, and p-values.
//' 
//' @param \code{respv} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{predm} A \emph{time series} or a \emph{matrix} of predictor
//'   data.
//' 
//' @return A named list with three elements: a \emph{matrix} of coefficients
//'   (named \emph{"coefficients"}), the \emph{z-score} of the last residual
//'   (named \emph{"zscore"}), and a \emph{vector} with the R-squared and
//'   F-statistic (named \emph{"stats"}). The numeric \emph{matrix} of
//'   coefficients named \emph{"coefficients"} contains the alpha and beta
//'   coefficients, and their \emph{t-values} and \emph{p-values}.
//'
//' @details
//'   The function \code{calc_lm()} performs the same calculations as the
//'   function \code{lm()} from package \emph{stats}. 
//'   It uses \code{RcppArmadillo} \code{C++} code so it's several times faster
//'   than \code{lm()}. The code was inspired by this article (but it's not
//'   identical to it):
//'   http://gallery.rcpp.org/articles/fast-linear-model-with-armadillo/
//'
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' retp <- na.omit(rutils::etfenv$returns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' respv <- retp[, 1]
//' # Predictor matrix equals VTI and IEF returns
//' predm <- retp[, -1]
//' # Perform multivariate regression using lm()
//' regmod <- lm(respv ~ predm)
//' regsum <- summary(regmod)
//' # Add unit intercept column to the predictor matrix
//' predm <- cbind(rep(1, NROW(predm)), predm)
//' # Perform multivariate regression using calc_lm()
//' regarma <- HighFreq::calc_lm(respv=respv, predm=predm)
//' # Compare the outputs of both functions
//' all.equal(regarma$coefficients[, "coeff"], unname(coef(regmod)))
//' all.equal(unname(regarma$coefficients), unname(regsum$coefficients))
//' all.equal(unname(regarma$stats), c(regsum$r.squared, unname(regsum$fstatistic[1])))
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_lm(respv=respv, predm=predm),
//'   Rcode=lm(respv ~ predm),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List calc_lm(const arma::vec& respv,  // Response vector
                   const arma::mat& predm) { // Predictor matrix
  
  // Add column for intercept to the predictor matrix
  arma::uword nrows = predm.n_rows;
  // arma::mat predm = arma::join_rows(ones(nrows), predm);
  arma::uword ncols = predm.n_cols;
  arma::uword degf = (nrows - ncols);
  
  // Calculate alpha and beta coefficients for the model response ~ predictor
  arma::colvec coeff = arma::solve(predm, respv);
  // Calculate residuals
  arma::colvec residuals = respv - predm*coeff;
  
  // Calculate TSS, RSS, and ESS
  double tot_sumsq = (nrows-1)*arma::var(respv);
  double res_sumsq = arma::dot(residuals, residuals);
  double exp_sumsq = tot_sumsq - res_sumsq;
  
  // Calculate R-squared and F-statistic
  double rsquared = exp_sumsq/tot_sumsq;
  double fstat = (exp_sumsq*degf)/(res_sumsq*(ncols-1));
  // arma::rowvec stats=join_horiz(rsquared, fstat);
  Rcpp::NumericVector stats(2);
  stats(0) = rsquared;
  stats(1) = fstat;
  stats.attr("names") = Rcpp::CharacterVector::create("R-squared", "F-statistic");
  
  // Calculate standard errors of beta coefficients
  arma::colvec stderrv = arma::sqrt(res_sumsq/degf*arma::diagvec(arma::pinv(arma::trans(predm)*predm)));
  // Calculate t-values and p-values of beta coefficients
  arma::colvec tvals = coeff/stderrv;
  arma::colvec pvals = 2*Rcpp::pt(-Rcpp::abs(Rcpp::wrap(tvals)), degf);
  Rcpp::NumericMatrix coeffmat = Rcpp::wrap(arma::join_rows(arma::join_rows(arma::join_rows(coeff, stderrv), tvals), pvals));
  Rcpp::colnames(coeffmat) = Rcpp::CharacterVector::create("coeff", "stderr", "tvals", "pvals");
  
  return Rcpp::List::create(Rcpp::Named("coefficients") = coeffmat,
                            // Named("residuals") = residuals,
                            Rcpp::Named("zscore") = residuals(nrows-1)/arma::stddev(residuals),
                            Rcpp::Named("stats") = stats);
  
}  // end calc_lm



////////////////////////////////////////////////////////////
//' Perform multivariate regression using different methods, and return a vector
//' of regression coefficients, their t-values, and the last residual z-score.
//' 
//' @param \code{respv} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{predm} A \emph{time series} or a \emph{matrix} of predictor
//'   data.
//' 
//' @param \code{controlv} A \emph{list} of model parameters (see Details).
//'
//' @return A single-row matrix with the regression coefficients, their
//'   t-values, and the last residual z-score.
//'
//' @details
//'   The function \code{calc_reg()} performs multivariate regression using
//'   different methods, and returns a vector of regression coefficients, their
//'   t-values, and the last residual z-score.
//'   
//'   The function \code{calc_reg()} accepts a list of regression model
//'   parameters through the argument \code{controlv}.
//'   The list of model parameters can be created using the function
//'   \code{param_reg()}.  Below is a description of the model parameters.
//'
//'   If \code{regmod = "least_squares"} (the default) then it performs the
//'   standard least squares regression, the same as the function
//'   \code{calc_lm()}, and the function \code{lm()} from the \code{R} package
//'   \emph{stats}.
//'   But it uses \code{RcppArmadillo} \code{C++} code so it's several times
//'   faster than \code{lm()}.
//'
//'   If \code{regmod = "regular"} then it performs shrinkage regression.  It
//'   calculates the \emph{reduced inverse} of the predictor matrix from its
//'   singular value decomposition.  It performs regularization by selecting
//'   only the largest \emph{singular values} equal in number to \code{dimax}.
//'   
//'   If \code{regmod = "quantile"} then it performs quantile regression (not
//'   implemented yet).
//' 
//'   The length of the return vector depends on the number of columns of the
//'   predictor matrix (including the intercept column, if it's been added in
//'   \code{R}).
//'   The number of regression coefficients is equal to the number of columns of
//'   the predictor matrix.
//'   The length of the return vector is equal to the number of regression
//'   coefficients, plus their t-values, plus the z-score.
//'   The number of t-values is equal to the number of coefficients.
//' 
//'   For example, if the number of columns of the predictor matrix is equal to
//'   \code{n}, then \code{calc_reg()} returns a vector with \code{2n+1}
//'   elements: \code{n} regression coefficients, \code{n} corresponding
//'   t-values, and \code{1} z-score value.
//'
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' retp <- na.omit(rutils::etfenv$returns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' respv <- retp[, 1]
//' # Predictor matrix equals VTI and IEF returns
//' predm <- retp[, -1]
//' # Perform multivariate regression using lm()
//' regmod <- lm(respv ~ predm)
//' regsum <- summary(regmod)
//' coeff <- regsum$coefficients
//' # Create a default list of regression parameters
//' controlv <- HighFreq::param_reg()
//' # Add unit intercept column to the predictor matrix
//' predm <- cbind(rep(1, NROW(predm)), predm)
//' # Perform multivariate regression using calc_reg()
//' regarma <- drop(HighFreq::calc_reg(respv=respv, predm=predm, controlv=controlv))
//' # Compare the outputs of both functions
//' all.equal(regarma[1:(2*NCOL(predm))], 
//'   c(coeff[, "Estimate"], coeff[, "t value"]), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_reg(respv=respv, predm=predm, controlv=controlv),
//'   Rcode=lm(respv ~ predm),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_reg(const arma::mat& respv,  // Response vector
                   const arma::mat& predm, // Predictor matrix
                   Rcpp::List controlv) { // List of regression model parameters
  
  // Unpack the control list
  // Type of regression model
  std::string regmod = Rcpp::as<std::string>(controlv["regmod"]);
  // Add unit intercept column to the predictor matrix?
  // bool intercept = Rcpp::as<int>(controlv["intercept"]);
  // Threshold level for discarding small singular values
  double singmin = Rcpp::as<double>(controlv["singmin"]);
  // Apply dimension reduction
  arma::uword dimax = Rcpp::as<int>(controlv["dimax"]);
  // Confidence level for calculating the quantiles of returns
  // double confl = Rcpp::as<double>(controlv["confl"]);
  // Shrinkage intensity of returns
  // double alpha = Rcpp::as<double>(controlv["alpha"]);
  
  // Add column for intercept to the predictor matrix - no, add it in R
  arma::uword nrows = predm.n_rows;
  // arma::mat predm = predm; // Predictor matrix with intercept column
  // if (intercept)
  //   predm = arma::join_rows(ones(nrows), predm);
  
  arma::uword ncols = predm.n_cols;
  arma::uword degf = (nrows - ncols);
  arma::vec coeff;
  arma::vec tvals;
  arma::mat reg_data = arma::zeros(2*ncols+1, 1);
  
  // Apply different calculation regmods for the regression coefficients
  switch(calc_method(regmod)) {
  case methodenum::least_squares: {
    // Calculate regression coefficients for the model response ~ predictor
    coeff = arma::solve(predm, respv);
    break;
  }  // end least_squares
  case methodenum::regular: {
    // Calculate shrinkage regression coefficients
    coeff = calc_inv(predm, dimax, singmin)*respv;
    break;
  }  // end regular
  case methodenum::quantile: {
    // Not implemented yet
    break;
  }  // end quantile
  default : {
    cout << "Warning: Invalid regmod parameter: " << regmod << endl;
    return reg_data;
  }  // end default
  }  // end switch
  
  // Calculate residuals
  arma::mat residuals = respv - predm*coeff;
  
  // Calculate TSS, RSS, and ESS
  // double tot_sumsq = (nrows-1)*arma::var(respv);
  double res_sumsq = arma::dot(residuals, residuals);
  // double exp_sumsq = tot_sumsq - res_sumsq;
  
  // Calculate standard errors of the beta coefficients
  arma::mat stderrv = arma::sqrt(res_sumsq/degf*arma::diagvec(arma::pinv(arma::trans(predm)*predm)));
  // Calculate t-values of the beta coefficients
  tvals = coeff/stderrv;
  
  // Calculate z-score
  arma::mat zscore = residuals(nrows-1, 0)/arma::stddev(residuals);
  
  // Combine regression data
  reg_data.rows(0, ncols-1) = coeff;
  reg_data.rows(ncols, 2*ncols-1) = tvals;
  reg_data.row(2*ncols) = zscore;
  
  return reg_data.t();
  
}  // end calc_reg




////////////////////////////////////////////////////////////
// Functions for rolling statistics
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of mean (location) estimates over a rolling
//' look-back interval attached at the end points of a \emph{time series} or a
//' \emph{matrix}.
//'
//' @param \code{lookb} The number of end points in the look-back interval
//'   (the default is \code{lookb = 1}).
//'   
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endd} An \emph{integer} vector of end points (the default is
//'   \code{endd = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @param \code{method} A \emph{character string} representing the type of mean 
//'   measure of (the default is \code{method = "moment"}).
//'
//' @return A \emph{matrix} of mean (location) estimates with the same number of
//'   columns as the input time series \code{tseries}, and the number of rows
//'   equal to the number of end points.
//'   
//' @details
//'   The function \code{roll_mean()} calculates a \emph{matrix} of mean
//'   (location) estimates over rolling look-back intervals attached at the end
//'   points of the \emph{time series} \code{tseries}.
//'   
//'   The function \code{roll_mean()} performs a loop over the end points, and at
//'   each end point it subsets the time series \code{tseries} over a look-back
//'   interval equal to \code{lookb} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_mean()}, which
//'   calculates the mean (location).
//'   See the function \code{calc_mean()} for a description of the mean methods.
//'   
//'   If the arguments \code{endd} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{tseries}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//' 
//'   For example, the rolling mean at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{lookb = 3}.
//'
//'   The function \code{roll_mean()} with the parameter \code{step = 1}
//'   performs the same calculation as the function \code{roll_mean()} from
//'   package
//'   \href{https://cran.r-project.org/web/packages/RcppRoll/index.html}{RcppRoll},
//'   but it's several times faster because it uses \code{C++}
//'   \code{RcppArmadillo} code.
//'
//'   The function \code{roll_mean()} is implemented in \code{RcppArmadillo}
//'   \code{RcppArmadillo} \code{C++} code, which makes it several times faster
//'   than \code{R} code.
//'
//'   If only a simple rolling mean is required (not the median) then other
//'   functions like \code{roll_sum()} or \code{roll_vec()} may be even faster.
//'
//' @examples
//' \dontrun{
//' # Define time series of returns using package rutils
//' retp <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the rolling means at 25 day end points, with a 75 day look-back
//' meanv <- HighFreq::roll_mean(retp, lookb=3, step=25)
//' # Compare the mean estimates over 11-period look-back intervals
//' all.equal(HighFreq::roll_mean(retp, lookb=11)[-(1:10), ], 
//'   drop(RcppRoll::roll_mean(retp, n=11)), check.attributes=FALSE)
//' # Define end points and start points
//' endd <- HighFreq::calc_endpoints(NROW(retp), step=25)
//' startp <- HighFreq::calc_startpoints(endd, lookb=3)
//' # Calculate the rolling means using RcppArmadillo
//' meanv <- HighFreq::roll_mean(retp, startp=startp, endd=endd)
//' # Calculate the rolling medians using RcppArmadillo
//' medianscpp <- HighFreq::roll_mean(retp, startp=startp, endd=endd, method="nonparametric")
//' # Calculate the rolling medians using R
//' medians = sapply(1:NROW(endd), function(i) {
//'   median(retp[startp[i]:endd[i] + 1])
//' })  # end sapply
//' all.equal(medians, drop(medianscpp))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_mean(retp, startp=startp, endd=endd, method="nonparametric"),
//'   Rcode=sapply(1:NROW(endd), function(i) {median(retp[startp[i]:endd[i] + 1])}),
//'   times=10))[, c(1, 4, 5)]
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_mean(const arma::mat& tseries, 
                    arma::uword lookb = 1, 
                    arma::uvec startp = 0, 
                    arma::uvec endd = 0, 
                    arma::uword step = 1, 
                    arma::uword stub = 0,
                    std::string method = "moment", 
                    double confl = 0.75) {
  
  // Allocate end points
  arma::uword nrows = tseries.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endd) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endd;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by lookb
    startpts = calc_startpoints(endpts, lookb);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  
  // Allocate mean matrix
  arma::uword numpts = endpts.n_elem;
  arma::mat meanm = arma::zeros(numpts, tseries.n_cols);
  meanm.row(0) = tseries.row(0);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate means
    if (endpts(ep) > startpts(ep)) {
      meanm.row(ep) = calc_mean(tseries.rows(startpts(ep), endpts(ep)), method, confl);
    }  // end if
  }  // end for
  
  return meanm;
  
}  // end roll_mean




////////////////////////////////////////////////////////////
//' Calculate a \emph{vector} of variance estimates over a rolling look-back
//' interval for a single-column \emph{time series} or a single-column
//' \emph{matrix}, using \code{RcppArmadillo}.
//'
//' @param \code{tseries} A single-column \emph{time series} or a single-column
//'   \emph{matrix}.
//' 
//' @param \code{lookb} The length of the look-back interval, equal to the
//'   number of \emph{vector} elements used for calculating a single variance
//'   estimate (the default is \code{lookb = 1}).
//'
//' @return A single-column \emph{matrix} with the same number of elements as
//'   the input argument \code{tseries}.
//'
//' @details
//'   The function \code{roll_varvec()} calculates a \emph{vector} of variance
//'   estimates over a rolling look-back interval for a single-column \emph{time
//'   series} or a single-column \emph{matrix}, using \code{RcppArmadillo} \code{C++}
//'   code.
//'   
//'   The function \code{roll_varvec()} uses an expanding look-back interval in
//'   the initial warmup period, to calculate the same number of elements as the
//'   input argument \code{tseries}.
//'
//'   The function \code{roll_varvec()} performs the same calculation as the
//'   function \code{roll_var()} from package
//'   \href{https://cran.r-project.org/web/packages/RcppRoll/index.html}{RcppRoll},
//'   but it's several times faster because it uses \code{RcppArmadillo} \code{C++}
//'   code.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random returns
//' retp <- rnorm(1e6)
//' # Compare the variance estimates over 11-period look-back intervals
//' all.equal(drop(HighFreq::roll_varvec(retp, lookb=11))[-(1:10)], 
//'   RcppRoll::roll_var(retp, n=11))
//' # Compare the speed of RcppArmadillo with RcppRoll
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_varvec(retp, lookb=11),
//'   RcppRoll=RcppRoll::roll_var(retp, n=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::vec roll_varvec(const arma::vec& tseries, arma::uword lookb = 1) {
  
  arma::uword length = tseries.n_elem;
  arma::vec vars = arma::zeros(length);
  
  // Warmup period
  for (arma::uword it = 1; it < lookb; it++) {
    vars(it) = arma::var(tseries.subvec(0, it));
  }  // end for
  
  // Remaining period
  for (arma::uword it = lookb; it < length; it++) {
    vars(it) = arma::var(tseries.subvec(it-lookb+1, it));
  }  // end for
  
  return vars;
  
}  // end roll_varvec




////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of dispersion (variance) estimates over a rolling
//' look-back interval attached at the end points of a \emph{time series} or a
//' \emph{matrix}.
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'
//' @param \code{lookb} The number of end points in the look-back interval
//'   (the default is \code{lookb = 1}).
//'   
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endd} An \emph{integer} vector of end points (the default is
//'   \code{endd = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @param \code{method} A \emph{character string} representing the type of the
//'   measure of dispersion (the default is \code{method = "moment"}).
//'
//' @return A \emph{matrix} dispersion (variance) estimates with the same number
//'   of columns as the input time series \code{tseries}, and the number of rows
//'   equal to the number of end points.
//'   
//' @details
//'   The function \code{roll_var()} calculates a \emph{matrix} of dispersion
//'   (variance) estimates over rolling look-back intervals attached at the end
//'   points of the \emph{time series} \code{tseries}.
//'   
//'   The function \code{roll_var()} performs a loop over the end points, and at
//'   each end point it subsets the time series \code{tseries} over a look-back
//'   interval equal to \code{lookb} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_var()}, which
//'   calculates the dispersion.
//'   See the function \code{calc_var()} for a description of the dispersion
//'   methods.
//'   
//'   If the arguments \code{endd} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{tseries}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//' 
//'   For example, the rolling variance at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{lookb = 3}.
//'
//'   The function \code{roll_var()} with the parameter \code{step = 1}
//'   performs the same calculation as the function \code{roll_var()} from
//'   package
//'   \href{https://cran.r-project.org/web/packages/RcppRoll/index.html}{RcppRoll},
//'   but it's several times faster because it uses \code{RcppArmadillo} \code{C++}
//'   code.
//'
//'   The function \code{roll_var()} is implemented in \code{RcppArmadillo}
//'   \code{RcppArmadillo} \code{C++} code, which makes it several times faster
//'   than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Define time series of returns using package rutils
//' retp <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the rolling variance at 25 day end points, with a 75 day look-back
//' varv <- HighFreq::roll_var(retp, lookb=3, step=25)
//' # Compare the variance estimates over 11-period look-back intervals
//' all.equal(HighFreq::roll_var(retp, lookb=11)[-(1:10), ], 
//'   drop(RcppRoll::roll_var(retp, n=11)), check.attributes=FALSE)
//' # Compare the speed of HighFreq::roll_var() with RcppRoll::roll_var()
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_var(retp, lookb=11),
//'   RcppRoll=RcppRoll::roll_var(retp, n=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Compare the speed of HighFreq::roll_var() with TTR::runMAD()
//' summary(microbenchmark(
//'     Rcpp=HighFreq::roll_var(retp, lookb=11, method="quantile"),
//'     TTR=TTR::runMAD(retp, n = 11),
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_var(const arma::mat& tseries, 
                   arma::uword lookb = 1, 
                   arma::uvec startp = 0, 
                   arma::uvec endd = 0, 
                   arma::uword step = 1, 
                   arma::uword stub = 0,
                   std::string method = "moment", 
                   double confl = 0.75) {
  
  // Allocate end points
  arma::uword nrows = tseries.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endd) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endd;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by lookb
    startpts = calc_startpoints(endpts, lookb);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  
  // Allocate variance matrix
  arma::uword numpts = endpts.n_elem;
  arma::mat varm = arma::zeros(numpts, tseries.n_cols);
  varm.row(0) = arma::square(tseries.row(0));
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate variance
    if (endpts(ep) > startpts(ep)) {
      varm.row(ep) = calc_var(tseries.rows(startpts(ep), endpts(ep)), method, confl);
    }  // end if
  }  // end for
  varm.row(1) = arma::square(tseries.row(1));
  
  return varm;
  
}  // end roll_var




////////////////////////////////////////////////////////////
//' Calculate a \emph{vector} of variance estimates over a rolling look-back
//' interval attached at the end points of a \emph{time series} or a
//' \emph{matrix} with \emph{OHLC} price data.
//' 
//' @param \code{ohlc} A \emph{time series} or a \emph{matrix} with \emph{OHLC}
//'   price data.
//'   
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endd} An \emph{integer} vector of end points (the default is
//'   \code{endd = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{lookb} The number of end points in the look-back interval
//'   (the default is \code{lookb = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @param \code{method} A \emph{character string} representing the price range
//'   estimator for calculating the variance.  The estimators include:
//'   \itemize{
//'     \item "close" close-to-close estimator,
//'     \item "rogers_satchell" Rogers-Satchell estimator,
//'     \item "garman_klass" Garman-Klass estimator,
//'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
//'     \item "yang_zhang" Yang-Zhang estimator,
//'    }
//'    (The default is the \emph{"yang_zhang"} estimator.)
//'    
//' @param \code{scale} \emph{Boolean} argument: Should the returns be divided
//'   by the time index, the number of seconds in each period?  (The default is
//'   \code{scale = TRUE}.)
//'   
//' @param \code{index} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument (the default is \code{index=0}).
//'
//' @return A column \emph{vector} of variance estimates, with the number of
//'   rows equal to the number of end points.
//'
//' @details
//'   The function \code{roll_var_ohlc()} calculates a \emph{vector} of variance
//'   estimates over a rolling look-back interval attached at the end points of
//'   the \emph{time series} \code{ohlc}.
//'   
//'   The input \emph{OHLC time series} \code{ohlc} is assumed to contain the
//'   log prices.
//'
//'   The function \code{roll_var_ohlc()} performs a loop over the end points,
//'   subsets the previous (past) rows of \code{ohlc}, and passes them into the
//'   function \code{calc_var_ohlc()}.
//' 
//'   At each end point, the variance is calculated over a look-back interval
//'   equal to \code{lookb} number of end points.
//'   In the initial warmup period, the variance is calculated over an expanding
//'   look-back interval.
//'   
//'   If the arguments \code{endd} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{ohlc}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//' 
//'   For example, the rolling variance at daily end points with an \code{11}
//'   day look-back, can be calculated using the parameters \code{step = 1} and
//'   \code{lookb = 1} (Assuming the \code{ohlc} data has daily
//'   frequency.)
//' 
//'   Similarly, the rolling variance at \code{25} day end points with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{lookb = 3} (because \code{3*25 = 75}).
//' 
//'   The function \code{roll_var_ohlc()} calculates the variance from all the
//'   different intra-day and day-over-day returns (defined as the differences
//'   between \emph{OHLC} prices), using several different variance estimation
//'   methods.
//'   
//'   The default \code{method} is \emph{"yang_zhang"}, which theoretically
//'   has the lowest standard error among unbiased estimators.
//'   The methods \emph{"close"}, \emph{"garman_klass_yz"}, and
//'   \emph{"yang_zhang"} do account for \emph{close-to-open} price jumps, while
//'   the methods \emph{"garman_klass"} and \emph{"rogers_satchell"} do not
//'   account for \emph{close-to-open} price jumps.
//'
//'   If \code{scale} is \code{TRUE} (the default), then the returns are
//'   divided by the differences of the time index (which scales the variance to
//'   the units of variance per second squared.) This is useful when calculating
//'   the variance from minutes bar data, because dividing returns by the
//'   number of seconds decreases the effect of overnight price jumps. If the
//'   time index is in days, then the variance is equal to the variance per day
//'   squared.
//'   
//'   The optional argument \code{index} is the time index of the \emph{time
//'   series} \code{ohlc}. If the time index is in seconds, then the
//'   differences of the index are equal to the number of seconds in each time
//'   period.  If the time index is in days, then the differences are equal to
//'   the number of days in each time period.
//'   
//'   The function \code{roll_var_ohlc()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, which makes it several times faster than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Extract the log OHLC prices of SPY
//' ohlc <- log(HighFreq::SPY)
//' # Extract the time index of SPY prices
//' indeks <- c(1, diff(xts::.index(ohlc)))
//' # Rolling variance at minutes end points, with a 21 minute look-back
//' varoll <- HighFreq::roll_var_ohlc(ohlc, 
//'                               step=1, lookb=21, 
//'                               method="yang_zhang", 
//'                               index=indeks, scale=TRUE)
//' # Daily OHLC prices
//' ohlc <- rutils::etfenv$VTI
//' indeks <- c(1, diff(xts::.index(ohlc)))
//' # Rolling variance at 5 day end points, with a 20 day look-back (20=4*5)
//' varoll <- HighFreq::roll_var_ohlc(ohlc, 
//'                               step=5, lookb=4, 
//'                               method="yang_zhang", 
//'                               index=indeks, scale=TRUE)
//' # Same calculation in R
//' nrows <- NROW(ohlc)
//' closel = HighFreq::lagit(ohlc[, 4])
//' endd <- drop(HighFreq::calc_endpoints(nrows, 3)) + 1
//' startp <- drop(HighFreq::calc_startpoints(endd, 2))
//' npts <- NROW(endd)
//' varollr <- sapply(2:npts, function(it) {
//'   rangev <- startp[it]:endd[it]
//'   sub_ohlc = ohlc[rangev, ]
//'   sub_close = closel[rangev]
//'   sub_index = indeks[rangev]
//'   HighFreq::calc_var_ohlc(sub_ohlc, closel=sub_close, scale=TRUE, index=sub_index)
//' })  # end sapply
//' varollr <- c(0, varollr)
//' all.equal(drop(var_rolling), varollr)
//' }
//' @export
// [[Rcpp::export]]
arma::vec roll_var_ohlc(const arma::mat& ohlc, 
                        arma::uvec startp = 0, 
                        arma::uvec endd = 0, 
                        arma::uword step = 1, 
                        arma::uword lookb = 1, 
                        arma::uword stub = 0,
                        std::string method = "yang_zhang", 
                        bool scale = true, // Divide the returns by time index
                        arma::colvec index = 0) {
  
  // Allocate end points
  arma::uword nrows = ohlc.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endd) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endd;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by lookb
    startpts = calc_startpoints(endpts, lookb);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  // Allocate variance matrix
  arma::uword numpts = endpts.n_elem;
  arma::vec varm = arma::zeros(numpts);
  
  // Extract OHLC close prices
  arma::colvec closep = ohlc.col(3);
  arma::colvec closel = lagit(closep, 1, false);
  
  // Set the time index to 1 if scale = FALSE
  if (!scale || (index.n_rows == 1)) {
    index = arma::ones(nrows);
  }  // end if
  
  // Define data subsets over look-back intervals
  arma::mat sub_ohlc;
  arma::colvec sub_close;
  arma::colvec sub_index;
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < numpts; ep++) {
    if (endpts(ep) > startpts(ep)) {
      sub_ohlc = ohlc.rows(startpts(ep), endpts(ep));
      sub_close = closel.rows(startpts(ep), endpts(ep));
      sub_index = index.subvec(startpts(ep), endpts(ep));
      // Calculate variance
      varm.row(ep) = calc_var_ohlc(sub_ohlc, method, sub_close, scale, sub_index);
    }  // end if
  }  // end for
  
  // Old code below
  
  // Warmup period
  // for (arma::uword it = 1; it < lookb; it++) {
  //   arma::mat sub_ohlc = ohlc.rows(0, it);
  //   arma::colvec sub_close = closel.rows(0, it);
  //   arma::colvec sub_index = index.subvec(0, it);
  //   variance(it) = calc_var_ohlc(sub_ohlc, method, sub_close, scale, sub_index);
  // }  // end for
  
  // Remaining period
  // for (arma::uword it = lookb; it < nrows; it++) {
  //   arma::mat sub_ohlc = ohlc.rows(it-lookb+1, it);
  //   arma::colvec sub_close = closel.rows(it-lookb+1, it);
  //   arma::colvec sub_index = index.subvec(it-lookb+1, it);
  //   variance(it) = calc_var_ohlc(sub_ohlc, method, sub_close, scale, sub_index);
  // }  // end for
  
  return varm;
  
}  // end roll_var_ohlc




////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of skewness estimates over a rolling look-back
//' interval attached at the end points of a \emph{time series} or a
//' \emph{matrix}.
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'    
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endd} An \emph{integer} vector of end points (the default is 
//'   \code{endd = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{lookb} The number of end points in the look-back interval
//'   (the default is \code{lookb = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @param \code{method} A \emph{character string} specifying the type of the
//'   skewness model (the default is \code{method = "moment"} - see Details).
//'
//' @param \code{confl} The confidence level for calculating the quantiles of
//'   returns (the default is \code{confl = 0.75}).
//'
//' @return A \emph{matrix} of skewness estimates with the same number of
//'   columns as the input time series \code{tseries}, and the number of rows
//'   equal to the number of end points.
//'   
//' @details
//'   The function \code{roll_skew()} calculates a \emph{matrix} of skewness
//'   estimates over rolling look-back intervals attached at the end points of
//'   the \emph{time series} \code{tseries}.
//'   
//'   The function \code{roll_skew()} performs a loop over the end points, and
//'   at each end point it subsets the time series \code{tseries} over a
//'   look-back interval equal to \code{lookb} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_skew()}, which
//'   calculates the skewness.
//'   See the function \code{calc_skew()} for a description of the skewness
//'   methods.
//'   
//'   If the arguments \code{endd} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{tseries}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//' 
//'   For example, the rolling skewness at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{lookb = 3}.
//'
//'   The function \code{roll_skew()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, which makes it several times faster than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Define time series of returns using package rutils
//' retp <- na.omit(rutils::etfenv$returns$VTI)
//' # Define end points and start points
//' endd <- 1 + HighFreq::calc_endpoints(NROW(retp), step=25)
//' startp <- HighFreq::calc_startpoints(endd, lookb=3)
//' # Calculate the rolling skewness at 25 day end points, with a 75 day look-back
//' skewv <- HighFreq::roll_skew(retp, step=25, lookb=3)
//' # Calculate the rolling skewness using R code
//' skewr <- sapply(1:NROW(endd), function(it) {
//'   HighFreq::calc_skew(retp[startp[it]:endd[it], ])
//' })  # end sapply
//' # Compare the skewness estimates
//' all.equal(drop(skewv), skewr, check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_skew(retp, step=25, lookb=3),
//'   Rcode=sapply(1:NROW(endd), function(it) {
//'     HighFreq::calc_skew(retp[startp[it]:endd[it], ])
//'   }),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_skew(const arma::mat& tseries, 
                    arma::uvec startp = 0, 
                    arma::uvec endd = 0, 
                    arma::uword step = 1, 
                    arma::uword lookb = 1, 
                    arma::uword stub = 0,
                    std::string method = "moment", 
                    double confl = 0.75) {
  
  // Allocate end points
  arma::uword nrows = tseries.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endd) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endd;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by lookb
    startpts = calc_startpoints(endpts, lookb);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  // Allocate skewness matrix
  arma::uword numpts = endpts.n_elem;
  arma::mat skewv = arma::zeros(numpts, tseries.n_cols);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate skewness
    if (endpts(ep) > startpts(ep)) {
      skewv.row(ep) = calc_skew(tseries.rows(startpts(ep), endpts(ep)), method, confl);
    }  // end if
  }  // end for
  
  return skewv;
  
}  // end roll_skew



////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of kurtosis estimates over a rolling look-back
//' interval attached at the end points of a \emph{time series} or a
//' \emph{matrix}.
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'    
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endd} An \emph{integer} vector of end points (the default is 
//'   \code{endd = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{lookb} The number of end points in the look-back interval
//'   (the default is \code{lookb = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @param \code{method} A \emph{character string} specifying the type of the
//'   kurtosis model (the default is \code{method = "moment"} - see Details).
//'
//' @param \code{confl} The confidence level for calculating the quantiles of
//'   returns (the default is \code{confl = 0.75}).
//'
//' @return A \emph{matrix} of kurtosis estimates with the same number of
//'   columns as the input time series \code{tseries}, and the number of rows
//'   equal to the number of end points.
//'   
//' @details
//'   The function \code{roll_kurtosis()} calculates a \emph{matrix} of kurtosis
//'   estimates over rolling look-back intervals attached at the end points of
//'   the \emph{time series} \code{tseries}.
//'   
//'   The function \code{roll_kurtosis()} performs a loop over the end points,
//'   and at each end point it subsets the time series \code{tseries} over a
//'   look-back interval equal to \code{lookb} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_kurtosis()},
//'   which calculates the kurtosis. See the function \code{calc_kurtosis()} for
//'   a description of the kurtosis methods.
//'   
//'   If the arguments \code{endd} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{tseries}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//' 
//'   For example, the rolling kurtosis at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{lookb = 3}.
//'
//'   The function \code{roll_kurtosis()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, which makes it several times faster than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Define time series of returns using package rutils
//' retp <- na.omit(rutils::etfenv$returns$VTI)
//' # Define end points and start points
//' endd <- 1 + HighFreq::calc_endpoints(NROW(retp), step=25)
//' startp <- HighFreq::calc_startpoints(endd, lookb=3)
//' # Calculate the rolling kurtosis at 25 day end points, with a 75 day look-back
//' kurtosisv <- HighFreq::roll_kurtosis(retp, step=25, lookb=3)
//' # Calculate the rolling kurtosis using R code
//' kurt_r <- sapply(1:NROW(endd), function(it) {
//'   HighFreq::calc_kurtosis(retp[startp[it]:endd[it], ])
//' })  # end sapply
//' # Compare the kurtosis estimates
//' all.equal(drop(kurtosisv), kurt_r, check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_kurtosis(retp, step=25, lookb=3),
//'   Rcode=sapply(1:NROW(endd), function(it) {
//'     HighFreq::calc_kurtosis(retp[startp[it]:endd[it], ])
//'   }),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_kurtosis(const arma::mat& tseries, 
                        arma::uvec startp = 0, 
                        arma::uvec endd = 0, 
                        arma::uword step = 1, 
                        arma::uword lookb = 1, 
                        arma::uword stub = 0,
                        std::string method = "moment", 
                        double confl = 0.75) {
  
  // Allocate end points
  arma::uword nrows = tseries.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endd) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endd;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by lookb
    startpts = calc_startpoints(endpts, lookb);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  // Allocate kurtosis matrix
  arma::uword numpts = endpts.n_elem;
  arma::mat kurtosisv = arma::zeros(numpts, tseries.n_cols);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate kurtosis
    if (endpts(ep) > startpts(ep)) {
      kurtosisv.row(ep) = calc_kurtosis(tseries.rows(startpts(ep), endpts(ep)), method, confl);
    }  // end if
  }  // end for
  
  return kurtosisv;
  
}  // end roll_kurtosis



////////////////////////////////////////////////////////////
//' Perform a rolling regression and calculate a matrix of regression
//' coefficients, their t-values, and z-scores.
//' 
//' @param \code{respv} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{predm} A \emph{time series} or a \emph{matrix} of predictor
//'   data.
//'   
//' @param \code{controlv} A \emph{list} of model parameters (see Details).
//'
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endd} An \emph{integer} vector of end points (the default is 
//'   \code{endd = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{lookb} The number of end points in the look-back interval
//'   (the default is \code{lookb = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @return A \emph{matrix} with the regression coefficients, their t-values,
//'   and z-scores, and with the same number of rows as \code{predm} a
//'   number of columns equal to \code{2n+1}, where \code{n} is the number of
//'   columns of \code{predm}.
//'
//' @details
//'   The function \code{roll_reg()} performs a rolling regression over the end
//'   points of the predictor matrix, and calculates a \emph{matrix} of
//'   regression coefficients, their t-values, and z-scores.
//'   
//'   The function \code{roll_reg()} performs a loop over the end points, and at
//'   each end point it subsets the time series \code{predm} over a look-back
//'   interval equal to \code{lookb} number of end points.
//'   
//'   If the arguments \code{endd} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{predm}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//'   
//'   For example, the rolling regression at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{lookb = 3}.
//'
//'   It passes the subset time series to the function \code{calc_reg()}, which
//'   calculates the regression coefficients, their t-values, and the z-score.
//'   The function \code{roll_reg()} accepts a list of model parameters
//'   through the argument \code{controlv}, and passes it to the function
//'   \code{calc_reg()}.
//'   The list of model parameters can be created using the function
//'   \code{param_reg()}.  See the function \code{param_reg()} for a
//'   description of the model parameters.
//'   
//'   The number of columns of the return matrix depends on the number of
//'   columns of the predictor matrix (including the intercept column, if it's
//'   been added in \code{R}).
//'   The number of regression coefficients is equal to the number of columns of
//'   the predictor matrix.
//'   If the predictor matrix contains an intercept column then the first
//'   regression coefficient is equal to the intercept value \eqn{\alpha}.
//'   
//'   The number of columns of the return matrix is equal to the number of
//'   regression coefficients, plus their t-values, plus the z-score column.
//'   The number of t-values is equal to the number of coefficients.
//'   If the number of columns of the predictor matrix is equal to \code{n},
//'   then \code{roll_reg()} returns a matrix with \code{2n+1} columns: \code{n}
//'   regression coefficients, \code{n} corresponding t-values, and \code{1}
//'   z-score column.
//' 
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' predm <- na.omit(rutils::etfenv$returns[, c("XLP", "VTI")])
//' # Add unit intercept column to the predictor matrix
//' predm <- cbind(rep(1, NROW(predm)), predm)
//' # Define monthly end points and start points
//' endd <- xts::endpoints(predm, on="months")[-1]
//' lookb <- 12
//' startp <- c(rep(1, lookb), endd[1:(NROW(endd)-lookb)])
//' # Create a default list of regression parameters
//' controlv <- HighFreq::param_reg()
//' # Calculate rolling betas using RcppArmadillo
//' regroll <- HighFreq::roll_reg(respv=predm[, 2], predm=predm[, -2], endd=(endd-1), startp=(startp-1), controlv=controlv)
//' betas <- regroll[, 2]
//' # Calculate rolling betas in R
//' betar <- sapply(1:NROW(endd), FUN=function(ep) {
//'   datav <- predm[startp[ep]:endd[ep], ]
//'   # HighFreq::calc_reg(datav[, 2], datav[, -2], controlv)
//'   drop(cov(datav[, 2], datav[, 3])/var(datav[, 3]))
//' })  # end sapply
//' # Compare the outputs of both functions
//' all.equal(betas, betar, check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_reg(const arma::mat& respv, // Response vector
                   const arma::mat& predm, // Predictor matrix
                   Rcpp::List controlv, 
                   arma::uvec startp = 0, 
                   arma::uvec endd = 0, 
                   arma::uword step = 1, 
                   arma::uword lookb = 1, 
                   arma::uword stub = 0) {
  
  // Allocate end points
  arma::uword nrows = predm.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endd) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endd;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by lookb
    startpts = calc_startpoints(endpts, lookb);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  
  // Allocate regression matrix
  arma::mat responsi;
  arma::mat predicti;
  arma::uword numpts = endpts.n_elem;
  arma::uword ncols = predm.n_cols;
  // Add unit intercept column to the predictor matrix?
  // bool intercept = Rcpp::as<int>(controlv["intercept"]);
  // if (intercept) ncols += 1;
  arma::mat regroll(numpts, (2*ncols + 1), fill::zeros);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate regression coefficients
    if (endpts(ep) > startpts(ep)) {
      // cout << "ep: " << ep << endl;
      responsi = respv.rows(startpts(ep), endpts(ep));
      predicti = predm.rows(startpts(ep), endpts(ep));
      regroll.row(ep) = calc_reg(responsi, predicti, controlv);
    }  // end if
  }  // end for
  
  // Warmup period
  // regroll.rows(0, ncols+1) = zeros(ncols+2, (ncols + 1));
  // for (arma::uword it = (ncols+2); it < lookb; it++) {
  //   responsi = respv.rows(0, it);
  //   predicti = predm.rows(0, it);
  //   reg_data = calc_reg(responsi, predicti);
  //   regroll.row(it) = arma::conv_to<rowvec>::from(reg_data);
  // }  // end for
  
  // Remaining periods
  // for (arma::uword it = lookb; it < nrows; it++) {
  //   responsi = respv.rows(it-lookb+1, it);
  //   predicti = predm.rows(it-lookb+1, it);
  //   reg_data = calc_reg(responsi, predicti, method, singmin, dimax, confl, alpha);
  //   regroll.row(it) = arma::conv_to<rowvec>::from(reg_data);
  // }  // end for
  
  return regroll;
  
}  // end roll_reg



////////////////////////////////////////////////////////////
//' Perform a rolling standardization (centering and scaling) of the columns of
//' a \emph{time series} of data using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A \emph{time series} or \emph{matrix} of data.
//' 
//' @param \code{lookb} The length of the look-back interval, equal to the
//'   number of rows of data used in the scaling.
//'   
//' @param \code{center} A \emph{Boolean} argument: if \code{TRUE} then center
//'   the columns so that they have zero mean or median (the default is
//'   \code{TRUE}).
//' 
//' @param \code{scale} A \emph{Boolean} argument: if \code{TRUE} then scale the
//'   columns so that they have unit standard deviation or MAD (the default is
//'   \code{TRUE}).
//' 
//' @param \code{use_median} A \emph{Boolean} argument: if \code{TRUE} then the 
//'   centrality (central tendency) is calculated as the \emph{median} and the 
//'   dispersion is calculated as the \emph{median absolute deviation}
//'   (\emph{MAD}) (the default is \code{FALSE}).
//'   If \code{use_median = FALSE} then the centrality is calculated as the
//'   \emph{mean} and the dispersion is calculated as the \emph{standard
//'   deviation}.
//'
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{roll_scale()} performs a rolling standardization
//'   (centering and scaling) of the columns of the \code{tseries} argument
//'   using \code{RcppArmadillo}.
//'   The function \code{roll_scale()} performs a loop over the rows of
//'   \code{tseries}, subsets a number of previous (past) rows equal to
//'   \code{lookb}, and standardizes the subset matrix by calling the
//'   function \code{calc_scale()}.  It assigns the last row of the standardized
//'   subset \emph{matrix} to the return matrix.
//'   
//'   If the arguments \code{center} and \code{scale} are both \code{TRUE} and
//'   \code{use_median} is \code{FALSE} (the defaults), then
//'   \code{calc_scale()} performs the same calculation as the function
//'   \code{roll::roll_scale()}.
//'   
//'   If the arguments \code{center} and \code{scale} are both \code{TRUE} (the
//'   defaults), then \code{calc_scale()} standardizes the data.
//'   If the argument \code{center} is \code{FALSE} then \code{calc_scale()}
//'   only scales the data (divides it by the standard deviations).
//'   If the argument \code{scale} is \code{FALSE} then \code{calc_scale()}
//'   only demeans the data (subtracts the means).
//'   
//'   If the argument \code{use_median} is \code{TRUE}, then it calculates the
//'   centrality as the \emph{median} and the dispersion as the \emph{median
//'   absolute deviation} (\emph{MAD}).
//'   
//' @examples
//' \dontrun{
//' # Calculate a time series of returns
//' retp <- zoo::coredata(na.omit(rutils::etfenv$returns[, c("IEF", "VTI")]))
//' lookb <- 11
//' rolled_scaled <- roll::roll_scale(retp, width=lookb, min_obs=1)
//' rolled_scaled2 <- HighFreq::roll_scale(retp, lookb=lookb)
//' all.equal(rolled_scaled[-(1:2), ], rolled_scaled2[-(1:2), ],
//'   check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_scale(const arma::mat& matrix, 
                     arma::uword lookb,
                     bool center = true, 
                     bool scale = true, 
                     bool use_median = false) {
  
  arma::uword nrows = matrix.n_rows;
  arma::mat scaledmat(nrows, matrix.n_cols);
  arma::mat sub_mat;
  
  // Warmup period
  scaledmat.row(0) = matrix.row(0);
  for (arma::uword it = 1; it < lookb; it++) {
    sub_mat = matrix.rows(0, it);
    calc_scale(sub_mat, center, scale, use_median);
    scaledmat.row(it) = sub_mat.row(sub_mat.n_rows-1);
  }  // end for
  
  // Perform loop over the remaining rows
  for (arma::uword it = lookb; it < nrows; it++) {
    sub_mat = matrix.rows(it-lookb+1, it);
    calc_scale(sub_mat, center, scale, use_median);
    scaledmat.row(it) = sub_mat.row(sub_mat.n_rows-1);
  }  // end for
  
  return scaledmat;
}  // end roll_scale



////////////////////////////////////////////////////////////
//' Standardize (center and scale) the columns of a \emph{time series} of data
//' over time and in place, without copying the data in memory, using
//' \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A \emph{time series} or \emph{matrix} of data.
//' 
//' @param \code{lambda} A decay factor which multiplies past estimates.
//' 
//' @param \code{center} A \emph{Boolean} argument: if \code{TRUE} then center
//'   the columns so that they have zero mean or median (the default is
//'   \code{TRUE}).
//' 
//' @param \code{scale} A \emph{Boolean} argument: if \code{TRUE} then scale the
//'   columns so that they have unit standard deviation or MAD (the default is
//'   \code{TRUE}).
//' 
//' @return Void (no return value - modifies the data in place).
//'
//' @details
//'   The function \code{run_scale()} performs a trailing standardization
//'   (centering and scaling) of the columns of the \code{tseries} argument
//'   using \code{RcppArmadillo}.
//' 
//'   The function \code{run_scale()} accepts a \emph{pointer} to the argument
//'   \code{tseries}, and it overwrites the old data with the standardized
//'   data. It performs the calculation in place, without copying the data in
//'   memory, which can significantly increase the computation speed for large
//'   time series.
//'
//'   The function \code{run_scale()} performs a loop over the rows of
//'   \code{tseries}, and standardizes the data using its trailing means and
//'   standard deviations.
//'
//'   The function \code{run_scale()} calculates the trailing mean and variance
//'   of streaming \emph{time series} data \eqn{r_t}, by recursively weighting
//'   the past estimates with the new data, using the decay factor \eqn{\lambda}:
//'   \deqn{
//'     \bar{r}_t = \lambda \bar{r}_{t-1} + (1-\lambda) r_t
//'   }
//'   \deqn{
//'     \sigma^2_t = \lambda \sigma^2_{t-1} + (1-\lambda) (r_t - \bar{r}_t)^2
//'   }
//'   Where \eqn{\bar{r}_t} is the trailing mean and \eqn{\sigma^2_t} is the
//'   trailing variance.
//'   
//'   It then calculates the standardized data as follows:
//'   \deqn{
//'     r^{\prime}_t = \frac{r_t - \bar{r}_t}{\sigma_t}
//'   }
//'
//'   If the arguments \code{center} and \code{scale} are both \code{TRUE} (the
//'   defaults), then \code{calc_scale()} standardizes the data.
//'   If the argument \code{center} is \code{FALSE} then \code{calc_scale()}
//'   only scales the data (divides it by the standard deviations).
//'   If the argument \code{scale} is \code{FALSE} then \code{calc_scale()}
//'   only demeans the data (subtracts the means).
//'   
//'   The value of the decay factor \eqn{\lambda} must be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the trailing variance values have a
//'   stronger dependence on past data.  This is equivalent to a long
//'   look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the trailing variance values have a
//'   weaker dependence on past data.  This is equivalent to a short look-back
//'   interval.
//' 
//'   The above online recursive formulas are convenient for processing live
//'   streaming data because they don't require maintaining a buffer of past
//'   data.
//'   The formulas are equivalent to a convolution with exponentially decaying
//'   weights, but they're much faster to calculate.
//'   Using exponentially decaying weights is more natural than using a sliding
//'   look-back interval, because it gradually "forgets" about the past data.
//' 
//'   The function \code{run_scale()} uses \code{RcppArmadillo} \code{C++} code,
//'   so it can be over \code{100} times faster than the equivalent \code{R}
//'   code.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' retp <- na.omit(rutils::etfenv$returns[, c("XLF", "VTI")])
//' # Calculate the trailing standardized returns using R code
//' lambdaf <- 0.9 # Decay factor
//' lambda1 <- 1 - lambdaf
//' scaled <- zoo::coredata(retp)
//' meanm <- scaled[1, ];
//' vars <- scaled[1, ]^2;
//' for (it in 2:NROW(retp)) {
//'   meanm <- lambdaf*meanm + lambda1*scaled[it, ];
//'   vars <- lambdaf*vars + lambda1*(scaled[it, ] - meanm)^2;
//'   scaled[it, ] <- (scaled[it, ] - meanm)/sqrt(vars)
//' }  # end for
//' # Calculate the trailing standardized returns using C++ code
//' HighFreq::run_scale(retp, lambda=lambdaf)
//' all.equal(zoo::coredata(retp), scaled, check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::run_scale(retp, lambda=lambdaf),
//'   Rcode={for (it in 2:NROW(retp)) {
//'    meanm <- lambdaf*meanm + lambda1*scaled[it, ];
//'    vars <- lambdaf*vars + lambda1*(scaled[it, ] - meanm)^2;
//'    scaled[it, ] <- (scaled[it, ] - meanm)/sqrt(vars)
//'   }},  # end for
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
void run_scale(arma::mat& tseries, 
               double lambda, // Decay factor which multiplies the past values
               bool center = true, 
               bool scale = true) {
  
  arma::uword nrows = tseries.n_rows;
  double lambda1 = 1-lambda;
  arma::mat meanm = tseries.row(0);
  arma::mat vars = arma::square(tseries.row(0));
  
  if (scale and center) {
    for (arma::uword it = 1; it < nrows; it++) {
      meanm = lambda*meanm + lambda1*tseries.row(it);
      vars = lambda*vars + lambda1*arma::square(tseries.row(it) - meanm);
      tseries.row(it) = (tseries.row(it) - meanm)/arma::sqrt(vars);
    }  // end for
  } else if (scale and (not center)) {
    for (arma::uword it = 1; it < nrows; it++) {
      meanm = lambda*meanm + lambda1*tseries.row(it);
      vars = lambda*vars + lambda1*arma::square(tseries.row(it) - meanm);
      tseries.row(it) = tseries.row(it)/arma::sqrt(vars);
    }  // end for
  } else if ((not scale) and center) {
    for (arma::uword it = 1; it < nrows; it++) {
      meanm = lambda*meanm + lambda1*tseries.row(it);
      tseries.row(it) = (tseries.row(it) - meanm);
    }  // end for
  }  // end if
  
}  // end run_scale



////////////////////////////////////////////////////////////
//' Calculate a \emph{vector} of z-scores of the residuals of rolling
//' regressions at the end points of the predictor matrix.
//' 
//' @param \code{respv} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{predm} A \emph{time series} or a \emph{matrix} of predictor
//'   data.
//'   
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endd} An \emph{integer} vector of end points (the default is
//'   \code{endd = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{lookb} The number of end points in the look-back interval
//'   (the default is \code{lookb = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @return A column \emph{vector} of the same length as the number of rows of
//'   \code{predm}.
//'
//' @details
//'   The function \code{roll_zscores()} calculates a \emph{vector} of z-scores
//'   of the residuals of rolling regressions at the end points of the
//'   \emph{time series} \code{predm}.
//'   
//'   The function \code{roll_zscores()} performs a loop over the end points,
//'   and at each end point it subsets the time series \code{predm} over a
//'   look-back interval equal to \code{lookb} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_lm()}, which
//'   calculates the regression data.
//'   
//'   If the arguments \code{endd} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{predm}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//'   
//'   For example, the rolling variance at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{lookb = 3}.
//'
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' retp <- na.omit(rutils::etfenv$returns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' respv <- retp[, 1]
//' # Predictor matrix equals VTI and IEF returns
//' predm <- retp[, -1]
//' # Calculate Z-scores from rolling time series regression using RcppArmadillo
//' lookb <- 11
//' zscores <- HighFreq::roll_zscores(respv=respv, predm=predm, lookb)
//' # Calculate z-scores in R from rolling multivariate regression using lm()
//' zscoresr <- sapply(1:NROW(predm), function(ro_w) {
//'   if (ro_w == 1) return(0)
//'   startpoint <- max(1, ro_w-lookb+1)
//'   responsi <- response[startpoint:ro_w]
//'   predicti <- predictor[startpoint:ro_w, ]
//'   regmod <- lm(responsi ~ predicti)
//'   residuals <- regmod$residuals
//'   residuals[NROW(residuals)]/sd(residuals)
//' })  # end sapply
//' # Compare the outputs of both functions
//' all.equal(zscores[-(1:lookb)], zscoresr[-(1:lookb)], 
//'   check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec roll_zscores(const arma::mat& respv, // Response vector
                       const arma::mat& predm, // Predictor matrix
                       arma::uvec startp = 0, 
                       arma::uvec endd = 0, 
                       arma::uword step = 1, 
                       arma::uword lookb = 1,
                       arma::uword stub = 0) {
  
  // Allocate end points
  arma::uword nrows = predm.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endd) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endd;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by lookb
    startpts = calc_startpoints(endpts, lookb);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  // Allocate regression matrix
  arma::mat responsi;
  arma::mat predicti;
  arma::uword numpts = endpts.n_elem;
  arma::vec zscores(numpts, fill::zeros);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate z-scores
    if (endpts(ep) > startpts(ep)) {
      responsi = respv.rows(startpts(ep), endpts(ep));
      predicti = predm.rows(startpts(ep), endpts(ep));
      zscores(ep) = calc_lm(responsi, predicti)["zscore"];
    }  // end if
  }  // end for
  
  // Old code below
  // Warmup period
  // for (arma::uword it = 1; it < lookb; it++) {
  //   responsi = respv.rows(0, it);
  //   predicti = predm.rows(0, it);
  //   zscores(it) = calc_lm(responsi, predicti)["zscore"];
  // }  // end for
  
  // Remaining periods
  // for (arma::uword it = lookb; it < nrows; it++) {
  //   responsi = respv.rows(it-lookb+1, it);
  //   predicti = predm.rows(it-lookb+1, it);
  //   zscores(it) = calc_lm(responsi, predicti)["zscore"];
  // }  // end for
  
  return zscores;
  
}  // end roll_zscores



// Define type for pointer to calc_* function
typedef arma::mat (*momptr)(const arma::mat&, std::string, double);


////////////////////////////////////////////////////////////
//' Calculate a pointer to a moment function from the function name (string).
//' 
//' @param \code{funame} A \emph{character} \emph{string} specifying the
//'   function name.
//'   
//' @return A pointer to a moment function.
//'
//' @details
//'   The function \code{calc_momptr()} calculates a pointer to a moment
//'   function from the function name (string).
//'   The function pointer is used internally in the \code{C++} code, but it's
//'   not exported to \code{R}.
//'   A moment function takes three arguments:
//'   \itemize{
//'     \item A \emph{time series} or a \emph{matrix} of data,
//'     \item A \emph{character string} specifying the type of the moment,
//'     \item A number specifying the confidence level for calculating the
//'   quantiles of returns.
//'   The function name must be one of the following:
//'   \itemize{
//'     \item "calc_mean" for the estimator of the mean (location),
//'     \item "calc_var" for the estimator of the dispersion (variance),
//'     \item "calc_skew" for the estimator of the skewness,
//'     \item "calc_kurtosis" for the estimator of the kurtosis.
//'    }
//'    (The default is the \code{funame = "calc_mean"}.)
//'    }
//'
//' @export
momptr calc_momptr(std::string funame = "calc_mean") {
  if (funame == "calc_mean")
    return (&calc_mean);
  else if (funame == "calc_var")
    return (&calc_var);
  else if (funame == "calc_skew")
    return (&calc_skew);
  else if (funame == "calc_kurtosis")
    return (&calc_kurtosis);
  else
    throw std::invalid_argument("No such function!");
}  // end calc_momptr



////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of moment values over a rolling look-back
//' interval attached at the end points of a \emph{time series} or a
//' \emph{matrix}.
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'    
//' @param \code{funame} A \emph{character string} specifying the moment
//'   function (the default is \code{funame = "calc_mean"}).
//' 
//' @param \code{method} A \emph{character string} specifying the type of the
//'   model for the moment (the default is \code{method = "moment"}).
//'
//' @param \code{confl} The confidence level for calculating the quantiles of
//'   returns (the default is \code{confl = 0.75}).
//'
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endd} An \emph{integer} vector of end points (the default is 
//'   \code{endd = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{lookb} The number of end points in the look-back interval
//'   (the default is \code{lookb = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//'
//' @return A \emph{matrix} with the same number of columns as the input time
//'   series \code{tseries}, and the number of rows equal to the number of end
//'   points.
//'   
//' @details
//'   The function \code{roll_moment()} calculates a \emph{matrix} of moment
//'   values, over rolling look-back intervals attached at the end points of the
//'   \emph{time series} \code{tseries}.
//'   
//'   The function \code{roll_moment()} performs a loop over the end points, and
//'   at each end point it subsets the time series \code{tseries} over a
//'   look-back interval equal to \code{lookb} number of end points.
//'   
//'   It passes the subset time series to the function specified by the argument
//'   \code{funame}, which calculates the statistic.
//'   See the functions \code{calc_*()} for a description of the different
//'   moments.
//'   The function name must be one of the following:
//'   \itemize{
//'     \item "calc_mean" for the estimator of the mean (location),
//'     \item "calc_var" for the estimator of the dispersion (variance),
//'     \item "calc_skew" for the estimator of the skewness,
//'     \item "calc_kurtosis" for the estimator of the kurtosis.
//'    }
//'    (The default is the \code{funame = "calc_mean"}).
//'   
//'   If the arguments \code{endd} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{tseries}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//' 
//'   For example, the rolling variance at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{lookb = 3}.
//'
//'   The function \code{roll_moment()} calls the function \code{calc_momptr()}
//'   to calculate a pointer to a moment function from the function name
//'   \code{funame} (string). The function pointer is used internally in the
//'   \code{C++} code, but the function \code{calc_momptr()} is not exported to
//'   \code{R}.
//'   
//'   The function \code{roll_moment()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, which makes it several times faster than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Define time series of returns using package rutils
//' retp <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the rolling variance at 25 day end points, with a 75 day look-back
//' var_rollfun <- HighFreq::roll_moment(retp, fun="calc_var", step=25, lookb=3)
//' # Calculate the rolling variance using roll_var()
//' var_roll <- HighFreq::roll_var(retp, step=25, lookb=3)
//' # Compare the two methods
//' all.equal(var_rollfun, var_roll, check.attributes=FALSE)
//' # Define end points and start points
//' endd <- HighFreq::calc_endpoints(NROW(retp), step=25)
//' startp <- HighFreq::calc_startpoints(endd, lookb=3)
//' # Calculate the rolling variance using RcppArmadillo
//' var_rollfun <- HighFreq::roll_moment(retp, fun="calc_var", startp=startp, endd=endd)
//' # Calculate the rolling variance using R code
//' var_roll <- sapply(1:NROW(endd), function(it) {
//'   var(retp[startp[it]:endd[it]+1, ])
//' })  # end sapply
//' var_roll[1] <- 0
//' # Compare the two methods
//' all.equal(drop(var_rollfun), var_roll, check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_moment(retp, fun="calc_var", startp=startp, endd=endd),
//'   Rcode=sapply(1:NROW(endd), function(it) {
//'     var(retp[startp[it]:endd[it]+1, ])
//'   }),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_moment(const arma::mat& tseries, 
                      std::string funame = "calc_mean",
                      std::string method = "moment", 
                      double confl = 0.75, 
                      arma::uvec startp = 0, 
                      arma::uvec endd = 0, 
                      arma::uword step = 1, 
                      arma::uword lookb = 1, 
                      arma::uword stub = 0) {
  
  // Allocate end points
  arma::uword nrows = tseries.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate the end points if missing
  if (sum(endd) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endd;
  }  // end if
  
  // Calculate the start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by lookb
    startpts = calc_startpoints(endpts, lookb);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  // Allocate matrix of statistics
  arma::uword numpts = endpts.n_elem;
  arma::mat stats = arma::zeros(numpts, tseries.n_cols);
  
  // Calculate a function pointer from the function name (string)
  momptr momfun = calc_momptr(funame);
  // momptr momfun = *momptr;
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate the statistics at the end point
    if (endpts(ep) > startpts(ep)) {
      stats.row(ep) = momfun(tseries.rows(startpts(ep), endpts(ep)), method, confl);
    }  // end if
  }  // end for
  
  return stats;
  
}  // end roll_moment




////////////////////////////////////////////////////////////
// Functions for simulation
////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////
//' Simulate or estimate the rolling variance under a \emph{GARCH(1,1)} process
//' using \emph{Rcpp}.
//' 
//' @param \code{omega} Parameter proportional to the long-term average level
//'   of variance.
//' 
//' @param \code{alpha} The weight associated with recent realized variance
//'   updates.
//' 
//' @param \code{beta} The weight associated with the past variance estimates.
//' 
//' @param \code{innov} A single-column \emph{matrix} of innovations.
//' 
//' @param \code{is_random} \emph{Boolean} argument: Are the innovations random
//'   numbers or historical returns? (The default is \code{is_random = TRUE}.)
//'
//' @return A \emph{matrix} with two columns and with the same number of rows as
//'   the argument \code{innov}.  The first column are the simulated returns and
//'   the second column is the variance.
//'
//' @details
//'   The function \code{sim_garch()} simulates or estimates the rolling variance
//'   under a \emph{GARCH(1,1)} process using \emph{Rcpp}.
//'
//'   If \code{is_random = TRUE} (the default) then the innovations \code{innov}
//'   are treated as random numbers \eqn{\xi_i} and the \emph{GARCH(1,1)}
//'   process is given by:
//'   \deqn{
//'     r_i = \sigma_{i-1} \xi_i
//'   }
//'   \deqn{
//'     \sigma^2_i = \omega + \alpha r^2_i + \beta \sigma_{i-1}^2
//'   }
//'   Where \eqn{r_i} and \eqn{\sigma^2_i} are the simulated returns and
//'   variance, and \eqn{\omega}, \eqn{\alpha}, and \eqn{\beta} are the
//'   \emph{GARCH} parameters, and \eqn{\xi_i} are standard normal
//'   \emph{innovations}.
//'
//'   The long-term equilibrium level of the simulated variance is proportional
//'   to the parameter \eqn{\omega}:
//'   \deqn{
//'     \sigma^2 = \frac{\omega}{1 - \alpha - \beta}
//'   }
//'   So the sum of \eqn{\alpha} plus \eqn{\beta} should be less than \eqn{1},
//'   otherwise the volatility becomes explosive.
//'   
//'   If \code{is_random = FALSE} then the function \code{sim_garch()}
//'   \emph{estimates} the rolling variance from the historical returns. The
//'   innovations \code{innov} are equal to the historical returns \eqn{r_i} and
//'   the \emph{GARCH(1,1)} process is simply:
//'   \deqn{
//'     \sigma^2_i = \omega + \alpha r^2_i + \beta \sigma_{i-1}^2
//'   }
//'   Where \eqn{\sigma^2_i} is the rolling variance.
//'   
//'   The above should be viewed as a formula for \emph{estimating} the rolling
//'   variance from the historical returns, rather than simulating them. It
//'   represents exponential smoothing of the squared returns with a decay
//'   factor equal to \eqn{\beta}.
//'
//'   The function \code{sim_garch()} simulates the \emph{GARCH} process using
//'   fast \emph{Rcpp} \code{C++} code.
//'
//' @examples
//' \dontrun{
//' # Define the GARCH model parameters
//' alpha <- 0.79
//' betav <- 0.2
//' om_ega <- 1e-4*(1-alpha-betav)
//' # Calculate matrix of standard normal innovations
//' innov <- matrix(rnorm(1e3))
//' # Simulate the GARCH process using Rcpp
//' garch_data <- HighFreq::sim_garch(omega=om_ega, alpha=alpha,  beta=betav, innov=innov)
//' # Plot the GARCH rolling volatility and cumulative returns
//' plot(sqrt(garch_data[, 2]), t="l", main="Simulated GARCH Volatility", ylab="volatility")
//' plot(cumsum(garch_data[, 1]), t="l", main="Simulated GARCH Cumulative Returns", ylab="cumulative returns")
//' # Calculate historical VTI returns
//' retp <- na.omit(rutils::etfenv$returns$VTI)
//' # Estimate the GARCH volatility of VTI returns
//' garch_data <- HighFreq::sim_garch(omega=om_ega, alpha=alpha,  beta=betav, 
//'   innov=retp, is_random=FALSE)
//' # Plot dygraph of the estimated GARCH volatility
//' dygraphs::dygraph(xts::xts(sqrt(garch_data[, 2]), index(retp)), 
//'   main="Estimated GARCH Volatility of VTI")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat sim_garch(double omega, 
                    double alpha, 
                    double beta, 
                    arma::mat& innov,
                    bool is_random = true) {
  
  arma::uword nrows = innov.n_rows;
  
  if (is_random) {
    // The innovations are random numbers
    arma::mat varm(nrows, 1);
    varm(0) = omega/(1-alpha-beta);
    arma::mat returns(nrows, 1);
    returns(0) = std::sqrt(varm(0))*innov(0);
    
    for (arma::uword it = 1; it < nrows; it++) {
      returns(it) = std::sqrt(varm(it-1))*innov(it);
      varm(it) = omega + alpha*pow(returns(it), 2) + beta*varm(it-1);
    }  // end for
    return arma::join_rows(returns, varm);
  } else {
    // The innovations are historical returns
    arma::mat varm = arma::square(innov);
    varm(0) = omega/(1-alpha-beta);
    for (arma::uword it = 1; it < nrows; it++) {
      varm(it) = omega + alpha*varm(it) + beta*varm(it-1);
    }  // end for
    return arma::join_rows(innov, varm);
  }  // end if
  
}  // end sim_garch



////////////////////////////////////////////////////////////
//' Simulate an \emph{Ornstein-Uhlenbeck} process using \emph{Rcpp}.
//' 
//' @param \code{init_price} The initial price. 
//' 
//' @param \code{eq_price} The equilibrium price. 
//' 
//' @param \code{theta} The strength of mean reversion.
//' 
//' @param \code{innov} A single-column \emph{matrix} of innovations (random
//'   numbers).
//' 
//' @return A single-column \emph{matrix} of simulated prices, with the same
//'   number of rows as the argument \code{innov}.
//'
//' @details
//'   The function \code{sim_ou()} simulates the following
//'   \emph{Ornstein-Uhlenbeck} process:
//'   \deqn{
//'     r_i = p_i - p_{i-1} = \theta \, (\mu - p_{i-1}) + \xi_i
//'   }
//'   \deqn{
//'     p_i = p_{i-1} + r_i
//'   }
//'   Where \eqn{r_i} and \eqn{p_i} are the simulated returns and prices,
//'   \eqn{\theta}, \eqn{\mu}, and \eqn{\sigma} are the
//'   \emph{Ornstein-Uhlenbeck} parameters, and \eqn{\xi_i} are the standard
//'   \emph{innovations}.
//'   The recursion starts with: \eqn{r_1 = \xi_1} and \eqn{p_1 = init\_price}.
//'
//'   The function \code{sim_ou()} simulates the percentage returns as equal to
//'   the difference between the equilibrium price \eqn{\mu} minus the latest
//'   price \eqn{p_{i-1}}, times the mean reversion parameter \eqn{\theta}, plus
//'   a random normal innovation. The log prices are calculated as the sum of
//'   returns (not compounded), so they can become negative.
//'
//'   The function \code{sim_ou()} simulates the \emph{Ornstein-Uhlenbeck}
//'   process using fast \emph{Rcpp} \code{C++} code.
//'
//'   The function \code{sim_ou()} returns a single-column \emph{matrix}
//'   representing the \emph{time series} of simulated prices.
//'
//' @examples
//' \dontrun{
//' # Define the Ornstein-Uhlenbeck model parameters
//' init_price <- 0.0
//' eq_price <- 1.0
//' sigmav <- 0.01
//' thetav <- 0.01
//' innov <- matrix(rnorm(1e3))
//' # Simulate Ornstein-Uhlenbeck process using Rcpp
//' prices <- HighFreq::sim_ou(init_price=init_price, eq_price=eq_price, volat=sigmav, theta=thetav, innov=innov)
//' plot(prices, t="l", main="Simulated Ornstein-Uhlenbeck Prices", ylab="prices")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat sim_ou(double init_price, 
                 double eq_price,
                 double theta, 
                 arma::mat& innov) {
  
  arma::uword nrows = innov.n_rows;
  arma::mat pricev = arma::zeros(nrows, 1);
  arma::mat retv = arma::zeros(nrows, 1);
  
  retv.row(0) = innov.row(0);
  pricev.row(0) = init_price;
  for (arma::uword it = 1; it < nrows; it++) {
    retv.row(it) = theta*(eq_price - pricev.row(it-1)) + innov.row(it);
    pricev.row(it) = pricev.row(it-1) + retv.row(it);
  }  // end for
  
  return pricev;
  
}  // end sim_ou



////////////////////////////////////////////////////////////
//' Simulate a \emph{Schwartz} process using \emph{Rcpp}.
//' 
//' @param \code{init_price} The initial price. 
//' 
//' @param \code{eq_price} The equilibrium price. 
//' 
//' @param \code{theta} The strength of mean reversion.
//' 
//' @param \code{innov} A single-column \emph{matrix} of innovations (random
//'   numbers).
//' 
//' @return A single-column \emph{matrix} of simulated prices, with the same
//'   number of rows as the argument \code{innov}.
//'
//' @details
//'   The function \code{sim_schwartz()} simulates a \emph{Schwartz} process
//'   using fast \emph{Rcpp} \code{C++} code.
//'   
//'   The \emph{Schwartz} process is the exponential of the
//'   \emph{Ornstein-Uhlenbeck} process, and similar comments apply to it.
//'   The prices are calculated as the exponentially compounded returns, so they
//'   are never negative. The log prices can be obtained by taking the logarithm
//'   of the prices.
//'   
//'   The function \code{sim_schwartz()} simulates the percentage returns as
//'   equal to the difference between the equilibrium price \eqn{\mu} minus the
//'   latest price \eqn{p_{i-1}}, times the mean reversion parameter
//'   \eqn{\theta}, plus a random normal innovation.
//'
//'   The function \code{sim_schwartz()} returns a single-column \emph{matrix}
//'   representing the \emph{time series} of simulated prices.
//'
//' @examples
//' \dontrun{
//' # Define the Schwartz model parameters
//' init_price <- 1.0
//' eq_price <- 2.0
//' thetav <- 0.01
//' innov <- matrix(rnorm(1e3, sd=0.01))
//' # Simulate Schwartz process using Rcpp
//' prices <- HighFreq::sim_schwartz(init_price=init_price, eq_price=eq_price, theta=thetav, innov=innov)
//' plot(prices, t="l", main="Simulated Schwartz Prices", ylab="prices")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat sim_schwartz(double init_price, 
                       double eq_price, 
                       double theta, 
                       arma::mat& innov) {
  
  arma::uword nrows = innov.n_rows;
  arma::mat prices = arma::zeros(nrows, 1);
  arma::mat returns = arma::zeros(nrows, 1);
  
  returns.row(0) = innov.row(0);
  prices.row(0) = init_price;
  for (arma::uword it = 1; it < nrows; it++) {
    returns.row(it) = theta*(eq_price - prices.row(it-1)) + innov.row(it);
    prices.row(it) = prices.row(it-1) * exp(returns.row(it));
  }  // end for
  
  return prices;
  
}  // end sim_schwartz



////////////////////////////////////////////////////////////
//' Simulate \emph{autoregressive} returns by recursively filtering a
//' \emph{matrix} of innovations through a \emph{matrix} of
//' \emph{autoregressive} coefficients.
//' 
//' @param \code{innov} A single-column \emph{matrix} of innovations.
//' 
//' @param \code{coeff} A single-column \emph{matrix} of \emph{autoregressive}
//'   coefficients.
//'
//' @return A single-column \emph{matrix} of simulated returns, with the same
//'   number of rows as the argument \code{innov}.
//'
//' @details
//'   The function \code{sim_ar()} recursively filters the \emph{matrix} of
//'   innovations \code{innov} through the \emph{matrix} of
//'   \emph{autoregressive} coefficients \code{coeff}, using fast
//'   \code{RcppArmadillo} \code{C++} code.
//'
//'   The function \code{sim_ar()} simulates an \emph{autoregressive} process
//'   \eqn{AR(n)} of order \eqn{n}:
//'   \deqn{
//'     r_i = \varphi_1 r_{i-1} + \varphi_2 r_{i-2} + \ldots + \varphi_n r_{i-n} + \xi_i
//'   }
//'   Where \eqn{r_i} is the simulated output time series, \eqn{\varphi_i} are
//'   the \emph{autoregressive} coefficients, and \eqn{\xi_i} are the standard
//'   normal \emph{innovations}.
//'
//'   The order \eqn{n} of the \emph{autoregressive} process \eqn{AR(n)}, is
//'   equal to the number of rows of the \emph{autoregressive} coefficients
//'   \code{coeff}.
//'
//'   The function \code{sim_ar()} performs the same calculation as the standard
//'   \code{R} function \cr\code{filter(x=innov, filter=coeff,
//'   method="recursive")}, but it's several times faster.
//'   
//' @examples
//' \dontrun{
//' # Define AR coefficients
//' coeff <- matrix(c(0.1, 0.3, 0.5))
//' # Calculate matrix of innovations
//' innov <- matrix(rnorm(1e4, sd=0.01))
//' # Calculate recursive filter using filter()
//' innof <- filter(innov, filter=coeff, method="recursive")
//' # Calculate recursive filter using RcppArmadillo
//' retp <- HighFreq::sim_ar(coeff, innov)
//' # Compare the two methods
//' all.equal(as.numeric(retp), as.numeric(innof))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::sim_ar(coeff, innov),
//'   Rcode=filter(innov, filter=coeff, method="recursive"),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat sim_ar(arma::mat& coeff, const arma::mat& innov) {
  
  arma::uword nrows = innov.n_rows;
  arma::uword ncoeff = coeff.n_rows;
  arma::mat coeffr = arma::reverse(coeff);
  arma::mat returns = arma::zeros(nrows, 1);

  // Warmup period
  returns.row(0) = innov.row(0);
  returns.row(1) = innov.row(1) + coeffr.row(ncoeff-1) * returns.row(0);
  for (arma::uword it = 2; it < ncoeff; it++) {
    returns.row(it) = innov.row(it) + arma::dot(coeffr.rows(ncoeff-it, ncoeff-1), returns.rows(0, it-1));
  }  // end for
  
  // Perform loop over the remaining rows
  for (arma::uword it = ncoeff; it < nrows; it++) {
    returns.row(it) = innov.row(it) + arma::dot(coeffr, returns.rows(it-ncoeff, it-1));
  }  // end for
  
  return returns;
  
}  // end sim_ar



////////////////////////////////////////////////////////////
//' Simulate a \emph{Dickey-Fuller} process using \emph{Rcpp}.
//' 
//' @param \code{init_price} The initial price. 
//' 
//' @param \code{eq_price} The equilibrium price. 
//' 
//' @param \code{theta} The strength of mean reversion.
//' 
//' @param \code{coeff} A single-column \emph{matrix} of \emph{autoregressive}
//'   coefficients.
//'
//' @param \code{innov} A single-column \emph{matrix} of innovations (random
//'   numbers).
//' 
//' @return A single-column \emph{matrix} of simulated prices, with the same
//'   number of rows as the argument \code{innov}.
//'
//' @details
//'   The function \code{sim_df()} simulates the following \emph{Dickey-Fuller}
//'   process:
//'   \deqn{
//'     r_i = \theta \, (\mu - p_{i-1}) + \varphi_1 r_{i-1} + \ldots + \varphi_n r_{i-n} + \xi_i
//'   }
//'   \deqn{
//'     p_i = p_{i-1} + r_i
//'   }
//'   Where \eqn{r_i} and \eqn{p_i} are the simulated returns and prices,
//'   \eqn{\theta} and \eqn{\mu} are the \emph{Ornstein-Uhlenbeck} parameters,
//'   \eqn{\varphi_i} are the \emph{autoregressive} coefficients, and
//'   \eqn{\xi_i} are the normal \emph{innovations}.
//'   The recursion starts with: \eqn{r_1 = \xi_1} and \eqn{p_1 = init\_price}.
//'
//'   The \emph{Dickey-Fuller} process is a combination of an
//'   \emph{Ornstein-Uhlenbeck} process and an \emph{autoregressive} process.
//'   The order \eqn{n} of the \emph{autoregressive} process \eqn{AR(n)}, is
//'   equal to the number of rows of the \emph{autoregressive} coefficients
//'   \code{coeff}.
//'
//'   The function \code{sim_df()} simulates the \emph{Dickey-Fuller}
//'   process using fast \emph{Rcpp} \code{C++} code.
//'
//'   The function \code{sim_df()} returns a single-column \emph{matrix}
//'   representing the \emph{time series} of prices.
//'
//' @examples
//' \dontrun{
//' # Define the Ornstein-Uhlenbeck model parameters
//' init_price <- 1.0
//' eq_price <- 2.0
//' thetav <- 0.01
//' # Define AR coefficients
//' coeff <- matrix(c(0.1, 0.3, 0.5))
//' # Calculate matrix of standard normal innovations
//' innov <- matrix(rnorm(1e3, sd=0.01))
//' # Simulate Dickey-Fuller process using Rcpp
//' prices <- HighFreq::sim_df(init_price=init_price, eq_price=eq_price, theta=thetav, coeff=coeff, innov=innov)
//' plot(prices, t="l", main="Simulated Dickey-Fuller Prices")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat sim_df(double init_price, 
                 double eq_price, 
                 double theta, 
                 arma::mat& coeff, 
                 arma::mat& innov) {
  
  arma::uword nrows = innov.n_rows;
  arma::uword ncoeff = coeff.n_rows;
  arma::mat coeffr = arma::reverse(coeff);
  arma::mat prices = arma::zeros(nrows, 1);
  arma::mat returns = arma::zeros(nrows, 1);

  // Warmup period
  returns.row(0) = innov.row(0);
  prices.row(0) = init_price;
  returns.row(1) = theta*(eq_price - prices.row(0)) + coeffr.row(ncoeff-1) * returns.row(0) + innov.row(1);
  prices.row(1) = prices.row(0) + returns.row(1);
  for (arma::uword it = 2; it < ncoeff; it++) {
    returns.row(it) = theta*(eq_price - prices.row(it-1)) + arma::dot(coeffr.rows(ncoeff-it, ncoeff-1), returns.rows(0, it-1)) + innov.row(it);
    prices.row(it) = prices.row(it-1) + returns.row(it);
  }  // end for
  
  // Perform loop over the remaining rows
  for (arma::uword it = ncoeff; it < nrows; it++) {
    returns.row(it) = theta*(eq_price - prices.row(it-1)) + arma::dot(coeffr, returns.rows(it-ncoeff, it-1)) + innov.row(it);
    prices.row(it) = prices.row(it-1) + returns.row(it);
  }  // end for
  
  return prices;
  
}  // end sim_df



////////////////////////////////////////////////////////////
//' Calculate the log-likelihood of a time series of returns assuming a
//' \emph{GARCH(1,1)} process.
//' 
//' @param \code{omega} Parameter proportional to the long-term average level
//'   of variance.
//' 
//' @param \code{alpha} The weight associated with recent realized variance
//'   updates.
//' 
//' @param \code{beta} The weight associated with the past variance estimates.
//' 
//' @param \code{returns} A single-column \emph{matrix} of returns.
//' 
//' @param \code{minval} The floor value applied to the variance, to avoid zero
//'   values. (The default is \code{minval = 0.000001}.)
//' 
//' @return The log-likelihood value.
//'
//' @details
//'   The function \code{lik_garch()} calculates the log-likelihood of a time
//'   series of returns assuming a \emph{GARCH(1,1)} process.
//'   
//'   It first estimates the rolling variance of the \code{returns} argument
//'   using function \code{sim_garch()}:
//'   \deqn{
//'     \sigma^2_i = \omega + \alpha r^2_i + \beta \sigma_{i-1}^2
//'   }
//'   Where \eqn{r_i} is the time series of returns, and \eqn{\sigma^2_i} is the
//'   estimated rolling variance.
//'   And \eqn{\omega}, \eqn{\alpha}, and \eqn{\beta} are the \emph{GARCH}
//'   parameters.
//'   It applies the floor value \code{minval} to the variance, to avoid zero
//'   values.  So the minimum value of the variance is equal to \code{minval}.
//'
//'   The function \code{lik_garch()} calculates the log-likelihood assuming a
//'   normal distribution of returns conditional on the variance
//'   \eqn{\sigma^2_{i-1}} in the previous period, as follows:
//'   \deqn{
//'     likelihood = - \sum_{i=1}^n (\frac{r^2_i}{\sigma^2_{i-1}} + \log(\sigma^2_{i-1}))
//'   }
//'
//' @examples
//' \dontrun{
//' # Define the GARCH model parameters
//' alpha <- 0.79
//' betav <- 0.2
//' om_ega <- 1e-4*(1-alpha-betav)
//' # Calculate historical VTI returns
//' retp <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the log-likelihood of VTI returns assuming GARCH(1,1)
//' HighFreq::lik_garch(omega=om_ega, alpha=alpha,  beta=betav, returns=retp)
//' }
//' 
//' @export
// [[Rcpp::export]]
double lik_garch(double omega, 
                 double alpha, 
                 double beta,
                 arma::mat& returns, 
                 double minval = 0.000001) {
  
  // Calculate the rolling variance of returns using function sim_garch()
  arma::mat garch_data = sim_garch(omega, alpha,  beta, returns, false);
  // Select the second column containing the variance of returns
  arma::mat varm = garch_data.col(1);
  // Apply floor to variance
  varm.transform([&minval](double x) {return max(x, minval);});
  // Lag the variance
  varm = lagit(varm, 1, false);
  // Calculate the log-likelihood
  double likelihood = -arma::conv_to<double>::from(arma::sum(pow(returns, 2)/varm + log(varm)));
  
  return likelihood;
  
}  // end lik_garch



////////////////////////////////////////////////////////////
// Functions for backtests
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
//' Simulate a portfolio optimization strategy using online (recursive) updating
//' of the covariance matrix.
//' 
//' @param \code{rets} A \emph{time series} or \emph{matrix} of asset returns.
//' 
//' @param \code{dimax} An \emph{integer} equal to the number of \emph{eigen
//'   values} used for calculating the reduced inverse of the covariance
//'   matrix (the default is \code{dimax = 0} - standard matrix inverse using
//'   all the \emph{eigen values}).
//'   
//' @param \code{lambda} A decay factor which multiplies the past asset returns.
//'   
//' @param \code{lambdacov} A decay factor which multiplies the past covariance.
//'   
//' @param \code{lambdaw} A decay factor which multiplies the past portfolio
//'   weights.
//'   
//' @return A \emph{matrix} of strategy returns and the portfolio weights, with
//'   the same number of rows as the argument \code{rets}.
//'   
//' @details
//'   The function \code{sim_portfoptim()} simulates a portfolio optimization
//'   strategy. The strategy calculates the maximum Sharpe portfolio weights
//'   \emph{in-sample} at every point in time, and applies them in the
//'   \emph{out-of-sample} time interval.  It updates the trailing covariance
//'   matrix recursively, instead of using past batches of data. The function
//'   \code{sim_portfoptim()} uses three different decay factors for averaging past
//'   values, to reduce the variance of its forecasts.
//'   
//'   The function \code{sim_portfoptim()} first scales the returns by their
//'   trailing volatilities:
//'   \deqn{
//'     r^s_t = \frac{r_t}{\sigma_{t-1}}
//'   }
//'   Returns scaled by their volatility are more stationary so they're easier
//'   to model.
//'   
//'   Then at every point in time, the function \code{sim_portfoptim()} calls
//'   the function \code{HighFreq::push_covar()} to update the trailing
//'   covariance matrix of the returns:
//'   \deqn{
//'     \bar{r}_t = \lambda_c \bar{r}_{t-1} + (1-\lambda_c) r^s_t
//'   }
//'   \deqn{
//'     \hat{r}_t = r^s_t - \bar{r}_t
//'   }
//'   \deqn{
//'     {cov}_t = \lambda_c {cov}_{t-1} + (1-\lambda_c) \hat{r}^T_t \hat{r}_t
//'   }
//'   Where \eqn{\lambda_c} is the decay factor which multiplies the past mean
//'   and covariance.
//'   
//'   It then calls the function \code{HighFreq::calc_inv()} to calculate the
//'   \emph{reduced inverse} of the covariance matrix using its eigen
//'   decomposition:
//'   \deqn{
//'     \strong{C}^{-1} = \strong{O}_{dimax} \, \Sigma^{-1}_{dimax} \, \strong{O}^T_{dimax}
//'   }
//'   See the function \code{HighFreq::calc_inv()} for details.
//'   
//'   It then calculates the \emph{in-sample} weights of the maximum Sharpe
//'   portfolio, by multiplying the inverse covariance matrix times the trailing
//'   means of the asset returns:
//'   \deqn{
//'     \bar{r}_t = \lambda \bar{r}_{t-1} + (1-\lambda) r^s_t
//'   }
//'   \deqn{
//'     \strong{w}_t = \strong{C}^{-1} \bar{r}_t
//'   }
//'   Note that the decay factor \eqn{\lambda} is different from the decay
//'   factor \eqn{\lambda_c} used for updating the trailing covariance
//'   matrix.
//'   
//'   It then scales the weights so their sum of squares is equal to one:
//'   \deqn{
//'     \strong{w}_t = \frac{\strong{w}_t}{\sqrt{\sum{\strong{w}^2_t}}}
//'   }
//'   
//'   It then calculates the trailing mean of the weights:
//'   \deqn{
//'     \bar{\strong{w}}_t = \lambda_w \bar{\strong{w}}_{t-1} + (1-\lambda_w) \strong{w}_t
//'   }
//'   Note that the decay factor \eqn{\lambda_w} is different from the decay
//'   factor \eqn{\lambda} used for updating the trailing means.
//'   
//'   It finally calculates the \emph{out-of-sample} portfolio returns by
//'   multiplying the trailing mean weights times the scaled asset returns:
//'   \deqn{
//'     r^p_t = \bar{\strong{w}}_{t-1} r^s_t
//'   }
//'   Applying weights to scaled returns means trading stock amounts with unit
//'   dollar volatility.  So if the weight is equal to \code{2} then we should
//'   purchase an amount of stock with dollar volatility equal to \code{2}
//'   dollars.  Trading stock amounts with unit dollar volatility improves
//'   portfolio diversification.
//'   
//'   The function \code{sim_portfoptim()} uses three different decay factors
//'   for averaging past values, to reduce the variance of its forecasts. The
//'   value of the decay factor \eqn{\lambda} must be in the range between
//'   \code{0} and \code{1}.
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, so the trailing values have a greater
//'   dependence on past data.  This is equivalent to a long look-back
//'   interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, so the trailing values have a weaker
//'   dependence on past data.  This is equivalent to a short look-back
//'   interval.
//' 
//'   The function \code{sim_portfoptim()} returns multiple columns of data,
//'   with the same number of rows as the input argument \code{rets}. The first
//'   column contains the strategy returns and the remaining columns contain the
//'   portfolio weights.
//'   
//' @examples
//' \dontrun{
//' # Load ETF returns
//' retp <- rutils::etfenv$returns[, c("VTI", "TLT", "DBC", "USO", "XLF", "XLK")]
//' retp <- na.omit(retp)
//' datev <- zoo::index(retp) # dates
//' # Simulate a portfolio optimization strategy
//' dimax <- 6
//' lambdaf <- 0.978 # Decay factor
//' lambdacov <- 0.995
//' lambdaw <- 0.9
//' pnls <- HighFreq::sim_portfoptim(retp, dimax, lambdaf, lambdacov, lambdaw)
//' colnames(pnls) <- c("pnls", "VTI", "TLT", "DBC", "USO", "XLF", "XLK")
//' pnls <- xts::xts(pnls, order.by=datev)
//' # Plot dygraph of strategy
//' wealthv <- cbind(retp$VTI, pnls$pnls*sd(retp$VTI)/sd(pnls$pnls))
//' colnames(wealthv) <- c("VTI", "Strategy")
//' endd <- rutils::calc_endpoints(wealthv, interval="weeks")
//' dygraphs::dygraph(cumsum(wealthv)[endd], main="Portfolio Optimization Strategy Returns") %>%
//'  dyOptions(colors=c("blue", "red"), strokeWidth=2) %>%
//'  dyLegend(width=300)
//' # Plot dygraph of weights
//' symbolv <- "VTI"
//' stockweights <- cbind(cumsum(get(symbolv, retp)), get(symbolv, pnls))
//' colnames(stockweights)[2] <- "Weight"
//' colnamev <- colnames(stockweights)
//' endd <- rutils::calc_endpoints(pnls, interval="weeks")
//' dygraphs::dygraph(stockweights[endd], main="Returns and Weight") %>%
//'   dyAxis("y", label=colnamev[1], independentTicks=TRUE) %>%
//'   dyAxis("y2", label=colnamev[2], independentTicks=TRUE) %>%
//'   dySeries(axis="y", label=colnamev[1], strokeWidth=2, col="blue") %>%
//'   dySeries(axis="y2", label=colnamev[2], strokeWidth=2, col="red")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat sim_portfoptim(const arma::mat& rets, // Asset returns
                         const arma::uword& dimax, // Number of eigen vectors for dimension reduction
                         const double& lambda, // Returns decay factor
                         const double& lambdacov, // Covariance decay factor
                         const double& lambdaw) { // Weight decay factor
  
  arma::uword nrows = rets.n_rows;
  arma::uword ncols = rets.n_cols;
  arma::rowvec retsc; // Row of scaled asset returns
  arma::rowvec retm(ncols, fill::zeros); // Trailing average of asset returns
  arma::rowvec varv; // Variance of asset returns
  arma::mat invmat(ncols, ncols, fill::ones); // Inverse covariance matrix
  arma::mat weightv(nrows, ncols, fill::zeros); // Portfolio weights
  arma::colvec stratret(nrows, fill::zeros); // Strategy returns
  double weightd; // Sum of squared weights
  double lambda1 = 1-lambda;
  double lambdaw1 = 1-lambdaw;
  
  // Initialize the covariance matrix with first row of data
  arma::rowvec meanv = rets.row(0);
  arma::mat covmat = arma::trans(meanv)*meanv;

  // Perform loop over the rows (time)
  for (arma::uword it = 1; it < nrows; it++) {
    // Scale the returns by the trailing volatility
    // retsc = rets.row(it);
    varv = arma::trans(covmat.diag());
    varv.replace(0, 1);
    retsc = rets.row(it)/arma::sqrt(varv);
    // Calculate the strategy returns - the products of the lagged weights times the asset returns
    stratret(it) = arma::dot(retsc, weightv.row(it-1));
    // Update the covariance matrix with new row of asset returns
    // if (scalit) {
    push_covar(retsc, covmat, meanv, lambdacov);
    // } else {
    //   push_covar(rets.row(it), covmat, meanv, lambdacov);
    // }  // end if
    // Calculate the trailing mean of asset returns
    retm = lambda*retm + lambda1*rets.row(it);
    // Calculate weights using reduced inverse
    invmat = calc_inv(covmat, dimax);
    // Calculate weights using recursive inverse
    // calc_invrec(covmat, invmat);
    // weightv.row(it) = retm*invmat;
    weightv.row(it) = lambdaw*weightv.row(it-1) + lambdaw1*retm*invmat;
    // Scale the weights so their sum of squares is equal to one
    weightd = std::sqrt(arma::sum(arma::square(weightv.row(it))));
    if (weightd == 0) weightd = 1;
    weightv.row(it) = weightv.row(it)/weightd;
  }  // end for
  
  // Return the strategy returns and the portfolio weights
  return arma::join_rows(stratret, weightv);
  
}  // end sim_portfoptim




////////////////////////////////////////////////////////////
//' Calculate the optimal portfolio weights using a variety of different
//' objective functions.
//' 
//' @param \code{returns} A \emph{time series} or a \emph{matrix} of returns
//'   data (the returns in excess of the risk-free rate).
//'   
//' @param \code{controlv} A \emph{list} of portfolio optimization model
//'   parameters (see Details).
//'
//'   
//' @return A column \emph{vector} of the same length as the number of columns
//'   of \code{returns}.
//'
//' @details
//'   The function \code{calc_weights()} calculates the optimal portfolio
//'   weights using a variety of different objective functions.
//' 
//'   The function \code{calc_weights()} accepts a list of portfolio
//'   optimization parameters through the argument \code{controlv}.
//'   The list of portfolio optimization parameters can be created using
//'   the function \code{param_portf()}.  Below is a description of the
//'   parameters.
//'
//'   If \code{method = "maxsharpe"} (the default) then \code{calc_weights()}
//'   calculates the weights of the maximum Sharpe portfolio, by multiplying the
//'   \emph{reduced inverse} of the \emph{covariance matrix}
//'   \eqn{\strong{C}^{-1}} times the mean column returns \eqn{\bar{r}}:
//'   \deqn{
//'     \strong{w} = \strong{C}^{-1} \bar{r}
//'   }
//'   
//'   If \code{method = "maxsharpemed"} then \code{calc_weights()} uses the
//'   medians instead of the means.
//'   
//'   If \code{method = "minvarlin"} then it calculates the weights of the
//'   minimum variance portfolio under linear constraint, by multiplying the
//'   \emph{reduced inverse} of the \emph{covariance matrix} times the unit
//'   vector:
//'   \deqn{
//'     \strong{w} = \strong{C}^{-1} \strong{1}
//'   }
//'   
//'   If \code{method = "minvarquad"} then it calculates the weights of the
//'   minimum variance portfolio under quadratic constraint (which is the
//'   highest order principal component).
//'
//'   If \code{method = "sharpem"} then it calculates the momentum weights equal
//'   to the Sharpe ratios (the \code{returns} divided by their standard
//'   deviations):
//'   \deqn{
//'     \strong{w} = \frac{\bar{r}}{\sigma}
//'   }
//'
//'   If \code{method = "kellym"} then it calculates the momentum weights equal
//'   to the Kelly ratios (the \code{returns} divided by their variance):
//'   \deqn{
//'     \strong{w} = \frac{\bar{r}}{\sigma^2}
//'   }
//'
//'   \code{calc_weights()} calls the function \code{calc_inv()} to calculate
//'   the \emph{reduced inverse} of the \emph{covariance matrix} of
//'   \code{returns}. It performs regularization by selecting only the largest
//'   eigenvalues equal in number to \code{dimax}.
//'   
//'   In addition, \code{calc_weights()} applies shrinkage to the columns of
//'   \code{returns}, by shrinking their means to their common mean value:
//'   \deqn{
//'     r^{\prime}_i = (1 - \alpha) \, \bar{r}_i + \alpha \, \mu
//'   }
//'   Where \eqn{\bar{r}_i} is the mean of column \eqn{i} and \eqn{\mu} is the
//'   average of all the column means.
//'   The shrinkage intensity \code{alpha} determines the amount of shrinkage
//'   that is applied, with \code{alpha = 0} representing no shrinkage (with the
//'   column means \eqn{\bar{r}_i} unchanged), and \code{alpha = 1} representing
//'   complete shrinkage (with the column means all equal to the single mean of
//'   all the columns: \eqn{\bar{r}_i = \mu}).
//'
//'   After the weights are calculated, they are scaled, depending on several
//'   arguments.
//'
//'   If \code{rankw = TRUE} then the weights are converted into their ranks.
//'   The default is \code{rankw = FALSE}.
//'
//'   If \code{centerw = TRUE} then the weights are centered so that their sum
//'   is equal to \code{0}.  The default is \code{centerw = FALSE}.
//'
//'   If \code{scalew = "voltarget"} (the default) then the weights are
//'   scaled (multiplied by a factor) so that the weighted portfolio has an
//'   in-sample volatility equal to \code{voltarget}.
//'   
//'   If \code{scalew = "voleqw"} then the weights are scaled so that the
//'   weighted portfolio has the same volatility as the equal weight portfolio.
//'   
//'   If \code{scalew = "sumone"} then the weights are scaled so that their
//'   sum is equal to \code{1}.
//'   If \code{scalew = "sumsq"} then the weights are scaled so that their
//'   sum of squares is equal to \code{1}.
//'   If \code{scalew = "none"} then the weights are not scaled.
//' 
//'   The function \code{calc_weights()} is written in \code{C++}
//'   \code{RcppArmadillo} code.
//'   
//' @examples
//' \dontrun{
//' # Calculate covariance matrix and eigen decomposition of ETF returns
//' retp <- na.omit(rutils::etfenv$returns[, 1:16])
//' ncols <- NCOL(retp)
//' eigend <- eigen(cov(retp))
//' # Calculate reduced inverse of covariance matrix
//' dimax <- 3
//' eigenvec <- eigend$vectors[, 1:dimax]
//' eigenval <- eigend$values[1:dimax]
//' invmat <- eigenvec %*% (t(eigenvec) / eigenval)
//' # Define shrinkage intensity and apply shrinkage to the mean returns
//' alpha <- 0.5
//' colmeans <- colMeans(retp)
//' colmeans <- ((1-alpha)*colmeans + alpha*mean(colmeans))
//' # Calculate weights using R
//' weightr <- drop(invmat %*% colmeans)
//' # Apply weights scaling
//' weightr <- weightr*sd(rowMeans(retp))/sd(retp %*% weightr)
//' weightr <- 0.01*weightr/sd(retp %*% weightr)
//' weightr <- weightr/sqrt(sum(weightr^2))
//' # Create a list of portfolio optimization parameters
//' controlv <- HighFreq::param_portf(method="maxsharpe", dimax=dimax, alpha=alpha, scalew="sumsq")
//' # Calculate weights using RcppArmadillo
//' weightcpp <- drop(HighFreq::calc_weights(retp, controlv=controlv))
//' all.equal(weightcpp, weightr)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec calc_weights(const arma::mat& returns, // Asset returns
                       Rcpp::List controlv) { // List of portfolio optimization parameters
  
  // Unpack the control list of portfolio optimization parameters
  // Type of portfolio optimization model
  std::string method = Rcpp::as<std::string>(controlv["method"]);
  // Threshold level for discarding small singular values
  double singmin = Rcpp::as<double>(controlv["singmin"]);
  // Dimension reduction
  arma::uword dimax = Rcpp::as<int>(controlv["dimax"]);
  // Confidence level for calculating the quantiles of returns
  double confl = Rcpp::as<double>(controlv["confl"]);
  // Shrinkage intensity of returns
  double alpha = Rcpp::as<double>(controlv["alpha"]);
  // Should the weights be ranked?
  bool rankw = Rcpp::as<int>(controlv["rankw"]);
  // Should the weights be centered?
  bool centerw = Rcpp::as<int>(controlv["centerw"]);
  // Method for scaling the weights
  std::string scalew = Rcpp::as<std::string>(controlv["scalew"]);
  // Volatility target for scaling the weights
  double voltarget = Rcpp::as<double>(controlv["voltarget"]);

  // Initialize the variables
  arma::uword ncols = returns.n_cols;
  arma::vec weightv(ncols, fill::zeros);
  // If no regularization then set dimax to ncols
  if (dimax == 0)  dimax = ncols;
  // Calculate the covariance matrix
  arma::mat covmat = calc_covar(returns);
  
  // Apply different calculation methods for the weights
  switch(calc_method(method)) {
  case methodenum::maxsharpe: {
    // Mean returns of columns
    arma::vec colmeans = arma::trans(arma::mean(returns, 0));
    // Shrink colmeans to the mean of returns
    colmeans = ((1-alpha)*colmeans + alpha*arma::mean(colmeans));
    // Calculate weights using reduced inverse
    weightv = calc_inv(covmat, dimax, singmin)*colmeans;
    break;
  }  // end maxsharpe
  case methodenum::maxsharpemed: {
    // Median returns of columns
    arma::vec colmeans = arma::trans(arma::median(returns, 0));
    // Shrink colmeans to the median of returns
    colmeans = ((1-alpha)*colmeans + alpha*arma::median(colmeans));
    // Calculate weights using reduced inverse
    weightv = calc_inv(covmat, dimax, singmin)*colmeans;
    break;
  }  // end maxsharpemed
  case methodenum::minvarlin: {
    // Minimum variance weights under linear constraint
    // Multiply reduced inverse times unit vector
    weightv = calc_inv(covmat, dimax, singmin)*arma::ones(ncols);
    break;
  }  // end minvarlin
  case methodenum::minvarquad: {
    // Minimum variance weights under quadratic constraint
    // Calculate highest order principal component
    arma::vec eigenval;
    arma::mat eigenvec;
    arma::eig_sym(eigenval, eigenvec, covmat);
    weightv = eigenvec.col(ncols-1);
    break;
  }  // end minvarquad
  case methodenum::sharpem: {
    // Momentum weights equal to Sharpe ratios
    // Mean returns of columns
    arma::vec colmeans = arma::trans(arma::mean(returns, 0));
    // Standard deviation of columns
    arma::vec colsd = arma::sqrt(covmat.diag());
    colsd.replace(0, 1);
    // Momentum weights equal to Sharpe ratios
    weightv = colmeans/colsd;
    break;
  }  // end sharpem
  case methodenum::kellym: {
    // Momentum weights equal to Kelly ratios
    // Mean returns of columns
    arma::vec colmeans = arma::trans(arma::mean(returns, 0));
    // Variance of columns
    arma::vec colvar = covmat.diag();
    colvar.replace(0, 1);
    // Momentum weights equal to Kelly ratios
    weightv = colmeans/colvar;
    break;
  }  // end kellym
  case methodenum::robustm: {
    // Momentum weights equal to robust Sharpe ratios
    // Median returns of columns
    arma::vec colmeans = arma::trans(arma::median(returns, 0));
    // Standard deviation of columns
    arma::vec colsd = arma::sqrt(covmat.diag());
    colsd.replace(0, 1);
    // Momentum weights equal to robust Sharpe ratios
    colmeans = colmeans/colsd;
    break;
  }  // end robustm
  case methodenum::quantile: {
    // Momentum weights equal to sum of quantiles for columns
    arma::vec levels = {confl, 1-confl};
    weightv = arma::conv_to<vec>::from(arma::sum(arma::quantile(returns, levels, 0), 0));
    break;
  }  // end quantile
  default : {
    cout << "Warning: Invalid method parameter: " << method << endl;
    return arma::ones(ncols);
  }  // end default
  }  // end switch
  
  if (rankw == TRUE) {
    // Convert the weights to their ranks
    weightv = arma::conv_to<vec>::from(calc_ranks_stl(weightv));
  }  // end if
  
  if (centerw == TRUE) {
    // Center the weights so their sum is equal to zero
    weightv = (weightv - arma::mean(weightv));
  }  // end if
  
  // Apply different scaling methods for the weights
  switch(calc_method(scalew)) {
  case methodenum::voltarget: {
    // Scale the weights so the portfolio has the volatility equal to voltarget
    weightv = weightv*voltarget/arma::stddev(returns*weightv);
    break;
  }  // end voltarget
  case methodenum::voleqw: {
    // Scale the weights to the volatility of the equal weight portfolio
    weightv = weightv*arma::stddev(arma::mean(returns, 1))/arma::stddev(returns*weightv);
    break;
  }  // end voleqw
  case methodenum::sumone: {
    // Scale the weights so their sum of squares is equal to one
    weightv = weightv/arma::sum(weightv*arma::ones(ncols));
    break;
  }  // end sumone
  case methodenum::sumsq: {
    // Scale the weights so their sum of squares is equal to one
    weightv = weightv/std::sqrt(arma::sum(arma::square(weightv)));
    break;
  }  // end sumsq
  default : {
    // No scaling
    break;
  }  // end default
  }  // end switch
  
  return weightv;
  
}  // end calc_weights



////////////////////////////////////////////////////////////
//' Simulate (backtest) a rolling portfolio optimization strategy, using
//' \code{RcppArmadillo}.
//' 
//' @param \code{retp} A \emph{time series} or a \emph{matrix} of asset
//'   returns data.
//'   
//' @param \code{retx} A \emph{time series} or a \emph{matrix} of excess
//'   returns data (the returns in excess of the risk-free rate).
//'   
//' @param \code{controlv} A \emph{list} of portfolio optimization model
//'   parameters (see Details).
//'
//' @param \code{startp} An \emph{integer vector} of start points.
//' 
//' @param \code{endd} An \emph{integer vector} of end points.
//' 
//' @param \code{lambda} A decay factor which multiplies the past portfolio
//'   weights.  (The default is \code{lambda = 0} - no memory.)
//'   
//' @param \code{coeff} A \emph{numeric} multiplier of the weights.  (The
//'   default is \code{1})
//'   
//' @param \code{bidask} A \emph{numeric} bid-ask spread (the default is
//'   \code{0})
//'   
//'   
//' @return A column \emph{vector} of strategy returns, with the same length as
//'   the number of rows of \code{retp}.
//'
//' @details
//'   The function \code{back_test()} performs a backtest simulation of a
//'   rolling portfolio optimization strategy over a \emph{vector} of the end
//'   points \code{endd}.
//'   
//'   It performs a loop over the end points \code{endd}, and subsets the
//'   \emph{matrix} of the excess asset returns \code{retx} along its rows,
//'   between the corresponding \emph{start point} and the \emph{end point}. 
//'   
//'   The function \code{back_test()} passes the subset matrix of excess returns
//'   into the function \code{calc_weights()}, which calculates the optimal
//'   portfolio weights at each \emph{end point}.
//'   It also passes to \code{calc_weights()} the argument \code{controlv},
//'   which is the list of portfolio optimization parameters.
//'   See the function \code{calc_weights()} for more details.
//'   The list of portfolio optimization parameters can be created using the
//'   function \code{param_portf()}.
//'   
//'   The function \code{back_test()} then recursively averages the weights
//'   \eqn{w_i} at the \emph{end point = i} with the weights \eqn{w_{i-1}} from
//'   the previous \emph{end point = (i-1)}, using the decay factor \code{lambda
//'   = \eqn{\lambda}}:
//'   \deqn{
//'     w_i = (1-\lambda) w_i + \lambda w_{i-1}
//'   }
//'   The purpose of averaging the weights is to reduce their variance, and
//'   improve their out-of-sample performance.  It is equivalent to extending
//'   the portfolio holding period beyond the time interval between neighboring
//'   \emph{end points}.
//'   
//'   The function \code{back_test()} then calculates the out-of-sample strategy
//'   returns by multiplying the average weights times the future asset returns.
//'   
//'   The function \code{back_test()} multiplies the out-of-sample strategy
//'   returns by the coefficient \code{coeff} (with default equal to \code{1}),
//'   which allows simulating either a trending strategy (if \code{coeff = 1}),
//'   or a reverting strategy (if \code{coeff = -1}).
//'   
//'   The function \code{back_test()} calculates the transaction costs by
//'   multiplying the bid-ask spread \code{bidask} times the absolute
//'   difference between the current weights minus the weights from the previous
//'   period. Then it subtracts the transaction costs from the out-of-sample
//'   strategy returns.
//'   
//'   The function \code{back_test()} returns a \emph{time series} (column
//'   \emph{vector}) of strategy returns, of the same length as the number of
//'   rows of \code{retp}.
//'
//' @examples
//' \dontrun{
//' # Calculate the ETF daily excess returns
//' retp <- na.omit(rutils::etfenv$returns[, 1:16])
//' # riskf is the daily risk-free rate
//' riskf <- 0.03/260
//' retx <- retp - riskf
//' # Define monthly end points without initial warmup period
//' endd <- rutils::calc_endpoints(retp, interval="months")
//' endd <- endd[endd > 0]
//' nrows <- NROW(endd)
//' # Define 12-month look-back interval and start points over sliding window
//' lookb <- 12
//' startp <- c(rep_len(1, lookb-1), endd[1:(nrows-lookb+1)])
//' # Define return shrinkage and dimension reduction
//' alpha <- 0.5
//' dimax <- 3
//' # Create a list of portfolio optimization parameters
//' controlv <- HighFreq::param_portf(method="maxsharpe", dimax=dimax, alpha=alpha, scalew="sumsq")
//' # Simulate a monthly rolling portfolio optimization strategy
//' pnls <- HighFreq::back_test(retx, retp, controlv=controlv, startp=(startp-1), endd=(endd-1))
//' pnls <- xts::xts(pnls, index(retp))
//' colnames(pnls) <- "strategy"
//' # Plot dygraph of strategy
//' dygraphs::dygraph(cumsum(pnls), 
//'   main="Cumulative Returns of Max Sharpe Portfolio Strategy")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat back_test(const arma::mat& retx, // Asset excess returns
                    const arma::mat& retp, // Asset returns
                    Rcpp::List controlv, // List of portfolio optimization model parameters
                    arma::uvec startp, // Start points
                    arma::uvec endd, // End points
                    double lambda = 0.0, // Decay factor for averaging the portfolio weights
                    double coeff = 1.0, // Multiplier of strategy returns
                    double bidask = 0.0) { // The bid-ask spread
  
  double lambda1 = 1-lambda;
  arma::uword nweights = retp.n_cols;
  arma::vec weightv(nweights, fill::zeros);
  arma::vec weights_past = arma::ones(nweights)/std::sqrt(nweights);
  arma::mat pnls = arma::zeros(retp.n_rows, 1);

  // Perform loop over the end points
  for (arma::uword it = 1; it < endd.size(); it++) {
    // cout << "it: " << it << endl;
    // Calculate the portfolio weights
    weightv = coeff*calc_weights(retx.rows(startp(it-1), endd(it-1)), controlv);
    // Calculate the weights as the weighted sum with past weights
    weightv = lambda1*weightv + lambda*weights_past;
    // Calculate out-of-sample returns
    pnls.rows(endd(it-1)+1, endd(it)) = retp.rows(endd(it-1)+1, endd(it))*weightv;
    // Add transaction costs
    pnls.row(endd(it-1)+1) -= bidask*sum(abs(weightv - weights_past))/2;
    // Copy the weights
    weights_past = weightv;
  }  // end for
  
  // Return the strategy pnls
  return pnls;
  
}  // end back_test

