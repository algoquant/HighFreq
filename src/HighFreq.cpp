// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
using namespace std;
using namespace Rcpp;
using namespace arma;

////////////////////////////////////////////////////////////
// Rcpp and RcppArmadillo functions for package HighFreq
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
// Functions for matrix algebra
////////////////////////////////////////////////////////////


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
//' returns <- rnorm(1e6)
//' # Compare lag_vec() with rutils::lagit()
//' all.equal(drop(HighFreq::lag_vec(returns)), 
//'   rutils::lagit(returns))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::lag_vec(returns),
//'   Rcode=rutils::lagit(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec lag_vec(const arma::vec& tseries, 
                  arma::sword lagg = 1, 
                  bool pad_zeros = true) {
  
  arma::uword nrows = (tseries.n_elem-1);
  
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
//' returns <- matrix(rnorm(5e6), nc=5)
//' # Compare lagit() with rutils::lagit()
//' all.equal(HighFreq::lagit(returns), rutils::lagit(returns))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::lagit(returns),
//'   Rcode=rutils::lagit(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat lagit(const arma::mat& tseries, 
                 arma::sword lagg = 1, 
                 bool pad_zeros = true) {
  
  arma::uword nrows = (tseries.n_rows-1);
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
//'   because it requires the copying of data.
//'   
//'   The function \code{diff_vec()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, which makes it several times faster than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random returns
//' returns <- rnorm(1e6)
//' # Compare diff_vec() with rutils::diffit()
//' all.equal(drop(HighFreq::diff_vec(returns, lagg=3, pad=TRUE)),
//'   rutils::diffit(returns, lagg=3))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::diff_vec(returns, lagg=3, pad=TRUE),
//'   Rcode=rutils::diffit(returns, lagg=3),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec diff_vec(const arma::vec& tseries, arma::uword lagg = 1, bool pad_zeros = true) {
  
  arma::uword length = (tseries.n_elem-1);
  
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
//'   \emph{matrix} be padded (extended) with zeros, in order to return a
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
//'   should be padded (extended) with the rows of the initial (warmup) period
//'   at the front, in order to return a \emph{matrix} with the same number of
//'   rows as the input \code{tseries}.  The default is \code{pad_zeros = TRUE}.
//'   The padding operation can be time-consuming, because it requires the
//'   copying of data.
//'   
//'   The function \code{diffit()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, which makes it much faster than \code{R} code.
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
arma::mat diffit(const arma::mat& tseries, 
                  arma::sword lagg = 1, 
                  bool pad_zeros = true) {
  
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
//' Calculate a vector of end points that divides a vector into equal intervals.
//'
//' @param \code{length} An \emph{integer} equal to the length of the vector to
//'   be divided into equal intervals.
//'   
//' @param \code{step} The number of elements in each interval between
//'   neighboring end points.
//' 
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points.
//'
//' @return A vector of equally spaced \emph{integers} representing the end
//'   points.
//'
//' @details
//'   The end points are a vector of integers which divide a vector of length
//'   equal to \code{length} into equally spaced intervals. If a whole number of
//'   intervals doesn't fit over the vector, then \code{calc_endpoints()} adds a
//'   stub interval at the end.
//'
//'   The first end point is equal to the argument \code{step}, unless the
//'   argument \code{stub} is provided, and then it becomes the first end point.
//'
//'   For example, consider the end points for a vector of length \code{20}
//'   divided into intervals of length \code{step=5}: \code{0, 5, 10, 15, 20}.
//'   In order for all the differences between neighboring end points to be
//'   equal to \code{5}, the first end point is set equal to \code{0}. But
//'   \code{0} doesn't correspond to any vector element, so
//'   \code{calc_endpoints()} doesn't include it and it only retains the
//'   non-zero end points equal to: \code{5, 10, 15, 20}. 
//'
//'   Since indexing in \code{C++} code starts at \code{0}, then
//'   \code{calc_endpoints()} shifts the end points by \code{-1} and returns the
//'   vector equal to \code{4, 9, 14, 19}.
//'
//'   If \code{stub = 1} then the first end point is equal to \code{1} and the
//'   end points are equal to: \code{1, 6, 11, 16, 20}.
//'   The extra stub interval at the end is equal to \code{4 = 20 - 16}.
//'   And \code{calc_endpoints()} returns \code{0, 5, 10, 15, 19}. The first
//'   value is equal to \code{0} which is the index of the first element in
//'   \code{C++} code.
//'
//'   If \code{stub = 2} then the first end point is equal to \code{2}, with an
//'   extra stub interval at the end, and the end points are equal to: \code{2,
//'   7, 12, 17, 20}.
//'   And \code{calc_endpoints()} returns \code{1, 6, 11, 16, 19}.
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
//'   This works in \code{R} code because the vector element corresponding to
//'   index \code{0} is empty.  For example, the \code{R} code: \code{(4:1)[c(0,
//'   1)]} produces \code{4}.  So in \code{R} we can select vector elements
//'   using the end points starting at zero.
//'   
//'   In \code{C++} the end points must be shifted by \code{-1} compared to
//'   \code{R} code, because indexing starts at \code{0}: \code{-1, 4, 9, 14,
//'   19}.  But there is no vector element corresponding to index \code{-1}. So
//'   in \code{C++} we cannot select vector elements using the end points
//'   starting at \code{-1}. The solution is to drop the first placeholder end
//'   point.
//'   
//' @examples
//' # Calculate end points without a stub interval
//' HighFreq::calc_endpoints(length=20, step=5)
//' # Calculate end points with a final stub interval
//' HighFreq::calc_endpoints(length=23, step=5)
//' # Calculate end points with initial and final stub intervals
//' HighFreq::calc_endpoints(length=20, step=5, stub=2)
//' # Calculate end points with initial and final stub intervals
//' HighFreq::calc_endpoints(length=20, step=5, stub=24)
//'
//' @export
// [[Rcpp::export]]
arma::uvec calc_endpoints(arma::uword length, arma::uword step = 1, arma::uword stub = 0) {
  
  arma::uword extra = length % step;
  arma::uvec endp;
  
  if ((stub == 0) & (extra == 0)) {
    // No stub interval
    endp = arma::regspace<uvec>(step, step, length);
  } else if ((stub == 0) & (extra > 0)) {
    // Add stub interval at end
    endp = arma::regspace<uvec>(step, step, length + step);
    endp.back() = length;
  } else if ((stub > 0) & (extra == 0)) {
    // Add initial stub interval equal to stub
    endp = arma::regspace<uvec>(stub, step, length + step);
    endp.back() = length;
  } else if ((stub > 0) & (extra > 0) & (stub == extra)) {
    // Add initial stub interval equal to stub without stub at end
    endp = arma::regspace<uvec>(stub, step, length);
  } else {
    // Add initial stub interval equal to stub and with extra stub at end
    endp = arma::regspace<uvec>(stub, step, length + step);
    endp.back() = length;
  }  // end if
  
  // Subtract 1 from endp because C++ indexing starts at 0
  endp = endp - 1;
  return endp;
  
}  // end calc_endpoints




////////////////////////////////////////////////////////////
//' Calculate a vector of start points by lagging (shifting) a vector of end
//' points.
//'
//' @param \code{endp} An \emph{integer} vector of end points.
//'   
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   lag (shift) applied to the end points.
//'   
//' @return An \emph{integer} vector with the same number of elements as the
//'   vector \code{endp}.
//'
//' @details
//'   The start points are equal to the values of the vector \code{endp} lagged
//'   (shifted) by an amount equal to \code{look_back}.  In addition, an extra
//'   value of \code{1} is added to them, to avoid data overlaps.  The lag
//'   operation requires appending a beginning warmup interval containing zeros,
//'   so that the vector of start points has the same length as the \code{endp}.
//'   
//'   For example, consider the end points for a vector of length \code{25}
//'   divided into equal intervals of length \code{5}: \code{4, 9, 14, 19, 24}.
//'   (In \code{C++} the vector indexing starts at \code{0} not \code{1}, so
//'   it's shifted by \code{-1}.)
//'   Then the start points for \code{look_back = 2} are equal to: \code{0, 0, 
//'   5, 10, 15}.  The differences between the end points minus the
//'   corresponding start points are equal to \code{9}, except for the warmup
//'   interval.
//'   
//' @examples
//' # Calculate end points
//' endp <- HighFreq::calc_endpoints(25, 5)
//' # Calculate start points corresponding to the end points
//' startp <- HighFreq::calc_startpoints(endp, 2)
//'
//' @export
// [[Rcpp::export]]
arma::uvec calc_startpoints(arma::uvec endp, arma::uword look_back) {
  
  arma::uword numpts = endp.n_elem;
  arma::uvec startp = arma::join_cols(arma::zeros<uvec>(look_back), 
                                      endp.subvec(0, numpts - look_back - 1) + 1);
  
  return startp;
  
}  // end calc_startpoints



////////////////////////////////////////////////////////////
//' Multiply in place (without copying) the columns or rows of a \emph{matrix}
//' times a \emph{vector}, element-wise.
//' 
//' @param \code{vector} A \emph{vector}.
//' 
//' @param \code{matrix} A \emph{matrix}.
//' 
//' @param \code{by_col} A \emph{Boolean} argument: if \code{TRUE} then multiply
//'   the columns, otherwise multiply the rows (the default is
//'   \code{by_col = TRUE}.)
//' 
//' @return A single \emph{integer} value, equal to either the number of
//'   \emph{matrix} columns or the number of rows.
//' 
//' @details
//'   The function \code{mult_vec_mat()} multiplies the columns or rows of a
//'   \emph{matrix} times a \emph{vector}, element-wise.
//'
//'   If the number of \emph{vector} elements is equal to the number of matrix
//'   columns, then it multiplies the columns by the \emph{vector}, and returns
//'   the number of columns. If the number of \emph{vector} elements is equal to
//'   the number of rows, then it multiplies the rows, and returns the number of
//'   rows.
//'
//'   If the \emph{matrix} is square and if \code{by_col} is \code{TRUE} then it
//'   multiplies the columns, otherwise it multiplies the rows.
//'   
//'   It accepts \emph{pointers} to the \emph{matrix} and \emph{vector}, and
//'   replaces the old \emph{matrix} values with the new values.
//'   It performs the calculation in place, without copying the \emph{matrix} in
//'   memory (which greatly increases the computation speed).
//'   It performs an implicit loop over the \emph{matrix} rows and columns using
//'   the \emph{Armadillo} operators \code{each_row()} and \code{each_col()},
//'   instead of performing explicit \code{for()} loops (both methods are
//'   equally fast).
//'
//'   The function \code{mult_vec_mat()} uses \code{RcppArmadillo} \code{C++}
//'   code, so when multiplying large \emph{matrix} columns it's several times
//'   faster than vectorized \code{R} code, and it's even much faster compared
//'   to \code{R} when multiplying the \emph{matrix} rows.
//'   
//' @examples
//' \dontrun{
//' # Multiply matrix columns using R
//' matrixv <- matrix(round(runif(25e4), 2), nc=5e2)
//' vectorv <- round(runif(5e2), 2)
//' prod_uct <- vectorv*matrixv
//' # Multiply the matrix in place
//' HighFreq::mult_vec_mat(vectorv, matrixv)
//' all.equal(prod_uct, matrixv)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_vec_mat(vectorv, matrixv),
//'     Rcode=vectorv*matrixv,
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' 
//' # Multiply matrix rows using R
//' matrixv <- matrix(round(runif(25e4), 2), nc=5e2)
//' vectorv <- round(runif(5e2), 2)
//' prod_uct <- t(vectorv*t(matrixv))
//' # Multiply the matrix in place
//' HighFreq::mult_vec_mat(vectorv, matrixv, by_col=FALSE)
//' all.equal(prod_uct, matrixv)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_vec_mat(vectorv, matrixv, by_col=FALSE),
//'     Rcode=t(vectorv*t(matrixv)),
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::uword mult_vec_mat(arma::vec& vector,
                         arma::mat& matrix,
                         bool by_col = true) {
  
  arma::uword numelem = vector.n_elem;
  arma::uword nrows = matrix.n_rows;
  arma::uword ncols = matrix.n_cols;
  
  if ((ncols == nrows) && (numelem == nrows)) {
    if (by_col) {
      // Multiply each column of matrix by vector
      matrix.each_col() %= vector;
      return nrows;
    } else {
      // Multiply each row of matrix by vector
      matrix.each_row() %= conv_to<rowvec>::from(vector);
      return ncols;
    }
  } else if (numelem == nrows) {
    // Multiply each column of matrix by vector
    matrix.each_col() %= vector;
    return nrows;
  } else if (numelem == ncols) {
    // Multiply each row of matrix by vector
    matrix.each_row() %= conv_to<rowvec>::from(vector);
    return ncols;
  } else 
    stop("Error: Vector length is neither equal to the number of columns nor rows of the matrix!");
  // Return NA_INTEGER;
  
}  // end mult_vec_mat




////////////////////////////////////////////////////////////
//' Calculate the eigen decomposition of the covariance \emph{matrix} of returns
//' data using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A \emph{time series} or \emph{matrix} of returns
//'   data.
//'
//' @return A list with two elements: a \emph{vector} of eigenvalues 
//'   (named "values"), and a \emph{matrix} of eigenvectors (named
//'   "vectors").
//'
//' @details
//'   The function \code{calc_eigen()} first calculates the covariance
//'   \emph{matrix} of \code{tseries}, and then calculates the eigen
//'   decomposition of the covariance \emph{matrix}.
//'
//' @examples
//' \dontrun{
//' # Create matrix of random data
//' datav <- matrix(rnorm(5e6), nc=5)
//' # Calculate eigen decomposition
//' eigend <- HighFreq::calc_eigen(scale(datav, scale=FALSE))
//' # Calculate PCA
//' pcad <- prcomp(datav)
//' # Compare PCA with eigen decomposition
//' all.equal(pcad$sdev^2, drop(eigend$values))
//' all.equal(abs(unname(pcad$rotation)), abs(eigend$vectors))
//' # Compare the speed of Rcpp with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_eigen(datav),
//'   Rcode=prcomp(datav),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List calc_eigen(const arma::mat& tseries) {
  
  arma::mat eigenvec;
  arma::vec eigenval;
  arma::eig_sym(eigenval, eigenvec, arma::cov(tseries));
  
  // Reverse the order of elements from largest eigenvalue to smallest, similar to R
  return Rcpp::List::create(Named("values") = arma::flipud(eigenval),
                            Named("vectors") = arma::fliplr(eigenvec));
  
}  // end calc_eigen




////////////////////////////////////////////////////////////
//' Calculate the shrinkage inverse of a \emph{matrix} of data using Singular
//' Value Decomposition (\emph{SVD}).
//' 
//' @param \code{tseries} A \emph{time series} or \emph{matrix} of data.
//' 
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small singular values in order to regularize the inverse of the
//'   matrix \code{tseries} (the default is \code{0.01}).
//'   
//' @param \code{eigen_max} An \emph{integer} equal to the number of singular
//'   values used for calculating the shrinkage inverse of the matrix
//'   \code{tseries} (the default is \code{eigen_max = 0} - equivalent to
//'   \code{eigen_max} equal to the number of columns of \code{tseries}).
//'
//' @return A \emph{matrix} equal to the shrinkage inverse of the matrix
//'   \code{tseries}.
//'
//' @details
//'   The function \code{calc_inv()} calculates the shrinkage inverse of the
//'   matrix \code{tseries} using Singular Value Decomposition (\emph{SVD}).
//'   
//'   The function \code{calc_inv()} first performs Singular Value Decomposition
//'   (\emph{SVD}) of the matrix \code{tseries}.  
//'   The \emph{SVD} of a matrix \eqn{\strong{A}} is defined as the
//'   factorization:
//'   \deqn{
//'     \strong{A} = \strong{U}  \, \Sigma  \, \strong{V}^T
//'   }
//'   Where \eqn{\strong{U}} and \eqn{\strong{V}} are the left and right
//'   \emph{singular matrices}, and \eqn{\Sigma} is a diagonal matrix of
//'   \emph{singular values} \eqn{\Sigma = \{\sigma_i\}}.
//'   
//'   The inverse \eqn{\strong{A}^{-1}} of the matrix \eqn{\strong{A}} can be
//'   calculated from the \emph{SVD} matrices as:
//'   \deqn{
//'     \strong{A}^{-1} = \strong{V} \, \Sigma^{-1} \, \strong{U}^T
//'   }
//'   
//'   The \emph{regularized inverse} of the matrix \eqn{\strong{A}} is given by:
//'   \deqn{
//'     \strong{A}^{-1} = \strong{V}_n \, \Sigma_n^{-1} \, \strong{U}_n^T
//'   }
//'   Where \eqn{\strong{U}_n}, \eqn{\strong{V}_n} and \eqn{\Sigma_n} are the
//'   \emph{SVD} matrices with the rows and columns corresponding to zero
//'   \emph{singular values} removed.
//'   
//'   The function \code{calc_inv()} applies regularization by discarding the
//'   smallest singular values \eqn{\sigma_i} that are less than the threshold
//'   level \code{eigen_thresh} times the sum of all the singular values:
//'   \deqn{\sigma_i < eigen\_thresh \cdot (\sum{\sigma_i})}
//'   
//'   It then discards additional singular values so that only the largest
//'   \code{eigen_max} singular values remain.  
//'   It calculates the shrinkage inverse from the \emph{SVD} matrices using
//'   only the largest singular values up to \code{eigen_max}.  For example, if
//'   \code{eigen_max = 3} then it only uses the \code{3} largest singular
//'   values. This has the effect of dimension shrinkage.
//'   
//'   If the matrix \code{tseries} has a large number of small singular values,
//'   then the number of remaining singular values may be less than
//'   \code{eigen_max}.
//'   
//' @examples
//' \dontrun{
//' # Calculate ETF returns
//' returns <- na.omit(rutils::etfenv$returns)
//' # Calculate covariance matrix
//' covmat <- cov(returns)
//' # Calculate shrinkage inverse using RcppArmadillo
//' inverse <- HighFreq::calc_inv(covmat, eigen_max=3)
//' # Calculate shrinkage inverse from SVD in R
//' svdec <- svd(covmat)
//' eigen_max <- 1:3
//' inverser <-  svdec$v[, eigen_max] %*% (t(svdec$u[, eigen_max]) / svdec$d[eigen_max])
//' # Compare RcppArmadillo with R
//' all.equal(inverse, inverser)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& tseries,
                   double eigen_thresh = 0.01, 
                   arma::uword eigen_max = 0) {
  
  // Allocate SVD variables
  arma::vec svdval;  // Singular values
  arma::mat svdu, svdv;  // Singular matrices
  // Calculate the SVD
  arma::svd(svdu, svdval, svdv, tseries);
  // Calculate the number of non-small singular values
  arma::uword svdnum = arma::sum(svdval > eigen_thresh*arma::sum(svdval));
  
  if (eigen_max == 0) {
    // Set eigen_max
    eigen_max = svdnum - 1;
  } else {
    // Adjust eigen_max
    eigen_max = min(eigen_max - 1, svdnum - 1);
  }  // end if
  
  // Remove all small singular values
  svdval = svdval.subvec(0, eigen_max);
  svdu = svdu.cols(0, eigen_max);
  svdv = svdv.cols(0, eigen_max);
  
  // Calculate the shrinkage inverse from the SVD decomposition
  return svdv*arma::diagmat(1/svdval)*svdu.t();
  
}  // end calc_inv



////////////////////////////////////////////////////////////
//' Scale (standardize) the columns of a \emph{matrix} of data using
//' \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A \emph{time series} or \emph{matrix} of data.
//' 
//' @param \code{use_median} A \emph{Boolean} argument: if \code{TRUE} then the 
//'   centrality (central tendency) is calculated as the \emph{median} and the 
//'   dispersion is calculated as the \emph{median absolute deviation}
//'   (\emph{MAD}).
//'   If \code{use_median = FALSE} then the centrality is calculated as the
//'   \emph{mean} and the dispersion is calculated as the \emph{standard
//'   deviation} (the default is \code{FALSE})
//'
//' @return A \emph{matrix} with the same dimensions as the input
//'   argument \code{tseries}.
//'
//' @details
//'   The function \code{calc_scaled()} scales (standardizes) the columns of the
//'   \code{tseries} argument using \code{RcppArmadillo}.
//'
//'   If the argument \code{use_median} is \code{FALSE} (the default), then it
//'   performs the same calculation as the standard \code{R} function
//'   \code{scale()}, and it calculates the centrality (central tendency) as the
//'   \emph{mean} and the dispersion as the \emph{standard deviation}.
//'
//'   If the argument \code{use_median} is \code{TRUE}, then it calculates the
//'   centrality as the \emph{median} and the dispersion as the \emph{median
//'   absolute deviation} (\emph{MAD}).
//'
//'   If the number of rows of \code{tseries} is less than \code{3} then it
//'   returns \code{tseries} unscaled.
//'   
//'   The function \code{calc_scaled()} uses \code{RcppArmadillo} \code{C++}
//'   code and is about \emph{5} times faster than function \code{scale()}, for
//'   a \emph{matrix} with \emph{1,000} rows and \emph{20} columns.
//'   
//' @examples
//' \dontrun{
//' # Create a matrix of random data
//' returns <- matrix(rnorm(20000), nc=20)
//' scaled <- calc_scaled(tseries=returns, use_median=FALSE)
//' scaled2 <- scale(returns)
//' all.equal(scaled, scaled2, check.attributes=FALSE)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=calc_scaled(tseries=returns, use_median=FALSE),
//'   Rcode=scale(returns),
//'   times=100))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_scaled(const arma::mat& tseries, bool use_median=false) {
  
  arma::uword nrows = tseries.n_rows;
  arma::uword ncols = tseries.n_cols;
  arma::mat scaledmat(nrows, ncols, fill::zeros);
  arma::vec scaled(nrows, fill::zeros);
  double center;
  
  // Return zeros if not enough data
  if (nrows < 3) {
    return tseries;
  }  // end if
  
  // Perform a loop over the columns
  for (arma::uword it=0; it < ncols; it++) {
    if (use_median) {
      center = arma::median(tseries.col(it));
      scaled = (tseries.col(it) - center);
      scaled = scaled/arma::median(arma::abs(scaled));
      scaledmat.col(it) = scaled;
    } else {
      center = arma::mean(tseries.col(it));
      scaled = (tseries.col(it) - center);
      scaled = scaled/arma::stddev(scaled);
      scaledmat.col(it) = scaled;
    }  // end if
  }  // end for
  
  return scaledmat;
  
}  // end calc_scaled




////////////////////////////////////////////////////////////
//' Calculate the ranks of the elements of a single-column \emph{time series} or
//' a \emph{vector} using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A single-column \emph{time series} or a \emph{vector}.
//'
//' @return An \emph{integer vector} with the ranks of the elements of the
//'   \code{tseries}.
//'
//' @details
//'   The function \code{calc_ranks()} calculates the ranks of the elements of a
//'   single-column \emph{time series} or a \emph{vector}. It uses the
//'   \code{RcppArmadillo} function \code{arma::sort_index()}. The function
//'   \code{arma::sort_index()} calculates the permutation index to sort a given
//'   vector into ascending order.
//'   
//'   Applying the function \code{arma::sort_index()} twice:
//'   \code{arma::sort_index(arma::sort_index())}, calculates the \emph{reverse}
//'   permutation index to sort the vector from ascending order back into its
//'   original unsorted order.
//'   The permutation index produced by:
//'   \code{arma::sort_index(arma::sort_index())} is the \emph{reverse} of the
//'   permutation index produced by: \code{arma::sort_index()}.
//'   
//'   The ranks of the elements are equal to the \emph{reverse} permutation
//'   index.
//'   The function \code{calc_ranks()} calculates the \emph{reverse} permutation
//'   index.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random data
//' datav <- round(runif(7), 2)
//' # Calculate the ranks of the elements in two ways
//' all.equal(rank(datav), drop(HighFreq::calc_ranks(datav)))
//' # Create a time series of random data
//' datav <- xts::xts(runif(7), seq.Date(Sys.Date(), by=1, length.out=7))
//' # Calculate the ranks of the elements in two ways
//' all.equal(rank(coredata(datav)), drop(HighFreq::calc_ranks(datav)))
//' # Compare the speed of RcppArmadillo with R code
//' datav <- runif(7)
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=calc_ranks(datav),
//'   Rcode=rank(datav),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::uvec calc_ranks(const arma::vec& tseries) {
  
  return (arma::sort_index(arma::sort_index(tseries)) + 1);
  
}  // end calc_ranks




////////////////////////////////////////////////////////////
// Functions for rolling aggregations
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
//' Aggregate a time series of data into a single bar of \emph{OHLC} data.
//'
//' @export
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
  count_true[0] = tseries[0];
  // Loop over tseries
  for (arma::uword it = 1; it < length; it++) {
    if (tseries[it])
      // Add count number
      count_true[it] = count_true[it-1] + 1;
    else
      // Reset count to zero
      count_true[it] = tseries[it];
  }  // end for
  
  return count_true;
  
}  // end roll_count




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
//' @param \emph{endp} An \emph{integer vector} of end points.
//'
//' @return A \emph{matrix} with \emph{OHLC} data, with the number of rows equal
//'   to the number of \emph{endp} minus one.
//'   
//' @details
//'   The function \code{roll_ohlc()} performs a loop over the end points
//'   \emph{endp}, along the rows of the data \code{tseries}. At each end point,
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
//' endp <- HighFreq::calc_endpoints(NROW(ohlc), step=25)
//' # Aggregate over endp:
//' ohlcagg <- HighFreq::roll_ohlc(tseries=ohlc, endp=endp)
//' # Compare with xts::to.period()
//' ohlcagg_xts <- .Call("toPeriod", ohlc, as.integer(endp+1), TRUE, NCOL(ohlc), FALSE, FALSE, colnames(ohlc), PACKAGE="xts")
//' all.equal(ohlcagg, coredata(ohlcagg_xts), check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_ohlc(const arma::mat& tseries, arma::uvec endp) {
  
  // arma::uword nrows = tseries.n_rows;
  arma::uword ncols = tseries.n_cols;
  arma::uword numpts = endp.size();
  arma::mat ohlcagg(numpts-1, ncols, fill::zeros);
  
  // Perform loop over the end points
  for (arma::uword it = 1; it < numpts; it++) {
    // cout << "it: " << it << endl;
    // Aggregate the OHLC
    ohlcagg.row(it-1) = agg_ohlc(tseries.rows(endp(it-1)+1, endp(it)));
  }  // end for
  
  // Return the aggregations
  return ohlcagg;
  
}  // end roll_ohlc




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
//' returns <- zoo::coredata(na.omit(rutils::etfenv$returns$VTI))
//' # Calculate rolling sums over 11-period look-back intervals
//' sum_rolling <- HighFreq::roll_vec(returns, look_back=11)
//' # Compare HighFreq::roll_vec() with rutils::roll_sum()
//' all.equal(HighFreq::roll_vec(returns, look_back=11), 
//'          rutils::roll_sum(returns, look_back=11), 
//'          check.attributes=FALSE)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_vec(returns, look_back=11),
//'   Rcode=rutils::roll_sum(returns, look_back=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_vec(const arma::mat& tseries, arma::uword look_back) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat rolling_sum(nrows, 1);
  
  // Warmup period
  rolling_sum[0] = tseries[0];
  for (arma::uword it = 1; it < look_back; it++) {
    rolling_sum[it] = rolling_sum[it-1] + tseries[it];
  }  // end for
  
  // Remaining period
  for (arma::uword it = look_back; it < nrows; it++) {
    rolling_sum[it] = rolling_sum[it-1] + tseries[it] - tseries[it-look_back];
  }  // end for
  
  return rolling_sum;
  
}  // end roll_vec




////////////////////////////////////////////////////////////
//' Calculate the rolling weighted sums over a single-column \emph{time series}
//' or a single-column \emph{matrix} using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A single-column \emph{time series} or a single-column
//'   \emph{matrix}.
//' 
//' @param \code{weights} A single-column \emph{matrix} of weights.
//'
//' @return A single-column \emph{matrix} of the same length as the argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{roll_vecw()} calculates the rolling weighted sums of a
//'   single-column \emph{matrix} over its past values (a convolution with the
//'   single-column \emph{matrix} of weights), using \code{RcppArmadillo}. It
//'   performs a similar calculation as the standard \code{R} function
//'   \cr\code{stats::filter(x=series, filter=weights, method="convolution",
//'   sides=1)}, but it's over \code{6} times faster, and it doesn't produce any
//'   \code{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Define a single-column matrix of returns
//' returns <- zoo::coredata(na.omit(rutils::etfenv$returns$VTI))
//' # Create simple weights
//' weights <- matrix(c(1, rep(0, 10)))
//' # Calculate rolling weighted sums
//' weighted <- HighFreq::roll_vecw(tseries=returns, weights=weights)
//' # Compare with original
//' all.equal(zoo::coredata(returns), weighted, check.attributes=FALSE)
//' # Second example
//' # Create exponentially decaying weights
//' weights <- matrix(exp(-0.2*1:11))
//' weights <- weights/sum(weights)
//' # Calculate rolling weighted sums
//' weighted <- HighFreq::roll_vecw(tseries=returns, weights=weights)
//' # Calculate rolling weighted sums using filter()
//' filtered <- stats::filter(x=returns, filter=weights, method="convolution", sides=1)
//' # Compare both methods
//' all.equal(filtered[-(1:11)], weighted[-(1:11)], check.attributes=FALSE)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_vecw(tseries=returns, weights=weights),
//'   Rcode=stats::filter(x=returns, filter=weights, method="convolution", sides=1),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_vecw(const arma::mat& tseries, arma::mat& weights) {
  
  arma::uword nrows = tseries.n_rows;
  arma::uword look_back = weights.n_rows;
  arma::mat rolling_sum(nrows, 1);
  arma::mat weightsr = arma::reverse(weights);
  // arma::mat weightsr = weights;
  
  // Warmup period
  rolling_sum.rows(0, look_back-2) = tseries.rows(0, look_back-2);
  
  // Remaining periods
  for (arma::uword it = look_back-1; it < nrows; it++) {
    rolling_sum(it) = arma::dot(weightsr, tseries.rows(it-look_back+1, it));
  }  // end for
  
  return rolling_sum;
  
}  // end roll_vecw




////////////////////////////////////////////////////////////
//' Calculate the rolling convolutions (weighted sums) of a \emph{time series}
//' with a single-column \emph{matrix} of weights.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//' 
//' @param \code{weights} A single-column \emph{matrix} of weights.
//'
//' @return A \emph{matrix} with the same dimensions as the input
//'   argument \code{tseries}.
//'
//' @details
//'   The function \code{roll_conv()} calculates the convolutions of the
//'   \emph{matrix} columns with a single-column \emph{matrix} of weights.  It
//'   performs a loop over the \emph{matrix} rows and multiplies the past
//'   (higher) values by the weights.  It calculates the rolling weighted sums
//'   of the past values.
//'   
//'   The function \code{roll_conv()} uses the \code{RcppArmadillo} function
//'   \code{arma::conv2()}. It performs a similar calculation to the standard
//'   \code{R} function \cr\code{filter(x=tseries, filter=weights,
//'   method="convolution", sides=1)}, but it's over \code{6} times faster, and
//'   it doesn't produce any leading \code{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Calculate a time series of returns
//' returns <- na.omit(rutils::etfenv$returns[, c("IEF", "VTI")])
//' # Create simple weights equal to a 1 value plus zeros
//' weights <- matrix(c(1, rep(0, 10)), nc=1)
//' # Calculate rolling weighted sums
//' weighted <- HighFreq::roll_conv(returns, weights)
//' # Compare with original
//' all.equal(coredata(returns), weighted, check.attributes=FALSE)
//' # Second example
//' # Calculate exponentially decaying weights
//' weights <- exp(-0.2*(1:11))
//' weights <- matrix(weights/sum(weights), nc=1)
//' # Calculate rolling weighted sums
//' weighted <- HighFreq::roll_conv(returns, weights)
//' # Calculate rolling weighted sums using filter()
//' filtered <- filter(x=returns, filter=weights, method="convolution", sides=1)
//' # Compare both methods
//' all.equal(filtered[-(1:11), ], weighted[-(1:11), ], check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_conv(const arma::mat& tseries, const arma::mat& weights) {
  
  arma::uword look_back = weights.n_rows-2;
  arma::uword nrows = tseries.n_rows-1;
  
  // Calculate the convolutions
  arma::mat convmat = arma::conv2(tseries, weights, "full");
  
  // Copy the warmup period
  convmat.rows(0, look_back) = tseries.rows(0, look_back);
  
  return convmat.rows(0, nrows);
  
}  // end roll_conv




////////////////////////////////////////////////////////////
//' Calculate the rolling sums over a \emph{time series} or a \emph{matrix}
//' using \emph{Rcpp}.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of data points included in calculating the rolling sum (the default
//'   is \code{look_back = 1}).
//'   
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{roll_sum()} calculates the rolling sums over the
//'   columns of the data \code{tseries}.
//'   
//'   The function \code{roll_sum()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{tseries}.
//' 
//'   The function \code{roll_sum()} uses the fast \code{RcppArmadillo} function
//'   \code{arma::cumsum()}, without explicit loops.
//'   The function \code{roll_sum()} is several times faster than
//'   \code{rutils::roll_sum()} which uses vectorized \code{R} code.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' returns <- na.omit(rutils::etfenv$returns[, c("VTI", "IEF")])
//' # Define parameters
//' look_back <- 22
//' # Calculate rolling sums and compare with rutils::roll_sum()
//' c_sum <- HighFreq::roll_sum(returns, look_back)
//' r_sum <- rutils::roll_sum(returns, look_back)
//' all.equal(c_sum, coredata(r_sum), check.attributes=FALSE)
//' # Calculate rolling sums using R code
//' r_sum <- apply(zoo::coredata(returns), 2, cumsum)
//' lag_sum <- rbind(matrix(numeric(2*look_back), nc=2), r_sum[1:(NROW(r_sum) - look_back), ])
//' r_sum <- (r_sum - lag_sum)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_sum(const arma::mat& tseries, arma::uword look_back = 1) {
  
  // Calculate the cumulative sum
  arma::mat cumsumv = arma::cumsum(tseries, 0);
  
  // Return the differences of the cumulative sum
  return diffit(cumsumv, look_back, true);
  
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
//' @param \code{endp} An \emph{integer} vector of end points (the default is
//'   \code{endp = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{look_back} The number of end points in the look-back interval
//'   (the default is \code{look_back = 1}).
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
//' returns <- na.omit(rutils::etfenv$returns[, c("VTI", "IEF")])
//' # Define end points at 25 day intervals
//' endp <- HighFreq::calc_endpoints(NROW(returns), step=25)
//' # Define start points as 75 day lag of end points
//' startp <- HighFreq::calc_startpoints(endp, look_back=3)
//' # Calculate rolling sums using Rcpp
//' c_sum <- HighFreq::roll_sumep(returns, startp=startp, endp=endp)
//' # Calculate rolling sums using R code
//' r_sum <- sapply(1:NROW(endp), function(ep) {
//' colSums(returns[(startp[ep]+1):(endp[ep]+1), ])
//'   })  # end sapply
//' r_sum <- t(r_sum)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_sumep(const arma::mat& tseries, 
                     arma::uvec startp = 0, 
                     arma::uvec endp = 0, 
                     arma::uword step = 1, 
                     arma::uword look_back = 1, 
                     arma::uword stub = 0) {
  
  // Allocate end points
  arma::uword nrows = tseries.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    startpts = calc_startpoints(endpts, look_back);
    // Old code for startpts
    // Start points equal to end points lagged by look_back - without adding +1
    // arma::uword numpts = endpts.n_elem;
    // arma::uvec startpts = arma::join_cols(arma::zeros<uvec>(look_back), 
    //                                        endpts.subvec(0, numpts - look_back - 1));
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  // Allocate sums matrix
  arma::uword numpts = endpts.n_elem;
  arma::mat sums = arma::zeros<mat>(numpts, tseries.n_cols);
  
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
//' @param \code{endp} An \emph{integer} vector of end points (the default is
//'   \code{endp = NULL}).
//'   
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of data points included in calculating the rolling sum (the default
//'   is \code{look_back = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = NULL}).
//' 
//' @param \code{weights} A single-column \emph{matrix} of weights (the default
//'   is \code{weights = NULL}).
//'
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{roll_wsum()} calculates the rolling weighted sums over
//'   the columns of the data \code{tseries}.
//' 
//'   The function \code{roll_wsum()} calculates the rolling weighted sums as
//'   convolutions of the columns of \code{tseries} with the \emph{column
//'   vector} of weights using the \code{RcppArmadillo} function
//'   \code{arma::conv2()}.  It performs a similar calculation to the standard
//'   \code{R} function \cr\code{stats::filter(x=returns, filter=weights,
//'   method="convolution", sides=1)}, but it can be many times faster, and it
//'   doesn't produce any leading \code{NA} values.
//'   
//'   The function \code{roll_wsum()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{tseries}.
//' 
//'   The arguments \code{weights}, \code{endp}, and \code{stub} are
//'   optional.
//'   
//'   If the argument \code{weights} is not supplied, then simple sums are
//'   calculated, not weighted sums.
//'   
//'   If either the \code{stub} or \code{endp} arguments are supplied,
//'   then the rolling sums are calculated at the end points. 
//'   
//'   If only the argument \code{stub} is supplied, then the end points are
//'   calculated from the \code{stub} and \code{look_back} arguments. The first
//'   end point is equal to \code{stub} and the end points are spaced
//'   \code{look_back} periods apart.
//'   
//'   If the arguments \code{weights}, \code{endp}, and \code{stub} are
//'   not supplied, then the sums are calculated over a number of data points
//'   equal to \code{look_back}.
//'   
//'   The function \code{roll_wsum()} is also several times faster than
//'   \code{rutils::roll_sum()} which uses vectorized \code{R} code.
//'   
//'   Technical note:
//'   The function \code{roll_wsum()} has arguments with default values equal to
//'   \code{NULL}, which are implemented in \code{Rcpp} code.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Calculate historical returns
//' returns <- na.omit(rutils::etfenv$returns[, c("VTI", "IEF")])
//' # Define parameters
//' look_back <- 22
//' # Calculate rolling sums and compare with rutils::roll_sum()
//' c_sum <- HighFreq::roll_sum(returns, look_back)
//' r_sum <- rutils::roll_sum(returns, look_back)
//' all.equal(c_sum, coredata(r_sum), check.attributes=FALSE)
//' # Calculate rolling sums using R code
//' r_sum <- apply(zoo::coredata(returns), 2, cumsum)
//' lag_sum <- rbind(matrix(numeric(2*look_back), nc=2), r_sum[1:(NROW(r_sum) - look_back), ])
//' r_sum <- (r_sum - lag_sum)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' 
//' # Calculate rolling sums at end points
//' stu_b <- 21
//' c_sum <- HighFreq::roll_wsum(returns, look_back, stub=stu_b)
//' endp <- (stu_b + look_back*(0:(NROW(returns) %/% look_back)))
//' endp <- endp[endp < NROW(returns)]
//' r_sum <- apply(zoo::coredata(returns), 2, cumsum)
//' r_sum <- r_sum[endp+1, ]
//' lag_sum <- rbind(numeric(2), r_sum[1:(NROW(r_sum) - 1), ])
//' r_sum <- (r_sum - lag_sum)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' 
//' # Calculate rolling sums at end points - pass in endp
//' c_sum <- HighFreq::roll_wsum(returns, endp=endp)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' 
//' # Create exponentially decaying weights
//' weights <- exp(-0.2*(1:11))
//' weights <- matrix(weights/sum(weights), nc=1)
//' # Calculate rolling weighted sum
//' c_sum <- HighFreq::roll_wsum(returns, weights=weights)
//' # Calculate rolling weighted sum using filter()
//' filtered <- filter(x=returns, filter=weights, method="convolution", sides=1)
//' all.equal(c_sum[-(1:11), ], filtered[-(1:11), ], check.attributes=FALSE)
//' 
//' # Calculate rolling weighted sums at end points
//' c_sum <- HighFreq::roll_wsum(returns, endp=endp, weights=weights)
//' all.equal(c_sum, filtered[endp+1, ], check.attributes=FALSE)
//' 
//' # Create simple weights equal to a 1 value plus zeros
//' weights <- matrix(c(1, rep(0, 10)), nc=1)
//' # Calculate rolling weighted sum
//' weighted <- HighFreq::roll_wsum(returns, weights=weights)
//' # Compare with original
//' all.equal(coredata(returns), weighted, check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_wsum(const arma::mat& tseries,
                    Rcpp::Nullable<Rcpp::IntegerVector> endp = R_NilValue, 
                    arma::uword look_back = 1,
                    Rcpp::Nullable<int> stub = R_NilValue, 
                    Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat cumsumv;
  
  // Calculate the rolling sums
  if (weights.isNotNull()) {
    // Coerce weights from Rcpp to Armadillo vector
    arma::vec weights_vec = Rcpp::as<vec>(weights);
    arma::uword nweights = weights_vec.n_elem;
    // Calculate the weighted averages as convolutions
    cumsumv = arma::conv2(tseries, weights_vec, "full");
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
  if (endp.isNotNull()) {
    // Copy endp
    endpts = Rcpp::as<uvec>(endp);
  } else if (stub.isNotNull()) {
    // Calculate end points with stub
    endpts = arma::regspace<uvec>(Rcpp::as<uword>(stub), look_back, nrows + look_back);
    endpts = endpts.elem(find(endpts < nrows));
  }  // end if
  
  
  // Subset the rolling sums according the end points
  if (endpts.is_empty() && weights.isNotNull()) {
    // Do nothing
    // Return the weighted averages (convolutions) at each point
    // return cumsumv;
  } else if (endpts.is_empty() && !weights.isNotNull()) {
    // Return unweighted rolling sums at each point
    cumsumv = diffit(cumsumv, look_back, true);
  } else if (!endpts.is_empty() && weights.isNotNull()) {
    // Return the weighted averages (convolutions) at end points
    cumsumv = cumsumv.rows(endpts);
  } else if (!endpts.is_empty() && !weights.isNotNull()) {
    // Return the unweighted rolling sums at end points
    cumsumv = cumsumv.rows(endpts);
    cumsumv = diffit(cumsumv, 1, true);
  }  // end if
  
  return cumsumv;
  
}  // end roll_wsum




////////////////////////////////////////////////////////////
// Functions for rolling aggregations of streaming data
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
//' Calculate the running weighted means of streaming \emph{time series} data.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lambda} A \emph{numeric} decay factor to multiply past
//'   estimates.
//'   
//' @param \code{weights} A single-column \emph{matrix} of weights.
//'
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{run_mean()} calculates the running weighted means of
//'   the streaming \emph{time series} data \eqn{p_t} by recursively weighing
//'   present and past values using the decay factor \eqn{\lambda}:
//'   \deqn{
//'     \mu^w_t = (1-\lambda) w_t + \lambda \mu^w_{t-1}
//'   }
//'   \deqn{
//'     \mu^p_t = (1-\lambda) w_t p_t + \lambda \mu^p_{t-1}
//'   }
//'   Where \eqn{p_t} is the streaming data, \eqn{w_t} are the streaming
//'   weights, \eqn{\mu^w_t} are the running mean weights, and \eqn{\mu^p_t} are
//'   the running mean products of the data and the weights. 
//'   
//'   The running mean weighted value \eqn{\mu_t} is equal to the ratio of the
//'   data and weights products, divided by the mean weights:
//'   \deqn{
//'     \mu_t = \frac{\mu^p_t}{\mu^w_t}
//'   }
//' 
//'   If the \code{weights} argument is omitted, then the function
//'   \code{run_mean()} simply calculates the running means of \eqn{p_t}:
//'   \deqn{
//'     \mu_t = (1-\lambda) p_t + \lambda \mu_{t-1}
//'   }
//'   
//'   The above recursive formulas are convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//'   The formula is equivalent to a convolution with exponentially decaying
//'   weights, but it's faster.
//'   
//'   The value of the decay factor \eqn{\lambda} should be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the running mean values have a stronger
//'   dependence on past values.  This is equivalent to a long look-back
//'   interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the running mean values have a
//'   weaker dependence on past values.  This is equivalent to a short look-back
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
//' # Calculate the running means
//' lambda <- 0.95
//' means <- HighFreq::run_mean(closep, lambda=lambda)
//' # Calculate running means using R code
//' filtered <- (1-lambda)*filter(prices, 
//'   filter=lambda, init=as.numeric(prices[1, 1])/(1-lambda), 
//'   method="recursive")
//' all.equal(drop(means), unclass(filtered), check.attributes=FALSE)
//' 
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::run_mean(prices, lambda=lambda),
//'   Rcode=filter(prices, filter=lambda, init=as.numeric(prices[1, 1])/(1-lambda), method="recursive"),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//'   
//' # Create weights equal to the trading volumes
//' weights <- quantmod::Vo(ohlc)
//' # Calculate the running weighted means
//' meanw <- HighFreq::run_mean(prices, lambda=lambda, weights=weights)
//' # Plot dygraph of the running weighted means
//' datav <- xts(cbind(means, meanw), zoo::index(ohlc))
//' colnames(datav) <- c("means running", "means weighted")
//' dygraphs::dygraph(datav, main="Running Means") %>%
//'   dyOptions(colors=c("blue", "red"), strokeWidth=1) %>%
//'   dyLegend(show="always", width=500)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_mean(const arma::mat& tseries, 
                   double lambda, 
                   const arma::colvec& weights = 0) {
  
  arma::uword nrows = tseries.n_rows;
  arma::uword weights_rows = weights.n_elem;
  arma::mat means = arma::zeros<mat>(nrows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  if (weights_rows == 1) {
    means.row(0) = tseries.row(0);
    // Calculate means without weights
    for (arma::uword it = 1; it < nrows; it++) {
      // Calculate the means as the weighted sums
      means.row(it) = lambda1*tseries.row(it) + lambda*means.row(it-1);
    }  // end for
  } else if (weights_rows == nrows) {
    // Calculate means with weights
    arma::mat meanw = arma::zeros<mat>(nrows, 1);
    means.row(0) = weights.row(0)*tseries.row(0);
    meanw.row(0) = weights.row(0);
    for (arma::uword it = 1; it < nrows; it++) {
      // Calculate the means as the weighted sums
      meanw.row(it) = lambda1*weights.row(it) + lambda*meanw.row(it-1);
      means.row(it) = lambda1*weights.row(it)*tseries.row(it) + lambda*means.row(it-1);
      // Divide by the mean weight
      means.row(it-1) = means.row(it-1)/meanw.row(it-1);
    }  // end for
    means.row(nrows-1) = means.row(nrows-1)/meanw.row(nrows-1);
  }  // end if
  
  return means;
  
}  // end run_mean




////////////////////////////////////////////////////////////
//' Calculate the running maximum of streaming \emph{time series} data.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lambda} A \emph{numeric} decay factor to multiply past
//'   estimates.
//'   
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{run_max()} calculates the running maximum of streaming
//'   \emph{time series} data by recursively weighing present and past values
//'   using the decay factor \eqn{\lambda}.
//'
//'   It first calculates the running mean of streaming data:
//'   \deqn{
//'     \mu_t = (1-\lambda) p_t + \lambda \mu_{t-1}
//'   }
//'   Where \eqn{\mu_t} is the mean value at time \eqn{t}, and \eqn{p_t} is the
//'   streaming data.
//'
//'   It then calculates the running maximums of streaming data, \eqn{p^{max}_t}:
//'   \deqn{
//'     p^{max}_t = max(p_t, p^{max}_{t-1}) + (1-\lambda) (\mu_{t-1} - p^{max}_{t-1})
//'   }
//' 
//'   The second term pulls the maximum value down to the mean value, allowing
//'   it to gradually "forget" the maximum value from the more distant past.
//' 
//'   The above recursive formulas are convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//'   
//'   The value of the decay factor \eqn{\lambda} should be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the running maximum values have a stronger
//'   dependence on past values.  This is equivalent to a long look-back
//'   interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the running maximum values have a
//'   weaker dependence on past values.  This is equivalent to a short look-back
//'   interval.
//' 
//'   The function \code{run_max()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{tseries}.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical prices
//' prices <- zoo::coredata(quantmod::Cl(rutils::etfenv$VTI))
//' # Calculate the running maximums
//' lambda <- 0.9
//' maxs <- HighFreq::run_max(prices, lambda=lambda)
//' # Plot dygraph of VTI prices and running maximums
//' datav <- cbind(quantmod::Cl(rutils::etfenv$VTI), maxs)
//' colnames(datav) <- c("prices", "max")
//' colnamev <- colnames(datav)
//' dygraphs::dygraph(datav, main="VTI Prices and Running Maximums") %>%
//'   dySeries(name=colnamev[1], label=colnamev[1], strokeWidth=2, col="blue") %>%
//'   dySeries(name=colnamev[2], label=colnamev[2], strokeWidth=2, col="red")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_max(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat maxs = arma::zeros<mat>(nrows, tseries.n_cols);
  arma::mat means = arma::zeros<mat>(nrows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  // Perform loop over the rows
  means.row(0) = tseries.row(0);
  maxs.row(0) = tseries.row(0);
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the mean as a weighted sum
    means.row(it) = lambda1*tseries.row(it) + lambda*means.row(it-1);
    // Calculate the max from a weighted sum
    maxs.row(it) = arma::max(tseries.row(it), maxs.row(it-1) + lambda1*(means.row(it-1) - maxs.row(it-1)));
  }  // end for
  
  return maxs;
  
}  // end run_max




////////////////////////////////////////////////////////////
//' Calculate the running minimum of streaming \emph{time series} data.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lambda} A \emph{numeric} decay factor to multiply past
//'   estimates.
//'   
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{run_min()} calculates the running minimum of streaming
//'   \emph{time series} data by recursively weighing present and past values
//'   using the decay factor \eqn{\lambda}.
//'
//'   It first calculates the running mean of streaming data:
//'   \deqn{
//'     \mu_t = (1-\lambda) p_t + \lambda \mu_{t-1}
//'   }
//'   Where \eqn{\mu_t} is the mean value at time \eqn{t}, and \eqn{p_t} is the
//'   streaming data.
//'
//'   It then calculates the running minimums of streaming data, \eqn{p^{min}_t}:
//'   \deqn{
//'     p^{min}_t = min(p_t, p^{min}_{t-1}) + (1-\lambda) (\mu_{t-1} - p^{min}_{t-1})
//'   }
//' 
//'   The second term pulls the minimum value up to the mean value, allowing
//'   it to gradually "forget" the minimum value from the more distant past.
//' 
//'   The above recursive formula is convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//' 
//'   The value of the decay factor \eqn{\lambda} should be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the running minimum values have a stronger
//'   dependence on past values.  This is equivalent to a long look-back
//'   interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the running minimum values have a
//'   weaker dependence on past values.  This is equivalent to a short look-back
//'   interval.
//' 
//'   The function \code{run_min()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{tseries}.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical prices
//' prices <- zoo::coredata(quantmod::Cl(rutils::etfenv$VTI))
//' # Calculate the running minimums
//' lambda <- 0.9
//' mins <- HighFreq::run_min(prices, lambda=lambda)
//' # Plot dygraph of VTI prices and running minimums
//' datav <- cbind(quantmod::Cl(rutils::etfenv$VTI), mins)
//' colnames(datav) <- c("prices", "min")
//' colnamev <- colnames(datav)
//' dygraphs::dygraph(datav, main="VTI Prices and Running Minimums") %>%
//'   dySeries(name=colnamev[1], label=colnamev[1], strokeWidth=2, col="blue") %>%
//'   dySeries(name=colnamev[2], label=colnamev[2], strokeWidth=2, col="red")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_min(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat mins = arma::zeros<mat>(nrows, tseries.n_cols);
  arma::mat means = arma::zeros<mat>(nrows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  // Perform loop over the rows
  means.row(0) = tseries.row(0);
  mins.row(0) = tseries.row(0);
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the mean as a weighted sum
    means.row(it) = lambda1*tseries.row(it) + lambda*means.row(it-1);
    // Calculate the min from a weighted sum
    mins.row(it) = arma::min(tseries.row(it), mins.row(it-1) + lambda1*(means.row(it-1) - mins.row(it-1)));
  }  // end for
  
  return mins;
  
}  // end run_min




////////////////////////////////////////////////////////////
//' Calculate the running variance of streaming \emph{time series} of returns.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of returns.
//' 
//' @param \code{lambda} A \emph{numeric} decay factor to multiply past
//'   estimates.
//'   
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details
//'   The function \code{run_var()} calculates the running variance of a
//'   streaming \emph{time series} of returns, by recursively weighing the
//'   squared returns \eqn{r^2_t} minus the squared means \eqn{\mu^2_t}, with
//'   the past variance estimates \eqn{\sigma^2_{t-1}}, using the decay factor
//'   \eqn{\lambda}:
//'   \deqn{
//'     \mu_t = (1-\lambda) r_t + \lambda \mu_{t-1}
//'   }
//'   \deqn{
//'     \sigma^2_t = (1-\lambda) (r^2_t - \mu^2_t) + \lambda \sigma^2_{t-1}
//'   }
//'   Where \eqn{\sigma^2_t} is the variance estimate at time \eqn{t}, and
//'   \eqn{r_t} are the streaming returns data.
//' 
//'   The above recursive formula is convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//'   The formula is equivalent to a convolution with exponentially decaying
//'   weights, but it's faster.
//' 
//'   The value of the decay factor \eqn{\lambda} should be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the running variance values have a
//'   stronger dependence on past values.  This is equivalent to a long
//'   look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the running variance values have a
//'   weaker dependence on past values.  This is equivalent to a short look-back
//'   interval.
//' 
//'   The function \code{run_var()} performs the same calculation as the
//'   standard \code{R} function\cr\code{stats::filter(x=series,
//'   filter=weights, method="recursive")}, but it's several times faster.
//' 
//'   The function \code{run_var()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{tseries}.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' returns <- zoo::coredata(na.omit(rutils::etfenv$returns$VTI))
//' # Calculate the running variance
//' lambda <- 0.9
//' vars <- HighFreq::run_var(returns, lambda=lambda)
//' # Calculate running variance using R code
//' filtered <- (1-lambda)*filter(returns^2, filter=lambda, 
//'   init=as.numeric(returns[1, 1])^2/(1-lambda), 
//'   method="recursive")
//' all.equal(vars, unclass(filtered), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::run_var(returns, lambda=lambda),
//'   Rcode=filter(returns^2, filter=lambda, init=as.numeric(returns[1, 1])^2/(1-lambda), method="recursive"),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_var(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat vars = arma::square(tseries);
  arma::mat means = arma::zeros<mat>(nrows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  // Perform loop over the rows
  means.row(0) = tseries.row(0);
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the mean as the weighted sum
    means.row(it) = lambda1*tseries.row(it) + lambda*means.row(it-1);
    // Calculate the variance as the weighted sum of squared returns minus the squared means
    vars.row(it) = lambda1*(vars.row(it) - arma::square(means.row(it))) + lambda*vars.row(it-1);
  }  // end for
  
  return vars;
  
}  // end run_var




////////////////////////////////////////////////////////////
//' Calculate the running variance of streaming \emph{OHLC} price data.
//' 
//' @param \code{ohlc} A \emph{time series} or a \emph{matrix} with \emph{OHLC}
//'   price data.
//'   
//' @param \code{lambda} A \emph{numeric} decay factor to multiply past
//'   estimates.
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
//'   
//'   The function \code{run_var_ohlc()} does not calculate the logarithm of
//'   the prices.
//'   So if the argument \code{ohlc} contains dollar prices then
//'   \code{run_var_ohlc()} calculates the dollar variance.
//'   If the argument \code{ohlc} contains the log prices then
//'   \code{run_var_ohlc()} calculates the percentage variance.
//'   
//'   The function \code{run_var_ohlc()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, so it's many times faster than the equivalent \code{R}
//'   code.
//'
//' @examples
//' \dontrun{
//' # Extract the log OHLC prices of VTI
//' ohlc <- log(rutils::etfenv$VTI)
//' # Calculate the running variance
//' var_running <- HighFreq::run_var_ohlc(ohlc, lambda=0.8)
//' # Calculate the rolling variance
//' var_rolling <- HighFreq::roll_var_ohlc(ohlc, look_back=5, method="yang_zhang", scale=FALSE)
//' datav <- cbind(var_running, var_rolling)
//' colnames(datav) <- c("running", "rolling")
//' colnamev <- colnames(datav)
//' datav <- xts::xts(datav, index(ohlc))
//' # dygraph plot of VTI running versus rolling volatility
//' dygraphs::dygraph(sqrt(datav[-(1:111), ]), main="Running and Rolling Volatility of VTI") %>%
//'   dyOptions(colors=c("red", "blue"), strokeWidth=1) %>%
//'   dyLegend(show="always", width=500)
//' # Compare the speed of running versus rolling volatility
//' library(microbenchmark)
//' summary(microbenchmark(
//'   running=HighFreq::run_var_ohlc(ohlc, lambda=0.8),
//'   rolling=HighFreq::roll_var_ohlc(ohlc, look_back=5, method="yang_zhang", scale=FALSE),
//'   times=10))[, c(1, 4, 5)]
//' }
//' @export
// [[Rcpp::export]]
arma::mat run_var_ohlc(const arma::mat& ohlc, 
                       double lambda) {
  
  // Allocate variance matrix
  arma::uword nrows = ohlc.n_rows;
  arma::mat vars = arma::zeros<mat>(nrows, 1);
  double lambda1 = 1-lambda;
  double coeff = 0.134;

  // Calculate all the different intra-day and day-over-day returns 
  // (differences of OHLC prices)
  arma::mat closep = ohlc.col(3);
  arma::mat open_close(closep.n_rows, 1);
  open_close = (ohlc.col(0) - lagit(closep, 1, false));
  arma::mat close_open = (closep - ohlc.col(0));
  arma::mat close_high = (closep - ohlc.col(1));
  arma::mat close_low = (closep - ohlc.col(2));
  arma::mat high_low = (ohlc.col(1) - ohlc.col(2));
  arma::mat high_open = (ohlc.col(1) - ohlc.col(0));
  arma::mat low_open = (ohlc.col(2) - ohlc.col(0));
  
  // Perform loop over the rows
  vars.row(0) = arma::square(open_close.row(0)) + coeff*arma::square(close_open.row(0)) +
    (coeff-1)*(close_high.row(0)*high_open.row(0) + close_low.row(0)*low_open.row(0));
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the variance as the weighted sum of squared returns minus the squared means
    vars.row(it) = lambda1*(arma::square(open_close.row(it)) + coeff*arma::square(close_open.row(it)) +
      (coeff-1)*(close_high.row(it)*high_open.row(it) + close_low.row(it)*low_open.row(it))) + lambda*vars.row(it-1);
  }  // end for
  
  return vars;
  
}  // end run_var_ohlc




////////////////////////////////////////////////////////////
//' Calculate the running covariance of two streaming \emph{time series} of
//' returns.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} with two
//'   columns of returns data.
//' 
//' @param \code{lambda} A \emph{numeric} decay factor to multiply past
//'   estimates.
//'   
//' @return A \emph{matrix} with three columns of data: the covariance and the
//'   variances of the two columns of the argument \code{tseries}.
//'
//' @details
//'   The function \code{run_covar()} calculates the running covariance of two
//'   streaming \emph{time series} of returns, by recursively weighing the
//'   products of their returns minus their means, with past covariance
//'   estimates \eqn{\sigma^{cov}_{t-1}}, using the decay factor \eqn{\lambda}:
//'   \deqn{
//'     \mu^1_t = (1-\lambda) r^1_t + \lambda \mu^1_{t-1}
//'   }
//'   \deqn{
//'     \mu^2_t = (1-\lambda) r^2_t + \lambda \mu^2_{t-1}
//'   }
//'   \deqn{
//'     \sigma^{cov}_t = (1-\lambda) (r^1_t - \mu^1_t) (r^2_t - \mu^2_t) + \lambda \sigma^{cov}_{t-1}
//'   }
//'   Where \eqn{\sigma^{cov}_t} is the covariance estimate at time \eqn{t},
//'   \eqn{r^1_t} and \eqn{r^2_t} are the two streaming returns data, and
//'   \eqn{\mu^1_t} and \eqn{\mu^2_t} are the means of the returns.
//' 
//'   The above recursive formula is convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//'   The formula is equivalent to a convolution with exponentially decaying
//'   weights, but it's faster.
//' 
//'   The value of the decay factor \eqn{\lambda} should be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the running covariance values have a
//'   stronger dependence on past values.  This is equivalent to a long
//'   look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the running covariance values have
//'   a weaker dependence on past values.  This is equivalent to a short
//'   look-back interval.
//' 
//'   The function \code{run_covar()} returns three columns of data: the
//'   covariance and the variances of the two columns of the argument
//'   \code{tseries}.  This allows calculating the running correlation.
//' 
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' returns <- zoo::coredata(na.omit(rutils::etfenv$returns[, c("IEF", "VTI")]))
//' # Calculate the running covariance
//' lambda <- 0.9
//' covars <- HighFreq::run_covar(returns, lambda=lambda)
//' # Calculate running covariance using R code
//' filtered <- (1-lambda)*filter(returns[, 1]*returns[, 2], 
//'   filter=lambda, init=as.numeric(returns[1, 1]*returns[1, 2])/(1-lambda), 
//'   method="recursive")
//' all.equal(covars[, 1], unclass(filtered), check.attributes=FALSE)
//' # Calculate the running correlation
//' correl <- covars[, 1]/sqrt(covars[, 2]*covars[, 3])
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_covar(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::mat vars = arma::square(tseries);
  arma::mat covar = arma::zeros<mat>(nrows, 1);
  arma::mat means = arma::zeros<mat>(nrows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  // Perform loop over the rows
  means.row(0) = tseries.row(0);
  covar(0) = tseries(0, 0)*tseries(0, 1);
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the mean as the weighted sum
    means.row(it) = lambda1*tseries.row(it) + lambda*means.row(it-1);
    // Calculate the covariance as the weighted sum of products of returns
    vars.row(it) = lambda1*(vars.row(it) - arma::square(means.row(it))) + lambda*vars.row(it-1);
    covar.row(it) = lambda1*((tseries(it, 0)-means(it, 0))*(tseries(it, 1)-means(it, 1))) + lambda*covar.row(it-1);
  }  // end for
  
  return arma::join_rows(covar, vars);
  
  // Old code version below - produces same output and it's slightly faster, but less elegant
  // arma::uword nrows = tseries.n_rows;
  // arma::mat var1 = arma::square(tseries.col(0));
  // arma::mat var2 = arma::square(tseries.col(1));
  // arma::mat covar = arma::zeros<mat>(nrows, 1);
  // arma::mat means = arma::zeros<mat>(nrows, tseries.n_cols);
  // double lambda1 = 1-lambda;
  // 
  // // Perform loop over the rows
  // means.row(0) = tseries.row(0);
  // covar(0) = tseries(0, 0)*tseries(0, 1);
  // for (arma::uword it = 1; it < nrows; it++) {
  //   // Calculate the mean as the weighted sum
  //   means.row(it) = lambda1*tseries.row(it) + lambda*means.row(it-1);
  //   // Calculate the covariance as the weighted sum of products of returns
  //   var1(it) = lambda1*(var1(it) - pow(means(it, 0), 2)) + lambda*var1(it-1);
  //   var2(it) = lambda1*(var2(it) - pow(means(it, 1), 2)) + lambda*var2(it-1);
  //   covar.row(it) = lambda1*((tseries(it, 0)-means(it, 0))*(tseries(it, 1)-means(it, 1))) + lambda*covar.row(it-1);
  // }  // end for
  // 
  // return arma::join_rows(covar, var1, var2);
  
}  // end run_covar




////////////////////////////////////////////////////////////
//' Perform running regressions of streaming \emph{time series} of response and
//' predictor data, and calculate the alphas, betas, and the residuals.
//' 
//' @param \code{response} A single-column \emph{time series} or a single-column
//'   \emph{matrix} of response data.
//' 
//' @param \code{predictor} A \emph{time series} or a \emph{matrix} of predictor
//'   data.
//' 
//' @param \code{lambda} A \emph{numeric} decay factor to multiply past
//'   estimates.
//'   
//' @param \code{method} A \emph{string} specifying the method for scaling the
//'   residuals (see Details) (the default is \code{method = "none"} - no
//'   scaling)
//'   
//' @return A \emph{matrix} with the regression alphas, betas, and residuals.
//'
//' @details
//'   The function \code{run_reg()} calculates the vectors of \emph{alphas}
//'   \eqn{\alpha_t}, \emph{betas} \eqn{\beta_t}, and the \emph{residuals}
//'   \eqn{\epsilon_t} of running regressions, by recursively weighing the
//'   current estimates with past estimates, using the decay factor
//'   \eqn{\lambda}:
//'   \deqn{
//'     \mu^r_t = (1-\lambda) r^r_t + \lambda \mu^r_{t-1}
//'   }
//'   \deqn{
//'     \mu^p_t = (1-\lambda) r^p_t + \lambda \mu^p_{t-1}
//'   }
//'   \deqn{
//'     \sigma^2_t = (1-\lambda) ({r^p_t}^2 - {\mu^p_t}^2) + \lambda \sigma^2_{t-1}
//'   }
//'   \deqn{
//'     \sigma^{cov}_t = (1-\lambda) (r^r_t - \mu^r_t) (r^p_t - \mu^p_t) + \lambda \sigma^{cov}_{t-1}
//'   }
//'   \deqn{
//'     \beta_t = (1-\lambda) \frac{\sigma^{cov}_t}{\sigma^2_t} + \lambda \beta_{t-1}
//'   }
//'   \deqn{
//'     \epsilon_t = (1-\lambda) (r^r_t - \beta_t r^p_t) + \lambda \epsilon_{t-1}
//'   }
//'   Where \eqn{\sigma^{cov}_t} are the covariances between the response and
//'   predictor data at time \eqn{t};
//'   \eqn{\sigma^2_t} is the vector of predictor variances,
//'   and \eqn{r^r_t} and \eqn{r^p_t} are the streaming data of the response
//'   and predictor data.
//' 
//'   The matrices \eqn{\sigma^2}, \eqn{\sigma^{cov}}, and \eqn{\beta} have the
//'   same dimensions as the input argument \code{predictor}.
//'
//'   The above recursive formulas are convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//'   The formula is equivalent to a convolution with exponentially decaying
//'   weights, but it's faster to calculate.
//'
//'   The value of the decay factor \eqn{\lambda} should be in the range between
//'   \code{0} and \code{1}.
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the running \emph{z-score} values have a
//'   stronger dependence on past values.  This is equivalent to a long
//'   look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the running \emph{z-score} values
//'   have a weaker dependence on past values.  This is equivalent to a short
//'   look-back interval.
//' 
//'   The \emph{residuals} may be scaled by their volatilities. The default is
//'   \code{method = "none"} - no scaling.
//'   If the argument \code{method = "scale"} then the \emph{residuals}
//'   \eqn{\epsilon_t} are divided by their volatilities \eqn{\sigma^{\epsilon}}
//'   without subtracting their means:
//'   \deqn{
//'     \epsilon_t = \frac{\epsilon_t}{\sigma^{\epsilon}}
//'   }
//'   If the argument \code{method = "standardize"} then the means
//'   \eqn{\mu_{\epsilon}} are subtracted from the \emph{residuals}, and then
//'   they are divided by their volatilities \eqn{\sigma^{\epsilon}}:
//'   \deqn{
//'     \epsilon_t = \frac{\epsilon_t - \mu_{\epsilon}}{\sigma^{\epsilon}}
//'   }
//' 
//'   The function \code{run_reg()} returns multiple columns of data. If the
//'   matrix \code{predictor} has \code{n} columns then \code{run_reg()} returns
//'   a matrix with \code{n+2} columns.  The first column contains the
//'   \emph{residuals}, the second the \emph{alphas}, and the last columns
//'   contain the \emph{betas}.
//' 
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' returns <- na.omit(rutils::etfenv$returns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' response <- returns[, 1]
//' # Predictor matrix equals VTI and IEF returns
//' predictor <- returns[, -1]
//' # Calculate the running regressions
//' lambda <- 0.9
//' regs <- HighFreq::run_reg(response=response, predictor=predictor, lambda=lambda)
//' # Plot the running alphas
//' datav <- cbind(cumsum(response), regs[, 1])
//' colnames(datav) <- c("XLF", "alphas")
//' colnamev <- colnames(datav)
//' dygraphs::dygraph(datav, main="Alphas of XLF Versus VTI and IEF") %>%
//'   dyAxis("y", label=colnamev[1], independentTicks=TRUE) %>%
//'   dyAxis("y2", label=colnamev[2], independentTicks=TRUE) %>%
//'   dySeries(name=colnamev[1], axis="y", label=colnamev[1], strokeWidth=1, col="blue") %>%
//'   dySeries(name=colnamev[2], axis="y2", label=colnamev[2], strokeWidth=1, col="red")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_reg(const arma::mat& response, 
                  const arma::mat& predictor,
                  double lambda, 
                  std::string method = "none") {
  
  arma::uword nrows = predictor.n_rows;
  arma::uword ncols = predictor.n_cols;
  arma::mat means_resp = arma::zeros<mat>(nrows, 1);
  arma::mat means_pred = arma::zeros<mat>(nrows, ncols);
  arma::mat vars = arma::square(predictor);
  arma::mat covars = arma::zeros<mat>(nrows, ncols);
  arma::mat betas = arma::zeros<mat>(nrows, ncols);
  arma::mat alphas = arma::zeros<mat>(nrows, 1);
  arma::mat resids = arma::zeros<mat>(nrows, 1);
  arma::mat varz = arma::ones<mat>(nrows, 1);
  arma::mat meanz = arma::zeros<mat>(nrows, 1);
  double lambda1 = 1-lambda;
  
  // Perform loop over the rows
  means_resp.row(0) = response.row(0);
  means_pred.row(0) = predictor.row(0);
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the mean as the weighted sum
    means_resp.row(it) = lambda1*response.row(it) + lambda*means_resp.row(it-1);
    means_pred.row(it) = lambda1*predictor.row(it) + lambda*means_pred.row(it-1);
    // cout << "Calculating vars: " << it << endl;
    vars.row(it) = lambda1*arma::square(predictor.row(it)-means_pred.row(it)) + lambda*vars.row(it-1);
    // cout << "Calculating covars: " << it << endl;
    covars.row(it) = lambda1*((response.row(it)-means_resp.row(it))*(predictor.row(it)-means_pred.row(it))) + lambda*covars.row(it-1);
    // cout << "Calculating betas: " << it << endl;
    // Calculate the alphas and betas.
    betas.row(it) = lambda1*covars.row(it)/vars.row(it) + lambda*betas.row(it-1);
    alphas.row(it) = lambda1*(means_resp.row(it) - arma::dot(betas.row(it), means_pred.row(it))) + lambda*alphas.row(it-1);
    // cout << "Calculating resids: " << it << endl;
    // Calculate the residuals.
    resids.row(it) = lambda1*(response.row(it) - arma::dot(betas.row(it), predictor.row(it))) + lambda*resids.row(it-1);
    // Calculate the mean and variance of the residuals.
    meanz.row(it) = lambda1*resids.row(it) + lambda*meanz.row(it-1);
    varz.row(it) = lambda1*arma::square(resids.row(it) - meanz.row(it)) + lambda*varz.row(it-1);
  }  // end for
  
  if (method == "scale") {
    // Divide the residuals by their volatility
    resids = resids/sqrt(varz);
  } else if (method == "standardize") {
    // De-mean the residuals and divide them by their volatility
    resids = (resids - meanz)/sqrt(varz);
  }  // end if
  
  return join_rows(resids, alphas, betas);

}  // end run_reg




////////////////////////////////////////////////////////////
//' Calculate the z-scores of running regressions of streaming \emph{time
//' series} of returns.
//' 
//' @param \code{response} A single-column \emph{time series} or a single-column
//'   \emph{matrix} of response data.
//' 
//' @param \code{predictor} A \emph{time series} or a \emph{matrix} of predictor
//'   data.
//' 
//' @param \code{lambda} A \emph{numeric} decay factor to multiply past
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
//'   \eqn{\beta_t} and the residuals \eqn{\epsilon_t} of running regressions by
//'   recursively weighing the current estimates with past estimates, using the
//'   decay factor \eqn{\lambda}:
//'   \deqn{
//'     \sigma^2_t = (1-\lambda) {r^p_t}^2 + \lambda \sigma^2_{t-1}
//'   }
//'   \deqn{
//'     \sigma^{cov}_t = (1-\lambda) r^r_t r^p_t + \lambda \sigma^{cov}_{t-1}
//'   }
//'   \deqn{
//'     \beta_t = (1-\lambda) \frac{\sigma^{cov}_t}{\sigma^2_t} + \lambda \beta_{t-1}
//'   }
//'   \deqn{
//'     \epsilon_t = (1-\lambda) (r^r_t - \beta_t r^p_t) + \lambda \epsilon_{t-1}
//'   }
//'   Where \eqn{\sigma^{cov}_t} is the vector of covariances between the
//'   response and predictor returns, at time \eqn{t};
//'   \eqn{\sigma^2_t} is the vector of predictor variances,
//'   and \eqn{r^r_t} and \eqn{r^p_t} are the streaming returns of the response
//'   and predictor data.
//' 
//'   The above formulas for \eqn{\sigma^2} and \eqn{\sigma^{cov}} are
//'   approximate because they don't subtract the means before squaring the
//'   returns.  But they're very good approximations for daily returns.
//' 
//'   The matrices \eqn{\sigma^2}, \eqn{\sigma^{cov}}, \eqn{\beta} have the same
//'   dimensions as the input argument \code{predictor}.
//'
//'   If the argument \code{demean = TRUE} (the default) then the
//'   \emph{z-scores} \eqn{z_t} are calculated as equal to the residuals
//'   \eqn{\epsilon_t} minus their means \eqn{\mu_{\epsilon}}, divided by their
//'   volatilities \eqn{\sigma^{\epsilon}}:
//'   \deqn{
//'     z_t = \frac{\epsilon_t - \mu_{\epsilon}}{\sigma^{\epsilon}}
//'   }
//'   If the argument \code{demean = FALSE} then the \emph{z-scores} are
//'   only divided by their volatilities without subtracting their means:
//'   \deqn{
//'     z_t = \frac{\epsilon_t}{\sigma^{\epsilon}}
//'   }
//' 
//'   The above recursive formulas are convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//'   The formula is equivalent to a convolution with exponentially decaying
//'   weights, but it's faster to calculate.
//' 
//'   The value of the decay factor \eqn{\lambda} should be in the range between
//'   \code{0} and \code{1}.
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the running \emph{z-score} values have a
//'   stronger dependence on past values.  This is equivalent to a long
//'   look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the running \emph{z-score} values
//'   have a weaker dependence on past values.  This is equivalent to a short
//'   look-back interval.
//' 
//'   The function \code{run_zscores()} returns multiple columns of data. 
//'   If the matrix \code{predictor} has \code{n} columns then \code{run_zscores()}
//'   returns a matrix with \code{2n+1} columns.  The first column contains the
//'   \emph{z-scores}, and the remaining columns contain the \emph{betas} and
//'   the \emph{variances} of the predictor data.
//' 
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' returns <- na.omit(rutils::etfenv$returns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' response <- returns[, 1]
//' # Predictor matrix equals VTI and IEF returns
//' predictor <- returns[, -1]
//' # Calculate the running z-scores
//' lambda <- 0.9
//' zscores <- HighFreq::run_zscores(response=response, predictor=predictor, lambda=lambda)
//' # Plot the running z-scores
//' datav <- cbind(cumsum(response), zscores[, 1])
//' colnames(datav) <- c("XLF", "zscores")
//' colnamev <- colnames(datav)
//' dygraphs::dygraph(datav, main="Z-Scores of XLF Versus VTI and IEF") %>%
//'   dyAxis("y", label=colnamev[1], independentTicks=TRUE) %>%
//'   dyAxis("y2", label=colnamev[2], independentTicks=TRUE) %>%
//'   dySeries(name=colnamev[1], axis="y", label=colnamev[1], strokeWidth=1, col="blue") %>%
//'   dySeries(name=colnamev[2], axis="y2", label=colnamev[2], strokeWidth=1, col="red")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_zscores(const arma::mat& response, 
                      const arma::mat& predictor,
                      double lambda, 
                      bool demean = true) {
  
  arma::uword nrows = predictor.n_rows;
  arma::uword ncols = predictor.n_cols;
  // arma::mat var1 = arma::square(tseries.col(0));
  arma::mat vars = arma::square(predictor);
  arma::mat betas = arma::ones<mat>(nrows, ncols);
  arma::mat zscores = arma::ones<mat>(nrows, 1);
  arma::mat varz = arma::ones<mat>(nrows, 1);
  arma::mat meanz = arma::zeros<mat>(nrows, 1);
  double lambda1 = 1-lambda;
  
  // Multiply each column of predictor by the response.
  arma::mat covars = predictor;
  covars.each_col() %= response;
  
  // Perform loop over the rows
  for (arma::uword it = 1; it < nrows; it++) {
    // Calculate the z-score as the weighted sum of products of returns.
    // cout << "Calculating vars: " << it << endl;
    vars.row(it) = lambda1*vars.row(it) + lambda*vars.row(it-1);
    // cout << "Calculating covars: " << it << endl;
    covars.row(it) = lambda1*covars.row(it) + lambda*covars.row(it-1);
    // cout << "Calculating betas: " << it << endl;
    betas.row(it) = lambda1*covars.row(it)/vars.row(it) + lambda*betas.row(it-1);
    // cout << "Calculating zscores: " << it << endl;
    zscores.row(it) = lambda1*(response.row(it) - arma::dot(betas.row(it), predictor.row(it))) + lambda*zscores.row(it-1);
    // Calculate the mean and variance of the z-scores.
    meanz.row(it) = lambda1*zscores.row(it) + lambda*meanz.row(it-1);
    varz.row(it) = lambda1*arma::square(zscores.row(it) - zscores.row(it-1)) + lambda*varz.row(it-1);
  }  // end for
  
  if (demean)
    return join_rows((zscores - meanz)/sqrt(varz), betas, vars);
  else
    return join_rows(zscores/sqrt(varz), betas, vars);

}  // end run_zscores




////////////////////////////////////////////////////////////
// Functions for statistics
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
// Define C++ enum type for the different methods for regularization,
// calculating variance, skewness, kurtosis, covariance, regression, 
// and matrix inverse.
enum methodenum {moment, least_squares, quantile, nonparametric, regular, ranksharpe, 
              max_sharpe, max_sharpe_median, min_var, min_varpca, rank, rankrob};

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
  else if (method == "ranksharpe")
    return methodenum::ranksharpe;
  else if (method == "max_sharpe")
    return methodenum::max_sharpe;
  else if (method == "max_sharpe_median")
    return methodenum::max_sharpe_median;
  else if (method == "min_var")
    return methodenum::min_var;
  else if (method == "min_varpca")
    return methodenum::min_varpca;
  else if (method == "rank")
    return methodenum::rank;
  else if (method == "rankrob")
    return methodenum::rankrob;
  else 
    return methodenum::moment;
}  // end calc_method




////////////////////////////////////////////////////////////
//' Calculate the mean (location) of the columns of a \emph{time series} or a
//' \emph{matrix} using \code{RcppArmadillo}.
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'
//' @param \code{method} A \emph{string} specifying the type of the mean
//'   (location) model (the default is \code{method = "moment"} - see Details).
//'
//' @param \code{conf_lev} The confidence level for calculating the
//'   quantiles (the default is \code{conf_lev = 0.75}).
//'
//' @return A single-row matrix with the mean (location) of the columns of
//'   \code{tseries}.
//'
//' @details
//'   The function \code{calc_mean()} calculates the mean (location) values of
//'   the columns of the \emph{time series} \code{tseries} using
//'   \code{RcppArmadillo} \code{C++} code.
//'
//'   If \code{method = "moment"} (the default) then \code{calc_mean()}
//'   calculates the location as the mean - the first moment of the data.
//'
//'   If \code{method = "quantile"} then it calculates the location \eqn{\mu} as
//'   the sum of the quantiles as follows:
//'   \deqn{
//'     \mu = q_{\alpha} + q_{1-\alpha}
//'   }
//'   Where \eqn{\alpha} is the confidence level for calculating the quantiles.
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
//' returns <- na.omit(rutils::etfenv$returns[, c("XLP", "VTI")])
//' # Calculate the column means in RcppArmadillo
//' HighFreq::calc_mean(returns)
//' # Calculate the column means in R
//' sapply(returns, mean)
//' # Compare the values
//' all.equal(drop(HighFreq::calc_mean(returns)), 
//'   sapply(returns, mean), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_mean(returns),
//'   Rcode=sapply(returns, mean),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the quantile mean (location)
//' HighFreq::calc_mean(returns, method="quantile", conf_lev=0.9)
//' # Calculate the quantile mean (location) in R
//' colSums(sapply(returns, quantile, c(0.9, 0.1), type=5))
//' # Compare the values
//' all.equal(drop(HighFreq::calc_mean(returns, method="quantile", conf_lev=0.9)), 
//'   colSums(sapply(returns, quantile, c(0.9, 0.1), type=5)), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_mean(returns, method="quantile", conf_lev=0.9),
//'   Rcode=colSums(sapply(returns, quantile, c(0.9, 0.1), type=5)),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the column medians in RcppArmadillo
//' HighFreq::calc_mean(returns, method="nonparametric")
//' # Calculate the column medians in R
//' sapply(returns, median)
//' # Compare the values
//' all.equal(drop(HighFreq::calc_mean(returns, method="nonparametric")), 
//'   sapply(returns, median), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_mean(returns, method="nonparametric"),
//'   Rcode=sapply(returns, median),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_mean(const arma::mat& tseries,
                    std::string method = "moment", 
                    double conf_lev = 0.75) {
  
  // Switch for the different methods of location
  switch(calc_method(method)) {
  case methodenum::moment: {  // moment
    return arma::mean(tseries);
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-conf_lev, conf_lev};
    arma::mat quantiles = arma::quantile(tseries, levels);
    return (quantiles.row(0) + quantiles.row(1));
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
//' returns <- rnorm(1e6)
//' # Compare calc_varvec() with standard var()
//' all.equal(HighFreq::calc_varvec(returns), 
//'   var(returns))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_varvec(returns),
//'   Rcode=var(returns),
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
//' @param \code{method} A \emph{string} specifying the type of the dispersion
//'   model (the default is \code{method = "moment"} - see Details).
//'    
//' @param \code{conf_lev} The confidence level for calculating the
//'   quantiles (the default is \code{conf_lev = 0.75}).
//'
//' @return A row vector equal to the dispersion of the columns of the matrix
//'   \code{tseries}.
//'
//' @details
//'   The dispersion is a measure of the variability of the data.  Examples of
//'   dispersion are the variance and the Median Absolute Deviation (\emph{MAD}).
//'
//'   The function \code{calc_var()} calculates the dispersion of the
//'   columns of a \emph{time series} or a \emph{matrix} of data using
//'   \code{RcppArmadillo} \code{C++} code.
//'   
//'   If \code{method = "moment"} (the default) then \code{calc_var()}
//'   calculates the dispersion as the second moment of the data \eqn{\sigma^2}
//'   (the variance).
//'
//'   If \code{method = "moment"} then \code{calc_var()} performs the same
//'   calculation as the function \code{colVars()} from package
//'   \href{https://cran.r-project.org/web/packages/matrixStats/index.html}{matrixStats},
//'   but it's much faster because it uses \code{RcppArmadillo} \code{C++} code.
//'
//'   If \code{method = "quantile"} then it calculates the dispersion as the
//'   difference between the quantiles as follows:
//'   \deqn{
//'     \mu = q_{\alpha} - q_{1-\alpha}
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
//' returns <- na.omit(rutils::etfenv$returns[, c("VTI", "XLF")])
//' # Compare HighFreq::calc_var() with standard var()
//' all.equal(drop(HighFreq::calc_var(returns)), 
//'   apply(returns, 2, var), check.attributes=FALSE)
//' # Compare HighFreq::calc_var() with matrixStats
//' all.equal(drop(HighFreq::calc_var(returns)), 
//'   matrixStats::colVars(returns), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with matrixStats and with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_var(returns),
//'   matrixStats=matrixStats::colVars(returns),
//'   Rcode=apply(returns, 2, var),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Compare HighFreq::calc_var() with stats::mad()
//' all.equal(drop(HighFreq::calc_var(returns, method="nonparametric")), 
//'   sapply(returns, mad), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with stats::mad()
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_var(returns, method="nonparametric"),
//'   Rcode=sapply(returns, mad),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_var(const arma::mat& tseries,
                   std::string method = "moment", 
                   double conf_lev = 0.75) {
  
  // Return zeros if not enough data
  if (tseries.n_rows < 3) {
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end if
  
  // Switch for the different methods of dispersion
  switch(calc_method(method)) {
  case methodenum::moment: {  // moment
    return arma::var(tseries);
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-conf_lev, conf_lev};
    arma::mat quantiles = arma::quantile(tseries, levels);
    return (quantiles.row(1) - quantiles.row(0));
  }  // end quantile
  case methodenum::nonparametric: {  // MAD
    double ncols = tseries.n_cols;
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
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end default
  }  // end switch
  
}  // end calc_var




////////////////////////////////////////////////////////////
//' Calculate the variance of returns aggregated over end points. 
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of prices.
//'
//' @param \code{step} The number of periods in each interval between
//'   neighboring end points.
//' 
//' @return The variance of aggregated returns.
//'
//' @details
//'   The function \code{calc_var_ag()} calculates the variance of returns
//'   aggregated over end points.
//'
//'   It first calculates the end points spaced apart by the number of periods
//'   equal to the argument \code{step}.  Then it calculates the aggregated
//'   returns by differencing the prices \code{tseries} calculated at the end
//'   points. Finally it calculates the variance of the returns.
//'
//'   If there are extra periods that don't fit over the length of
//'   \code{tseries}, then \code{calc_var_ag()} loops over all possible stub
//'   intervals, then it calculates all the corresponding variance values, and
//'   averages them.
//'
//'   For example, if the number of rows of \code{tseries} is equal to
//'   \code{20}, and \code{step=3} then \code{6} end points fit over the length
//'   of \code{tseries}, and there are \code{2} extra periods that must fit into
//'   stubs, either at the beginning or at the end (or both).
//' 
//'   The aggregated volatility \eqn{\sigma_t} scales (increases) with the
//'   length of the aggregation interval \eqn{\Delta t} raised to the power of
//'   the \emph{Hurst exponent} \eqn{H}:
//'     \deqn{
//'       \sigma_t = \sigma {\Delta t}^H
//'     }
//'   Where \eqn{\sigma} is the daily return volatility.
//' 
//'   The function \code{calc_var_ag()} can therefore be used to calculate the
//'   \emph{Hurst exponent} from the volatility ratio.
//'
//' @examples
//' \dontrun{
//' # Calculate the log prices
//' prices <- na.omit(rutils::etfenv$prices[, c("XLP", "VTI")])
//' prices <- log(prices)
//' # Calculate the daily variance of percentage returns
//' calc_var_ag(prices, step=1)
//' # Calculate the daily variance using R
//' sapply(rutils::diffit(prices), var)
//' # Calculate the variance of returns aggregated over 21 days
//' calc_var_ag(prices, step=21)
//' # The variance over 21 days is approximately 21 times the daily variance
//' 21*calc_var_ag(prices, step=1)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_var_ag(const arma::mat& tseries, 
                      arma::uword step = 1) {
  
  if (step == 1)
    // Calculate the variance without aggregations.
    return arma::var(diffit(tseries, 1, false));
  else {
    // Calculate the number of extra periods that don't fit over nrows.
    arma::uword nrows = tseries.n_rows;
    arma::uword remainder = nrows % step;
    // Allocate aggregations, end points, and variance.
    arma::mat aggs;
    arma::uvec endp;
    // The number of rows is (remainder+1) so that it works for remainder=0
    arma::mat vars(remainder+1, tseries.n_cols);
    // Perform loop over the stubs
    for (arma::uword stub = 0; stub <= remainder; stub++) {
      endp = calc_endpoints(tseries.n_rows, step, stub);
      // endp = arma::regspace<uvec>(stub, step, nrows + step);
      // endp = endp.elem(find(endp < nrows));
      aggs = tseries.rows(endp);
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
//' @param \code{method} A \emph{character} string representing the price range
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
//' @param \code{close_lag} A \emph{vector} with the lagged \emph{close} prices
//'   of the \emph{OHLC time series}.  This is an optional argument. (The
//'   default is \code{close_lag = 0}).
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
//'   the variance from minutely bar data, because dividing returns by the
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
//'   The optional argument \code{close_lag} are the lagged \emph{close} prices
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
//' close_lag <- HighFreq::lagit(ohlc[, 4])
//' all.equal(HighFreq::calc_var_ohlc(ohlc), 
//'   HighFreq::calc_var_ohlc(ohlc, close_lag=close_lag))
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
                     arma::colvec close_lag = 0, 
                     bool scale = true, 
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
  arma::colvec closep = ohlc.col(3);
  arma::colvec open_close(closep.n_rows);
  if (close_lag.n_rows == 1) {
    open_close = arma::join_cols(closep.subvec(0, 0), closep.subvec(0, closep.n_elem-2));
    open_close = (ohlc.col(0) - open_close)/index;
  } else {
    open_close = (ohlc.col(0) - close_lag)/index;
  }  // end if
  arma::colvec close_open = (closep - ohlc.col(0))/index;
  arma::colvec close_high = (closep - ohlc.col(1))/index;
  arma::colvec close_low = (closep - ohlc.col(2))/index;
  arma::colvec high_low = (ohlc.col(1) - ohlc.col(2))/index;
  arma::colvec high_open = (ohlc.col(1) - ohlc.col(0))/index;
  arma::colvec low_open = (ohlc.col(2) - ohlc.col(0))/index;
  
  if (method == "close") {
    // cout << "Calc method is Close" << endl;
    return arma::var(arma::diff(closep));
  } else if (method == "rogers_satchell") {
    // cout << "Calc method is Rogers-Satchell" << endl;
    return -(arma::dot(close_high, high_open) +
             arma::dot(close_low, low_open))/nrows;
  } else if (method == "garman_klass") {
    // cout << "Calc method is Garman-Klass" << endl;
    return (0.5*arma::dot(high_low, high_low) -
            (2*log(2)-1)*arma::dot(close_open, close_open))/nrows;
  } else if (method == "garman_klass_yz") {
    // cout << "Calc method is Garman-Klass-YZ" << endl;
    return (0.5*arma::dot(high_low, high_low) -
            (2*log(2)-1)*arma::dot(close_open, close_open))/nrows + 
            arma::var(open_close);
  } else if (method == "yang_zhang") {
    // cout << "Calc method is Yang-Zhang" << endl;
    return arma::var(open_close) + coeff*arma::var(close_open) +
      (coeff-1)*(arma::dot(close_high, high_open) + 
      arma::dot(close_low, low_open))/nrows;
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
//' @param \code{step} The number of periods in each interval between
//'   neighboring end points.
//' 
//' @param \code{method} A \emph{character} string representing the price range
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
//' @param \code{close_lag} A \emph{vector} with the lagged \emph{close} prices
//'   of the \emph{OHLC time series}.  This is an optional argument. (The
//'   default is \code{close_lag = 0}).
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
//'   \emph{OHLC} prices aggregated over end points.
//'
//'   It first calculates the end points spaced apart by the number of periods
//'   equal to the argument \code{step}.  Then it aggregates the \emph{OHLC}
//'   prices to the end points. Finally it calculates the variance of the
//'   aggregated \emph{OHLC} prices.
//'
//'   If there are extra periods that don't fit over the length of \code{ohlc},
//'   then \code{calc_var_ohlc_ag()} loops over all possible stub intervals,
//'   it calculates all the corresponding variance values, and it averages
//'   them.
//'
//'   For example, if the number of rows of \code{ohlc} is equal to \code{20},
//'   and \code{step=3} then \code{6} end points fit over the length of
//'   \code{ohlc}, and there are \code{2} extra periods that must fit into
//'   stubs, either at the beginning or at the end (or both).
//' 
//'   The aggregated volatility \eqn{\sigma_t} scales (increases) with the
//'   length of the aggregation interval \eqn{\Delta t} raised to the power of
//'   the \emph{Hurst exponent} \eqn{H}:
//'     \deqn{
//'       \sigma_t = \sigma {\Delta t}^H
//'     }
//'   Where \eqn{\sigma} is the daily return volatility.
//' 
//'   The function \code{calc_var_ohlc_ag()} can therefore be used to calculate the
//'   \emph{Hurst exponent} from the volatility ratio.
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
                        arma::uword step = 1, 
                        std::string method = "yang_zhang", 
                        arma::colvec close_lag = 0, 
                        bool scale = true, 
                        arma::colvec index = 0) {
  
  if (step == 1)
    // Calculate the variance without aggregations.
    return calc_var_ohlc(ohlc, method, close_lag, scale, index);
  else {
    // Calculate the number of extra periods that don't fit over nrows.
    arma::uword nrows = ohlc.n_rows;
    arma::uword remainder = nrows % step;
    // Allocate aggregations, end points, and variance.
    arma::mat aggs;
    arma::uvec endp;
    arma::mat vars(remainder, 1);
    // Perform loop over the stubs
    for (arma::uword stub = 0; stub < remainder; stub++) {
      endp = calc_endpoints(nrows, step, stub);
      // endp = arma::regspace<uvec>(stub, step, nrows + step);
      // endp = endp.elem(find(endp < nrows));
      // roll_ohlc
      aggs = roll_ohlc(ohlc, endp);
      vars.row(stub) = calc_var_ohlc(aggs, method, close_lag, scale, index);
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
//' @param \code{method} A \emph{string} specifying the type of the skewness
//'   model (the default is \code{method = "moment"} - see Details).
//'
//' @param \code{conf_lev} The confidence level for calculating the
//'   quantiles (the default is \code{conf_lev = 0.75}).
//'
//' @return A single-row matrix with the skewness of the columns of
//'   \code{tseries}.
//'
//' @details
//'   The function \code{calc_skew()} calculates the skewness of the columns of
//'   a \emph{time series} or a \emph{matrix} of data using \code{RcppArmadillo}
//'   \code{C++} code.
//'
//'   If \code{method = "moment"} (the default) then \code{calc_skew()}
//'   calculates the skewness as the third moment of the data.
//'
//'   If \code{method = "quantile"} then it calculates the skewness
//'   \eqn{\varsigma} from the differences between the quantiles of the data as
//'   follows:
//'   \deqn{
//'     \varsigma = \frac{q_{\alpha} + q_{1-\alpha} - 2*q_{0.5}}{q_{\alpha} - q_{1-\alpha}}
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
//' returns <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the moment skewness
//' HighFreq::calc_skew(returns)
//' # Calculate the moment skewness in R
//' calc_skewr <- function(x) {
//'   x <- (x-mean(x))
//'   sum(x^3)/var(x)^1.5/NROW(x)
//' }  # end calc_skewr
//' all.equal(HighFreq::calc_skew(returns), 
//'   calc_skewr(returns), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(returns),
//'   Rcode=calc_skewr(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the quantile skewness
//' HighFreq::calc_skew(returns, method="quantile", conf_lev=0.9)
//' # Calculate the quantile skewness in R
//' calc_skewq <- function(x, a = 0.75) {
//'   	quantiles <- quantile(x, c(1-a, 0.5, a), type=5)
//'   	(quantiles[3] + quantiles[1] - 2*quantiles[2])/(quantiles[3] - quantiles[1])
//' }  # end calc_skewq
//' all.equal(drop(HighFreq::calc_skew(returns, method="quantile", conf_lev=0.9)), 
//'   calc_skewq(returns, a=0.9), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(returns, method="quantile"),
//'   Rcode=calc_skewq(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the nonparametric skewness
//' HighFreq::calc_skew(returns, method="nonparametric")
//' # Compare HighFreq::calc_skew() with R nonparametric skewness
//' all.equal(drop(HighFreq::calc_skew(returns, method="nonparametric")), 
//'   (mean(returns)-median(returns))/sd(returns), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(returns, method="nonparametric"),
//'   Rcode=(mean(returns)-median(returns))/sd(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_skew(const arma::mat& tseries,
                    std::string method = "moment", 
                    double conf_lev = 0.75) {
  
  // Return zeros if not enough data
  if (tseries.n_rows < 3) {
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end if
  
  // Switch for the different methods of skew
  switch(calc_method(method)) {
  case methodenum::moment: {  // moment
    double nrows = tseries.n_rows;
    double ncols = tseries.n_cols;
    arma::mat means = arma::mean(tseries);
    arma::mat vars = arma::var(tseries);
    arma::mat skewness(1, ncols);
    // De-mean the columns of tseries
    // tseries.each_row() -= means;
    // return arma::sum(arma::pow(tseries, 3))/arma::pow(vars, 1.5)/nrows;
    for (arma::uword it = 0; it < ncols; it++) {
      skewness.col(it) = arma::sum(arma::pow(tseries.col(it) - arma::as_scalar(means.col(it)), 3))/arma::pow(vars.col(it), 1.5)/nrows;
    }  // end for
    return skewness;
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-conf_lev, 0.5, conf_lev};
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
//' @param \code{method} A \emph{string} specifying the type of the kurtosis
//'   model (the default is \code{method = "moment"} - see Details).
//'
//' @param \code{conf_lev} The confidence level for calculating the
//'   quantiles (the default is \code{conf_lev = 0.75}).
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
//'   But it doesn't de-mean the columns of \code{tseries} because that requires
//'   copying the matrix \code{tseries}, so it's time-consuming.
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
//' returns <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the moment kurtosis
//' HighFreq::calc_kurtosis(returns)
//' # Calculate the moment kurtosis in R
//' calc_kurtr <- function(x) {
//'   x <- (x-mean(x))
//'   sum(x^4)/var(x)^2/NROW(x)
//' }  # end calc_kurtr
//' all.equal(HighFreq::calc_kurtosis(returns), 
//'   calc_kurtr(returns), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_kurtosis(returns),
//'   Rcode=calc_kurtr(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the quantile kurtosis
//' HighFreq::calc_kurtosis(returns, method="quantile", conf_lev=0.9)
//' # Calculate the quantile kurtosis in R
//' calc_kurtq <- function(x, a=0.9) {
//'   	quantiles <- quantile(x, c(1-a, 0.25, 0.75, a), type=5)
//'   	(quantiles[4] - quantiles[1])/(quantiles[3] - quantiles[2])
//' }  # end calc_kurtq
//' all.equal(drop(HighFreq::calc_kurtosis(returns, method="quantile", conf_lev=0.9)), 
//'   calc_kurtq(returns, a=0.9), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_kurtosis(returns, method="quantile"),
//'   Rcode=calc_kurtq(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the nonparametric kurtosis
//' HighFreq::calc_kurtosis(returns, method="nonparametric")
//' # Compare HighFreq::calc_kurtosis() with R nonparametric kurtosis
//' all.equal(drop(HighFreq::calc_kurtosis(returns, method="nonparametric")), 
//'   (mean(returns)-median(returns))/sd(returns), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_kurtosis(returns, method="nonparametric"),
//'   Rcode=(mean(returns)-median(returns))/sd(returns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_kurtosis(const arma::mat& tseries,
                        std::string method = "moment", 
                        double conf_lev = 0.75) {
  
  // Return zeros if not enough data
  if (tseries.n_rows < 3) {
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end if
  
  // Switch for the different methods of kurtosis
  switch(calc_method(method)) {
  case methodenum::moment: {  // Fourth moment
    double nrows = tseries.n_rows;
    double ncols = tseries.n_cols;
    arma::mat means = arma::mean(tseries);
    arma::mat vars = arma::var(tseries);
    arma::mat kurtosis(1, ncols);
    // Don't de-mean the columns of tseries because that requires copying the matrix of data, so it's time-consuming
    // Loop over columns of tseries
    for (arma::uword it = 0; it < ncols; it++) {
      kurtosis.col(it) = arma::sum(arma::pow(tseries.col(it) - arma::as_scalar(means.col(it)), 4))/arma::pow(vars.col(it), 2)/nrows;
    }  // end for
    // tseries.each_row() -= means;
    // return arma::sum(arma::pow(tseries, 4))/arma::pow(vars, 2)/nrows;
    return kurtosis;
  }  // end moment
  case methodenum::quantile: {  // quantile
    arma::vec levels = {1-conf_lev, 0.25, 0.75, conf_lev};
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
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of prices.
//'
//' @param \code{step} The number of periods in each interval between
//'   neighboring end points.
//' 
//' @return The Hurst exponent calculated from the variance of aggregated
//'   returns.
//'
//' @details
//'   The function \code{calc_hurst()} calculates the Hurst exponent from the
//'   ratios of the volatilities of aggregated returns.
//'
//'   The aggregated volatility \eqn{\sigma_t} scales (increases) with the
//'   length of the aggregation interval \eqn{\Delta t} raised to the power of
//'   the \emph{Hurst exponent} \eqn{H}:
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
//'   The function \code{calc_hurst()} calls the function \code{calc_var_ag()}
//'   to calculate the aggregated volatility \eqn{\sigma_t}.
//' 
//' @examples
//' \dontrun{
//' # Calculate the log prices
//' prices <- na.omit(rutils::etfenv$prices[, c("XLP", "VTI")])
//' prices <- log(prices)
//' # Calculate the Hurst exponent from 21 day aggregations
//' calc_hurst(prices, step=21)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_hurst(const arma::mat& tseries, 
                     arma::uword step = 1) {
  
  return 0.5*arma::log(calc_var_ag(tseries, step)/calc_var_ag(tseries, 1))/log(step);
  
}  // end calc_hurst




////////////////////////////////////////////////////////////
//' Calculate the Hurst exponent from the volatility ratio of aggregated
//' \emph{OHLC} prices.
//'
//' @param \code{ohlc} A \emph{time series} or a \emph{matrix} of \emph{OHLC}
//'   prices.
//'
//' @param \code{step} The number of periods in each interval between
//'   neighboring end points.
//' 
//' @param \code{method} A \emph{character} string representing the price range
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
//' @param \code{close_lag} A \emph{vector} with the lagged \emph{close} prices
//'   of the \emph{OHLC time series}.  This is an optional argument. (The
//'   default is \code{close_lag = 0}).
//'   
//' @param \code{scale} \emph{Boolean} argument: Should the returns be divided
//'   by the time index, the number of seconds in each period? (The default is
//'   \code{scale = TRUE}).
//'
//' @param \code{index} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument (the default is \code{index = 0}).
//'   
//' @return The Hurst exponent calculated from the variance ratio of aggregated
//' \emph{OHLC} prices.
//'
//' @details
//' The function \code{calc_hurst_ohlc()} calculates the Hurst exponent from the
//' ratios of the volatilities of aggregated \emph{OHLC} prices.
//'
//'   The aggregated volatility \eqn{\sigma_t} scales (increases) with the
//'   length of the aggregation interval \eqn{\Delta t} raised to the power of
//'   the \emph{Hurst exponent} \eqn{H}:
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
//'   \code{calc_var_ohlc_ag()} to calculate the aggregated volatility
//'   \eqn{\sigma_t}.
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
                       arma::uword step = 1, 
                       std::string method = "yang_zhang", 
                       arma::colvec close_lag = 0, 
                       bool scale = true, 
                       arma::colvec index = 0) {
  
  return 0.5*log(calc_var_ohlc_ag(ohlc, step, method, close_lag, scale, index)/
                 calc_var_ohlc_ag(ohlc, 1, method, close_lag, scale, index))/log(step);
  
}  // end calc_hurst_ohlc




////////////////////////////////////////////////////////////
//' Perform multivariate linear regression using least squares and return a
//' named list of regression coefficients, their t-values, and p-values.
//' 
//' @param \code{response} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{predictor} A \emph{time series} or a \emph{matrix} of predictor
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
//' returns <- na.omit(rutils::etfenv$returns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' response <- returns[, 1]
//' # Predictor matrix equals VTI and IEF returns
//' predictor <- returns[, -1]
//' # Perform multivariate regression using lm()
//' reg_model <- lm(response ~ predictor)
//' sum_mary <- summary(reg_model)
//' # Perform multivariate regression using calc_lm()
//' reg_arma <- HighFreq::calc_lm(response=response, predictor=predictor)
//' # Compare the outputs of both functions
//' all.equal(reg_arma$coefficients[, "coeff"], unname(coef(reg_model)))
//' all.equal(unname(reg_arma$coefficients), unname(sum_mary$coefficients))
//' all.equal(unname(reg_arma$stats), c(sum_mary$r.squared, unname(sum_mary$fstatistic[1])))
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_lm(response=response, predictor=predictor),
//'   Rcode=lm(response ~ predictor),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List calc_lm(const arma::vec& response, const arma::mat& predictor) {
  
  // Add column for intercept to predictor matrix
  arma::uword nrows = predictor.n_rows;
  arma::mat predictori = join_rows(ones(nrows), predictor);
  arma::uword ncols = predictori.n_cols;
  arma::uword deg_free = (nrows - ncols);
  
  // Calculate alpha and beta coefficients for the model response ~ predictor
  arma::colvec coeff = arma::solve(predictori, response);
  // Calculate residuals
  arma::colvec residuals = response - predictori*coeff;
  
  // Calculate TSS, RSS, and ESS
  double tot_sumsq = (nrows-1)*arma::var(response);
  double res_sumsq = arma::dot(residuals, residuals);
  double exp_sumsq = tot_sumsq - res_sumsq;
  
  // Calculate R-squared and F-statistic
  double rsquared = exp_sumsq/tot_sumsq;
  double fstat = (exp_sumsq*deg_free)/(res_sumsq*(ncols-1));
  // arma::rowvec stats=join_horiz(rsquared, fstat);
  Rcpp::NumericVector stats(2);
  stats(0) = rsquared;
  stats(1) = fstat;
  stats.attr("names") = Rcpp::CharacterVector::create("R-squared", "F-statistic");
  
  // Calculate standard errors of beta coefficients
  arma::colvec stderr = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(predictori)*predictori)));
  // Calculate t-values and p-values of beta coefficients
  arma::colvec tvals = coeff/stderr;
  arma::colvec pvals = 2*Rcpp::pt(-abs(wrap(tvals)), deg_free);
  Rcpp::NumericMatrix coeffmat = Rcpp::wrap(join_rows(join_rows(join_rows(coeff, stderr), tvals), pvals));
  Rcpp::colnames(coeffmat) = Rcpp::CharacterVector::create("coeff", "stderr", "tvals", "pvals");
  
  return Rcpp::List::create(Named("coefficients") = coeffmat,
                            // Named("residuals") = residuals,
                            Named("zscore") = residuals(nrows-1)/arma::stddev(residuals),
                            Named("stats") = stats);
  
}  // end calc_lm



////////////////////////////////////////////////////////////
//' Perform multivariate regression using different methods, and return a vector
//' of regression coefficients, their t-values, and the last residual z-score.
//' 
//' @param \code{response} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{predictor} A \emph{time series} or a \emph{matrix} of predictor
//'   data.
//' 
//' @param \code{intercept} A \emph{Boolean} specifying whether an intercept
//'   term should be added to the predictor (the default is \code{intercept =
//'   FALSE}).
//'
//' @param \code{method} A \emph{string} specifying the type of the regression
//'   model the default is \code{method = "least_squares"} - see Details).
//'   
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small singular values in order to regularize the inverse of the
//'   \code{predictor} matrix (the default is \code{1e-5}).
//'   
//' @param \code{eigen_max} An \emph{integer} equal to the number of singular
//'   values used for calculating the shrinkage inverse of the \code{predictor}
//'   matrix (the default is \code{0} - equivalent to \code{eigen_max} equal to
//'   the number of columns of \code{predictor}).
//'   
//' @param \code{conf_lev} The confidence level for calculating the
//'   quantiles (the default is \code{conf_lev = 0.75}).
//'
//' @param \code{alpha} The shrinkage intensity between \code{0} and \code{1}.
//'   (the default is \code{0}).
//' 
//' @return A single-row matrix with
//' A vector with the regression coefficients, their t-values, and the
//'   last residual z-score.
//'
//' @details
//'   The function \code{calc_reg()} performs multivariate regression using
//'   different methods, and returns a vector of regression coefficients, their
//'   t-values, and the last residual z-score.
//'
//'   If \code{method = "least_squares"} (the default) then it performs the
//'   standard least squares regression, the same as the function
//'   \code{calc_lm()}, and the function \code{lm()} from the \code{R} package
//'   \emph{stats}.
//'   But it uses \code{RcppArmadillo} \code{C++} code so it's several times
//'   faster than \code{lm()}.
//'
//'   If \code{method = "regular"} then it performs shrinkage regression.  It
//'   calculates the shrinkage inverse of the \code{predictor} matrix from its
//'   singular value decomposition.  It applies dimension regularization by
//'   selecting only the largest singular values equal in number to
//'   \code{eigen_max}.
//'   
//'   If \code{method = "quantile"} then it performs quantile regression (not
//'   implemented yet).
//' 
//'   If \code{intercept = TRUE} then an extra intercept column (unit column) is
//'   added to the predictor matrix (the default is \code{intercept = FALSE}).
//'   
//'   The length of the return vector depends on the number of columns of the
//'   \code{predictor} matrix (including the intercept column, if it's added).
//'   The length of the return vector is equal to the number of regression
//'   coefficients, plus their t-values, plus the z-score.
//'   The number of regression coefficients is equal to the number of columns of
//'   the \code{predictor} matrix (including the intercept column, if it's
//'   added).
//'   The number of t-values is equal to the number of coefficients.
//' 
//'   For example, if the number of columns of the \code{predictor} matrix is
//'   equal to \code{n}, and if \code{intercept = TRUE}, then \code{calc_reg()}
//'   returns a vector with \code{2n+3} elements: \code{n+1} regression
//'   coefficients (including the intercept coefficient), \code{n+1}
//'   corresponding t-values, and \code{1} z-score.
//'
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' returns <- na.omit(rutils::etfenv$returns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' response <- returns[, 1]
//' # Predictor matrix equals VTI and IEF returns
//' predictor <- returns[, -1]
//' # Perform multivariate regression using lm()
//' reg_model <- lm(response ~ predictor)
//' sum_mary <- summary(reg_model)
//' coeff <- sum_mary$coefficients
//' # Perform multivariate regression using calc_reg()
//' reg_arma <- drop(HighFreq::calc_reg(response=response, predictor=predictor))
//' # Compare the outputs of both functions
//' all.equal(reg_arma[1:(2*(1+NCOL(predictor)))], 
//'   c(coeff[, "Estimate"], coeff[, "t value"]), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_reg(response=response, predictor=predictor),
//'   Rcode=lm(response ~ predictor),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_reg(const arma::mat& response, 
                   const arma::mat& predictor,
                   bool intercept = false,
                   std::string method = "least_squares",
                   double eigen_thresh = 1e-5,
                   arma::uword eigen_max = 0,
                   double conf_lev = 0.1,
                   double alpha = 0.0) {
  
  // Add column for intercept to predictor matrix
  arma::uword nrows = predictor.n_rows;
  arma::mat predictori = predictor;
  if (intercept)
    predictori = join_rows(ones(nrows), predictor);

  arma::uword ncols = predictori.n_cols;
  arma::uword deg_free = (nrows - ncols);
  arma::vec coeff;
  arma::vec tvals;
  arma::mat reg_data = arma::zeros<mat>(2*ncols+1, 1);
  
  // Switch for the different methods for weights
  switch(calc_method(method)) {
  case methodenum::least_squares: {
    // Calculate regression coefficients for the model response ~ predictor
    coeff = arma::solve(predictori, response);
    break;
  }  // end least_squares
  case methodenum::regular: {
    // Calculate shrinkage regression coefficients
    coeff = calc_inv(predictori, eigen_thresh, eigen_max)*response;
    break;
  }  // end regular
  case methodenum::quantile: {
    // Not implemented yet
    break;
  }  // end quantile
  default : {
    cout << "Warning: Invalid method parameter: " << method << endl;
    return reg_data;
  }  // end default
  }  // end switch
  
  // Calculate residuals
  arma::mat residuals = response - predictori*coeff;
  
  // Calculate TSS, RSS, and ESS
  // double tot_sumsq = (nrows-1)*arma::var(response);
  double res_sumsq = arma::dot(residuals, residuals);
  // double exp_sumsq = tot_sumsq - res_sumsq;
  
  // Calculate standard errors of the beta coefficients
  arma::mat stderr = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(predictori)*predictori)));
  // Calculate t-values of the beta coefficients
  tvals = coeff/stderr;
  
  // Calculate z-score
  mat zscore = residuals(nrows-1, 0)/arma::stddev(residuals);
  
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
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endp} An \emph{integer} vector of end points (the default is
//'   \code{endp = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{look_back} The number of end points in the look-back interval
//'   (the default is \code{look_back = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @param \code{method} A \emph{character} string representing the type of mean 
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
//'   interval equal to \code{look_back} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_mean()}, which
//'   calculates the mean (location).
//'   See the function \code{calc_mean()} for a description of the mean methods.
//'   
//'   If the arguments \code{endp} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{tseries}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//' 
//'   For example, the rolling mean at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{look_back = 3}.
//'
//'   The function \code{roll_mean()} with the parameter \code{step = 1}
//'   performs the same calculation as the function \code{roll_mean()} from
//'   package
//'   \href{https://cran.r-project.org/web/packages/RcppRoll/index.html}{RcppRoll},
//'   but it's several times faster because it uses \code{RcppArmadillo}
//'   \code{C++} code.
//'
//'   The function \code{roll_mean()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, so it's many times faster than the equivalent \code{R}
//'   code.
//'
//'   If only a simple rolling mean is required (not the median) then other
//'   functions like \code{roll_sum()} or \code{roll_vec()} may be even faster.
//'
//' @examples
//' \dontrun{
//' # Define time series of returns using package rutils
//' returns <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the rolling means at 25 day end points, with a 75 day look-back
//' means <- HighFreq::roll_mean(returns, step=25, look_back=3)
//' # Compare the mean estimates over 11-period look-back intervals
//' all.equal(HighFreq::roll_mean(returns, look_back=11)[-(1:10), ], 
//'   drop(RcppRoll::roll_mean(returns, n=11)), check.attributes=FALSE)
//' # Define end points and start points
//' endp <- HighFreq::calc_endpoints(NROW(returns), step=25)
//' startp <- HighFreq::calc_startpoints(endp, look_back=3)
//' # Calculate the rolling means using RcppArmadillo
//' means <- HighFreq::roll_mean(returns, startp=startp, endp=endp)
//' # Calculate the rolling medians using RcppArmadillo
//' medianscpp <- HighFreq::roll_mean(returns, startp=startp, endp=endp, method="nonparametric")
//' # Calculate the rolling medians using R
//' medians = sapply(1:NROW(endp), function(i) {
//'   median(returns[startp[i]:endp[i] + 1])
//' })  # end sapply
//' all.equal(medians, drop(medianscpp))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_mean(returns, startp=startp, endp=endp, method="nonparametric"),
//'   Rcode=sapply(1:NROW(endp), function(i) {median(returns[startp[i]:endp[i] + 1])}),
//'   times=10))[, c(1, 4, 5)]
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_mean(const arma::mat& tseries, 
                    arma::uvec startp = 0, 
                    arma::uvec endp = 0, 
                    arma::uword step = 1, 
                    arma::uword look_back = 1, 
                    arma::uword stub = 0,
                    std::string method = "moment", 
                    double conf_lev = 0.75) {
  
  // Allocate end points
  arma::uword nrows = tseries.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    startpts = calc_startpoints(endpts, look_back);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  
  // Allocate mean matrix
  arma::uword numpts = endpts.n_elem;
  arma::mat means = arma::zeros<mat>(numpts, tseries.n_cols);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate means
    if (endpts(ep) > startpts(ep)) {
      means.row(ep) = calc_mean(tseries.rows(startpts(ep), endpts(ep)), method, conf_lev);
    }  // end if
  }  // end for
  
  return means;
  
}  // end roll_mean




////////////////////////////////////////////////////////////
//' Calculate a \emph{vector} of variance estimates over a rolling look-back
//' interval for a single-column \emph{time series} or a single-column
//' \emph{matrix}, using \code{RcppArmadillo}.
//'
//' @param \code{tseries} A single-column \emph{time series} or a single-column
//'   \emph{matrix}.
//' 
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of \emph{vector} elements used for calculating a single variance
//'   estimate (the default is \code{look_back = 1}).
//'
//' @return A single-column \emph{matrix} with the same number of elements as
//'   the input argument \code{tseries}.
//'
//' @details
//'   The function \code{roll_varvec()} calculates a \emph{vector} of variance
//'   estimates over a rolling look-back interval for a single-column \emph{time
//'   series} or a single-column \emph{matrix}, using \code{RcppArmadillo}
//'   \code{C++} code.
//'   
//'   The function \code{roll_varvec()} uses an expanding look-back interval in
//'   the initial warmup period, to calculate the same number of elements as the
//'   input argument \code{tseries}.
//'
//'   The function \code{roll_varvec()} performs the same calculation as the
//'   function \code{roll_var()} from package
//'   \href{https://cran.r-project.org/web/packages/RcppRoll/index.html}{RcppRoll},
//'   but it's several times faster because it uses \code{RcppArmadillo}
//'   \code{C++} code.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random returns
//' returns <- rnorm(1e6)
//' # Compare the variance estimates over 11-period look-back intervals
//' all.equal(drop(HighFreq::roll_varvec(returns, look_back=11))[-(1:10)], 
//'   RcppRoll::roll_var(returns, n=11))
//' # Compare the speed of RcppArmadillo with RcppRoll
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_varvec(returns, look_back=11),
//'   RcppRoll=RcppRoll::roll_var(returns, n=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::vec roll_varvec(const arma::vec& tseries, arma::uword look_back = 1) {
  
  arma::uword length = tseries.n_elem;
  arma::vec varvec = arma::zeros<vec>(length);
  
  // Warmup period
  for (arma::uword it = 1; it < look_back; it++) {
    varvec(it) = arma::var(tseries.subvec(0, it));
  }  // end for
  
  // Remaining period
  for (arma::uword it = look_back; it < length; it++) {
    varvec(it) = arma::var(tseries.subvec(it-look_back+1, it));
  }  // end for
  
  return varvec;
  
}  // end roll_varvec




////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of dispersion (variance) estimates over a rolling
//' look-back interval attached at the end points of a \emph{time series} or a
//' \emph{matrix}.
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endp} An \emph{integer} vector of end points (the default is
//'   \code{endp = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{look_back} The number of end points in the look-back interval
//'   (the default is \code{look_back = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @param \code{method} A \emph{character} string representing the type of the
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
//'   interval equal to \code{look_back} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_var()}, which
//'   calculates the dispersion.
//'   See the function \code{calc_var()} for a description of the dispersion
//'   methods.
//'   
//'   If the arguments \code{endp} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{tseries}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//' 
//'   For example, the rolling variance at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{look_back = 3}.
//'
//'   The function \code{roll_var()} with the parameter \code{step = 1}
//'   performs the same calculation as the function \code{roll_var()} from
//'   package
//'   \href{https://cran.r-project.org/web/packages/RcppRoll/index.html}{RcppRoll},
//'   but it's several times faster because it uses \code{RcppArmadillo}
//'   \code{C++} code.
//'
//'   The function \code{roll_var()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, so it's many times faster than the equivalent \code{R}
//'   code.
//'
//' @examples
//' \dontrun{
//' # Define time series of returns using package rutils
//' returns <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the rolling variance at 25 day end points, with a 75 day look-back
//' variance <- HighFreq::roll_var(returns, step=25, look_back=3)
//' # Compare the variance estimates over 11-period look-back intervals
//' all.equal(HighFreq::roll_var(returns, look_back=11)[-(1:10), ], 
//'   drop(RcppRoll::roll_var(returns, n=11)), check.attributes=FALSE)
//' # Compare the speed of HighFreq::roll_var() with RcppRoll::roll_var()
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_var(returns, look_back=11),
//'   RcppRoll=RcppRoll::roll_var(returns, n=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Compare the speed of HighFreq::roll_var() with TTR::runMAD()
//' summary(microbenchmark(
//'     Rcpp=HighFreq::roll_var(returns, look_back=11, method="quantile"),
//'     TTR=TTR::runMAD(returns, n = 11),
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_var(const arma::mat& tseries, 
                   arma::uvec startp = 0, 
                   arma::uvec endp = 0, 
                   arma::uword step = 1, 
                   arma::uword look_back = 1, 
                   arma::uword stub = 0,
                   std::string method = "moment", 
                   double conf_lev = 0.75) {
  
  // Allocate end points
  arma::uword nrows = tseries.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    startpts = calc_startpoints(endpts, look_back);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  
  // Allocate variance matrix
  arma::uword numpts = endpts.n_elem;
  arma::mat variance = arma::zeros<mat>(numpts, tseries.n_cols);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate variance
    if (endpts(ep) > startpts(ep)) {
      variance.row(ep) = calc_var(tseries.rows(startpts(ep), endpts(ep)), method, conf_lev);
    }  // end if
  }  // end for
  
  return variance;
  
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
//' @param \code{endp} An \emph{integer} vector of end points (the default is
//'   \code{endp = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{look_back} The number of end points in the look-back interval
//'   (the default is \code{look_back = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @param \code{method} A \emph{character} string representing the price range
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
//'   equal to \code{look_back} number of end points.
//'   In the initial warmup period, the variance is calculated over an expanding
//'   look-back interval.
//'   
//'   If the arguments \code{endp} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{ohlc}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//' 
//'   For example, the rolling variance at daily end points with an \code{11}
//'   day look-back, can be calculated using the parameters \code{step = 1} and
//'   \code{look_back = 1} (Assuming the \code{ohlc} data has daily
//'   frequency.)
//' 
//'   Similarly, the rolling variance at \code{25} day end points with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{look_back = 3} (because \code{3*25 = 75}).
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
//'   the variance from minutely bar data, because dividing returns by the
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
//'   \code{C++} code, so it's many times faster than the equivalent \code{R}
//'   code.
//'
//' @examples
//' \dontrun{
//' # Extract the log OHLC prices of SPY
//' ohlc <- log(HighFreq::SPY)
//' # Extract the time index of SPY prices
//' indeks <- c(1, diff(xts::.index(ohlc)))
//' # Rolling variance at minutely end points, with a 21 minute look-back
//' var_rolling <- HighFreq::roll_var_ohlc(ohlc, 
//'                               step=1, look_back=21, 
//'                               method="yang_zhang", 
//'                               index=indeks, scale=TRUE)
//' # Daily OHLC prices
//' ohlc <- rutils::etfenv$VTI
//' indeks <- c(1, diff(xts::.index(ohlc)))
//' # Rolling variance at 5 day end points, with a 20 day look-back (20=4*5)
//' var_rolling <- HighFreq::roll_var_ohlc(ohlc, 
//'                               step=5, look_back=4, 
//'                               method="yang_zhang", 
//'                               index=indeks, scale=TRUE)
//' # Same calculation in R
//' nrows <- NROW(ohlc)
//' close_lag = HighFreq::lagit(ohlc[, 4])
//' endp <- drop(HighFreq::calc_endpoints(nrows, 3)) + 1
//' startp <- drop(HighFreq::calc_startpoints(endp, 2))
//' n_pts <- NROW(endp)
//' var_rollingr <- sapply(2:n_pts, function(it) {
//'   rangev <- startp[it]:endp[it]
//'   sub_ohlc = ohlc[rangev, ]
//'   sub_close = close_lag[rangev]
//'   sub_index = indeks[rangev]
//'   HighFreq::calc_var_ohlc(sub_ohlc, close_lag=sub_close, scale=TRUE, index=sub_index)
//' })  # end sapply
//' var_rollingr <- c(0, var_rollingr)
//' all.equal(drop(var_rolling), var_rollingr)
//' }
//' @export
// [[Rcpp::export]]
arma::vec roll_var_ohlc(const arma::mat& ohlc, 
                        arma::uvec startp = 0, 
                        arma::uvec endp = 0, 
                        arma::uword step = 1, 
                        arma::uword look_back = 1, 
                        arma::uword stub = 0,
                        std::string method = "yang_zhang", 
                        bool scale = true, 
                        arma::colvec index = 0) {
  
  // Allocate end points
  arma::uword nrows = ohlc.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    startpts = calc_startpoints(endpts, look_back);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  // Allocate variance matrix
  arma::uword numpts = endpts.n_elem;
  arma::vec variance = arma::zeros<vec>(numpts);
  
  // Extract OHLC close prices
  arma::colvec closep = ohlc.col(3);
  arma::colvec close_lag = lagit(closep, 1, false);
  
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
      sub_close = close_lag.rows(startpts(ep), endpts(ep));
      sub_index = index.subvec(startpts(ep), endpts(ep));
      // Calculate variance
      variance.row(ep) = calc_var_ohlc(sub_ohlc, method, sub_close, scale, sub_index);
    }  // end if
  }  // end for
  
  // Old code below
  
  // Warmup period
  // for (arma::uword it = 1; it < look_back; it++) {
  //   arma::mat sub_ohlc = ohlc.rows(0, it);
  //   arma::colvec sub_close = close_lag.rows(0, it);
  //   arma::colvec sub_index = index.subvec(0, it);
  //   variance(it) = calc_var_ohlc(sub_ohlc, method, sub_close, scale, sub_index);
  // }  // end for
  
  // Remaining period
  // for (arma::uword it = look_back; it < nrows; it++) {
  //   arma::mat sub_ohlc = ohlc.rows(it-look_back+1, it);
  //   arma::colvec sub_close = close_lag.rows(it-look_back+1, it);
  //   arma::colvec sub_index = index.subvec(it-look_back+1, it);
  //   variance(it) = calc_var_ohlc(sub_ohlc, method, sub_close, scale, sub_index);
  // }  // end for
  
  return variance;
  
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
//' @param \code{endp} An \emph{integer} vector of end points (the default is 
//'   \code{endp = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{look_back} The number of end points in the look-back interval
//'   (the default is \code{look_back = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @param \code{method} A \emph{string} specifying the type of the skewness
//'   model (the default is \code{method = "moment"} - see Details).
//'
//' @param \code{conf_lev} The confidence level for calculating the
//'   quantiles (the default is \code{conf_lev = 0.75}).
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
//'   look-back interval equal to \code{look_back} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_skew()}, which
//'   calculates the skewness.
//'   See the function \code{calc_skew()} for a description of the skewness
//'   methods.
//'   
//'   If the arguments \code{endp} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{tseries}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//' 
//'   For example, the rolling skewness at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{look_back = 3}.
//'
//'   The function \code{roll_skew()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, so it's many times faster than the equivalent \code{R}
//'   code.
//'
//' @examples
//' \dontrun{
//' # Define time series of returns using package rutils
//' returns <- na.omit(rutils::etfenv$returns$VTI)
//' # Define end points and start points
//' endp <- 1 + HighFreq::calc_endpoints(NROW(returns), step=25)
//' startp <- HighFreq::calc_startpoints(endp, look_back=3)
//' # Calculate the rolling skewness at 25 day end points, with a 75 day look-back
//' skew_ness <- HighFreq::roll_skew(returns, step=25, look_back=3)
//' # Calculate the rolling skewness using R code
//' skew_r <- sapply(1:NROW(endp), function(it) {
//'   HighFreq::calc_skew(returns[startp[it]:endp[it], ])
//' })  # end sapply
//' # Compare the skewness estimates
//' all.equal(drop(skew_ness), skew_r, check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_skew(returns, step=25, look_back=3),
//'   Rcode=sapply(1:NROW(endp), function(it) {
//'     HighFreq::calc_skew(returns[startp[it]:endp[it], ])
//'   }),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_skew(const arma::mat& tseries, 
                    arma::uvec startp = 0, 
                    arma::uvec endp = 0, 
                    arma::uword step = 1, 
                    arma::uword look_back = 1, 
                    arma::uword stub = 0,
                    std::string method = "moment", 
                    double conf_lev = 0.75) {
  
  // Allocate end points
  arma::uword nrows = tseries.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    startpts = calc_startpoints(endpts, look_back);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  // Allocate skewness matrix
  arma::uword numpts = endpts.n_elem;
  arma::mat skew_ness = arma::zeros<mat>(numpts, tseries.n_cols);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate skewness
    if (endpts(ep) > startpts(ep)) {
      skew_ness.row(ep) = calc_skew(tseries.rows(startpts(ep), endpts(ep)), method, conf_lev);
    }  // end if
  }  // end for
  
  return skew_ness;
  
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
//' @param \code{endp} An \emph{integer} vector of end points (the default is 
//'   \code{endp = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{look_back} The number of end points in the look-back interval
//'   (the default is \code{look_back = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @param \code{method} A \emph{string} specifying the type of the kurtosis
//'   model (the default is \code{method = "moment"} - see Details).
//'
//' @param \code{conf_lev} The confidence level for calculating the
//'   quantiles (the default is \code{conf_lev = 0.75}).
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
//'   look-back interval equal to \code{look_back} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_kurtosis()},
//'   which calculates the kurtosis. See the function \code{calc_kurtosis()} for
//'   a description of the kurtosis methods.
//'   
//'   If the arguments \code{endp} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{tseries}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//' 
//'   For example, the rolling kurtosis at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{look_back = 3}.
//'
//'   The function \code{roll_kurtosis()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, so it's many times faster than the equivalent \code{R}
//'   code.
//'
//' @examples
//' \dontrun{
//' # Define time series of returns using package rutils
//' returns <- na.omit(rutils::etfenv$returns$VTI)
//' # Define end points and start points
//' endp <- 1 + HighFreq::calc_endpoints(NROW(returns), step=25)
//' startp <- HighFreq::calc_startpoints(endp, look_back=3)
//' # Calculate the rolling kurtosis at 25 day end points, with a 75 day look-back
//' kurto_sis <- HighFreq::roll_kurtosis(returns, step=25, look_back=3)
//' # Calculate the rolling kurtosis using R code
//' kurt_r <- sapply(1:NROW(endp), function(it) {
//'   HighFreq::calc_kurtosis(returns[startp[it]:endp[it], ])
//' })  # end sapply
//' # Compare the kurtosis estimates
//' all.equal(drop(kurto_sis), kurt_r, check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_kurtosis(returns, step=25, look_back=3),
//'   Rcode=sapply(1:NROW(endp), function(it) {
//'     HighFreq::calc_kurtosis(returns[startp[it]:endp[it], ])
//'   }),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_kurtosis(const arma::mat& tseries, 
                        arma::uvec startp = 0, 
                        arma::uvec endp = 0, 
                        arma::uword step = 1, 
                        arma::uword look_back = 1, 
                        arma::uword stub = 0,
                        std::string method = "moment", 
                        double conf_lev = 0.75) {
  
  // Allocate end points
  arma::uword nrows = tseries.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    startpts = calc_startpoints(endpts, look_back);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  // Allocate kurtosis matrix
  arma::uword numpts = endpts.n_elem;
  arma::mat kurto_sis = arma::zeros<mat>(numpts, tseries.n_cols);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate kurtosis
    if (endpts(ep) > startpts(ep)) {
      kurto_sis.row(ep) = calc_kurtosis(tseries.rows(startpts(ep), endpts(ep)), method, conf_lev);
    }  // end if
  }  // end for
  
  return kurto_sis;
  
}  // end roll_kurtosis



////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of regression coefficients, their t-values, and
//' z-scores, at the end points of the predictor matrix.
//' 
//' @param \code{response} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{predictor} A \emph{time series} or a \emph{matrix} of predictor
//'   data.
//'   
//' @param \code{intercept} A \emph{Boolean} specifying whether an intercept
//'   term should be added to the predictor (the default is \code{intercept =
//'   FALSE}).
//'
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endp} An \emph{integer} vector of end points (the default is 
//'   \code{endp = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{look_back} The number of end points in the look-back interval
//'   (the default is \code{look_back = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @param \code{intercept} A \emph{Boolean} specifying whether an intercept
//'   term should be added to the predictor (the default is \code{intercept =
//'   FALSE}).
//'
//' @param \code{method} A \emph{string} specifying the type of the regression
//'   model the default is \code{method = "least_squares"} - see Details).
//'   
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small singular values in order to regularize the inverse of the
//'   \code{predictor} matrix (the default is \code{1e-5}).
//'   
//' @param \code{eigen_max} An \emph{integer} equal to the number of singular
//'   values used for calculating the shrinkage inverse of the \code{predictor}
//'   matrix (the default is \code{0} - equivalent to \code{eigen_max} equal to
//'   the number of columns of \code{predictor}).
//'   
//' @param \code{conf_lev} The confidence level for calculating the
//'   quantiles (the default is \code{conf_lev = 0.75}).
//'
//' @param \code{alpha} The shrinkage intensity between \code{0} and \code{1}.
//'   (the default is \code{0}).
//' 
//' @return A \emph{matrix} with the regression coefficients, their t-values, and
//' z-scores, and with the 
//' same number of rows as \code{predictor}
//' a
//'   number of columns equal to \code{2n+3}, where \code{n} is the number of
//'   columns of \code{predictor}.
//'
//' @details
//'   The function \code{roll_reg()} calculates a \emph{matrix} of regression
//'   coefficients, their t-values, and z-scores at the end points of the predictor
//'   matrix.
//'   
//'   The function \code{roll_reg()} performs a loop over the end points, and at
//'   each end point it subsets the time series \code{predictor} over a look-back
//'   interval equal to \code{look_back} number of end points.
//'   
//'   If the arguments \code{endp} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{predictor}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//'   
//'   For example, the rolling regression at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{look_back = 3}.
//'
//'   It passes the subset time series to the function \code{calc_reg()}, which
//'   calculates the regression coefficients, their t-values, and the z-score.
//'   
//'   If \code{intercept = TRUE} then an extra intercept column (unit column) is
//'   added to the predictor matrix (the default is \code{intercept = FALSE}).
//'   
//'   The number of columns of the return matrix depends on the number of
//'   columns of the \code{predictor} matrix (including the intercept column, if
//'   it's added).
//'   The number of columns of the return matrix is equal to the number of
//'   regression coefficients, plus their t-values, plus the z-score column.
//'   The number of regression coefficients is equal to the number of columns of
//'   the \code{predictor} matrix (including the intercept column, if it's
//'   added).
//'   The number of t-values is equal to the number of coefficients.
//'   For example, if the number of columns of the \code{predictor} matrix is
//'   equal to \code{n}, and if \code{intercept = TRUE}, then \code{roll_reg()}
//'   returns a matrix with \code{2n+3} columns: \code{n+1} regression
//'   coefficients (including the intercept coefficient), \code{n+1}
//'   corresponding t-values, and \code{1} z-score column.
//' 
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' returns <- na.omit(rutils::etfenv$returns[, c("XLP", "VTI")])
//' # Define monthly end points and start points
//' endp <- xts::endpoints(returns, on="months")[-1]
//' look_back <- 12
//' startp <- c(rep(1, look_back), endp[1:(NROW(endp)-look_back)])
//' # Calculate rolling betas using RcppArmadillo
//' reg_stats <- HighFreq::roll_reg(response=returns[, 1], predictor=returns[, 2], endp=(endp-1), startp=(startp-1))
//' betas <- reg_stats[, 2]
//' # Calculate rolling betas in R
//' betas_r <- sapply(1:NROW(endp), FUN=function(ep) {
//'   datav <- returns[startp[ep]:endp[ep], ]
//'   drop(cov(datav[, 1], datav[, 2])/var(datav[, 2]))
//' })  # end sapply
//' # Compare the outputs of both functions
//' all.equal(betas, betas_r, check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_reg(const arma::mat& response, 
                   const arma::mat& predictor, 
                   bool intercept = false,
                   arma::uvec startp = 0, 
                   arma::uvec endp = 0, 
                   arma::uword step = 1, 
                   arma::uword look_back = 1, 
                   arma::uword stub = 0,
                   std::string method = "least_squares",
                   double eigen_thresh = 1e-5,
                   arma::uword eigen_max = 0,
                   double conf_lev = 0.1,
                   double alpha = 0.0) {
  
  // Allocate end points
  arma::uword nrows = predictor.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    startpts = calc_startpoints(endpts, look_back);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  
  // Allocate regression matrix
  arma::mat responsi;
  arma::mat predicti;
  arma::uword numpts = endpts.n_elem;
  arma::uword ncols = predictor.n_cols;
  if (intercept) ncols += 1;
  arma::mat reg_stats(numpts, (2*ncols + 1), fill::zeros);

  // Perform loop over the endpts
  for (arma::uword ep = 0; ep < numpts; ep++) {
    // Calculate regression coefficients
    if (endpts(ep) > startpts(ep)) {
      responsi = response.rows(startpts(ep), endpts(ep));
      predicti = predictor.rows(startpts(ep), endpts(ep));
      reg_stats.row(ep) = calc_reg(responsi, predicti, intercept, method, eigen_thresh, eigen_max, conf_lev, alpha);
    }  // end if
  }  // end for
  
  // Warmup period
  // reg_stats.rows(0, ncols+1) = zeros(ncols+2, (ncols + 1));
  // for (arma::uword it = (ncols+2); it < look_back; it++) {
  //   responsi = response.rows(0, it);
  //   predicti = predictor.rows(0, it);
  //   reg_data = calc_reg(responsi, predicti);
  //   reg_stats.row(it) = conv_to<rowvec>::from(reg_data);
  // }  // end for
  
  // Remaining periods
  // for (arma::uword it = look_back; it < nrows; it++) {
  //   responsi = response.rows(it-look_back+1, it);
  //   predicti = predictor.rows(it-look_back+1, it);
  //   reg_data = calc_reg(responsi, predicti, method, eigen_thresh, eigen_max, conf_lev, alpha);
  //   reg_stats.row(it) = conv_to<rowvec>::from(reg_data);
  // }  // end for
  
  return reg_stats;
  
}  // end roll_reg



////////////////////////////////////////////////////////////
//' Perform a rolling scaling (standardization) of the columns of a
//' \emph{matrix} of data using \code{RcppArmadillo}.
//' 
//' @param \code{matrix} A \emph{matrix} of data.
//' 
//' @param \code{look_back} The length of the look-back interval, equal to the number 
//'   of rows of data used in the scaling.
//'   
//' @param use_median A \emph{Boolean} argument: if \code{TRUE} then the 
//'   centrality (central tendency) is calculated as the \emph{median} and the 
//'   dispersion is calculated as the \emph{median absolute deviation}
//'   (\emph{MAD}).
//'   If \code{use_median} is \code{FALSE} then the centrality is calculated as 
//'   the \emph{mean} and the dispersion is calculated as the \emph{standard
//'   deviation} (the default is \code{use_median = FALSE})
//'
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{matrix}.
//'
//' @details
//'   The function \code{roll_scale()} performs a rolling scaling
//'   (standardization) of the columns of the \code{matrix} argument using
//'   \code{RcppArmadillo}.
//'   The function \code{roll_scale()} performs a loop over the rows of 
//'   \code{matrix}, subsets a number of previous (past) rows equal to 
//'   \code{look_back}, and scales the subset matrix.  It assigns the last row
//'   of the scaled subset \emph{matrix} to the return matrix.
//'   
//'   If the argument \code{use_median} is \code{FALSE} (the default), then it
//'   performs the same calculation as the function \code{roll::roll_scale()}.
//'   If the argument \code{use_median} is \code{TRUE}, then it calculates the
//'   centrality as the \emph{median} and the dispersion as the \emph{median
//'   absolute deviation} (\emph{MAD}).
//'   
//' @examples
//' \dontrun{
//' matrixv <- matrix(rnorm(20000), nc=2)
//' look_back <- 11
//' rolled_scaled <- roll::roll_scale(data=matrixv, width = look_back, min_obs=1)
//' rolled_scaled2 <- roll_scale(matrix=matrixv, look_back = look_back, use_median=FALSE)
//' all.equal(rolled_scaled[-1, ], rolled_scaled2[-1, ])
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_scale(const arma::mat& matrix, 
                     arma::uword look_back,
                     bool use_median=false) {
  
  arma::uword nrows = matrix.n_rows;
  arma::mat scaledmat(nrows, matrix.n_cols);
  arma::mat sub_mat;
  
  // Warmup period
  scaledmat.row(0) = matrix.row(0);
  for (arma::uword it = 1; it < look_back; it++) {
    sub_mat = matrix.rows(0, it);
    sub_mat = calc_scaled(sub_mat, use_median);
    scaledmat.row(it) = sub_mat.row(sub_mat.n_rows-1);
  }  // end for
  
  // Remaining periods
  for (arma::uword it = look_back; it < nrows; it++) {
    sub_mat = matrix.rows(it-look_back+1, it);
    sub_mat = calc_scaled(sub_mat, use_median);
    scaledmat.row(it) = sub_mat.row(sub_mat.n_rows-1);
  }  // end for
  
  return scaledmat;
}  // end roll_scale



////////////////////////////////////////////////////////////
//' Calculate a \emph{vector} of z-scores of the residuals of rolling
//' regressions at the end points of the predictor matrix.
//' 
//' @param \code{response} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{predictor} A \emph{time series} or a \emph{matrix} of predictor
//'   data.
//'   
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endp} An \emph{integer} vector of end points (the default is
//'   \code{endp = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{look_back} The number of end points in the look-back interval
//'   (the default is \code{look_back = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @return A column \emph{vector} of the same length as the number of rows of
//'   \code{predictor}.
//'
//' @details
//'   The function \code{roll_zscores()} calculates a \emph{vector} of z-scores
//'   of the residuals of rolling regressions at the end points of the
//'   \emph{time series} \code{predictor}.
//'   
//'   The function \code{roll_zscores()} performs a loop over the end points,
//'   and at each end point it subsets the time series \code{predictor} over a
//'   look-back interval equal to \code{look_back} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_lm()}, which
//'   calculates the regression data.
//'   
//'   If the arguments \code{endp} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{predictor}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//'   
//'   For example, the rolling variance at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{look_back = 3}.
//'
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' returns <- na.omit(rutils::etfenv$returns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' response <- returns[, 1]
//' # Predictor matrix equals VTI and IEF returns
//' predictor <- returns[, -1]
//' # Calculate Z-scores from rolling time series regression using RcppArmadillo
//' look_back <- 11
//' z_scores <- HighFreq::roll_zscores(response=response, predictor=predictor, look_back)
//' # Calculate z-scores in R from rolling multivariate regression using lm()
//' z_scoresr <- sapply(1:NROW(predictor), function(ro_w) {
//'   if (ro_w == 1) return(0)
//'   startpoint <- max(1, ro_w-look_back+1)
//'   responsi <- response[startpoint:ro_w]
//'   predicti <- predictor[startpoint:ro_w, ]
//'   reg_model <- lm(responsi ~ predicti)
//'   residuals <- reg_model$residuals
//'   residuals[NROW(residuals)]/sd(residuals)
//' })  # end sapply
//' # Compare the outputs of both functions
//' all.equal(z_scores[-(1:look_back)], z_scoresr[-(1:look_back)], 
//'   check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec roll_zscores(const arma::mat& response, 
                       const arma::mat& predictor, 
                       arma::uvec startp = 0, 
                       arma::uvec endp = 0, 
                       arma::uword step = 1, 
                       arma::uword look_back = 1,
                       arma::uword stub = 0) {
  
  // Allocate end points
  arma::uword nrows = predictor.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    startpts = calc_startpoints(endpts, look_back);
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
      responsi = response.rows(startpts(ep), endpts(ep));
      predicti = predictor.rows(startpts(ep), endpts(ep));
      zscores(ep) = calc_lm(responsi, predicti)["zscore"];
    }  // end if
  }  // end for
  
  // Old code below
  // Warmup period
  // for (arma::uword it = 1; it < look_back; it++) {
  //   responsi = response.rows(0, it);
  //   predicti = predictor.rows(0, it);
  //   zscores(it) = calc_lm(responsi, predicti)["zscore"];
  // }  // end for
  
  // Remaining periods
  // for (arma::uword it = look_back; it < nrows; it++) {
  //   responsi = response.rows(it-look_back+1, it);
  //   predicti = predictor.rows(it-look_back+1, it);
  //   zscores(it) = calc_lm(responsi, predicti)["zscore"];
  // }  // end for
  
  return zscores;
  
}  // end roll_zscores




////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of estimator values over a rolling look-back
//' interval attached at the end points of a \emph{time series} or a
//' \emph{matrix}.
//'
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'    
//' @param \code{fun} A \emph{string} specifying the estimator function (the
//'   default is \code{fun = "calc_var"}.)
//'
//' @param \code{startp} An \emph{integer} vector of start points (the default
//'   is \code{startp = 0}).
//' 
//' @param \code{endp} An \emph{integer} vector of end points (the default is 
//'   \code{endp = 0}).
//' 
//' @param \code{step} The number of time periods between the end points (the
//'   default is \code{step = 1}).
//'
//' @param \code{look_back} The number of end points in the look-back interval
//'   (the default is \code{look_back = 1}).
//'   
//' @param \code{stub} An \emph{integer} value equal to the first end point for
//'   calculating the end points (the default is \code{stub = 0}).
//' 
//' @param \code{method} A \emph{string} specifying the type of the model for the
//'   estimator (the default is \code{method = "moment"}.)
//'
//' @param \code{conf_lev} The confidence level for calculating the
//'   quantiles (the default is \code{conf_lev = 0.75}).
//'
//' @return A \emph{matrix} with the same number of columns as the input time
//'   series \code{tseries}, and the number of rows equal to the number of end
//'   points.
//'   
//' @details
//'   The function \code{roll_fun()} calculates a \emph{matrix} of estimator
//'   values, over rolling look-back intervals attached at the end points of the
//'   \emph{time series} \code{tseries}.
//'   
//'   The function \code{roll_fun()} performs a loop over the end points, and at
//'   each end point it subsets the time series \code{tseries} over a look-back
//'   interval equal to \code{look_back} number of end points.
//'   
//'   It passes the subset time series to the function specified by the argument
//'   \code{fun}, which calculates the statistic.
//'   See the functions \code{calc_*()} for a description of the different
//'   estimators.
//'   
//'   If the arguments \code{endp} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{tseries}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//' 
//'   For example, the rolling variance at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{look_back = 3}.
//'
//'   The function \code{roll_fun()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, so it's many times faster than the equivalent \code{R}
//'   code.
//'
//' @examples
//' \dontrun{
//' # Define time series of returns using package rutils
//' returns <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the rolling variance at 25 day end points, with a 75 day look-back
//' var_rollfun <- HighFreq::roll_fun(returns, fun="calc_var", step=25, look_back=3)
//' # Calculate the rolling variance using roll_var()
//' var_roll <- HighFreq::roll_var(returns, step=25, look_back=3)
//' # Compare the two methods
//' all.equal(var_rollfun, var_roll, check.attributes=FALSE)
//' # Define end points and start points
//' endp <- HighFreq::calc_endpoints(NROW(returns), step=25)
//' startp <- HighFreq::calc_startpoints(endp, look_back=3)
//' # Calculate the rolling variance using RcppArmadillo
//' var_rollfun <- HighFreq::roll_fun(returns, fun="calc_var", startp=startp, endp=endp)
//' # Calculate the rolling variance using R code
//' var_roll <- sapply(1:NROW(endp), function(it) {
//'   var(returns[startp[it]:endp[it]+1, ])
//' })  # end sapply
//' # Compare the two methods
//' all.equal(drop(var_rollfun), var_roll, check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_fun(returns, fun="calc_var", startp=startp, endp=endp),
//'   Rcode=sapply(1:NROW(endp), function(it) {
//'     var(returns[startp[it]:endp[it]+1, ])
//'   }),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_fun(const arma::mat& tseries, 
                   std::string fun = "calc_var", 
                   arma::uvec startp = 0, 
                   arma::uvec endp = 0, 
                   arma::uword step = 1, 
                   arma::uword look_back = 1, 
                   arma::uword stub = 0,
                   std::string method = "moment", 
                   double conf_lev = 0.75) {
  
  // Allocate end points
  arma::uword nrows = tseries.n_rows;
  arma::uvec endpts;
  arma::uvec startpts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    endpts = calc_endpoints(nrows, step, stub);
  } else {
    // Copy end points
    endpts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    startpts = calc_startpoints(endpts, look_back);
  } else {
    // Copy start points
    startpts = startp;
  }  // end if
  
  // Allocate matrix of statistics
  arma::uword numpts = endpts.n_elem;
  arma::mat stats = arma::zeros<mat>(numpts, tseries.n_cols);
  
  // Perform loop over the end points
  if (fun == "calc_mean") {
    // Calculate the dispersion (variance)
    for (arma::uword ep = 0; ep < numpts; ep++) {
      // Calculate kurtosis
      if (endpts(ep) > startpts(ep)) {
        stats.row(ep) = calc_mean(tseries.rows(startpts(ep), endpts(ep)), method, conf_lev);
      }  // end if
    }  // end for
  } else if (fun == "calc_var") {
    // Calculate the dispersion (variance)
    for (arma::uword ep = 0; ep < numpts; ep++) {
      // Calculate kurtosis
      if (endpts(ep) > startpts(ep)) {
        stats.row(ep) = calc_var(tseries.rows(startpts(ep), endpts(ep)), method, conf_lev);
      }  // end if
    }  // end for
  } else if (fun == "calc_skew") {
    // Perform loop over the end points
    for (arma::uword ep = 0; ep < numpts; ep++) {
      // Calculate kurtosis
      if (endpts(ep) > startpts(ep)) {
        stats.row(ep) = calc_skew(tseries.rows(startpts(ep), endpts(ep)), method, conf_lev);
      }  // end if
    }  // end for
  } else if (fun == "calc_kurtosis") {
    // Perform loop over the end points
    for (arma::uword ep = 0; ep < numpts; ep++) {
      // Calculate kurtosis
      if (endpts(ep) > startpts(ep)) {
        stats.row(ep) = calc_kurtosis(tseries.rows(startpts(ep), endpts(ep)), method, conf_lev);
      }  // end if
    }  // end for
  } else {
    cout << "Wrong calc method!" << endl;
    return stats;
  }  // end if
  
  return stats;
  
}  // end roll_fun




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
//' returns <- na.omit(rutils::etfenv$returns$VTI)
//' # Estimate the GARCH volatility of VTI returns
//' garch_data <- HighFreq::sim_garch(omega=om_ega, alpha=alpha,  beta=betav, 
//'   innov=returns, is_random=FALSE)
//' # Plot dygraph of the estimated GARCH volatility
//' dygraphs::dygraph(xts::xts(sqrt(garch_data[, 2]), index(returns)), 
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
    arma::mat variance(nrows, 1);
    variance[0] = omega/(1-alpha-beta);
    arma::mat returns(nrows, 1);
    returns[0] = std::sqrt(variance[0])*innov[0];
    
    for (arma::uword it = 1; it < nrows; it++) {
      returns[it] = std::sqrt(variance[it-1])*innov[it];
      variance[it] = omega + alpha*pow(returns[it], 2) + beta*variance[it-1];
    }  // end for
    return join_rows(returns, variance);
  } else {
    // The innovations are historical returns
    arma::mat variance = arma::square(innov);
    variance[0] = omega/(1-alpha-beta);
    for (arma::uword it = 1; it < nrows; it++) {
      variance[it] = omega + alpha*variance[it] + beta*variance[it-1];
    }  // end for
    return join_rows(innov, variance);
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
  arma::mat prices = arma::zeros<mat>(nrows, 1);
  arma::mat returns = arma::zeros<mat>(nrows, 1);
  
  returns.row(0) = innov.row(0);
  prices.row(0) = init_price;
  for (arma::uword it = 1; it < nrows; it++) {
    returns.row(it) = theta*(eq_price - prices.row(it-1)) + innov.row(it);
    prices.row(it) = prices.row(it-1) + returns.row(it);
  }  // end for
  
  return prices;
  
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
  arma::mat prices = arma::zeros<mat>(nrows, 1);
  arma::mat returns = arma::zeros<mat>(nrows, 1);
  
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
//' filtered <- filter(innov, filter=coeff, method="recursive")
//' # Calculate recursive filter using RcppArmadillo
//' returns <- HighFreq::sim_ar(coeff, innov)
//' # Compare the two methods
//' all.equal(as.numeric(returns), as.numeric(filtered))
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
  arma::mat returns = arma::zeros<mat>(nrows, 1);

  // Warmup period
  returns.row(0) = innov.row(0);
  returns.row(1) = innov.row(1) + coeffr.row(ncoeff-1) * returns.row(0);
  for (arma::uword it = 2; it < ncoeff; it++) {
    returns.row(it) = innov.row(it) + arma::dot(coeffr.rows(ncoeff-it, ncoeff-1), returns.rows(0, it-1));
  }  // end for
  
  // Remaining periods
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
  arma::mat prices = arma::zeros<mat>(nrows, 1);
  arma::mat returns = arma::zeros<mat>(nrows, 1);

  // Warmup period
  returns.row(0) = innov.row(0);
  prices.row(0) = init_price;
  returns.row(1) = theta*(eq_price - prices.row(0)) + coeffr.row(ncoeff-1) * returns.row(0) + innov.row(1);
  prices.row(1) = prices.row(0) + returns.row(1);
  for (arma::uword it = 2; it < ncoeff; it++) {
    returns.row(it) = theta*(eq_price - prices.row(it-1)) + arma::dot(coeffr.rows(ncoeff-it, ncoeff-1), returns.rows(0, it-1)) + innov.row(it);
    prices.row(it) = prices.row(it-1) + returns.row(it);
  }  // end for
  
  // Remaining periods
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
//' returns <- na.omit(rutils::etfenv$returns$VTI)
//' # Calculate the log-likelihood of VTI returns assuming GARCH(1,1)
//' HighFreq::lik_garch(omega=om_ega, alpha=alpha,  beta=betav, returns=returns)
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
  arma::mat variance = garch_data.col(1);
  // Apply floor to variance
  variance.transform([&minval](double x) {return max(x, minval);});
  // Lag the variance
  variance = lagit(variance, 1, false);
  // Calculate the log-likelihood
  double likelihood = -conv_to<double>::from(arma::sum(pow(returns, 2)/variance + log(variance)));
  
  return likelihood;
  
}  // end lik_garch



////////////////////////////////////////////////////////////
// Functions for backtests
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
//' Calculate the optimal portfolio weights for different types of objective
//' functions.
//' 
//' @param \code{returns} A \emph{time series} or a \emph{matrix} of returns
//'   data (the returns in excess of the risk-free rate).
//'   
//' @param \code{method} A \emph{string} specifying the method for
//'   calculating the weights (see Details) (the default is \code{method =
//'   "ranksharpe"})
//'   
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small singular values in order to regularize the inverse of the
//'   \code{returns} matrix (the default is \code{1e-5}).
//'   
//' @param \code{eigen_max} An \emph{integer} equal to the number of singular
//'   values used for calculating the shrinkage inverse of the \code{returns}
//'   matrix (the default is \code{0} - equivalent to \code{eigen_max} equal to
//'   the number of columns of \code{returns}).
//'   
//' @param \code{conf_lev} The confidence level for calculating the
//'   quantiles (the default is \code{conf_lev = 0.75}).
//'
//' @param \code{alpha} The shrinkage intensity between \code{0} and \code{1}.
//'   (the default is \code{0}).
//' 
//' @param \code{scale} A \emph{Boolean} specifying whether the weights should
//'   be scaled (the default is \code{scale = TRUE}).
//'
//' @param \code{vol_target} A \emph{numeric} volatility target for scaling the
//'   weights (the default is \code{0.001})
//'   
//' @return A column \emph{vector} of the same length as the number of columns
//'   of \code{returns}.
//'
//' @details
//'   The function \code{calc_weights()} calculates the optimal portfolio
//'   weights for different types of methods, using \code{RcppArmadillo}
//'   \code{C++} code.
//' 
//'   If \code{method = "ranksharpe"} (the default) then it calculates the
//'   weights as the ranks (order index) of the trailing Sharpe ratios of the
//'   asset \code{returns}.
//'
//'   If \code{method = "rank"} then it calculates the weights as the ranks
//'   (order index) of the last row of the \code{returns}.
//'
//'   If \code{method = "max_sharpe"} then \code{calc_weights()} calculates
//'   the weights of the maximum Sharpe portfolio, by multiplying the inverse of
//'   the covariance \emph{matrix} times the mean column returns.
//'
//'   If \code{method = "min_var"} then it calculates the weights of the
//'   minimum variance portfolio under linear constraints.
//'
//'   If \code{method = "min_varpca"} then it calculates the weights of the
//'   minimum variance portfolio under quadratic constraints (which is the
//'   highest order principal component).
//'
//'   If \code{scale = TRUE} (the default) then the weights are scaled so that
//'   the resulting portfolio has a volatility equal to \code{vol_target}.
//'
//'   \code{calc_weights()} calculates the shrinkage inverse of the covariance
//'   \emph{matrix} of \code{returns} from its eigen decomposition.  It applies
//'   dimension regularization by selecting only the largest eigenvalues equal
//'   in number to \code{eigen_max}. 
//'   
//'   In addition, \code{calc_weights()} applies shrinkage to the columns of
//'   \code{returns}, by shrinking their means to their common mean value. The
//'   shrinkage intensity \code{alpha} determines the amount of shrinkage that
//'   is applied, with \code{alpha = 0} representing no shrinkage (with the
//'   column means of \code{returns} unchanged), and \code{alpha = 1}
//'   representing complete shrinkage (with the column means of \code{returns}
//'   all equal to the single mean of all the columns).
//' 
//' @examples
//' \dontrun{
//' # Calculate covariance matrix of ETF returns
//' returns <- na.omit(rutils::etfenv$returns[, 1:16])
//' eigend <- eigen(cov(returns))
//' # Calculate shrinkage inverse of covariance matrix
//' eigen_max <- 3
//' eigenvec <- eigend$vectors[, 1:eigen_max]
//' eigenval <- eigend$values[1:eigen_max]
//' inverse <- eigenvec %*% (t(eigenvec) / eigenval)
//' # Define shrinkage intensity and apply shrinkage to the mean returns
//' alpha <- 0.5
//' colmeans <- colMeans(returns)
//' colmeans <- ((1-alpha)*colmeans + alpha*mean(colmeans))
//' # Calculate weights using R
//' weights <- inverse %*% colmeans
//' n_col <- NCOL(returns)
//' weightsr <- weightsr*sd(returns %*% rep(1/n_col, n_col))/sd(returns %*% weightsr)
//' # Calculate weights using RcppArmadillo
//' weights <- drop(HighFreq::calc_weights(returns, eigen_max, alpha=alpha))
//' all.equal(weights, weightsr)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec calc_weights(const arma::mat& returns, // Portfolio returns
                       std::string method = "ranksharpe",
                       double eigen_thresh = 1e-5,
                       arma::uword eigen_max = 0,
                       double conf_lev = 0.1,
                       double alpha = 0.0,
                       bool scale = true,
                       double vol_target = 0.01) {
  // Initialize
  arma::vec weights(returns.n_cols, fill::zeros);
  if (eigen_max == 0)  eigen_max = returns.n_cols;
  
  // Switch for the different methods for weights
  switch(calc_method(method)) {
  case methodenum::ranksharpe: {
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to<vec>::from(arma::sort_index(arma::sort_index(meancols)));
    weights = (weights - arma::mean(weights));
    break;
  }  // end ranksharpe
  case methodenum::max_sharpe: {
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Shrink meancols to the mean of returns
    meancols = ((1-alpha)*meancols + alpha*arma::mean(meancols));
    // Apply shrinkage inverse
    // arma::mat inverse = calc_inv(cov(returns), eigen_max);
    // weights = calc_inv(cov(returns), eigen_max)*meancols;
    weights = calc_inv(cov(returns), eigen_thresh, eigen_max)*meancols;
    break;
  }  // end max_sharpe
  case methodenum::max_sharpe_median: {
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::median(returns, 0));
    // Shrink meancols to the mean of returns
    meancols = ((1-alpha)*meancols + alpha*arma::median(meancols));
    // Apply shrinkage inverse
    // arma::mat inverse = calc_inv(cov(returns), eigen_max);
    weights = calc_inv(cov(returns), eigen_thresh, eigen_max)*meancols;
    break;
  }  // end max_sharpe_median
  case methodenum::min_var: {
    // Apply shrinkage inverse to unit vector
    weights = calc_inv(cov(returns), eigen_thresh, eigen_max)*arma::ones(returns.n_cols);
    break;
  }  // end min_var
  case methodenum::min_varpca: {
    // Calculate highest order principal component
    arma::vec eigenval;
    arma::mat eigenvec;
    arma::eig_sym(eigenval, eigenvec, arma::cov(returns));
    weights = eigenvec.col(0);
    break;
  }  // end min_varpca
  case methodenum::rank: {
    // Mean returns by columns
    arma::vec meancols = arma::trans(arma::mean(returns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to<vec>::from(arma::sort_index(arma::sort_index(meancols)));
    weights = (weights - arma::mean(weights));
    break;
  }  // end rank
  case methodenum::rankrob: {
    // Median returns by columns
    arma::vec meancols = arma::trans(arma::median(returns, 0));
    // meancols = ((1-alpha)*meancols + alpha*arma::mean(meancols));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    meancols = meancols/sd_cols;
    // Apply shrinkage inverse
    // arma::mat inverse = calc_inv(cov(returns), eigen_max);
    // weights = calc_inv(cov(returns), eigen_max)*meancols;
    // weights = calc_inv(cov(returns), eigen_max)*meancols;
    // // Standard deviation by columns
    // arma::vec sd_cols = meancols;
    // for (arma::uword it=0; it < returns.n_cols; it++) {
    //   sd_cols(it) = arma::median(arma::abs((returns.col(it) - sd_cols)));
    // }  // end for
    // sd_cols.replace(0, 1);
    // meancols = meancols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to<vec>::from(arma::sort_index(arma::sort_index(meancols)));
    // level;
    weights = (weights - arma::mean(weights));
    break;
  }  // end rankrob
  case methodenum::quantile: {
    // Sum of quantiles for columns
    arma::vec levels = {conf_lev, 1-conf_lev};
    weights = conv_to<vec>::from(arma::sum(arma::quantile(returns, levels, 0), 0));
    // Weights equal to ranks
    weights = conv_to<vec>::from(arma::sort_index(arma::sort_index(weights)));
    weights = (weights - arma::mean(weights));
    break;
  }  // end quantile
  default : {
    cout << "Warning: Invalid method parameter: " << method << endl;
    return arma::ones(returns.n_cols);
  }  // end default
  }  // end switch
  
  if (scale == TRUE) {
    // return weights/std::sqrt(sum(square(weights)));
    // return weights/sum(weights);
    // Returns of equally weighted portfolio
    // arma::vec meanrows = arma::mean(returns, 1);
    // Returns of weighted portfolio
    // arma::vec returns_portf = returns*weights;
    // Scale weights to equally weighted portfolio and return them
    // return weights*arma::stddev(arma::mean(returns, 1))/arma::stddev(returns*weights);
    // Scale weights so the resulting portfolio has a volatility equal to vol_target
    return weights*vol_target/arma::stddev(returns*weights);
  }  // end if
  
  return weights;
  
}  // end calc_weights



////////////////////////////////////////////////////////////
//' Simulate (backtest) a rolling portfolio optimization strategy, using
//' \code{RcppArmadillo}.
//' 
//' @param \code{returns} A \emph{time series} or a \emph{matrix} of returns
//'   data (the returns in excess of the risk-free rate).
//'   
//' @param \code{excess} A \emph{time series} or a \emph{matrix} of excess
//'   returns data (the returns in excess of the risk-free rate).
//'   
//' @param \code{startp} An \emph{integer vector} of start points.
//' 
//' @param \code{endp} An \emph{integer vector} of end points.
//' 
//' @param \code{lambda} A \emph{numeric} decay factor to multiply the past
//'   portfolio weights.  (The default is \code{lambda = 0} - no memory.)
//'   
//' @param \code{coeff} A \emph{numeric} multiplier of the weights.  (The
//'   default is \code{1})
//'   
//' @param \code{bid_offer} A \emph{numeric} bid-offer spread (the default is
//'   \code{0})
//'
//' @param \code{method} A \emph{string} specifying the method for calculating
//'   the weights (see Details) (the default is \code{method = "ranksharpe"})
//'   
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small singular values in order to regularize the inverse of the
//'   \code{returns} matrix (the default is \code{1e-5}).
//'   
//' @param \code{eigen_max} An \emph{integer} equal to the number of singular
//'   values used for calculating the shrinkage inverse of the \code{returns}
//'   matrix (the default is \code{0} - equivalent to \code{eigen_max} equal to
//'   the number of columns of \code{returns}).
//'   
//' @param \code{conf_lev} The confidence level for calculating the
//'   quantiles (the default is \code{conf_lev = 0.75}).
//'
//' @param \code{alpha} The shrinkage intensity between \code{0} and \code{1}.
//'   (the default is \code{0}).
//' 
//' @param \code{scale} A \emph{Boolean} specifying whether the weights should
//'   be scaled (the default is \code{scale = TRUE}).
//'
//' @param \code{vol_target} A \emph{numeric} volatility target for scaling the
//'   weights (the default is \code{1e-5})
//'   
//' @return A column \emph{vector} of strategy returns, with the same length as
//'   the number of rows of \code{returns}.
//'
//' @details
//'   The function \code{back_test()} performs a backtest simulation of a
//'   rolling portfolio optimization strategy over a \emph{vector} of the end
//'   points \code{endp}.
//'   
//'   It performs a loop over the end points \code{endp}, and subsets the
//'   \emph{matrix} of the excess asset returns \code{excess} along its rows,
//'   between the corresponding \emph{start point} and the \emph{end point}. It
//'   passes the subset matrix of excess returns into the function
//'   \code{calc_weights()}, which calculates the optimal portfolio weights at
//'   each \emph{end point}. The arguments \code{eigen_max}, \code{alpha},
//'   \code{method}, and \code{scale} are also passed to the function
//'   \code{calc_weights()}.
//'   
//'   It then recursively averages the weights \eqn{w_i} at the \emph{end point
//'   = i} with the weights \eqn{w_{i-1}} from the previous \emph{end point =
//'   (i-1)}, using the decay factor \code{lambda = \eqn{\lambda}}:
//'   \deqn{
//'     w_i = (1-\lambda) w_i + \lambda w_{i-1}
//'   }
//'   The purpose of averaging the weights is to reduce their variance to
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
//'   multiplying the bid-offer spread \code{bid_offer} times the absolute
//'   difference between the current weights minus the weights from the previous
//'   period. Then it subtracts the transaction costs from the out-of-sample
//'   strategy returns.
//'   
//'   The function \code{back_test()} returns a \emph{time series} (column
//'   \emph{vector}) of strategy returns, of the same length as the number of
//'   rows of \code{returns}.
//'
//' @examples
//' \dontrun{
//' # Calculate the ETF daily excess returns
//' returns <- na.omit(rutils::etfenv$returns[, 1:16])
//' # riskf is the daily risk-free rate
//' riskf <- 0.03/260
//' excess <- returns - riskf
//' # Define monthly end points without initial warmpup period
//' endp <- rutils::calc_endpoints(returns, interval="months")
//' endp <- endp[endp > 0]
//' nrows <- NROW(endp)
//' # Define 12-month look-back interval and start points over sliding window
//' look_back <- 12
//' startp <- c(rep_len(1, look_back-1), endp[1:(nrows-look_back+1)])
//' # Define shrinkage and regularization intensities
//' alpha <- 0.5
//' eigen_max <- 3
//' # Simulate a monthly rolling portfolio optimization strategy
//' pnls <- HighFreq::back_test(excess, returns, 
//'                             startp-1, endp-1, 
//'                             eigen_max = eigen_max, 
//'                             alpha = alpha)
//' pnls <- xts::xts(pnls, index(returns))
//' colnames(pnls) <- "strat_rets"
//' # Plot dygraph of strategy
//' dygraphs::dygraph(cumsum(pnls), 
//'   main="Cumulative Returns of Max Sharpe Portfolio Strategy")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat back_test(const arma::mat& excess, // Asset excess returns
                    const arma::mat& returns, // Asset returns
                    arma::uvec startp, // Start points
                    arma::uvec endp, // End points
                    double lambda, // Decay factor for averaging the portfolio weights
                    std::string method = "ranksharpe",  // The method for calculating the weights
                    double eigen_thresh = 1e-5,
                    arma::uword eigen_max = 0,
                    double conf_lev = 0.1,
                    double alpha = 0.0,
                    bool scale = true,
                    double vol_target = 0.01,
                    double coeff = 1.0,
                    double bid_offer = 0.0) {
  
  double lambda1 = 1-lambda;
  arma::uword nweights = returns.n_cols;
  arma::vec weights(nweights, fill::zeros);
  arma::vec weights_past = ones(nweights)/sqrt(nweights);
  arma::mat pnls = zeros(returns.n_rows, 1);

  // Perform loop over the end points
  for (arma::uword it = 1; it < endp.size(); it++) {
    // cout << "it: " << it << endl;
    // Calculate portfolio weights
    weights = coeff*calc_weights(excess.rows(startp(it-1), endp(it-1)), method, eigen_thresh, eigen_max, conf_lev, alpha, scale, vol_target);
    // Calculate the weights as the weighted sum with past weights
    weights = lambda1*weights + lambda*weights_past;
    // Calculate out-of-sample returns
    pnls.rows(endp(it-1)+1, endp(it)) = returns.rows(endp(it-1)+1, endp(it))*weights;
    // Add transaction costs
    pnls.row(endp(it-1)+1) -= bid_offer*sum(abs(weights - weights_past))/2;
    // Copy the weights
    weights_past = weights;
  }  // end for
  
  // Return the strategy pnls
  return pnls;
  
}  // end back_test

