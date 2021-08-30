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
//' re_turns <- rnorm(1e6)
//' # Compare lag_vec() with rutils::lag_it()
//' all.equal(drop(HighFreq::lag_vec(re_turns)), 
//'   rutils::lag_it(re_turns))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::lag_vec(re_turns),
//'   Rcode=rutils::lag_it(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec lag_vec(const arma::vec& tseries, 
                  arma::sword lagg = 1, 
                  bool pad_zeros = true) {
  
  arma::uword num_rows = (tseries.n_elem-1);
  
  if (lagg > 0) {
    if (pad_zeros) {
      // Pad front with zeros
      return arma::join_cols(arma::zeros<vec>(lagg), 
                             tseries.subvec(0, num_rows-lagg));
    } else {
      // Pad front with first element of tseries
      return arma::join_cols(arma::repelem(tseries.subvec(0, 0), lagg, 1), 
                             tseries.subvec(0, num_rows-lagg));
    }  // end if
  } else {
    if (pad_zeros) {
      // Pad back with zeros
      return arma::join_cols(tseries.subvec(-lagg, num_rows), 
                             arma::zeros<vec>(-lagg));
    } else {
      // Pad back with last element of tseries
      return arma::join_cols(tseries.subvec(-lagg, num_rows), 
                             arma::repelem(tseries.subvec(num_rows, num_rows), -lagg, 1));
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
//'   The function \code{lag_it()} applies a lag to the input \emph{matrix} by
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
//' re_turns <- matrix(rnorm(5e6), nc=5)
//' # Compare lag_it() with rutils::lag_it()
//' all.equal(HighFreq::lag_it(re_turns), 
//'   rutils::lag_it(re_turns))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::lag_it(re_turns),
//'   Rcode=rutils::lag_it(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat lag_it(const arma::mat& tseries, 
                 arma::sword lagg = 1, 
                 bool pad_zeros = true) {
  
  arma::uword num_rows = (tseries.n_rows-1);
  arma::uword num_cols = tseries.n_cols;
  
  if (lagg > 0) {
    // Positive lag
    if (pad_zeros) {
      // Pad front with zeros
      return arma::join_cols(arma::zeros<mat>(lagg, num_cols), 
                             tseries.rows(0, num_rows-lagg));
    } else {
      // Pad front with first element of tseries
      return arma::join_cols(arma::repmat(tseries.rows(0, 0), lagg, 1), 
                             tseries.rows(0, num_rows-lagg));
    }  // end if
  } else {
    // Negative lag
    if (pad_zeros) {
      // Pad back with zeros
      return arma::join_cols(tseries.rows(-lagg, num_rows), 
                             arma::zeros<mat>(-lagg, num_cols));
    } else {
      // Pad back with last element of tseries
      return arma::join_cols(tseries.rows(-lagg, num_rows), 
                             arma::repmat(tseries.rows(num_rows, num_rows), -lagg, 1));
    }  // end if
  }  // end if
  
  // Old code below
  // if (lagg > 0)
  //   // Positive lag
  //   return arma::join_cols(arma::repelem(tseries.row(0), lagg, 1), 
  //                          tseries.rows(0, num_rows-lagg));
  // else
  //   // Negative lag
  //   return arma::join_cols(tseries.rows(-lagg, num_rows), 
  //                          arma::repelem(tseries.row(num_rows), -lagg, 1));
  
}  // end lag_it




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
//' re_turns <- rnorm(1e6)
//' # Compare diff_vec() with rutils::diff_it()
//' all.equal(drop(HighFreq::diff_vec(re_turns, lagg=3, pad=TRUE)),
//'   rutils::diff_it(re_turns, lagg=3))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::diff_vec(re_turns, lagg=3, pad=TRUE),
//'   Rcode=rutils::diff_it(re_turns, lagg=3),
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
//' Calculate the row differences of a a \emph{time series} or a \emph{matrix}
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
//'   The function \code{diff_it()} calculates the differences between the rows
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
//'   The function \code{diff_it()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, which makes it much faster than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Create a matrix of random data
//' da_ta <- matrix(sample(15), nc=3)
//' # Calculate differences with lagged rows
//' HighFreq::diff_it(da_ta, lagg=2)
//' # Calculate differences with advanced rows
//' HighFreq::diff_it(da_ta, lagg=-2)
//' # Compare HighFreq::diff_it() with rutils::diff_it()
//' all.equal(HighFreq::diff_it(da_ta, lagg=2), 
//'   rutils::diff_it(da_ta, lagg=2), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::diff_it(da_ta, lagg=2),
//'   Rcode=rutils::diff_it(da_ta, lagg=2),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat diff_it(const arma::mat& tseries, 
                  arma::sword lagg = 1, 
                  bool pad_zeros = true) {
  
  arma::uword num_rows = (tseries.n_rows-1);
  arma::uword num_cols = tseries.n_cols;
  
  if (lagg > 0) {
    // Positive lag
    // Matrix difference without padding
    arma::mat diff_mat = (tseries.rows(lagg, num_rows) - tseries.rows(0, num_rows - lagg));
    if (pad_zeros) {
      // Pad diff_mat with zeros at the front
      return arma::join_cols(arma::zeros<mat>(lagg, num_cols), diff_mat);
    } else {
      // Don't pad the output
      return diff_mat;
    }  // end if pad_zeros
  } else {
    // Negative lag
    // Matrix difference without padding
    arma::mat diff_mat = (tseries.rows(0, num_rows + lagg) - tseries.rows(-lagg, num_rows));
    if (pad_zeros) {
      // Pad diff_mat with zeros at the back
      return arma::join_cols(diff_mat, arma::zeros<mat>(-lagg, num_cols));
    } else {
      // Don't pad the output
      return diff_mat;
    }  // end if pad_zeros
  }  // end if lagg
  
}  // end diff_it




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
//' end_p <- HighFreq::calc_endpoints(25, 5)
//' # Calculate start points corresponding to the end points
//' start_p <- HighFreq::calc_startpoints(end_p, 2)
//'
//' @export
// [[Rcpp::export]]
arma::uvec calc_startpoints(arma::uvec endp, arma::uword look_back) {
  
  arma::uword num_pts = endp.n_elem;
  arma::uvec startp = arma::join_cols(arma::zeros<uvec>(look_back), 
                                      endp.subvec(0, num_pts - look_back - 1) + 1);
  
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
//' mat_rix <- matrix(round(runif(25e4), 2), nc=5e2)
//' vec_tor <- round(runif(5e2), 2)
//' prod_uct <- vec_tor*mat_rix
//' # Multiply the matrix in place
//' HighFreq::mult_vec_mat(vec_tor, mat_rix)
//' all.equal(prod_uct, mat_rix)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_vec_mat(vec_tor, mat_rix),
//'     Rcode=vec_tor*mat_rix,
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' 
//' # Multiply matrix rows using R
//' mat_rix <- matrix(round(runif(25e4), 2), nc=5e2)
//' vec_tor <- round(runif(5e2), 2)
//' prod_uct <- t(vec_tor*t(mat_rix))
//' # Multiply the matrix in place
//' HighFreq::mult_vec_mat(vec_tor, mat_rix, by_col=FALSE)
//' all.equal(prod_uct, mat_rix)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     Rcpp=HighFreq::mult_vec_mat(vec_tor, mat_rix, by_col=FALSE),
//'     Rcode=t(vec_tor*t(mat_rix)),
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::uword mult_vec_mat(arma::vec& vector,
                         arma::mat& matrix,
                         bool by_col = true) {
  
  arma::uword num_elem = vector.n_elem;
  arma::uword num_rows = matrix.n_rows;
  arma::uword num_cols = matrix.n_cols;
  
  if ((num_cols == num_rows) && (num_elem == num_rows)) {
    if (by_col) {
      // Multiply each column of matrix by vector
      matrix.each_col() %= vector;
      return num_rows;
    } else {
      // Multiply each row of matrix by vector
      matrix.each_row() %= conv_to<rowvec>::from(vector);
      return num_cols;
    }
  } else if (num_elem == num_rows) {
    // Multiply each column of matrix by vector
    matrix.each_col() %= vector;
    return num_rows;
  } else if (num_elem == num_cols) {
    // Multiply each row of matrix by vector
    matrix.each_row() %= conv_to<rowvec>::from(vector);
    return num_cols;
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
//' da_ta <- matrix(rnorm(5e6), nc=5)
//' # Calculate eigen decomposition
//' ei_gen <- HighFreq::calc_eigen(scale(da_ta, scale=FALSE))
//' # Calculate PCA
//' pc_a <- prcomp(da_ta)
//' # Compare PCA with eigen decomposition
//' all.equal(pc_a$sdev^2, drop(ei_gen$values))
//' all.equal(abs(unname(pc_a$rotation)), abs(ei_gen$vectors))
//' # Compare the speed of Rcpp with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_eigen(da_ta),
//'   Rcode=prcomp(da_ta),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List calc_eigen(const arma::mat& tseries) {
  
  arma::mat eigen_vec;
  arma::vec eigen_val;
  arma::eig_sym(eigen_val, eigen_vec, arma::cov(tseries));
  
  // Reverse the order of elements from largest eigenvalue to smallest, similar to R
  return Rcpp::List::create(Named("values") = arma::flipud(eigen_val),
                            Named("vectors") = arma::fliplr(eigen_vec));
  
}  // end calc_eigen



////////////////////////////////////////////////////////////
//' Calculate the regularized inverse of a rectangular \emph{matrix} of data
//' using Singular Value Decomposition (\emph{SVD}).
//' 
//' @param \code{tseries} A \emph{time series} or \emph{matrix} of returns data.
//' 
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small singular values in order to regularize the inverse of the
//'   matrix \code{tseries} (the default is \code{0.001}).
//'   
//' @param \code{eigen_max} An \emph{integer} equal to the number of singular
//'   values used for calculating the regularized inverse of the matrix
//'   \code{tseries} (the default is \code{0} - equivalent to \code{eigen_max}
//'   equal to the number of columns of \code{tseries}).
//'
//' @return A \emph{matrix} equal to the regularized inverse of the matrix
//'   \code{tseries}.
//'
//' @details 
//'   The function calc_inv() calculates the regularized inverse of
//'   \code{tseries} using Singular Value Decomposition (\emph{SVD}).
//'   
//'   If \code{eigen_max} is given, then it calculates the regularized inverse
//'   from the \emph{SVD} using the first \code{eigen_max} largest singular
//'   values.  For example, if \code{eigen_max = 3} then it only uses the
//'   \code{3} largest singular values.
//'   If \code{eigen_max} is set equal to the number of columns of
//'   \code{tseries} then it uses all the singular values without any
//'   regularization.
//'
//'   If \code{eigen_max} is not given then it calculates the regularized
//'   inverse using the function \code{arma::pinv()}. Then it discards small
//'   singular values that are less than the threshold level
//'   \code{eigen_thresh}.
//'   
//' @examples
//' \dontrun{
//' # Calculate ETF returns
//' re_turns <- na.omit(rutils::etf_env$re_turns)
//' # Calculate regularized inverse using RcppArmadillo
//' in_verse <- HighFreq::calc_inv(re_turns, eigen_max=3)
//' # Calculate regularized inverse from SVD in R
//' s_vd <- svd(re_turns)
//' eigen_max <- 1:3
//' inverse_r <-  s_vd$v[, eigen_max] %*% (t(s_vd$u[, eigen_max]) / s_vd$d[eigen_max])
//' # Compare RcppArmadillo with R
//' all.equal(in_verse, inverse_r)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& tseries,
                   double eigen_thresh = 0.001, 
                   arma::uword eigen_max = 0) {
  
  if (eigen_max == 0) {
    // Calculate the inverse using arma::pinv()
    return arma::pinv(tseries, eigen_thresh);
  } else {
    // Calculate the regularized inverse using SVD decomposition
    
    // Allocate SVD
    arma::vec svd_val;
    arma::mat svd_u, svd_v;
    
    // Calculate the SVD
    arma::svd(svd_u, svd_val, svd_v, tseries);
    
    // Subset the SVD
    eigen_max = eigen_max - 1;
    // For no regularization: eigen_max = tseries.n_cols
    svd_u = svd_u.cols(0, eigen_max);
    svd_v = svd_v.cols(0, eigen_max);
    svd_val = svd_val.subvec(0, eigen_max);
    
    // Calculate the inverse from the SVD
    return svd_v*arma::diagmat(1/svd_val)*svd_u.t();
    
  }  // end if
  
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
//' re_turns <- matrix(rnorm(20000), nc=20)
//' scale_d <- calc_scaled(tseries=re_turns, use_median=FALSE)
//' scale_d2 <- scale(re_turns)
//' all.equal(scale_d, scale_d2, check.attributes=FALSE)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=calc_scaled(tseries=re_turns, use_median=FALSE),
//'   Rcode=scale(re_turns),
//'   times=100))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_scaled(const arma::mat& tseries, bool use_median=false) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::uword num_cols = tseries.n_cols;
  arma::mat scaled_mat(num_rows, num_cols, fill::zeros);
  arma::vec scale_d(num_rows, fill::zeros);
  double cen_ter;
  
  // Return zeros if not enough data
  if (num_rows < 3) {
    return tseries;
  }  // end if
  
  // Perform a loop over the columns
  for (arma::uword it=0; it < num_cols; it++) {
    if (use_median) {
      cen_ter = arma::median(tseries.col(it));
      scale_d = (tseries.col(it) - cen_ter);
      scale_d = scale_d/arma::median(arma::abs(scale_d));
      scaled_mat.col(it) = scale_d;
    } else {
      cen_ter = arma::mean(tseries.col(it));
      scale_d = (tseries.col(it) - cen_ter);
      scale_d = scale_d/arma::stddev(scale_d);
      scaled_mat.col(it) = scale_d;
    }  // end if
  }  // end for
  
  return scaled_mat;
  
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
//' da_ta <- round(runif(7), 2)
//' # Calculate the ranks of the elements in two ways
//' all.equal(rank(da_ta), drop(HighFreq::calc_ranks(da_ta)))
//' # Create a time series of random data
//' da_ta <- xts::xts(runif(7), seq.Date(Sys.Date(), by=1, length.out=7))
//' # Calculate the ranks of the elements in two ways
//' all.equal(rank(coredata(da_ta)), drop(HighFreq::calc_ranks(da_ta)))
//' # Compare the speed of RcppArmadillo with R code
//' da_ta <- runif(7)
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=calc_ranks(da_ta),
//'   Rcode=rank(da_ta),
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
//' oh_lc <- coredata(rutils::etf_env$VTI[, 1:5])
//' # Aggregate to single row matrix
//' ohlc_agg <- HighFreq::agg_ohlc(oh_lc)
//' # Compare with calculation in R
//' all.equal(drop(ohlc_agg),
//'   c(oh_lc[1, 1], max(oh_lc[, 2]), min(oh_lc[, 3]), oh_lc[NROW(oh_lc), 4], sum(oh_lc[, 5])), 
//'   check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat agg_ohlc(const arma::mat& tseries) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::uword num_cols = tseries.n_cols;
  
  // Number of output columns
  arma::uword num_ohlc = num_cols;
  if (num_cols < 4)
    // Add volume column for non-OHLC data
    num_ohlc = 4 + num_cols - 1;
  // Allocate output matrix
  arma::mat ohlc(1, num_ohlc);
  
  if (num_cols < 4) {
    // Aggregate time series into a single bar of OHLC data.
    ohlc(0, 0) = tseries(0, 0);
    ohlc(0, 1) = arma::max(tseries.col(0));
    ohlc(0, 2) = arma::min(tseries.col(0));
    ohlc(0, 3) = tseries(num_rows-1, 0);
    if (num_cols == 2) {
      // Aggregate volume data.
      ohlc(0, 4) = arma::sum(tseries.col(1));
    }  // end if
  } else {
    // Aggregate OHLC time series into a single bar of OHLC data.
    ohlc(0, 0) = tseries(0, 0);
    ohlc(0, 1) = arma::max(tseries.col(1));
    ohlc(0, 2) = arma::min(tseries.col(2));
    ohlc(0, 3) = tseries(num_rows-1, 3);
    if (num_cols == 5) {
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
//' oh_lc <- rutils::etf_env$VTI[, 1:5]
//' # Define end points at 25 day intervals
//' end_p <- HighFreq::calc_endpoints(NROW(oh_lc), step=25)
//' # Aggregate over end_p:
//' ohlc_agg <- HighFreq::roll_ohlc(tseries=oh_lc, endp=end_p)
//' # Compare with xts::to.period()
//' ohlc_agg_xts <- .Call("toPeriod", oh_lc, as.integer(end_p+1), TRUE, NCOL(oh_lc), FALSE, FALSE, colnames(oh_lc), PACKAGE="xts")
//' all.equal(ohlc_agg, coredata(ohlc_agg_xts), check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_ohlc(const arma::mat& tseries, arma::uvec endp) {
  
  // arma::uword num_rows = tseries.n_rows;
  arma::uword num_cols = tseries.n_cols;
  arma::uword num_pts = endp.size();
  arma::mat ohlc_agg(num_pts-1, num_cols, fill::zeros);
  
  // Perform loop over the end points
  for (arma::uword it = 1; it < num_pts; it++) {
    // cout << "it: " << it << endl;
    // Aggregate the OHLC
    ohlc_agg.row(it-1) = agg_ohlc(tseries.rows(endp(it-1)+1, endp(it)));
  }  // end for
  
  // Return the aggregations
  return ohlc_agg;
  
}  // end roll_ohlc




////////////////////////////////////////////////////////////
//' Calculate the rolling sums over a single-column \emph{time series} or a
//' \emph{column vector} using \emph{Rcpp}.
//' 
//' @param \code{tseries} A single-column \emph{time series} or a \emph{column
//'   vector} (a single-column matrix).
//' 
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of elements of data used for calculating the sum.
//'
//' @return A \emph{column vector} of the same length as the argument
//'   \code{tseries}.
//'
//' @details 
//'   The function \code{roll_vec()} calculates a \emph{column vector} of
//'   rolling sums, over a \emph{column vector} of data, using fast \emph{Rcpp}
//'   \code{C++} code.  The function \code{roll_vec()} is several times faster
//'   than \code{rutils::roll_sum()} which uses vectorized \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Define a single-column matrix of returns
//' re_turns <- zoo::coredata(na.omit(rutils::etf_env$re_turns$VTI))
//' # Calculate rolling sums over 11-period look-back intervals
//' sum_rolling <- HighFreq::roll_vec(re_turns, look_back=11)
//' # Compare HighFreq::roll_vec() with rutils::roll_sum()
//' all.equal(HighFreq::roll_vec(re_turns, look_back=11), 
//'          rutils::roll_sum(re_turns, look_back=11), 
//'          check.attributes=FALSE)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_vec(re_turns, look_back=11),
//'   Rcode=rutils::roll_sum(re_turns, look_back=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_vec(const arma::mat& tseries, arma::uword look_back) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::mat rolling_sum(num_rows, 1);
  
  // Warmup period
  rolling_sum[0] = tseries[0];
  for (arma::uword it = 1; it < look_back; it++) {
    rolling_sum[it] = rolling_sum[it-1] + tseries[it];
  }  // end for
  
  // Remaining period
  for (arma::uword it = look_back; it < num_rows; it++) {
    rolling_sum[it] = rolling_sum[it-1] + tseries[it] - tseries[it-look_back];
  }  // end for
  
  return rolling_sum;
  
}  // end roll_vec




////////////////////////////////////////////////////////////
//' Calculate the rolling weighted sums over a single-column \emph{time series}
//' or a \emph{column vector} using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A single-column \emph{time series} or a \emph{column
//'   vector} (a single-column matrix).
//' 
//' @param \code{weights} A \emph{column vector} of weights.
//'
//' @return A \emph{column vector} of the same length as the argument
//'   \code{tseries}.
//'
//' @details 
//'   The function \code{roll_vecw()} calculates the rolling weighted sums of a
//'   \emph{column vector} over its past values (a convolution with the \emph{column vector}
//'   of weights), using \code{RcppArmadillo}. It performs a similar calculation
//'   as the standard \code{R} function \cr\code{stats::filter(x=series,
//'   filter=weight_s, method="convolution", sides=1)}, but it's over \code{6}
//'   times faster, and it doesn't produce any \code{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Define a single-column matrix of returns
//' re_turns <- zoo::coredata(na.omit(rutils::etf_env$re_turns$VTI))
//' # Create simple weights
//' weight_s <- matrix(c(1, rep(0, 10)))
//' # Calculate rolling weighted sums
//' weight_ed <- HighFreq::roll_vecw(tseries=re_turns, weights=weight_s)
//' # Compare with original
//' all.equal(zoo::coredata(re_turns), weight_ed, check.attributes=FALSE)
//' # Second example
//' # Create exponentially decaying weights
//' weight_s <- matrix(exp(-0.2*1:11))
//' weight_s <- weight_s/sum(weight_s)
//' # Calculate rolling weighted sums
//' weight_ed <- HighFreq::roll_vecw(tseries=re_turns, weights=weight_s)
//' # Calculate rolling weighted sums using filter()
//' filter_ed <- stats::filter(x=re_turns, filter=weight_s, method="convolution", sides=1)
//' # Compare both methods
//' all.equal(filter_ed[-(1:11)], weight_ed[-(1:11)], check.attributes=FALSE)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_vecw(tseries=re_turns, weights=weight_s),
//'   Rcode=stats::filter(x=re_turns, filter=weight_s, method="convolution", sides=1),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_vecw(const arma::mat& tseries, arma::mat& weights) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::uword look_back = weights.n_rows;
  arma::mat rolling_sum(num_rows, 1);
  arma::mat rev_weights = arma::reverse(weights);
  // arma::mat rev_weights = weights;
  
  // Warmup period
  rolling_sum.rows(0, look_back-2) = tseries.rows(0, look_back-2);
  
  // Remaining periods
  for (arma::uword it = look_back-1; it < num_rows; it++) {
    rolling_sum(it) = arma::dot(rev_weights, tseries.rows(it-look_back+1, it));
  }  // end for
  
  return rolling_sum;
  
}  // end roll_vecw




////////////////////////////////////////////////////////////
//' Calculate the rolling convolutions (weighted sums) of a \emph{time series}
//' with a \emph{column vector} of weights.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//' 
//' @param \code{weights} A \emph{column vector} of weights.
//'
//' @return A \emph{matrix} with the same dimensions as the input
//'   argument \code{tseries}.
//'
//' @details 
//'   The function \code{roll_conv()} calculates the convolutions of the
//'   \emph{matrix} columns with a \emph{column vector} of weights.  It performs
//'   a loop over the \emph{matrix} rows and multiplies the past (higher) values
//'   by the weights.  It calculates the rolling weighted sums of the past
//'   values.
//'   
//'   The function \code{roll_conv()} uses the \code{RcppArmadillo} function
//'   \code{arma::conv2()}. It performs a similar calculation to the standard
//'   \code{R} function \cr\code{filter(x=tseries, filter=weight_s,
//'   method="convolution", sides=1)}, but it's over \code{6} times faster, and
//'   it doesn't produce any leading \code{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Calculate a time series of returns
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("IEF", "VTI")])
//' # Create simple weights equal to a 1 value plus zeros
//' weight_s <- matrix(c(1, rep(0, 10)), nc=1)
//' # Calculate rolling weighted sums
//' weight_ed <- HighFreq::roll_conv(re_turns, weight_s)
//' # Compare with original
//' all.equal(coredata(re_turns), weight_ed, check.attributes=FALSE)
//' # Second example
//' # Calculate exponentially decaying weights
//' weight_s <- exp(-0.2*(1:11))
//' weight_s <- matrix(weight_s/sum(weight_s), nc=1)
//' # Calculate rolling weighted sums
//' weight_ed <- HighFreq::roll_conv(re_turns, weight_s)
//' # Calculate rolling weighted sums using filter()
//' filter_ed <- filter(x=re_turns, filter=weight_s, method="convolution", sides=1)
//' # Compare both methods
//' all.equal(filter_ed[-(1:11), ], weight_ed[-(1:11), ], check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_conv(const arma::mat& tseries, const arma::mat& weights) {
  
  arma::uword look_back = weights.n_rows-2;
  arma::uword num_rows = tseries.n_rows-1;
  
  // Calculate the convolutions
  arma::mat convmat = arma::conv2(tseries, weights, "full");
  
  // Copy the warmup period
  convmat.rows(0, look_back) = tseries.rows(0, look_back);
  
  return convmat.rows(0, num_rows);
  
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
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("VTI", "IEF")])
//' # Define parameters
//' look_back <- 22
//' # Calculate rolling sums and compare with rutils::roll_sum()
//' c_sum <- HighFreq::roll_sum(re_turns, look_back)
//' r_sum <- rutils::roll_sum(re_turns, look_back)
//' all.equal(c_sum, coredata(r_sum), check.attributes=FALSE)
//' # Calculate rolling sums using R code
//' r_sum <- apply(zoo::coredata(re_turns), 2, cumsum)
//' lag_sum <- rbind(matrix(numeric(2*look_back), nc=2), r_sum[1:(NROW(r_sum) - look_back), ])
//' r_sum <- (r_sum - lag_sum)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_sum(const arma::mat& tseries, arma::uword look_back = 1) {
  
  // Calculate the cumulative sum
  arma::mat cum_sum = arma::cumsum(tseries, 0);
  
  // Return the differences of the cumulative sum
  return diff_it(cum_sum, look_back, true);
  
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
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("VTI", "IEF")])
//' # Define end points at 25 day intervals
//' end_p <- HighFreq::calc_endpoints(NROW(re_turns), step=25)
//' # Define start points as 75 day lag of end points
//' start_p <- HighFreq::calc_startpoints(end_p, look_back=3)
//' # Calculate rolling sums using Rcpp
//' c_sum <- HighFreq::roll_sumep(re_turns, startp=start_p, endp=end_p)
//' # Calculate rolling sums using R code
//' r_sum <- sapply(1:NROW(end_p), function(ep) {
//' colSums(re_turns[(start_p[ep]+1):(end_p[ep]+1), ])
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
  arma::uword num_rows = tseries.n_rows;
  arma::uvec end_pts;
  arma::uvec start_pts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    end_pts = calc_endpoints(num_rows, step, stub);
  } else {
    // Copy end points
    end_pts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    start_pts = calc_startpoints(end_pts, look_back);
    // Old code for start_pts
    // Start points equal to end points lagged by look_back - without adding +1
    // arma::uword num_pts = end_pts.n_elem;
    // arma::uvec start_pts = arma::join_cols(arma::zeros<uvec>(look_back), 
    //                                        end_pts.subvec(0, num_pts - look_back - 1));
  } else {
    // Copy start points
    start_pts = startp;
  }  // end if
  
  // Allocate sums matrix
  arma::uword num_pts = end_pts.n_elem;
  arma::mat sums = arma::zeros<mat>(num_pts, tseries.n_cols);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < num_pts; ep++) {
    // Calculate sums
    if (end_pts(ep) > start_pts(ep)) {
      sums.row(ep) = arma::sum(tseries.rows(start_pts(ep), end_pts(ep)));
    }  // end if
  }  // end for
  
  return sums;
  
  // Old code using arma::cumsum() - sums exclude start_pts
  // Calculate cumulative sums at end points.
  // arma::mat cum_sum = arma::cumsum(tseries, 0);
  // arma::mat cum_start = cum_sum.rows(start_pts);
  // arma::mat cum_end = cum_sum.rows(end_pts);
  
  // Return the differences of the cumulative returns
  // return diff_it(cum_sum, 1, true);
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
//' @param \code{weights} A \emph{column vector} of weights (the default is
//'   \code{weights = NULL}).
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
//'   \code{R} function \cr\code{stats::filter(x=re_turns, filter=weight_s,
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
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("VTI", "IEF")])
//' # Define parameters
//' look_back <- 22
//' # Calculate rolling sums and compare with rutils::roll_sum()
//' c_sum <- HighFreq::roll_sum(re_turns, look_back)
//' r_sum <- rutils::roll_sum(re_turns, look_back)
//' all.equal(c_sum, coredata(r_sum), check.attributes=FALSE)
//' # Calculate rolling sums using R code
//' r_sum <- apply(zoo::coredata(re_turns), 2, cumsum)
//' lag_sum <- rbind(matrix(numeric(2*look_back), nc=2), r_sum[1:(NROW(r_sum) - look_back), ])
//' r_sum <- (r_sum - lag_sum)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' 
//' # Calculate rolling sums at end points
//' stu_b <- 21
//' c_sum <- HighFreq::roll_wsum(re_turns, look_back, stub=stu_b)
//' end_p <- (stu_b + look_back*(0:(NROW(re_turns) %/% look_back)))
//' end_p <- end_p[end_p < NROW(re_turns)]
//' r_sum <- apply(zoo::coredata(re_turns), 2, cumsum)
//' r_sum <- r_sum[end_p+1, ]
//' lag_sum <- rbind(numeric(2), r_sum[1:(NROW(r_sum) - 1), ])
//' r_sum <- (r_sum - lag_sum)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' 
//' # Calculate rolling sums at end points - pass in end_p
//' c_sum <- HighFreq::roll_wsum(re_turns, endp=end_p)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' 
//' # Create exponentially decaying weights
//' weight_s <- exp(-0.2*(1:11))
//' weight_s <- matrix(weight_s/sum(weight_s), nc=1)
//' # Calculate rolling weighted sum
//' c_sum <- HighFreq::roll_wsum(re_turns, weights=weight_s)
//' # Calculate rolling weighted sum using filter()
//' filter_ed <- filter(x=re_turns, filter=weight_s, method="convolution", sides=1)
//' all.equal(c_sum[-(1:11), ], filter_ed[-(1:11), ], check.attributes=FALSE)
//' 
//' # Calculate rolling weighted sums at end points
//' c_sum <- HighFreq::roll_wsum(re_turns, endp=end_p, weights=weight_s)
//' all.equal(c_sum, filter_ed[end_p+1, ], check.attributes=FALSE)
//' 
//' # Create simple weights equal to a 1 value plus zeros
//' weight_s <- matrix(c(1, rep(0, 10)), nc=1)
//' # Calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_wsum(re_turns, weights=weight_s)
//' # Compare with original
//' all.equal(coredata(re_turns), weight_ed, check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_wsum(const arma::mat& tseries,
                    Rcpp::Nullable<Rcpp::IntegerVector> endp = R_NilValue, 
                    arma::uword look_back = 1,
                    Rcpp::Nullable<int> stub = R_NilValue, 
                    Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::mat cum_sum;
  
  // Calculate the rolling sums
  if (weights.isNotNull()) {
    // Coerce weights from Rcpp to Armadillo vector
    arma::vec weights_vec = Rcpp::as<vec>(weights);
    arma::uword num_weights = weights_vec.n_elem;
    // Calculate the weighted averages as convolutions
    cum_sum = arma::conv2(tseries, weights_vec, "full");
    // Copy the warmup period
    // cout << "num_weights = " << num_weights << endl;
    cum_sum.rows(0, num_weights-2) = tseries.rows(0, num_weights-2);
    cum_sum = cum_sum.rows(0, num_rows-1);
    // cout << "cum_sum.n_rows = " << cum_sum.n_rows << endl;
  } else {
    // Calculate the cumulative sum
    cum_sum = arma::cumsum(tseries, 0);
  }  // end if
  
  
  // Declare empty end points
  arma::uvec end_pts;
  // Update end points
  if (endp.isNotNull()) {
    // Copy endp
    end_pts = Rcpp::as<uvec>(endp);
  } else if (stub.isNotNull()) {
    // Calculate end points with stub
    end_pts = arma::regspace<uvec>(Rcpp::as<uword>(stub), look_back, num_rows + look_back);
    end_pts = end_pts.elem(find(end_pts < num_rows));
  }  // end if
  
  
  // Subset the rolling sums according the end points
  if (end_pts.is_empty() && weights.isNotNull()) {
    // Do nothing
    // Return the weighted averages (convolutions) at each point
    // return cum_sum;
  } else if (end_pts.is_empty() && !weights.isNotNull()) {
    // Return unweighted rolling sums at each point
    cum_sum = diff_it(cum_sum, look_back, true);
  } else if (!end_pts.is_empty() && weights.isNotNull()) {
    // Return the weighted averages (convolutions) at end points
    cum_sum = cum_sum.rows(end_pts);
  } else if (!end_pts.is_empty() && !weights.isNotNull()) {
    // Return the unweighted rolling sums at end points
    cum_sum = cum_sum.rows(end_pts);
    cum_sum = diff_it(cum_sum, 1, true);
  }  // end if
  
  return cum_sum;
  
}  // end roll_wsum




////////////////////////////////////////////////////////////
// Functions for rolling aggregations of streaming data
////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////
//' Calculate the rolling mean of streaming \emph{time series} data.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lambda} A \emph{numeric} decay factor.
//'   
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details 
//'   The function \code{run_mean()} calculates the rolling mean of streaming
//'   \emph{time series} data by recursively weighing present and past values
//'   using the decay factor \eqn{\lambda}:
//'   \deqn{
//'     \mu_t = (1-\lambda) p_t + \lambda \mu_{t-1}
//'   }
//'   Where \eqn{\mu_t} is the mean value at time \eqn{t}, and \eqn{p_t} is the
//'   streaming data.
//' 
//'   The value of the decay factor \eqn{\lambda} should be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the rolling mean values have a stronger
//'   dependence on past values.  This is equivalent to a long look-back
//'   interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the rolling mean values have a
//'   weaker dependence on past values.  This is equivalent to a short look-back
//'   interval.
//' 
//'   The above recursive formula is convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//'   The formula is equivalent to a convolution with exponentially decaying
//'   weights, but it's faster.
//' 
//'   The function \code{run_mean()} performs the same calculation as the
//'   standard \code{R} function\cr\code{stats::filter(x=series, filter=lamb_da,
//'   method="recursive")}, but it's several times faster.
//' 
//'   The function \code{run_mean()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{tseries}.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical prices
//' price_s <- zoo::coredata(quantmod::Cl(rutils::etf_env$VTI))
//' # Calculate the rolling means
//' lamb_da <- 0.9
//' means <- HighFreq::run_mean(price_s, lambda=lamb_da)
//' # Calculate rolling means using R code
//' filter_ed <- (1-lamb_da)*filter(price_s, 
//'   filter=lamb_da, init=as.numeric(price_s[1, 1])/(1-lamb_da), 
//'   method="recursive")
//' all.equal(means, unclass(filter_ed), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::run_mean(price_s, lambda=lamb_da),
//'   Rcode=filter(price_s, filter=lamb_da, init=as.numeric(price_s[1, 1])/(1-lamb_da), method="recursive"),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_mean(const arma::mat& tseries, double lambda) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::mat means = arma::zeros<mat>(num_rows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  // Perform loop over rows
  means.row(0) = tseries.row(0);
  for (arma::uword it = 1; it < num_rows; it++) {
    // Calculate the mean as the weighted sum
    means.row(it) = lambda1*tseries.row(it) + lambda*means.row(it-1);
  }  // end for
  
  return means;
  
}  // end run_mean




////////////////////////////////////////////////////////////
//' Calculate the rolling maximum of streaming \emph{time series} data.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lambda} A \emph{numeric} decay factor.
//'   
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details 
//'   The function \code{run_max()} calculates the rolling maximum of streaming
//'   \emph{time series} data by recursively weighing present and past values
//'   using the decay factor \eqn{\lambda}.
//'
//'   It first calculates the rolling mean of streaming data:
//'   \deqn{
//'     \mu_t = (1-\lambda) p_t + \lambda \mu_{t-1}
//'   }
//'   Where \eqn{\mu_t} is the mean value at time \eqn{t}, and \eqn{p_t} is the
//'   streaming data.
//'
//'   It then calculates the rolling maximums of streaming data, \eqn{p^{max}_t}:
//'   \deqn{
//'     p^{max}_t = max(p_t, p^{max}_{t-1}) + (1-\lambda) (\mu_{t-1} - p^{max}_{t-1})
//'   }
//' 
//'   The second term pulls the maximum value down to the mean value, allowing
//'   it to gradually "forget" the maximum value from the more distant past.
//' 
//'   The value of the decay factor \eqn{\lambda} should be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the rolling maximum values have a stronger
//'   dependence on past values.  This is equivalent to a long look-back
//'   interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the rolling maximum values have a
//'   weaker dependence on past values.  This is equivalent to a short look-back
//'   interval.
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
//' price_s <- zoo::coredata(quantmod::Cl(rutils::etf_env$VTI))
//' # Calculate the rolling maximums
//' lamb_da <- 0.9
//' maxs <- HighFreq::run_max(price_s, lambda=lamb_da)
//' # Plot dygraph of VTI prices and rolling maximums
//' da_ta <- cbind(quantmod::Cl(rutils::etf_env$VTI), maxs)
//' colnames(da_ta) <- c("prices", "max")
//' col_names <- colnames(da_ta)
//' dygraphs::dygraph(da_ta, main="VTI Prices and Rolling Maximums") %>%
//'   dySeries(name=col_names[1], label=col_names[1], strokeWidth=2, col="blue") %>%
//'   dySeries(name=col_names[2], label=col_names[2], strokeWidth=2, col="red")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_max(const arma::mat& tseries, double lambda) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::mat maxs = arma::zeros<mat>(num_rows, tseries.n_cols);
  arma::mat means = arma::zeros<mat>(num_rows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  // Perform loop over rows
  means.row(0) = tseries.row(0);
  maxs.row(0) = tseries.row(0);
  for (arma::uword it = 1; it < num_rows; it++) {
    // Calculate the mean as a weighted sum
    means.row(it) = lambda1*tseries.row(it) + lambda*means.row(it-1);
    // Calculate the max from a weighted sum
    maxs.row(it) = arma::max(tseries.row(it), maxs.row(it-1) + lambda1*(means.row(it-1) - maxs.row(it-1)));
  }  // end for
  
  return maxs;
  
}  // end run_max




////////////////////////////////////////////////////////////
//' Calculate the rolling minimum of streaming \emph{time series} data.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lambda} A \emph{numeric} decay factor.
//'   
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details 
//'   The function \code{run_min()} calculates the rolling minimum of streaming
//'   \emph{time series} data by recursively weighing present and past values
//'   using the decay factor \eqn{\lambda}.
//'
//'   It first calculates the rolling mean of streaming data:
//'   \deqn{
//'     \mu_t = (1-\lambda) p_t + \lambda \mu_{t-1}
//'   }
//'   Where \eqn{\mu_t} is the mean value at time \eqn{t}, and \eqn{p_t} is the
//'   streaming data.
//'
//'   It then calculates the rolling minimums of streaming data, \eqn{p^{min}_t}:
//'   \deqn{
//'     p^{min}_t = min(p_t, p^{min}_{t-1}) + (1-\lambda) (\mu_{t-1} - p^{min}_{t-1})
//'   }
//' 
//'   The second term pulls the minimum value up to the mean value, allowing
//'   it to gradually "forget" the minimum value from the more distant past.
//' 
//'   The value of the decay factor \eqn{\lambda} should be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the rolling minimum values have a stronger
//'   dependence on past values.  This is equivalent to a long look-back
//'   interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the rolling minimum values have a
//'   weaker dependence on past values.  This is equivalent to a short look-back
//'   interval.
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
//' price_s <- zoo::coredata(quantmod::Cl(rutils::etf_env$VTI))
//' # Calculate the rolling minimums
//' lamb_da <- 0.9
//' mins <- HighFreq::run_min(price_s, lambda=lamb_da)
//' # Plot dygraph of VTI prices and rolling minimums
//' da_ta <- cbind(quantmod::Cl(rutils::etf_env$VTI), mins)
//' colnames(da_ta) <- c("prices", "min")
//' col_names <- colnames(da_ta)
//' dygraphs::dygraph(da_ta, main="VTI Prices and Rolling Minimums") %>%
//'   dySeries(name=col_names[1], label=col_names[1], strokeWidth=2, col="blue") %>%
//'   dySeries(name=col_names[2], label=col_names[2], strokeWidth=2, col="red")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_min(const arma::mat& tseries, double lambda) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::mat mins = arma::zeros<mat>(num_rows, tseries.n_cols);
  arma::mat means = arma::zeros<mat>(num_rows, tseries.n_cols);
  double lambda1 = 1-lambda;
  
  // Perform loop over rows
  means.row(0) = tseries.row(0);
  mins.row(0) = tseries.row(0);
  for (arma::uword it = 1; it < num_rows; it++) {
    // Calculate the mean as a weighted sum
    means.row(it) = lambda1*tseries.row(it) + lambda*means.row(it-1);
    // Calculate the min from a weighted sum
    mins.row(it) = arma::min(tseries.row(it), mins.row(it-1) + lambda1*(means.row(it-1) - mins.row(it-1)));
  }  // end for
  
  return mins;
  
}  // end run_min




////////////////////////////////////////////////////////////
//' Calculate the rolling variance of streaming \emph{time series} of returns.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of returns.
//' 
//' @param \code{lambda} A \emph{numeric} decay factor.
//'   
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{tseries}.
//'
//' @details 
//'   The function \code{run_var()} calculates the rolling variance of a
//'   streaming \emph{time series} of returns by recursively weighing the
//'   squared present returns with past variance estimates, using the decay
//'   factor \eqn{\lambda}:
//'   \deqn{
//'     \sigma^2_t = (1-\lambda) r^2_t + \lambda \sigma^2_{t-1}
//'   }
//'   Where \eqn{\sigma^2_t} is the variance estimate at time \eqn{t}, and
//'   \eqn{r_t} are the streaming returns data.
//' 
//'   The value of the decay factor \eqn{\lambda} should be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the rolling variance values have a
//'   stronger dependence on past values.  This is equivalent to a long
//'   look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the rolling variance values have a
//'   weaker dependence on past values.  This is equivalent to a short look-back
//'   interval.
//' 
//'   The above formula slightly overestimates the variance because it doesn't
//'   subtract the mean returns.
//' 
//'   The above recursive formula is convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//'   The formula is equivalent to a convolution with exponentially decaying
//'   weights, but it's faster.
//' 
//'   The function \code{run_var()} performs the same calculation as the
//'   standard \code{R} function\cr\code{stats::filter(x=series,
//'   filter=weight_s, method="recursive")}, but it's several times faster.
//' 
//'   The function \code{run_var()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{tseries}.
//'   
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' re_turns <- zoo::coredata(na.omit(rutils::etf_env$re_turns$VTI))
//' # Calculate the rolling variance
//' lamb_da <- 0.9
//' vars <- HighFreq::run_var(re_turns, lambda=lamb_da)
//' # Calculate rolling variance using R code
//' filter_ed <- (1-lamb_da)*filter(re_turns^2, filter=lamb_da, 
//'   init=as.numeric(re_turns[1, 1])^2/(1-lamb_da), 
//'   method="recursive")
//' all.equal(vars, unclass(filter_ed), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::run_var(re_turns, lambda=lamb_da),
//'   Rcode=filter(re_turns^2, filter=lamb_da, init=as.numeric(re_turns[1, 1])^2/(1-lamb_da), method="recursive"),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_var(const arma::mat& tseries, double lambda) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::mat vars = arma::square(tseries);
  double lambda1 = 1-lambda;
  
  // Perform loop over rows
  for (arma::uword it = 1; it < num_rows; it++) {
    // Calculate the variance as the weighted sum of squared returns
    vars.row(it) = lambda1*vars.row(it) + lambda*vars.row(it-1);
  }  // end for
  
  return vars;
  
}  // end run_var




////////////////////////////////////////////////////////////
//' Calculate the rolling covariance of two streaming \emph{time series} of
//' returns.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} with two
//'   columns of returns data.
//' 
//' @param \code{lambda} A \emph{numeric} decay factor.
//'   
//' @return A \emph{matrix} with three columns of data: the covariance and the
//'   variances of the two columns of the argument \code{tseries}.
//'
//' @details 
//'   The function \code{run_covar()} calculates the rolling covariance of 
//'   two streaming \emph{time series} of returns by recursively weighing the
//'   products of their present returns with past covariance estimates, using
//'   the decay factor \eqn{\lambda}:
//'   \deqn{
//'     \sigma^{12}_t = (1-\lambda) r^1_t r^2_t + \lambda \sigma^{12}_{t-1}
//'   }
//'   Where \eqn{\sigma^{12}_t} is the covariance estimate at time \eqn{t}, and
//'   \eqn{r^1_t} and \eqn{r^2_t} are the streaming returns data.
//' 
//'   The value of the decay factor \eqn{\lambda} should be in the range between
//'   \code{0} and \code{1}.  
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the rolling covariance values have a
//'   stronger dependence on past values.  This is equivalent to a long
//'   look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the rolling covariance values have
//'   a weaker dependence on past values.  This is equivalent to a short
//'   look-back interval.
//' 
//'   The above formula slightly overestimates the covariance because it doesn't
//'   subtract the mean returns.
//' 
//'   The above recursive formula is convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//'   The formula is equivalent to a convolution with exponentially decaying
//'   weights, but it's faster.
//' 
//'   The function \code{run_covar()} returns three columns of data: the
//'   covariance and the variances of the two columns of the argument
//'   \code{tseries}.  This allows calculating the rolling correlation.
//' 
//'   The function \code{run_covar()} performs the same calculation as the
//'   standard \code{R} function\cr\code{stats::filter(x=series,
//'   filter=weight_s, method="recursive")}, but it's several times faster.
//' 
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' re_turns <- zoo::coredata(na.omit(rutils::etf_env$re_turns[, c("IEF", "VTI")]))
//' # Calculate the rolling covariance
//' lamb_da <- 0.9
//' covars <- HighFreq::run_covar(re_turns, lambda=lamb_da)
//' # Calculate rolling covariance using R code
//' filter_ed <- (1-lamb_da)*filter(re_turns[, 1]*re_turns[, 2], 
//'   filter=lamb_da, init=as.numeric(re_turns[1, 1]*re_turns[1, 2])/(1-lamb_da), 
//'   method="recursive")
//' all.equal(covars[, 1], unclass(filter_ed), check.attributes=FALSE)
//' # Calculate the rolling correlation
//' correl <- covars[, 1]/sqrt(covars[, 2]*covars[, 3])
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_covar(const arma::mat& tseries, double lambda) {
  
  arma::uword num_rows = tseries.n_rows;
  arma::mat var1 = arma::square(tseries.col(0));
  arma::mat var2 = arma::square(tseries.col(1));
  arma::mat covar = tseries.col(0) % tseries.col(1);
  double lambda1 = 1-lambda;
  
  // Perform loop over rows
  for (arma::uword it = 1; it < num_rows; it++) {
    // Calculate the covariance as the weighted sum of products of returns
    var1.row(it) = lambda1*var1.row(it) + lambda*var1.row(it-1);
    var2.row(it) = lambda1*var2.row(it) + lambda*var2.row(it-1);
    covar.row(it) = lambda1*covar.row(it) + lambda*covar.row(it-1);
  }  // end for
  
  return arma::join_rows(covar, var1, var2);
  
}  // end run_covar




////////////////////////////////////////////////////////////
//' Calculate the z-scores of rolling regressions of streaming \emph{time
//' series} of returns.
//' 
//' @param \code{response} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{design} A \emph{time series} or a \emph{matrix} of design data
//'   (predictor or explanatory data).
//' 
//' @param \code{lambda} A \emph{numeric} decay factor.
//'   
//' @return A \emph{matrix} with the z-scores, betas, and variances of the
//'   design data.
//'
//' @details 
//'   The function \code{run_zscores()} calculates the vectors of \emph{betas}
//'   \eqn{\beta_t} and the residuals \eqn{\epsilon_t} of rolling regressions by
//'   recursively weighing the current estimates with past estimates, using the
//'   decay factor \eqn{\lambda}:
//'   \deqn{
//'     \beta_t = (1-\lambda) \frac{\sigma^{cov}_t}{\sigma^2_t} + \lambda \beta_{t-1}
//'   }
//'   \deqn{
//'     \epsilon_t = (1-\lambda) (r^r_t - \beta_t r^d_t) + \lambda \epsilon_{t-1}
//'   }
//'   Where \eqn{\sigma^{cov}_t} is the vector of covariances between the
//'   response and design returns, at time \eqn{t};
//'   \eqn{\sigma^2_t} is the vector of design variances,
//'   and \eqn{r^r_t} and \eqn{r^d_t} are the streaming returns of the response
//'   and design data.
//' 
//'   The matrices \eqn{\sigma^2}, \eqn{\sigma^{cov}}, \eqn{\beta} have the same
//'   dimensions as the input argument \code{design}.
//'
//'   The above formula is approximate because it doesn't subtract the mean
//'   returns.
//' 
//'   The \emph{z-score} \eqn{z_t} is equal to the residual \eqn{\epsilon_t} divided by
//'   its volatility \eqn{\sigma^{\epsilon}_t}: 
//'   \deqn{
//'     z_t = \frac{\epsilon_t}{\sigma^{\epsilon}_t}
//'   }
//' 
//'   The value of the decay factor \eqn{\lambda} should be in the range between
//'   \code{0} and \code{1}.
//'   If \eqn{\lambda} is close to \code{1} then the decay is weak and past
//'   values have a greater weight, and the rolling \emph{z-score} values have a
//'   stronger dependence on past values.  This is equivalent to a long
//'   look-back interval.
//'   If \eqn{\lambda} is much less than \code{1} then the decay is strong and
//'   past values have a smaller weight, and the rolling \emph{z-score} values have
//'   a weaker dependence on past values.  This is equivalent to a short
//'   look-back interval.
//' 
//'   The above recursive formula is convenient for processing live streaming
//'   data because it doesn't require maintaining a buffer of past data.
//'   The formula is equivalent to a convolution with exponentially decaying
//'   weights, but it's faster.
//' 
//'   The function \code{run_zscores()} returns multiple columns of data. 
//'   If the matrix \code{design} has \code{n} columns then \code{run_zscores()}
//'   returns a matrix with \code{2n+1} columns.  The first column contains the
//'   \emph{z-scores}, and the remaining columns contain the \emph{betas} and
//'   the \emph{variances} of the design data.
//' 
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' res_ponse <- re_turns[, 1]
//' # Design matrix equals VTI and IEF returns
//' de_sign <- re_turns[, -1]
//' # Calculate the running z-scores
//' lamb_da <- 0.9
//' zscores <- HighFreq::run_zscores(response=res_ponse, design=de_sign, lambda=lamb_da)
//' # Plot the running z-scores
//' da_ta <- cbind(cumsum(res_ponse), zscores[, 1])
//' colnames(da_ta) <- c("XLF", "zscores")
//' col_names <- colnames(da_ta)
//' dygraphs::dygraph(da_ta, main="Z-Scores of XLF Versus VTI and IEF") %>%
//'   dyAxis("y", label=col_names[1], independentTicks=TRUE) %>%
//'   dyAxis("y2", label=col_names[2], independentTicks=TRUE) %>%
//'   dySeries(name=col_names[1], axis="y", label=col_names[1], strokeWidth=1, col="blue") %>%
//'   dySeries(name=col_names[2], axis="y2", label=col_names[2], strokeWidth=1, col="red")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat run_zscores(const arma::mat& response, 
                      const arma::mat& design, 
                      double lambda) {
  
  arma::uword num_rows = design.n_rows;
  arma::uword num_cols = design.n_cols;
  // arma::mat var1 = arma::square(tseries.col(0));
  arma::mat vars = arma::square(design);
  arma::mat betas = arma::ones<mat>(num_rows, num_cols);
  arma::mat zscores = arma::ones<mat>(num_rows, 1);
  arma::mat varz = arma::ones<mat>(num_rows, 1);
  double lambda1 = 1-lambda;
  
  // Multiply each column of design by the response.
  arma::mat covars = design;
  covars.each_col() %= response;
  
  // Perform loop over rows
  for (arma::uword it = 1; it < num_rows; it++) {
    // Calculate the z-score as the weighted sum of products of returns.
    // cout << "Calculating vars: " << it << endl;
    vars.row(it) = lambda1*vars.row(it) + lambda*vars.row(it-1);
    // cout << "Calculating covars: " << it << endl;
    covars.row(it) = lambda1*covars.row(it) + lambda*covars.row(it-1);
    // cout << "Calculating betas: " << it << endl;
    betas.row(it) = lambda1*covars.row(it)/vars.row(it) + lambda*betas.row(it-1);
    // cout << "Calculating zscores: " << it << endl;
    zscores.row(it) = lambda1*(response.row(it) - arma::dot(betas.row(it), design.row(it))) + lambda*zscores.row(it-1);
    // Calculate the variance of the z-scores.
    varz.row(it) = lambda1*arma::square(zscores.row(it) - zscores.row(it-1)) + lambda*varz.row(it-1);
  }  // end for
  
  return join_rows(zscores/sqrt(varz), betas, vars);

}  // end run_zscores





////////////////////////////////////////////////////////////
// Functions for statistics
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
// Define C++ enum type for the different methods for regularization,
// calculating variance, skewness, kurtosis, covariance, regression, 
// and matrix inverse.
enum meth_od {moment, least_squares, quantile, nonparametric, regular, rank_sharpe, 
              max_sharpe, max_sharpe_median, min_var, min_varpca, rank, rankrob};

// Map string to C++ enum type for switch statement.
// This is needed because Rcpp can't map C++ enum type to R variable SEXP.
meth_od calc_method(std::string method) {
  if (method == "moment" || method == "m") 
    return meth_od::moment;
  else if (method == "least_squares" || method == "l")
    return meth_od::least_squares;
  else if (method == "quantile" || method == "q")
    return meth_od::quantile;
  else if (method == "nonparametric" || method == "n")
    return meth_od::nonparametric;
  else if (method == "regular")
    return meth_od::regular;
  else if (method == "rank_sharpe")
    return meth_od::rank_sharpe;
  else if (method == "max_sharpe")
    return meth_od::max_sharpe;
  else if (method == "max_sharpe_median")
    return meth_od::max_sharpe_median;
  else if (method == "min_var")
    return meth_od::min_var;
  else if (method == "min_varpca")
    return meth_od::min_varpca;
  else if (method == "rank")
    return meth_od::rank;
  else if (method == "rankrob")
    return meth_od::rankrob;
  else 
    return meth_od::moment;
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
//' @param \code{con_fi} The confidence level for calculating the
//'   quantiles (the default is \code{con_fi = 0.75}).
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
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("XLP", "VTI")])
//' # Calculate the column means in RcppArmadillo
//' HighFreq::calc_mean(re_turns)
//' # Calculate the column means in R
//' sapply(re_turns, mean)
//' # Compare the values
//' all.equal(drop(HighFreq::calc_mean(re_turns)), 
//'   sapply(re_turns, mean), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_mean(re_turns),
//'   Rcode=sapply(re_turns, mean),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the quantile mean (location)
//' HighFreq::calc_mean(re_turns, method="quantile", con_fi=0.9)
//' # Calculate the quantile mean (location) in R
//' colSums(sapply(re_turns, quantile, c(0.9, 0.1), type=5))
//' # Compare the values
//' all.equal(drop(HighFreq::calc_mean(re_turns, method="quantile", con_fi=0.9)), 
//'   colSums(sapply(re_turns, quantile, c(0.9, 0.1), type=5)), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_mean(re_turns, method="quantile", con_fi=0.9),
//'   Rcode=colSums(sapply(re_turns, quantile, c(0.9, 0.1), type=5)),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the column medians in RcppArmadillo
//' HighFreq::calc_mean(re_turns, method="nonparametric")
//' # Calculate the column medians in R
//' sapply(re_turns, median)
//' # Compare the values
//' all.equal(drop(HighFreq::calc_mean(re_turns, method="nonparametric")), 
//'   sapply(re_turns, median), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_mean(re_turns, method="nonparametric"),
//'   Rcode=sapply(re_turns, median),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_mean(const arma::mat& tseries,
                    std::string method = "moment", 
                    double con_fi = 0.75) {
  
  // Switch for the different methods of location
  switch(calc_method(method)) {
  case meth_od::moment: {  // moment
    return arma::mean(tseries);
  }  // end moment
  case meth_od::quantile: {  // quantile
    arma::vec level_s = {1-con_fi, con_fi};
    arma::mat quantile_s = arma::quantile(tseries, level_s);
    return (quantile_s.row(0) + quantile_s.row(1));
  }  // end quantile
  case meth_od::nonparametric: {  // nonparametric
    return arma::median(tseries);
  }  // end nonparametric
  default : {
    cout << "Warning: Invalid method parameter" << endl;
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end default
  }  // end switch
  
}  // end calc_mean




////////////////////////////////////////////////////////////
//' Calculate the variance of a a single-column \emph{time series} or a
//' \emph{vector} using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A single-column \emph{time series} or a \emph{vector}.
//'
//' @return A \emph{numeric} value equal to the variance of the \emph{vector}.
//'
//' @details 
//'   The function \code{calc_var_vec()} calculates the variance of a
//'   \emph{vector} using \code{RcppArmadillo} \code{C++} code, so it's
//'   significantly faster than the \code{R} function \code{var()}.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random returns
//' re_turns <- rnorm(1e6)
//' # Compare calc_var_vec() with standard var()
//' all.equal(HighFreq::calc_var_vec(re_turns), 
//'   var(re_turns))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_var_vec(re_turns),
//'   Rcode=var(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
double calc_var_vec(const arma::vec& tseries) {
  
  return arma::var(tseries);
  
}  // end calc_var_vec




////////////////////////////////////////////////////////////
//' Calculate the dispersion (variance) of the columns of a \emph{time series}
//' or a \emph{matrix} using \code{RcppArmadillo}.
//' 
//' @param \code{tseries} A \emph{time series} or a \emph{matrix} of data.
//'   
//' @param \code{method} A \emph{string} specifying the type of the dispersion model
//'   (the default is \code{method = "moment"} - see Details).
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
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("VTI", "XLF")])
//' # Compare HighFreq::calc_var() with standard var()
//' all.equal(drop(HighFreq::calc_var(re_turns)), 
//'   apply(re_turns, 2, var), check.attributes=FALSE)
//' # Compare HighFreq::calc_var() with matrixStats
//' all.equal(drop(HighFreq::calc_var(re_turns)), 
//'   matrixStats::colVars(re_turns), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with matrixStats and with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_var(re_turns),
//'   matrixStats=matrixStats::colVars(re_turns),
//'   Rcode=apply(re_turns, 2, var),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Compare HighFreq::calc_var() with stats::mad()
//' all.equal(drop(HighFreq::calc_var(re_turns, method="nonparametric")), 
//'   sapply(re_turns, mad), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with stats::mad()
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_var(re_turns, method="nonparametric"),
//'   Rcode=sapply(re_turns, mad),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_var(const arma::mat& tseries,
                   std::string method = "moment", 
                   double con_fi = 0.75) {
  
  // Return zeros if not enough data
  if (tseries.n_rows < 3) {
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end if
  
  // Switch for the different methods of dispersion
  switch(calc_method(method)) {
  case meth_od::moment: {  // moment
    return arma::var(tseries);
  }  // end moment
  case meth_od::quantile: {  // MAD
    double num_cols = tseries.n_cols;
    arma::mat medians = arma::median(tseries);
    arma::mat mads(1, num_cols);
    // Loop over columns of tseries
    for (arma::uword it = 0; it < num_cols; it++) {
      mads.col(it) = arma::median(arma::abs(tseries.col(it) - arma::as_scalar(medians.col(it))));
    }  // end for
    // tseries.each_row() -= arma::median(tseries, 0);
    // return 1.4826*arma::median(arma::abs(tseries), 0);
    return 1.4826*mads;
  }  // end quantile
  case meth_od::nonparametric: {  // nonparametric
    return (arma::mean(tseries) - arma::median(tseries))/arma::stddev(tseries);
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
//' price_s <- na.omit(rutils::etf_env$price_s[, c("XLP", "VTI")])
//' price_s <- log(price_s)
//' # Calculate the daily variance of percentage returns
//' calc_var_ag(price_s, step=1)
//' # Calculate the daily variance using R
//' sapply(rutils::diff_it(price_s), var)
//' # Calculate the variance of returns aggregated over 21 days
//' calc_var_ag(price_s, step=21)
//' # The variance over 21 days is approximately 21 times the daily variance
//' 21*calc_var_ag(price_s, step=1)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_var_ag(const arma::mat& tseries, 
                      arma::uword step = 1) {
  
  if (step == 1)
    // Calculate the variance without aggregations.
    return arma::var(diff_it(tseries, 1, false));
  else {
    // Calculate the number of extra periods that don't fit over num_rows.
    arma::uword num_rows = tseries.n_rows;
    arma::uword remainder = num_rows % step;
    // Allocate aggregations, end points, and variance.
    arma::mat aggs;
    arma::uvec end_p;
    arma::mat var_s(remainder, tseries.n_cols);
    // Perform loop over the stubs
    for (arma::uword stub = 0; stub < remainder; stub++) {
      end_p = calc_endpoints(tseries.n_rows, step, stub);
      // end_p = arma::regspace<uvec>(stub, step, num_rows + step);
      // end_p = end_p.elem(find(end_p < num_rows));
      aggs = tseries.rows(end_p);
      var_s.row(stub) = arma::var(diff_it(aggs, 1, false));
    }  // end for
    return mean(var_s);
  }  // end if
  
}  // end calc_var_ag




////////////////////////////////////////////////////////////
//' Calculate the variance of \emph{OHLC} prices using different price range
//' estimators.
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
//' @param \code{lag_close} A \emph{vector} with the lagged \emph{close} prices
//'   of the \emph{OHLC time series}.  This is an optional argument. (The
//'   default is \code{lag_close = 0}).
//'   
//' @param \code{scale} \emph{Boolean} argument: Should the returns be divided
//'   by the time index, the number of seconds in each period? (The default is
//'   \code{scale = TRUE}).
//'
//' @param \code{in_dex} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument (the default is \code{in_dex = 0}).
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
//'   The input \emph{OHLC time series} \code{ohlc} is assumed to be the log
//'   prices.
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
//'   The optional argument \code{in_dex} is the time index of the \emph{time
//'   series} \code{ohlc}. If the time index is in seconds, then the
//'   differences of the index are equal to the number of seconds in each time
//'   period.  If the time index is in days, then the differences are equal to
//'   the number of days in each time period.
//'   
//'   The optional argument \code{lag_close} are the lagged \emph{close} prices
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
//' oh_lc <- log(HighFreq::SPY)
//' # Extract the time index of SPY prices
//' in_dex <- c(1, diff(xts::.index(oh_lc)))
//' # Calculate the variance of SPY returns, with scaling of the returns
//' HighFreq::calc_var_ohlc(oh_lc, 
//'  method="yang_zhang", scale=TRUE, in_dex=in_dex)
//' # Calculate variance without accounting for overnight jumps
//' HighFreq::calc_var_ohlc(oh_lc, 
//'  method="rogers_satchell", scale=TRUE, in_dex=in_dex)
//' # Calculate the variance without scaling the returns
//' HighFreq::calc_var_ohlc(oh_lc, scale=FALSE)
//' # Calculate the variance by passing in the lagged close prices
//' lag_close <- HighFreq::lag_it(oh_lc[, 4])
//' all.equal(HighFreq::calc_var_ohlc(oh_lc), 
//'   HighFreq::calc_var_ohlc(oh_lc, lag_close=lag_close))
//' # Compare with HighFreq::calc_var_ohlc_r()
//' all.equal(HighFreq::calc_var_ohlc(oh_lc, in_dex=in_dex), 
//'   HighFreq::calc_var_ohlc_r(oh_lc))
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_var_ohlc(oh_lc),
//'   Rcode=HighFreq::calc_var_ohlc_r(oh_lc),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
double calc_var_ohlc(const arma::mat& ohlc, 
                     std::string method = "yang_zhang", 
                     arma::colvec lag_close = 0, 
                     bool scale = true, 
                     arma::colvec in_dex = 0) {
  
  arma::uword num_rows = ohlc.n_rows;
  double coeff = 0.34/(1.34 + (num_rows+1)/(num_rows-1));
  
  if (num_rows < 3) {
    // Return zero if not enough data
    return 0;
  }  // end if
  
  if (!scale || (in_dex.n_rows == 1)) {
    in_dex = arma::ones(num_rows);
    // cout << "ohlc.n_rows = " << num_rows << endl;
    // cout << "in_dex.n_rows = " << in_dex.n_rows << endl;
  }  // end if
  
  // Calculate all the different intra-day and day-over-day returns 
  // (differences of OHLC prices)
  arma::colvec clo_se = ohlc.col(3);
  arma::colvec open_close(clo_se.n_rows);
  if (lag_close.n_rows == 1) {
    open_close = arma::join_cols(clo_se.subvec(0, 0), clo_se.subvec(0, clo_se.n_elem-2));
    open_close = (ohlc.col(0) - open_close)/in_dex;
  } else {
    open_close = (ohlc.col(0) - lag_close)/in_dex;
  }  // end if
  arma::colvec close_open = (clo_se - ohlc.col(0))/in_dex;
  arma::colvec close_high = (clo_se - ohlc.col(1))/in_dex;
  arma::colvec close_low = (clo_se - ohlc.col(2))/in_dex;
  arma::colvec high_low = (ohlc.col(1) - ohlc.col(2))/in_dex;
  arma::colvec high_open = (ohlc.col(1) - ohlc.col(0))/in_dex;
  arma::colvec low_open = (ohlc.col(2) - ohlc.col(0))/in_dex;
  
  if (method == "close") {
    // cout << "Calc method is Close" << endl;
    return arma::var(arma::diff(clo_se));
  } else if (method == "rogers_satchell") {
    // cout << "Calc method is Rogers-Satchell" << endl;
    return -(arma::dot(close_high, high_open) +
             arma::dot(close_low, low_open))/num_rows;
  } else if (method == "garman_klass") {
    // cout << "Calc method is Garman-Klass" << endl;
    return (0.5*arma::dot(high_low, high_low) -
            (2*log(2)-1)*arma::dot(close_open, close_open))/num_rows;
  } else if (method == "garman_klass_yz") {
    // cout << "Calc method is Garman-Klass-YZ" << endl;
    return (0.5*arma::dot(high_low, high_low) -
            (2*log(2)-1)*arma::dot(close_open, close_open))/num_rows + 
            arma::var(open_close);
  } else if (method == "yang_zhang") {
    // cout << "Calc method is Yang-Zhang" << endl;
    return arma::var(open_close) + coeff*arma::var(close_open) +
      (coeff-1)*(arma::dot(close_high, high_open) + 
      arma::dot(close_low, low_open))/num_rows;
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
//' @param \code{lag_close} A \emph{vector} with the lagged \emph{close} prices
//'   of the \emph{OHLC time series}.  This is an optional argument. (The
//'   default is \code{lag_close = 0}).
//'   
//' @param \code{scale} \emph{Boolean} argument: Should the returns be divided
//'   by the time index, the number of seconds in each period? (The default is
//'   \code{scale = TRUE}).
//'
//' @param \code{in_dex} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument (the default is \code{in_dex = 0}).
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
//' oh_lc <- log(rutils::etf_env$VTI)
//' # Calculate the daily variance of percentage returns
//' calc_var_ohlc_ag(oh_lc, step=1)
//' # Calculate the variance of returns aggregated over 21 days
//' calc_var_ohlc_ag(oh_lc, step=21)
//' # The variance over 21 days is approximately 21 times the daily variance
//' 21*calc_var_ohlc_ag(oh_lc, step=1)
//' }
//' 
//' @export
// [[Rcpp::export]]
double calc_var_ohlc_ag(const arma::mat& ohlc,
                        arma::uword step = 1, 
                        std::string method = "yang_zhang", 
                        arma::colvec lag_close = 0, 
                        bool scale = true, 
                        arma::colvec in_dex = 0) {
  
  if (step == 1)
    // Calculate the variance without aggregations.
    return calc_var_ohlc(ohlc, method, lag_close, scale, in_dex);
  else {
    // Calculate the number of extra periods that don't fit over num_rows.
    arma::uword num_rows = ohlc.n_rows;
    arma::uword remainder = num_rows % step;
    // Allocate aggregations, end points, and variance.
    arma::mat aggs;
    arma::uvec end_p;
    arma::mat var_s(remainder, 1);
    // Perform loop over the stubs
    for (arma::uword stub = 0; stub < remainder; stub++) {
      end_p = calc_endpoints(num_rows, step, stub);
      // end_p = arma::regspace<uvec>(stub, step, num_rows + step);
      // end_p = end_p.elem(find(end_p < num_rows));
      // roll_ohlc
      aggs = roll_ohlc(ohlc, end_p);
      var_s.row(stub) = calc_var_ohlc(aggs, method, lag_close, scale, in_dex);
    }  // end for
    return arma::as_scalar(mean(var_s));
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
//' @param \code{con_fi} The confidence level for calculating the
//'   quantiles (the default is \code{con_fi = 0.75}).
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
//' re_turns <- na.omit(rutils::etf_env$re_turns$VTI)
//' # Calculate the moment skewness
//' HighFreq::calc_skew(re_turns)
//' # Calculate the moment skewness in R
//' calc_skewr <- function(x) {
//'   x <- (x-mean(x))
//'   sum(x^3)/var(x)^1.5/NROW(x)
//' }  # end calc_skewr
//' all.equal(HighFreq::calc_skew(re_turns), 
//'   calc_skewr(re_turns), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(re_turns),
//'   Rcode=calc_skewr(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the quantile skewness
//' HighFreq::calc_skew(re_turns, method="quantile", con_fi=0.9)
//' # Calculate the quantile skewness in R
//' calc_skewq <- function(x, a = 0.75) {
//'   	quantile_s <- quantile(x, c(1-a, 0.5, a), type=5)
//'   	(quantile_s[3] + quantile_s[1] - 2*quantile_s[2])/(quantile_s[3] - quantile_s[1])
//' }  # end calc_skewq
//' all.equal(drop(HighFreq::calc_skew(re_turns, method="quantile", con_fi=0.9)), 
//'   calc_skewq(re_turns, a=0.9), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(re_turns, method="quantile"),
//'   Rcode=calc_skewq(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the nonparametric skewness
//' HighFreq::calc_skew(re_turns, method="nonparametric")
//' # Compare HighFreq::calc_skew() with R nonparametric skewness
//' all.equal(drop(HighFreq::calc_skew(re_turns, method="nonparametric")), 
//'   (mean(re_turns)-median(re_turns))/sd(re_turns), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(re_turns, method="nonparametric"),
//'   Rcode=(mean(re_turns)-median(re_turns))/sd(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_skew(const arma::mat& tseries,
                    std::string method = "moment", 
                    double con_fi = 0.75) {
  
  // Return zeros if not enough data
  if (tseries.n_rows < 3) {
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end if
  
  // Switch for the different methods of skew
  switch(calc_method(method)) {
  case meth_od::moment: {  // moment
    double num_rows = tseries.n_rows;
    double num_cols = tseries.n_cols;
    arma::mat means = arma::mean(tseries);
    arma::mat var_s = arma::var(tseries);
    arma::mat skewness(1, num_cols);
    // De-mean the columns of tseries
    // tseries.each_row() -= means;
    // return arma::sum(arma::pow(tseries, 3))/arma::pow(var_s, 1.5)/num_rows;
    for (arma::uword it = 0; it < num_cols; it++) {
      skewness.col(it) = arma::sum(arma::pow(tseries.col(it) - arma::as_scalar(means.col(it)), 3))/arma::pow(var_s.col(it), 1.5)/num_rows;
    }  // end for
    return skewness;
  }  // end moment
  case meth_od::quantile: {  // quantile
    arma::vec level_s = {1-con_fi, 0.5, con_fi};
    arma::mat quantile_s = arma::quantile(tseries, level_s);
    return (quantile_s.row(2) + quantile_s.row(0) - 2*quantile_s.row(1))/(quantile_s.row(2) - quantile_s.row(0));
  }  // end quantile
  case meth_od::nonparametric: {  // nonparametric
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
//' @param \code{con_fi} The confidence level for calculating the
//'   quantiles (the default is \code{con_fi = 0.75}).
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
//' re_turns <- na.omit(rutils::etf_env$re_turns$VTI)
//' # Calculate the moment kurtosis
//' HighFreq::calc_kurtosis(re_turns)
//' # Calculate the moment kurtosis in R
//' calc_kurtr <- function(x) {
//'   x <- (x-mean(x))
//'   sum(x^4)/var(x)^2/NROW(x)
//' }  # end calc_kurtr
//' all.equal(HighFreq::calc_kurtosis(re_turns), 
//'   calc_kurtr(re_turns), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_kurtosis(re_turns),
//'   Rcode=calc_kurtr(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the quantile kurtosis
//' HighFreq::calc_kurtosis(re_turns, method="quantile", con_fi=0.9)
//' # Calculate the quantile kurtosis in R
//' calc_kurtq <- function(x, a=0.9) {
//'   	quantile_s <- quantile(x, c(1-a, 0.25, 0.75, a), type=5)
//'   	(quantile_s[4] - quantile_s[1])/(quantile_s[3] - quantile_s[2])
//' }  # end calc_kurtq
//' all.equal(drop(HighFreq::calc_kurtosis(re_turns, method="quantile", con_fi=0.9)), 
//'   calc_kurtq(re_turns, a=0.9), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_kurtosis(re_turns, method="quantile"),
//'   Rcode=calc_kurtq(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the nonparametric kurtosis
//' HighFreq::calc_kurtosis(re_turns, method="nonparametric")
//' # Compare HighFreq::calc_kurtosis() with R nonparametric kurtosis
//' all.equal(drop(HighFreq::calc_kurtosis(re_turns, method="nonparametric")), 
//'   (mean(re_turns)-median(re_turns))/sd(re_turns), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_kurtosis(re_turns, method="nonparametric"),
//'   Rcode=(mean(re_turns)-median(re_turns))/sd(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_kurtosis(const arma::mat& tseries,
                        std::string method = "moment", 
                        double con_fi = 0.75) {
  
  // Return zeros if not enough data
  if (tseries.n_rows < 3) {
    return arma::zeros<rowvec>(tseries.n_cols);
  }  // end if
  
  // Switch for the different methods of kurtosis
  switch(calc_method(method)) {
  case meth_od::moment: {  // Fourth moment
    double num_rows = tseries.n_rows;
    double num_cols = tseries.n_cols;
    arma::mat means = arma::mean(tseries);
    arma::mat var_s = arma::var(tseries);
    arma::mat kurtosis(1, num_cols);
    // Don't de-mean the columns of tseries because that requires copying the matrix of data, so it's time-consuming
    // Loop over columns of tseries
    for (arma::uword it = 0; it < num_cols; it++) {
      kurtosis.col(it) = arma::sum(arma::pow(tseries.col(it) - arma::as_scalar(means.col(it)), 4))/arma::pow(var_s.col(it), 2)/num_rows;
    }  // end for
    // tseries.each_row() -= means;
    // return arma::sum(arma::pow(tseries, 4))/arma::pow(var_s, 2)/num_rows;
    return kurtosis;
  }  // end moment
  case meth_od::quantile: {  // quantile
    arma::vec level_s = {1-con_fi, 0.25, 0.75, con_fi};
    arma::mat quantile_s = arma::quantile(tseries, level_s);
    return (quantile_s.row(3) - quantile_s.row(0))/(quantile_s.row(2) - quantile_s.row(1));
  }  // end quantile
  case meth_od::nonparametric: {  // nonparametric
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
//' price_s <- na.omit(rutils::etf_env$price_s[, c("XLP", "VTI")])
//' price_s <- log(price_s)
//' # Calculate the Hurst exponent from 21 day aggregations
//' calc_hurst(price_s, step=21)
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
//' @param \code{lag_close} A \emph{vector} with the lagged \emph{close} prices
//'   of the \emph{OHLC time series}.  This is an optional argument. (The
//'   default is \code{lag_close = 0}).
//'   
//' @param \code{scale} \emph{Boolean} argument: Should the returns be divided
//'   by the time index, the number of seconds in each period? (The default is
//'   \code{scale = TRUE}).
//'
//' @param \code{in_dex} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument (the default is \code{in_dex = 0}).
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
//' oh_lc <- log(rutils::etf_env$VTI)
//' # Calculate the Hurst exponent from 21 day aggregations
//' calc_hurst_ohlc(oh_lc, step=21)
//' }
//' 
//' @export
// [[Rcpp::export]]
double calc_hurst_ohlc(const arma::mat& ohlc,
                       arma::uword step = 1, 
                       std::string method = "yang_zhang", 
                       arma::colvec lag_close = 0, 
                       bool scale = true, 
                       arma::colvec in_dex = 0) {
  
  return 0.5*log(calc_var_ohlc_ag(ohlc, step, method, lag_close, scale, in_dex)/
                 calc_var_ohlc_ag(ohlc, 1, method, lag_close, scale, in_dex))/log(step);
  
}  // end calc_hurst_ohlc




////////////////////////////////////////////////////////////
//' Perform multivariate linear regression using least squares and return a
//' named list of regression coefficients, their t-values, and p-values.
//' 
//' @param \code{response} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{design} A \emph{time series} or a \emph{matrix} of design data
//'   (predictor or explanatory data).
//' 
//' @return A named list with three elements: a \emph{matrix} of coefficients
//'   (named \emph{"coefficients"}), the \emph{z-score} of the last residual
//'   (named \emph{"z_score"}), and a \emph{vector} with the R-squared and
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
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' res_ponse <- re_turns[, 1]
//' # Design matrix equals VTI and IEF returns
//' de_sign <- re_turns[, -1]
//' # Perform multivariate regression using lm()
//' reg_model <- lm(res_ponse ~ de_sign)
//' sum_mary <- summary(reg_model)
//' # Perform multivariate regression using calc_lm()
//' reg_arma <- HighFreq::calc_lm(response=res_ponse, design=de_sign)
//' # Compare the outputs of both functions
//' all.equal(reg_arma$coefficients[, "coeff"], unname(coef(reg_model)))
//' all.equal(unname(reg_arma$coefficients), unname(sum_mary$coefficients))
//' all.equal(unname(reg_arma$stats), c(sum_mary$r.squared, unname(sum_mary$fstatistic[1])))
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_lm(response=res_ponse, design=de_sign),
//'   Rcode=lm(res_ponse ~ de_sign),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List calc_lm(const arma::vec& response, const arma::mat& design) {
  
  // Add column for intercept to explanatory matrix
  arma::uword num_rows = design.n_rows;
  arma::mat design_p = join_rows(ones(num_rows), design);
  arma::uword num_cols = design_p.n_cols;
  arma::uword deg_free = (num_rows - num_cols);
  
  // Calculate alpha and beta coefficients for the model response ~ design
  arma::colvec coeff = arma::solve(design_p, response);
  // Calculate residuals
  arma::colvec resid_uals = response - design_p*coeff;
  
  // Calculate TSS, RSS, and ESS
  double tot_sumsq = (num_rows-1)*arma::var(response);
  double res_sumsq = arma::dot(resid_uals, resid_uals);
  double exp_sumsq = tot_sumsq - res_sumsq;
  
  // Calculate R-squared and F-statistic
  double r_squared = exp_sumsq/tot_sumsq;
  double f_stat = (exp_sumsq*deg_free)/(res_sumsq*(num_cols-1));
  // arma::rowvec stat_s=join_horiz(r_squared, f_stat);
  Rcpp::NumericVector stat_s(2);
  stat_s(0) = r_squared;
  stat_s(1) = f_stat;
  stat_s.attr("names") = Rcpp::CharacterVector::create("R-squared", "F-statistic");
  
  // Calculate standard errors of beta coefficients
  arma::colvec std_err = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(design_p)*design_p)));
  // Calculate t-values and p-values of beta coefficients
  arma::colvec t_vals = coeff/std_err;
  arma::colvec p_vals = 2*Rcpp::pt(-abs(wrap(t_vals)), deg_free);
  Rcpp::NumericMatrix coeff_mat = Rcpp::wrap(join_rows(join_rows(join_rows(coeff, std_err), t_vals), p_vals));
  Rcpp::colnames(coeff_mat) = Rcpp::CharacterVector::create("coeff", "std_err", "tvals", "pvals");
  
  return Rcpp::List::create(Named("coefficients") = coeff_mat,
                            // Named("residuals") = resid_uals,
                            Named("z_score") = resid_uals(num_rows-1)/arma::stddev(resid_uals),
                            Named("stats") = stat_s);
  
}  // end calc_lm



////////////////////////////////////////////////////////////
//' Perform multivariate regression using different methods, and return a vector
//' of regression coefficients, their t-values, and the last residual z-score.
//' 
//' @param \code{response} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{design} A \emph{time series} or a \emph{matrix} of design data
//'   (predictor or explanatory data).
//' 
//' @param \code{method} A \emph{string} specifying the type of the regression
//'   model the default is \code{method = "least_squares"} - see Details).
//'   
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small singular values in order to regularize the inverse of the
//'   \code{design} matrix (the default is \code{0.001}).
//'   
//' @param \code{eigen_max} An \emph{integer} equal to the number of singular
//'   values used for calculating the regularized inverse of the \code{design}
//'   matrix (the default is \code{0} - equivalent to \code{eigen_max} equal to
//'   the number of columns of \code{design}).
//'   
//' @param \code{con_fi} The confidence level for calculating the
//'   quantiles (the default is \code{con_fi = 0.75}).
//'
//' @param \code{alpha} The shrinkage intensity between \code{0} and \code{1}.
//'   (the default is \code{0}).
//' 
//' @return A vector with the regression coefficients, their t-values, and the
//'   last residual z-score.
//'
//' @details 
//'   The function \code{calc_reg()} performs multivariate regression using
//'   different methods, and returns a vector of regression coefficients, their
//'   t-values, and the last residual z-score.
//' 
//'   The length of the return vector depends on the number of columns of
//'   \code{design}.
//'   The number of regression coefficients is equal to the number of columns of
//'   \code{design} plus \code{1}.  The number of t-values is equal to the
//'   number of coefficients.  And there is only \code{1} z-score.
//'   So if the number of columns of \code{design} is equal to \code{n}, then
//'   the return vector will have \code{2n+3} elements.
//' 
//'   For example, if the design matrix has \code{2} columns of data, then
//'   \code{calc_reg()} returns a vector with \code{7} elements: \code{3}
//'   regression coefficients (including the intercept coefficient), \code{3}
//'   corresponding t-values, and \code{1} z-score.
//'
//'   If \code{method = "least_squares"} (the default) then it performs the
//'   standard least squares regression, the same as the function
//'   \code{calc_reg()}, and the function \code{lm()} from package \emph{stats}.
//'   It uses \code{RcppArmadillo} \code{C++} code so it's several times faster
//'   than \code{lm()}.
//'
//'   If \code{method = "regular"} then it performs regularized regression.  It
//'   calculates the regularized inverse of the \code{design} matrix from its
//'   singular value decomposition.  It applies dimension regularization by
//'   selecting only the largest singular values equal in number to
//'   \code{eigen_max}.
//'   
//'   If \code{method = "quantile"} then it performs quantile regression (not
//'   implemented yet).
//'
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' res_ponse <- re_turns[, 1]
//' # Design matrix equals VTI and IEF returns
//' de_sign <- re_turns[, -1]
//' # Perform multivariate regression using lm()
//' reg_model <- lm(res_ponse ~ de_sign)
//' sum_mary <- summary(reg_model)
//' co_eff <- sum_mary$coefficients
//' # Perform multivariate regression using calc_reg()
//' reg_arma <- drop(HighFreq::calc_reg(response=res_ponse, design=de_sign))
//' # Compare the outputs of both functions
//' all.equal(reg_arma[1:(2*(1+NCOL(de_sign)))], 
//'   c(co_eff[, "Estimate"], co_eff[, "t value"]), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_reg(response=res_ponse, design=de_sign),
//'   Rcode=lm(res_ponse ~ de_sign),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::colvec calc_reg(const arma::vec& response, 
                      const arma::mat& design,
                      std::string method = "least_squares",
                      double eigen_thresh = 0.001,
                      arma::uword eigen_max = 0,
                      double con_fi = 0.1,
                      double alpha = 0.0) {
  
  // Add column for intercept to explanatory matrix
  arma::uword num_rows = design.n_rows;
  arma::mat design_p = join_rows(ones(num_rows), design);
  arma::uword num_cols = design_p.n_cols;
  arma::uword deg_free = (num_rows - num_cols);
  arma::colvec coeff(num_cols, fill::zeros);
  arma::colvec reg_data(2*num_cols+1, fill::zeros);
  
  // Switch for the different methods for weights
  switch(calc_method(method)) {
  case meth_od::least_squares: {
    // Calculate regression coefficients for the model response ~ design
    coeff = arma::solve(design_p, response);
    break;
  }  // end least_squares
  case meth_od::regular: {
    // Calculate regularized regression coefficients
    coeff = calc_inv(design_p, eigen_thresh, eigen_max)*response;
    break;
  }  // end regular
  case meth_od::quantile: {
    // Not implemented yet
    break;
  }  // end quantile
  default : {
    cout << "Warning: Invalid method parameter: " << method << endl;
    return reg_data;
  }  // end default
  }  // end switch
  
  // Calculate residuals
  arma::colvec resid_uals = response - design_p*coeff;
  
  // Calculate TSS, RSS, and ESS
  // double tot_sumsq = (num_rows-1)*arma::var(response);
  double res_sumsq = arma::dot(resid_uals, resid_uals);
  // double exp_sumsq = tot_sumsq - res_sumsq;
  
  // Calculate standard errors of beta coefficients
  arma::colvec std_err = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(design_p)*design_p)));
  // Calculate t-values and p-values of beta coefficients
  arma::colvec t_vals = coeff/std_err;
  
  // Calculate z-score
  double z_score = resid_uals(num_rows-1)/arma::stddev(resid_uals);
  
  // Combine regression data
  reg_data.subvec(0, num_cols-1) = coeff;
  reg_data.subvec(num_cols, 2*num_cols-1) = t_vals;
  reg_data(2*num_cols) = z_score;
  
  return reg_data;
  
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
//' re_turns <- na.omit(rutils::etf_env$re_turns$VTI)
//' # Calculate the rolling means at 25 day end points, with a 75 day look-back
//' means <- HighFreq::roll_mean(re_turns, step=25, look_back=3)
//' # Compare the mean estimates over 11-period look-back intervals
//' all.equal(HighFreq::roll_mean(re_turns, look_back=11)[-(1:10), ], 
//'   drop(RcppRoll::roll_mean(re_turns, n=11)), check.attributes=FALSE)
//' # Define end points and start points
//' end_p <- HighFreq::calc_endpoints(NROW(re_turns), step=25)
//' start_p <- HighFreq::calc_startpoints(end_p, look_back=3)
//' # Calculate the rolling means using RcppArmadillo
//' means <- HighFreq::roll_mean(re_turns, startp=start_p, endp=end_p)
//' # Calculate the rolling medians using RcppArmadillo
//' medianscpp <- HighFreq::roll_mean(re_turns, startp=start_p, endp=end_p, method="nonparametric")
//' # Calculate the rolling medians using R
//' medians = sapply(1:NROW(end_p), function(i) {
//'   median(re_turns[start_p[i]:end_p[i] + 1])
//' })  # end sapply
//' all.equal(medians, drop(medianscpp))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_mean(re_turns, startp=start_p, endp=end_p, method="nonparametric"),
//'   Rcode=sapply(1:NROW(end_p), function(i) {median(re_turns[start_p[i]:end_p[i] + 1])}),
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
                    double con_fi = 0.75) {
  
  // Allocate end points
  arma::uword num_rows = tseries.n_rows;
  arma::uvec end_pts;
  arma::uvec start_pts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    end_pts = calc_endpoints(num_rows, step, stub);
  } else {
    // Copy end points
    end_pts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    start_pts = calc_startpoints(end_pts, look_back);
  } else {
    // Copy start points
    start_pts = startp;
  }  // end if
  
  
  // Allocate mean matrix
  arma::uword num_pts = end_pts.n_elem;
  arma::mat means = arma::zeros<mat>(num_pts, tseries.n_cols);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < num_pts; ep++) {
    // Calculate means
    if (end_pts(ep) > start_pts(ep)) {
      means.row(ep) = calc_mean(tseries.rows(start_pts(ep), end_pts(ep)), method, con_fi);
    }  // end if
  }  // end for
  
  return means;
  
}  // end roll_mean




////////////////////////////////////////////////////////////
//' Calculate a \emph{vector} of variance estimates over a rolling look-back
//' interval for a single-column \emph{time series} or a \emph{column vector}, using
//' \code{RcppArmadillo}.
//'
//' @param \code{tseries} A single-column \emph{time series} or a \emph{column vector}.
//' 
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of \emph{vector} elements used for calculating a single variance
//'   estimate (the default is \code{look_back = 1}).
//'
//' @return A \emph{column vector} with the same number of elements as the input
//'   argument \code{tseries}.
//'
//' @details 
//'   The function \code{roll_var_vec()} calculates a \emph{vector} of variance
//'   estimates over a rolling look-back interval for a single-column \emph{time
//'   series} or a \emph{column vector}, using \code{RcppArmadillo} \code{C++} code.
//'   
//'   The function \code{roll_var_vec()} uses an expanding look-back interval in
//'   the initial warmup period, to calculate the same number of elements as the
//'   input argument \code{tseries}.
//'
//'   The function \code{roll_var_vec()} performs the same calculation as the
//'   function \code{roll_var()} from package
//'   \href{https://cran.r-project.org/web/packages/RcppRoll/index.html}{RcppRoll},
//'   but it's several times faster because it uses \code{RcppArmadillo}
//'   \code{C++} code.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random returns
//' re_turns <- rnorm(1e6)
//' # Compare the variance estimates over 11-period look-back intervals
//' all.equal(drop(HighFreq::roll_var_vec(re_turns, look_back=11))[-(1:10)], 
//'   RcppRoll::roll_var(re_turns, n=11))
//' # Compare the speed of RcppArmadillo with RcppRoll
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_var_vec(re_turns, look_back=11),
//'   RcppRoll=RcppRoll::roll_var(re_turns, n=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::vec roll_var_vec(const arma::vec& tseries, arma::uword look_back = 1) {
  
  arma::uword length = tseries.n_elem;
  arma::vec var_vec = arma::zeros<vec>(length);
  
  // Warmup period
  for (arma::uword it = 1; it < look_back; it++) {
    var_vec(it) = arma::var(tseries.subvec(0, it));
  }  // end for
  
  // Remaining period
  for (arma::uword it = look_back; it < length; it++) {
    var_vec(it) = arma::var(tseries.subvec(it-look_back+1, it));
  }  // end for
  
  return var_vec;
  
}  // end roll_var_vec




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
//' re_turns <- na.omit(rutils::etf_env$re_turns$VTI)
//' # Calculate the rolling variance at 25 day end points, with a 75 day look-back
//' vari_ance <- HighFreq::roll_var(re_turns, step=25, look_back=3)
//' # Compare the variance estimates over 11-period look-back intervals
//' all.equal(HighFreq::roll_var(re_turns, look_back=11)[-(1:10), ], 
//'   drop(RcppRoll::roll_var(re_turns, n=11)), check.attributes=FALSE)
//' # Compare the speed of HighFreq::roll_var() with RcppRoll::roll_var()
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_var(re_turns, look_back=11),
//'   RcppRoll=RcppRoll::roll_var(re_turns, n=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Compare the speed of HighFreq::roll_var() with TTR::runMAD()
//' summary(microbenchmark(
//'     Rcpp=HighFreq::roll_var(re_turns, look_back=11, method="quantile"),
//'     TTR=TTR::runMAD(re_turns, n = 11),
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
                   double con_fi = 0.75) {
  
  // Allocate end points
  arma::uword num_rows = tseries.n_rows;
  arma::uvec end_pts;
  arma::uvec start_pts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    end_pts = calc_endpoints(num_rows, step, stub);
  } else {
    // Copy end points
    end_pts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    start_pts = calc_startpoints(end_pts, look_back);
  } else {
    // Copy start points
    start_pts = startp;
  }  // end if
  
  
  // Allocate variance matrix
  arma::uword num_pts = end_pts.n_elem;
  arma::mat variance = arma::zeros<mat>(num_pts, tseries.n_cols);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < num_pts; ep++) {
    // Calculate variance
    if (end_pts(ep) > start_pts(ep)) {
      variance.row(ep) = calc_var(tseries.rows(start_pts(ep), end_pts(ep)), method, con_fi);
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
//' @param \code{in_dex} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument (the default is \code{in_dex=0}).
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
//'   The optional argument \code{in_dex} is the time index of the \emph{time
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
//' oh_lc <- log(HighFreq::SPY)
//' # Extract the time index of SPY prices
//' in_dex <- c(1, diff(xts::.index(oh_lc)))
//' # Rolling variance at minutely end points, with a 21 minute look-back
//' var_rolling <- HighFreq::roll_var_ohlc(oh_lc, 
//'                               step=1, look_back=21, 
//'                               method="yang_zhang", 
//'                               in_dex=in_dex, scale=TRUE)
//' # Daily OHLC prices
//' oh_lc <- rutils::etf_env$VTI
//' in_dex <- c(1, diff(xts::.index(oh_lc)))
//' # Rolling variance at 5 day end points, with a 20 day look-back (20=4*5)
//' var_rolling <- HighFreq::roll_var_ohlc(oh_lc, 
//'                               step=5, look_back=4, 
//'                               method="yang_zhang", 
//'                               in_dex=in_dex, scale=TRUE)
//' # Same calculation in R
//' n_rows <- NROW(oh_lc)
//' lag_close = HighFreq::lag_it(oh_lc[, 4])
//' end_p <- drop(HighFreq::calc_endpoints(n_rows, 3)) + 1
//' start_p <- drop(HighFreq::calc_startpoints(end_p, 2))
//' n_pts <- NROW(end_p)
//' var_rollingr <- sapply(2:n_pts, function(it) {
//'   ran_ge <- start_p[it]:end_p[it]
//'   sub_ohlc = oh_lc[ran_ge, ]
//'   sub_close = lag_close[ran_ge]
//'   sub_index = in_dex[ran_ge]
//'   HighFreq::calc_var_ohlc(sub_ohlc, lag_close=sub_close, scale=TRUE, in_dex=sub_index)
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
                        arma::colvec in_dex = 0) {
  
  // Allocate end points
  arma::uword num_rows = ohlc.n_rows;
  arma::uvec end_pts;
  arma::uvec start_pts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    end_pts = calc_endpoints(num_rows, step, stub);
  } else {
    // Copy end points
    end_pts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    start_pts = calc_startpoints(end_pts, look_back);
  } else {
    // Copy start points
    start_pts = startp;
  }  // end if
  
  // Allocate variance matrix
  arma::uword num_pts = end_pts.n_elem;
  arma::vec variance = arma::zeros<vec>(num_pts);
  
  // Extract OHLC close prices
  arma::colvec clo_se = ohlc.col(3);
  arma::colvec lag_close = lag_it(clo_se);
  
  // Set the time index to 1 if scale = FALSE
  if (!scale || (in_dex.n_rows == 1)) {
    in_dex = arma::ones(num_rows);
  }  // end if
  
  // Define data subsets over look-back intervals
  arma::mat sub_ohlc;
  arma::colvec sub_close;
  arma::colvec sub_index;
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < num_pts; ep++) {
    if (end_pts(ep) > start_pts(ep)) {
      sub_ohlc = ohlc.rows(start_pts(ep), end_pts(ep));
      sub_close = lag_close.rows(start_pts(ep), end_pts(ep));
      sub_index = in_dex.subvec(start_pts(ep), end_pts(ep));
      // Calculate variance
      variance.row(ep) = calc_var_ohlc(sub_ohlc, method, sub_close, scale, sub_index);
    }  // end if
  }  // end for
  
  // Old code below
  
  // Warmup period
  // for (arma::uword it = 1; it < look_back; it++) {
  //   arma::mat sub_ohlc = ohlc.rows(0, it);
  //   arma::colvec sub_close = lag_close.rows(0, it);
  //   arma::colvec sub_index = in_dex.subvec(0, it);
  //   variance(it) = calc_var_ohlc(sub_ohlc, method, sub_close, scale, sub_index);
  // }  // end for
  
  // Remaining period
  // for (arma::uword it = look_back; it < num_rows; it++) {
  //   arma::mat sub_ohlc = ohlc.rows(it-look_back+1, it);
  //   arma::colvec sub_close = lag_close.rows(it-look_back+1, it);
  //   arma::colvec sub_index = in_dex.subvec(it-look_back+1, it);
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
//' @param \code{con_fi} The confidence level for calculating the
//'   quantiles (the default is \code{con_fi = 0.75}).
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
//' re_turns <- na.omit(rutils::etf_env$re_turns$VTI)
//' # Define end points and start points
//' end_p <- 1 + HighFreq::calc_endpoints(NROW(re_turns), step=25)
//' start_p <- HighFreq::calc_startpoints(end_p, look_back=3)
//' # Calculate the rolling skewness at 25 day end points, with a 75 day look-back
//' skew_ness <- HighFreq::roll_skew(re_turns, step=25, look_back=3)
//' # Calculate the rolling skewness using R code
//' skew_r <- sapply(1:NROW(end_p), function(it) {
//'   HighFreq::calc_skew(re_turns[start_p[it]:end_p[it], ])
//' })  # end sapply
//' # Compare the skewness estimates
//' all.equal(drop(skew_ness), skew_r, check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_skew(re_turns, step=25, look_back=3),
//'   Rcode=sapply(1:NROW(end_p), function(it) {
//'     HighFreq::calc_skew(re_turns[start_p[it]:end_p[it], ])
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
                    double con_fi = 0.75) {
  
  // Allocate end points
  arma::uword num_rows = tseries.n_rows;
  arma::uvec end_pts;
  arma::uvec start_pts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    end_pts = calc_endpoints(num_rows, step, stub);
  } else {
    // Copy end points
    end_pts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    start_pts = calc_startpoints(end_pts, look_back);
  } else {
    // Copy start points
    start_pts = startp;
  }  // end if
  
  // Allocate skewness matrix
  arma::uword num_pts = end_pts.n_elem;
  arma::mat skew_ness = arma::zeros<mat>(num_pts, tseries.n_cols);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < num_pts; ep++) {
    // Calculate skewness
    if (end_pts(ep) > start_pts(ep)) {
      skew_ness.row(ep) = calc_skew(tseries.rows(start_pts(ep), end_pts(ep)), method, con_fi);
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
//' @param \code{con_fi} The confidence level for calculating the
//'   quantiles (the default is \code{con_fi = 0.75}).
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
//'   periods. It calculates the end points along the rows of \code{design}
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
//' re_turns <- na.omit(rutils::etf_env$re_turns$VTI)
//' # Define end points and start points
//' end_p <- 1 + HighFreq::calc_endpoints(NROW(re_turns), step=25)
//' start_p <- HighFreq::calc_startpoints(end_p, look_back=3)
//' # Calculate the rolling kurtosis at 25 day end points, with a 75 day look-back
//' kurto_sis <- HighFreq::roll_kurtosis(re_turns, step=25, look_back=3)
//' # Calculate the rolling kurtosis using R code
//' kurt_r <- sapply(1:NROW(end_p), function(it) {
//'   HighFreq::calc_kurtosis(re_turns[start_p[it]:end_p[it], ])
//' })  # end sapply
//' # Compare the kurtosis estimates
//' all.equal(drop(kurto_sis), kurt_r, check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_kurtosis(re_turns, step=25, look_back=3),
//'   Rcode=sapply(1:NROW(end_p), function(it) {
//'     HighFreq::calc_kurtosis(re_turns[start_p[it]:end_p[it], ])
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
                        double con_fi = 0.75) {
  
  // Allocate end points
  arma::uword num_rows = tseries.n_rows;
  arma::uvec end_pts;
  arma::uvec start_pts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    end_pts = calc_endpoints(num_rows, step, stub);
  } else {
    // Copy end points
    end_pts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    start_pts = calc_startpoints(end_pts, look_back);
  } else {
    // Copy start points
    start_pts = startp;
  }  // end if
  
  // Allocate kurtosis matrix
  arma::uword num_pts = end_pts.n_elem;
  arma::mat kurto_sis = arma::zeros<mat>(num_pts, tseries.n_cols);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < num_pts; ep++) {
    // Calculate kurtosis
    if (end_pts(ep) > start_pts(ep)) {
      kurto_sis.row(ep) = calc_kurtosis(tseries.rows(start_pts(ep), end_pts(ep)), method, con_fi);
    }  // end if
  }  // end for
  
  return kurto_sis;
  
}  // end roll_kurtosis



////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of regression coefficients, their t-values, and
//' z-scores, at the end points of the design matrix.
//' 
//' @param \code{response} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{design} A \emph{time series} or a \emph{matrix} of design data
//'   (predictor or explanatory data).
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
//' @param \code{method} A \emph{string} specifying the type of the regression
//'   model the default is \code{method = "least_squares"} - see Details).
//'   
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small singular values in order to regularize the inverse of the
//'   \code{design} matrix (the default is \code{0.001}).
//'   
//' @param \code{eigen_max} An \emph{integer} equal to the number of singular
//'   values used for calculating the regularized inverse of the \code{design}
//'   matrix (the default is \code{0} - equivalent to \code{eigen_max} equal to
//'   the number of columns of \code{design}).
//'   
//' @param \code{con_fi} The confidence level for calculating the
//'   quantiles (the default is \code{con_fi = 0.75}).
//'
//' @param \code{alpha} The shrinkage intensity between \code{0} and \code{1}.
//'   (the default is \code{0}).
//' 
//' @return A \emph{matrix} with the same number of rows as \code{design}, and a
//'   number of columns equal to \code{2n+3}, where \code{n} is the number of
//'   columns of \code{design}.
//'
//' @details 
//'   The function \code{roll_reg()} calculates a \emph{matrix} of regression
//'   coefficients, their t-values, and z-scores at the end points of the design
//'   matrix.
//'   
//'   The function \code{roll_reg()} performs a loop over the end points, and at
//'   each end point it subsets the time series \code{design} over a look-back
//'   interval equal to \code{look_back} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_reg()}, which
//'   calculates the regression data.
//'   
//'   If the arguments \code{endp} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{design}
//'   using the function \code{calc_endpoints()}, with the number of time
//'   periods between the end points equal to \code{step} time periods.
//'   
//'   For example, the rolling regression at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{step = 25} and \code{look_back = 3}.
//'
//' @examples
//' \dontrun{
//' # Calculate historical returns
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("XLP", "VTI")])
//' # Define monthly end points and start points
//' end_p <- xts::endpoints(re_turns, on="months")[-1]
//' look_back <- 12
//' start_p <- c(rep(1, look_back), end_p[1:(NROW(end_p)-look_back)])
//' # Calculate rolling betas using RcppArmadillo
//' reg_stats <- HighFreq::roll_reg(response=re_turns[, 1], design=re_turns[, 2], endp=(end_p-1), startp=(start_p-1))
//' beta_s <- reg_stats[, 2]
//' # Calculate rolling betas in R
//' betas_r <- sapply(1:NROW(end_p), FUN=function(ep) {
//'   da_ta <- re_turns[start_p[ep]:end_p[ep], ]
//'   drop(cov(da_ta[, 1], da_ta[, 2])/var(da_ta[, 2]))
//' })  # end sapply
//' # Compare the outputs of both functions
//' all.equal(beta_s, betas_r, check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_reg(const arma::vec& response, 
                   const arma::mat& design, 
                   arma::uvec startp = 0, 
                   arma::uvec endp = 0, 
                   arma::uword step = 1, 
                   arma::uword look_back = 1, 
                   arma::uword stub = 0,
                   std::string method = "least_squares",
                   double eigen_thresh = 0.001,
                   arma::uword eigen_max = 0,
                   double con_fi = 0.1,
                   double alpha = 0.0) {
  
  // Allocate end points
  arma::uword num_rows = design.n_rows;
  arma::uvec end_pts;
  arma::uvec start_pts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    end_pts = calc_endpoints(num_rows, step, stub);
  } else {
    // Copy end points
    end_pts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    start_pts = calc_startpoints(end_pts, look_back);
  } else {
    // Copy start points
    start_pts = startp;
  }  // end if
  
  
  // Allocate regression matrix
  arma::vec sub_response;
  arma::mat sub_design;
  arma::colvec reg_data;
  arma::uword num_cols = design.n_cols;
  arma::uword num_pts = end_pts.n_elem;
  arma::mat reg_stats(num_pts, (2*(num_cols + 1) + 1), fill::zeros);
  
  
  // Perform loop over the end_pts
  for (arma::uword ep = 0; ep < num_pts; ep++) {
    // Calculate regression coefficients
    if (end_pts(ep) > start_pts(ep)) {
      sub_response = response.subvec(start_pts(ep), end_pts(ep));
      sub_design = design.rows(start_pts(ep), end_pts(ep));
      reg_data = calc_reg(sub_response, sub_design, method, eigen_thresh, eigen_max, con_fi, alpha);
      reg_stats.row(ep) = conv_to<rowvec>::from(reg_data);
    }  // end if
  }  // end for
  
  // Warmup period
  // reg_stats.rows(0, num_cols+1) = zeros(num_cols+2, (num_cols + 1));
  // for (arma::uword it = (num_cols+2); it < look_back; it++) {
  //   sub_response = response.subvec(0, it);
  //   sub_design = design.rows(0, it);
  //   reg_data = calc_reg(sub_response, sub_design);
  //   reg_stats.row(it) = conv_to<rowvec>::from(reg_data);
  // }  // end for
  
  // Remaining periods
  // for (arma::uword it = look_back; it < num_rows; it++) {
  //   sub_response = response.subvec(it-look_back+1, it);
  //   sub_design = design.rows(it-look_back+1, it);
  //   reg_data = calc_reg(sub_response, sub_design, method, eigen_thresh, eigen_max, con_fi, alpha);
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
//' mat_rix <- matrix(rnorm(20000), nc=2)
//' look_back <- 11
//' rolled_scaled <- roll::roll_scale(data=mat_rix, width = look_back, min_obs=1)
//' rolled_scaled2 <- roll_scale(matrix=mat_rix, look_back = look_back, use_median=FALSE)
//' all.equal(rolled_scaled[-1, ], rolled_scaled2[-1, ])
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_scale(const arma::mat& matrix, 
                     arma::uword look_back,
                     bool use_median=false) {
  
  arma::uword num_rows = matrix.n_rows;
  arma::mat scaled_mat(num_rows, matrix.n_cols);
  arma::mat sub_mat;
  
  // Warmup period
  scaled_mat.row(0) = matrix.row(0);
  for (arma::uword it = 1; it < look_back; it++) {
    sub_mat = matrix.rows(0, it);
    sub_mat = calc_scaled(sub_mat, use_median);
    scaled_mat.row(it) = sub_mat.row(sub_mat.n_rows-1);
  }  // end for
  
  // Remaining periods
  for (arma::uword it = look_back; it < num_rows; it++) {
    sub_mat = matrix.rows(it-look_back+1, it);
    sub_mat = calc_scaled(sub_mat, use_median);
    scaled_mat.row(it) = sub_mat.row(sub_mat.n_rows-1);
  }  // end for
  
  return scaled_mat;
}  // end roll_scale



////////////////////////////////////////////////////////////
//' Calculate a \emph{vector} of z-scores of the residuals of rolling
//' regressions at the end points of the design matrix.
//' 
//' @param \code{response} A single-column \emph{time series} or a \emph{vector}
//'   of response data.
//' 
//' @param \code{design} A \emph{time series} or a \emph{matrix} of design data
//'   (predictor or explanatory data).
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
//'   \code{design}.
//'
//' @details 
//'   The function \code{roll_zscores()} calculates a \emph{vector} of z-scores
//'   of the residuals of rolling regressions at the end points of the
//'   \emph{time series} \code{design}.
//'   
//'   The function \code{roll_zscores()} performs a loop over the end points,
//'   and at each end point it subsets the time series \code{design} over a
//'   look-back interval equal to \code{look_back} number of end points.
//'   
//'   It passes the subset time series to the function \code{calc_lm()}, which
//'   calculates the regression data.
//'   
//'   If the arguments \code{endp} and \code{startp} are not given then it
//'   first calculates a vector of end points separated by \code{step} time
//'   periods. It calculates the end points along the rows of \code{design}
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
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("XLF", "VTI", "IEF")])
//' # Response equals XLF returns
//' res_ponse <- re_turns[, 1]
//' # Design matrix equals VTI and IEF returns
//' de_sign <- re_turns[, -1]
//' # Calculate Z-scores from rolling time series regression using RcppArmadillo
//' look_back <- 11
//' z_scores <- HighFreq::roll_zscores(response=res_ponse, design=de_sign, look_back)
//' # Calculate z-scores in R from rolling multivariate regression using lm()
//' z_scoresr <- sapply(1:NROW(de_sign), function(ro_w) {
//'   if (ro_w == 1) return(0)
//'   start_point <- max(1, ro_w-look_back+1)
//'   sub_response <- res_ponse[start_point:ro_w]
//'   sub_design <- de_sign[start_point:ro_w, ]
//'   reg_model <- lm(sub_response ~ sub_design)
//'   resid_uals <- reg_model$residuals
//'   resid_uals[NROW(resid_uals)]/sd(resid_uals)
//' })  # end sapply
//' # Compare the outputs of both functions
//' all.equal(z_scores[-(1:look_back)], z_scoresr[-(1:look_back)], 
//'   check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec roll_zscores(const arma::vec& response, 
                       const arma::mat& design, 
                       arma::uvec startp = 0, 
                       arma::uvec endp = 0, 
                       arma::uword step = 1, 
                       arma::uword look_back = 1,
                       arma::uword stub = 0) {
  
  // Allocate end points
  arma::uword num_rows = design.n_rows;
  arma::uvec end_pts;
  arma::uvec start_pts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    end_pts = calc_endpoints(num_rows, step, stub);
  } else {
    // Copy end points
    end_pts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    start_pts = calc_startpoints(end_pts, look_back);
  } else {
    // Copy start points
    start_pts = startp;
  }  // end if
  
  // Allocate regression matrix
  arma::vec sub_response;
  arma::mat sub_design;
  arma::uword num_pts = end_pts.n_elem;
  arma::vec z_scores(num_pts, fill::zeros);
  
  // Perform loop over the end points
  for (arma::uword ep = 0; ep < num_pts; ep++) {
    // Calculate z-scores
    if (end_pts(ep) > start_pts(ep)) {
      sub_response = response.subvec(start_pts(ep), end_pts(ep));
      sub_design = design.rows(start_pts(ep), end_pts(ep));
      z_scores(ep) = calc_lm(sub_response, sub_design)["z_score"];
    }  // end if
  }  // end for
  
  // Old code below
  // Warmup period
  // for (arma::uword it = 1; it < look_back; it++) {
  //   sub_response = response.subvec(0, it);
  //   sub_design = design.rows(0, it);
  //   z_scores(it) = calc_lm(sub_response, sub_design)["z_score"];
  // }  // end for
  
  // Remaining periods
  // for (arma::uword it = look_back; it < num_rows; it++) {
  //   sub_response = response.subvec(it-look_back+1, it);
  //   sub_design = design.rows(it-look_back+1, it);
  //   z_scores(it) = calc_lm(sub_response, sub_design)["z_score"];
  // }  // end for
  
  return z_scores;
  
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
//' @param \code{con_fi} The confidence level for calculating the
//'   quantiles (the default is \code{con_fi = 0.75}).
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
//'   periods. It calculates the end points along the rows of \code{design}
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
//' re_turns <- na.omit(rutils::etf_env$re_turns$VTI)
//' # Calculate the rolling variance at 25 day end points, with a 75 day look-back
//' var_rollfun <- HighFreq::roll_fun(re_turns, fun="calc_var", step=25, look_back=3)
//' # Calculate the rolling variance using roll_var()
//' var_roll <- HighFreq::roll_var(re_turns, step=25, look_back=3)
//' # Compare the two methods
//' all.equal(var_rollfun, var_roll, check.attributes=FALSE)
//' # Define end points and start points
//' end_p <- HighFreq::calc_endpoints(NROW(re_turns), step=25)
//' start_p <- HighFreq::calc_startpoints(end_p, look_back=3)
//' # Calculate the rolling variance using RcppArmadillo
//' var_rollfun <- HighFreq::roll_fun(re_turns, fun="calc_var", startp=start_p, endp=end_p)
//' # Calculate the rolling variance using R code
//' var_roll <- sapply(1:NROW(end_p), function(it) {
//'   var(re_turns[start_p[it]:end_p[it]+1, ])
//' })  # end sapply
//' # Compare the two methods
//' all.equal(drop(var_rollfun), var_roll, check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::roll_fun(re_turns, fun="calc_var", startp=start_p, endp=end_p),
//'   Rcode=sapply(1:NROW(end_p), function(it) {
//'     var(re_turns[start_p[it]:end_p[it]+1, ])
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
                   double con_fi = 0.75) {
  
  // Allocate end points
  arma::uword num_rows = tseries.n_rows;
  arma::uvec end_pts;
  arma::uvec start_pts;
  
  // Calculate end points if missing
  if (sum(endp) == 0) {
    end_pts = calc_endpoints(num_rows, step, stub);
  } else {
    // Copy end points
    end_pts = endp;
  }  // end if
  
  // Calculate start points if missing
  if (sum(startp) == 0) {
    // Start points equal to end points lagged by look_back
    start_pts = calc_startpoints(end_pts, look_back);
  } else {
    // Copy start points
    start_pts = startp;
  }  // end if
  
  // Allocate matrix of statistics
  arma::uword num_pts = end_pts.n_elem;
  arma::mat stats = arma::zeros<mat>(num_pts, tseries.n_cols);
  
  // Perform loop over the end points
  if (fun == "calc_mean") {
    // Calculate the dispersion (variance)
    for (arma::uword ep = 0; ep < num_pts; ep++) {
      // Calculate kurtosis
      if (end_pts(ep) > start_pts(ep)) {
        stats.row(ep) = calc_mean(tseries.rows(start_pts(ep), end_pts(ep)), method, con_fi);
      }  // end if
    }  // end for
  } else if (fun == "calc_var") {
    // Calculate the dispersion (variance)
    for (arma::uword ep = 0; ep < num_pts; ep++) {
      // Calculate kurtosis
      if (end_pts(ep) > start_pts(ep)) {
        stats.row(ep) = calc_var(tseries.rows(start_pts(ep), end_pts(ep)), method, con_fi);
      }  // end if
    }  // end for
  } else if (fun == "calc_skew") {
    // Perform loop over the end points
    for (arma::uword ep = 0; ep < num_pts; ep++) {
      // Calculate kurtosis
      if (end_pts(ep) > start_pts(ep)) {
        stats.row(ep) = calc_skew(tseries.rows(start_pts(ep), end_pts(ep)), method, con_fi);
      }  // end if
    }  // end for
  } else if (fun == "calc_kurtosis") {
    // Perform loop over the end points
    for (arma::uword ep = 0; ep < num_pts; ep++) {
      // Calculate kurtosis
      if (end_pts(ep) > start_pts(ep)) {
        stats.row(ep) = calc_kurtosis(tseries.rows(start_pts(ep), end_pts(ep)), method, con_fi);
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
//'   rolling variance from the historical returns, rather than simulating them.
//'   It represents exponential smoothing of the squared returns with a decay
//'   factor equal to \eqn{\beta}.
//'
//'   The function \code{sim_garch()} simulates the \emph{GARCH} process using
//'   fast \emph{Rcpp} \code{C++} code.
//'
//' @examples
//' \dontrun{
//' # Define the GARCH model parameters
//' al_pha <- 0.79
//' be_ta <- 0.2
//' om_ega <- 1e-4*(1-al_pha-be_ta)
//' # Calculate matrix of standard normal innovations
//' in_nov <- matrix(rnorm(1e3))
//' # Simulate the GARCH process using Rcpp
//' garch_data <- HighFreq::sim_garch(omega=om_ega, alpha=al_pha,  beta=be_ta, innov=in_nov)
//' # Plot the GARCH rolling volatility and cumulative returns
//' plot(sqrt(garch_data[, 2]), t="l", main="Simulated GARCH Volatility", ylab="volatility")
//' plot(cumsum(garch_data[, 1]), t="l", main="Simulated GARCH Cumulative Returns", ylab="cumulative returns")
//' # Calculate historical VTI returns
//' re_turns <- na.omit(rutils::etf_env$re_turns$VTI)
//' # Estimate the volatility of VTI returns
//' garch_data <- HighFreq::sim_garch(omega=om_ega, alpha=al_pha,  beta=be_ta, 
//'   innov=re_turns, is_random=FALSE)
//' # Plot dygraph of the estimated GARCH volatility
//' dygraphs::dygraph(xts::xts(sqrt(garch_data[, 2]), index(re_turns)), 
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
  
  arma::uword num_rows = innov.n_rows;
  arma::vec variance(num_rows);
  
  if (is_random) {
    // The innovations are random numbers
    arma::vec returns(num_rows);
    variance[0] = omega/(1-alpha-beta);
    returns[0] = std::sqrt(variance[0])*innov[0];
    
    for (arma::uword it = 1; it < num_rows; it++) {
      returns[it] = std::sqrt(variance[it-1])*innov[it];
      variance[it] = omega + alpha*pow(returns[it], 2) + beta*variance[it-1];
    }  // end for
    return join_rows(returns, variance);
  } else {
    // The innovations are historical returns
    for (arma::uword it = 1; it < num_rows; it++) {
      variance[it] = omega + alpha*pow(innov[it], 2) + beta*variance[it-1];
    }  // end for
    return join_rows(innov, variance);
  }  // end if
  
}  // end sim_garch



////////////////////////////////////////////////////////////
//' Simulate an \emph{Ornstein-Uhlenbeck} process using \emph{Rcpp}.
//' 
//' @param \code{volat} The volatility of returns.
//' 
//' @param \code{eq_price} The equilibrium price. 
//' 
//' @param \code{theta} The strength of mean reversion.
//' 
//' @param \code{innov} A single-column \emph{matrix} of innovations (random
//'   numbers).
//' 
//' @return A single-column \emph{matrix} of simulated returns, with the same
//'   number of rows as the argument \code{innov}.
//'
//' @details 
//'   The function \code{sim_ou()} simulates the following
//'   \emph{Ornstein-Uhlenbeck} process:
//'   \deqn{
//'     r_i = p_i - p_{i-1} = \theta \, (\mu - p_{i-1}) + \sigma \, \xi_i
//'   }
//'   \deqn{
//'     p_i = p_{i-1} + r_i
//'   }
//'   Where \eqn{r_i} and \eqn{p_i} are the simulated returns and prices,
//'   \eqn{\theta}, \eqn{\mu}, and \eqn{\sigma} are the
//'   \emph{Ornstein-Uhlenbeck} parameters, and \eqn{\xi_i} are the standard
//'   normal \emph{innovations}.
//'   The recursion starts with: \eqn{p_1 = r_1 = \sigma \, \xi_1}.
//'
//'   The function \code{sim_ou()} simulates the percentage returns as equal to
//'   the difference between the equilibrium price \eqn{\mu} minus the latest
//'   price \eqn{p_{i-1}}, times the mean reversion parameter \eqn{\theta}, plus
//'   a random innovation proportional to the volatility \eqn{\sigma}. The log
//'   prices are calculated as the sum of returns (not compounded), so they can
//'   become negative.
//'
//'   The function \code{sim_ou()} simulates the \emph{Ornstein-Uhlenbeck}
//'   process using fast \emph{Rcpp} \code{C++} code.
//'
//'   The function \code{sim_ou()} returns a single-column \emph{matrix}
//'   representing the \emph{time series} of simulated returns.
//'
//' @examples
//' \dontrun{
//' # Define the Ornstein-Uhlenbeck model parameters
//' eq_price <- 1.0
//' sig_ma <- 0.01
//' the_ta <- 0.01
//' in_nov <- matrix(rnorm(1e3))
//' # Simulate Ornstein-Uhlenbeck process using Rcpp
//' re_turns <- HighFreq::sim_ou(eq_price=eq_price, volat=sig_ma, theta=the_ta, innov=in_nov)
//' plot(cumsum(re_turns), t="l", main="Simulated Ornstein-Uhlenbeck Prices", ylab="prices")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat sim_ou(double eq_price, 
                 double volat, 
                 double theta, 
                 arma::mat& innov) {
  
  arma::uword num_rows = innov.n_rows;
  arma::mat prices = arma::zeros<mat>(num_rows, 1);
  arma::mat returns = arma::zeros<mat>(num_rows, 1);
  
  prices.row(0) = volat*innov.row(0);
  for (arma::uword it = 1; it < num_rows; it++) {
    returns.row(it) = theta*(eq_price - prices.row(it-1)) + volat*innov.row(it);
    prices.row(it) = prices.row(it-1) + returns.row(it);
  }  // end for
  
  return returns;
  
}  // end sim_ou



////////////////////////////////////////////////////////////
//' Simulate a \emph{Schwartz} process using \emph{Rcpp}.
//' 
//' @param \code{volat} The volatility of returns.
//' 
//' @param \code{eq_price} The equilibrium price. 
//' 
//' @param \code{theta} The strength of mean reversion.
//' 
//' @param \code{innov} A single-column \emph{matrix} of innovations (random
//'   numbers).
//' 
//' @return A single-column \emph{matrix} of simulated returns, with the same
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
//'   \eqn{\theta}, plus a random innovation proportional to the volatility
//'   \eqn{\sigma}.
//'
//'   The function \code{sim_schwartz()} returns a single-column \emph{matrix}
//'   representing the \emph{time series} of simulated returns.
//'
//' @examples
//' \dontrun{
//' # Define the Schwartz model parameters
//' eq_price <- 2.0
//' sig_ma <- 0.01
//' the_ta <- 0.01
//' in_nov <- matrix(rnorm(1e3))
//' # Simulate Schwartz process using Rcpp
//' re_turns <- HighFreq::sim_schwartz(eq_price=eq_price, volat=sig_ma, theta=the_ta, innov=in_nov)
//' plot(exp(cumsum(re_turns)), t="l", main="Simulated Schwartz Prices", ylab="prices")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat sim_schwartz(double eq_price, 
                       double volat, 
                       double theta, 
                       arma::mat& innov) {
  
  arma::uword num_rows = innov.n_rows;
  arma::mat prices = arma::zeros<mat>(num_rows, 1);
  arma::mat returns = arma::zeros<mat>(num_rows, 1);
  
  prices.row(0) = exp(volat*innov.row(0));
  for (arma::uword it = 1; it < num_rows; it++) {
    returns.row(it) = theta*(eq_price - prices.row(it-1)) + volat*innov.row(it);
    prices.row(it) = prices.row(it-1) * exp(returns.row(it));
  }  // end for
  
  return returns;
  
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
//'   \eqn{AR(p)} of order \eqn{p}:
//'   \deqn{
//'     r_i = \varphi_1 r_{i-1} + \varphi_2 r_{i-2} + \ldots + \varphi_p r_{i-p} + \xi_i
//'   }
//'   Where \eqn{r_i} is the simulated output time series, \eqn{\varphi_i} are
//'   the \emph{autoregressive} coefficients, and \eqn{\xi_i} are the standard
//'   normal \emph{innovations}.
//'
//'   The order \eqn{p} of the \emph{autoregressive} process \eqn{AR(p)}, is
//'   equal to the number of rows of the \emph{autoregressive} coefficients
//'   \code{coeff}.
//'
//'   The function \code{sim_ar()} performs the same calculation as the standard
//'   \code{R} function \cr\code{filter(x=innov, filter=co_eff,
//'   method="recursive")}, but it's several times faster.
//'   
//' @examples
//' \dontrun{
//' # Define AR coefficients
//' co_eff <- matrix(c(0.2, 0.2))
//' # Calculate matrix of innovations
//' in_nov <- matrix(rnorm(1e4, sd=0.01))
//' # Calculate recursive filter using filter()
//' filter_ed <- filter(in_nov, filter=co_eff, method="recursive")
//' # Calculate recursive filter using RcppArmadillo
//' re_turns <- HighFreq::sim_ar(co_eff, in_nov)
//' # Compare the two methods
//' all.equal(as.numeric(re_turns), as.numeric(filter_ed))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::sim_ar(co_eff, in_nov),
//'   Rcode=filter(in_nov, filter=co_eff, method="recursive"),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat sim_ar(arma::mat& coeff, const arma::mat& innov) {
  
  arma::uword num_rows = innov.n_rows;
  arma::uword look_back = coeff.n_rows;
  arma::mat rev_coeff = arma::reverse(coeff);
  arma::mat returns = arma::zeros<mat>(num_rows, 1);

  // Warmup period
  returns.row(0) = innov.row(0);
  returns.row(1) = innov.row(1) + rev_coeff.row(look_back-1) * returns.row(0);
  for (arma::uword it=2; it < look_back-1; it++) {
    returns.row(it) = innov.row(it) + arma::dot(rev_coeff.rows(look_back-it, look_back-1), returns.rows(0, it-1));
  }  // end for
  
  // Remaining periods
  for (arma::uword it = look_back; it < num_rows; it++) {
    returns.row(it) = innov.row(it) + arma::dot(rev_coeff, returns.rows(it-look_back, it-1));
  }  // end for
  
  return returns;
  
}  // end sim_ar



////////////////////////////////////////////////////////////
//' Simulate a \emph{Dickey-Fuller} process using \emph{Rcpp}.
//' 
//' @param \code{volat} The volatility of returns.
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
//' @return A single-column \emph{matrix} of simulated returns, with the same
//'   number of rows as the argument \code{innov}.
//'
//' @details 
//'   The function \code{sim_df()} simulates the following \emph{Dickey-Fuller}
//'   process:
//'   \deqn{
//'     r_i = \theta \, (\mu - p_{i-1}) + \varphi_1 r_{i-1} + \ldots + \varphi_p r_{i-p} + \sigma \, \xi_i
//'   }
//'   \deqn{
//'     p_i = p_{i-1} + r_i
//'   }
//'   Where \eqn{r_i} and \eqn{p_i} are the simulated returns and prices,
//'   \eqn{\theta}, \eqn{\mu}, and \eqn{\sigma} are the
//'   \emph{Ornstein-Uhlenbeck} parameters, \eqn{\varphi_i} are the
//'   \emph{autoregressive} coefficients, and \eqn{\xi_i} are the standard
//'   normal \emph{innovations}.
//'   The recursion starts with: \eqn{p_1 = r_1 = \sigma \, \xi_1}.
//'
//'   The \emph{Dickey-Fuller} process is a combination of an
//'   \emph{Ornstein-Uhlenbeck} process and an \emph{autoregressive} process.
//'   The order \eqn{p} of the \emph{autoregressive} process \eqn{AR(p)}, is
//'   equal to the number of rows of the \emph{autoregressive} coefficients
//'   \code{coeff}.
//'
//'   The function \code{sim_df()} simulates the \emph{Dickey-Fuller}
//'   process using fast \emph{Rcpp} \code{C++} code.
//'
//'   The function \code{sim_df()} returns a single-column \emph{matrix}
//'   representing the \emph{time series} of returns.
//'
//' @examples
//' \dontrun{
//' # Define the Ornstein-Uhlenbeck model parameters
//' eq_price <- 1.0
//' sig_ma <- 0.01
//' the_ta <- 0.01
//' # Define AR coefficients
//' co_eff <- matrix(c(0.2, 0.2))
//' # Calculate matrix of standard normal innovations
//' in_nov <- matrix(rnorm(1e3))
//' # Simulate Dickey-Fuller process using Rcpp
//' re_turns <- HighFreq::sim_df(eq_price=eq_price, volat=sig_ma, theta=the_ta, co_eff, innov=in_nov)
//' plot(cumsum(re_turns), t="l", main="Simulated Dickey-Fuller Prices")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat sim_df(double eq_price, 
                 double volat, 
                 double theta, 
                 arma::mat& coeff, 
                 arma::mat& innov) {
  
  arma::uword num_rows = innov.n_rows;
  arma::uword look_back = coeff.n_rows;
  arma::mat rev_coeff = arma::reverse(coeff);
  arma::mat prices = arma::zeros<mat>(num_rows, 1);
  arma::mat returns = arma::zeros<mat>(num_rows, 1);

  // Warmup period
  returns.row(0) = volat*innov.row(0);
  returns.row(1) = volat*innov.row(1) + rev_coeff.row(look_back-1) * returns.row(0);
  prices.row(0) = volat*innov.row(0);
  for (arma::uword it=2; it < look_back-1; it++) {
    returns.row(it) = volat*innov.row(it) + theta*(eq_price - prices.row(it-1)) + arma::dot(rev_coeff.rows(look_back-it, look_back-1), returns.rows(0, it-1));
    prices.row(it) = prices.row(it-1) + returns.row(it);
  }  // end for
  
  // Remaining periods
  for (arma::uword it = look_back; it < num_rows; it++) {
    returns.row(it) = volat*innov.row(it) + theta*(eq_price - prices.row(it-1)) + arma::dot(rev_coeff, returns.rows(it-look_back, it-1));
    prices.row(it) = prices.row(it-1) + returns.row(it);
  }  // end for
  
  return returns;
  
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
//'   normal distribution of returns as follows:
//'   \deqn{
//'     likelihood = - \sum_{i=1}^n (\frac{r^2_i}{\sigma^2_i} + \log(\sigma^2_i))
//'   }
//'
//' @examples
//' \dontrun{
//' # Define the GARCH model parameters
//' al_pha <- 0.79
//' be_ta <- 0.2
//' om_ega <- 1e-4*(1-al_pha-be_ta)
//' # Calculate historical VTI returns
//' re_turns <- na.omit(rutils::etf_env$re_turns$VTI)
//' # Calculate the log-likelihood of VTI returns assuming GARCH(1,1)
//' HighFreq::lik_garch(omega=om_ega, alpha=al_pha,  beta=be_ta, returns=re_turns)
//' }
//' 
//' @export
// [[Rcpp::export]]
double lik_garch(double omega, 
                 double alpha, 
                 double beta,
                 arma::mat& returns, 
                 double minval = 0.000001) {
  
  // Calculate the rolling volatility of returns using function sim_garch()
  arma::mat garch_data = sim_garch(omega, alpha,  beta, returns, FALSE);
  // Select the second column containing the volatility of returns
  arma::mat variance = garch_data.col(1);
  // Apply floor to volatility
  variance.transform([&minval](double x) {return max(x, minval);});
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
//' @param \code{method} A \emph{string} specifying the objective function for
//'   calculating the weights (see Details) (the default is \code{method =
//'   "rank_sharpe"})
//'   
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small singular values in order to regularize the inverse of the
//'   \code{returns} matrix (the default is \code{0.001}).
//'   
//' @param \code{eigen_max} An \emph{integer} equal to the number of singular
//'   values used for calculating the regularized inverse of the \code{returns}
//'   matrix (the default is \code{0} - equivalent to \code{eigen_max} equal to
//'   the number of columns of \code{returns}).
//'   
//' @param \code{con_fi} The confidence level for calculating the
//'   quantiles (the default is \code{con_fi = 0.75}).
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
//'   weights for different types of objective functions, using
//'   \code{RcppArmadillo} \code{C++} code.
//' 
//'   If \code{method = "rank_sharpe"} (the default) then it calculates the
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
//'   \code{calc_weights()} calculates the regularized inverse of the covariance
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
//' re_turns <- na.omit(rutils::etf_env$re_turns[, 1:16])
//' ei_gen <- eigen(cov(re_turns))
//' # Calculate regularized inverse of covariance matrix
//' eigen_max <- 3
//' eigen_vec <- ei_gen$vectors[, 1:eigen_max]
//' eigen_val <- ei_gen$values[1:eigen_max]
//' in_verse <- eigen_vec %*% (t(eigen_vec) / eigen_val)
//' # Define shrinkage intensity and apply shrinkage to the mean returns
//' al_pha <- 0.5
//' col_means <- colMeans(re_turns)
//' col_means <- ((1-al_pha)*col_means + al_pha*mean(col_means))
//' # Calculate weights using R
//' weight_s <- in_verse %*% col_means
//' n_col <- NCOL(re_turns)
//' weights_r <- weights_r*sd(re_turns %*% rep(1/n_col, n_col))/sd(re_turns %*% weights_r)
//' # Calculate weights using RcppArmadillo
//' weight_s <- drop(HighFreq::calc_weights(re_turns, eigen_max, alpha=al_pha))
//' all.equal(weight_s, weights_r)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec calc_weights(const arma::mat& returns, // Portfolio returns
                       std::string method = "rank_sharpe",
                       double eigen_thresh = 0.001,
                       arma::uword eigen_max = 0,
                       double con_fi = 0.1,
                       double alpha = 0.0,
                       bool scale = true,
                       double vol_target = 0.01) {
  // Initialize
  arma::vec weights(returns.n_cols, fill::zeros);
  if (eigen_max == 0)  eigen_max = returns.n_cols;
  
  // Switch for the different methods for weights
  switch(calc_method(method)) {
  case meth_od::rank_sharpe: {
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(returns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to<vec>::from(arma::sort_index(arma::sort_index(mean_cols)));
    weights = (weights - arma::mean(weights));
    break;
  }  // end rank_sharpe
  case meth_od::max_sharpe: {
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(returns, 0));
    // Shrink mean_cols to the mean of returns
    mean_cols = ((1-alpha)*mean_cols + alpha*arma::mean(mean_cols));
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(cov(returns), eigen_max);
    // weights = calc_inv(cov(returns), eigen_max)*mean_cols;
    weights = calc_inv(cov(returns), eigen_thresh, eigen_max)*mean_cols;
    break;
  }  // end max_sharpe
  case meth_od::max_sharpe_median: {
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::median(returns, 0));
    // Shrink mean_cols to the mean of returns
    mean_cols = ((1-alpha)*mean_cols + alpha*arma::median(mean_cols));
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(cov(returns), eigen_max);
    weights = calc_inv(cov(returns), eigen_thresh, eigen_max)*mean_cols;
    break;
  }  // end max_sharpe_median
  case meth_od::min_var: {
    // Apply regularized inverse to unit vector
    weights = calc_inv(cov(returns), eigen_thresh, eigen_max)*arma::ones(returns.n_cols);
    break;
  }  // end min_var
  case meth_od::min_varpca: {
    // Calculate highest order principal component
    arma::vec eigen_val;
    arma::mat eigen_vec;
    arma::eig_sym(eigen_val, eigen_vec, arma::cov(returns));
    weights = eigen_vec.col(0);
    break;
  }  // end min_varpca
  case meth_od::rank: {
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(returns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to<vec>::from(arma::sort_index(arma::sort_index(mean_cols)));
    weights = (weights - arma::mean(weights));
    break;
  }  // end rank
  case meth_od::rankrob: {
    // Median returns by columns
    arma::vec mean_cols = arma::trans(arma::median(returns, 0));
    // mean_cols = ((1-alpha)*mean_cols + alpha*arma::mean(mean_cols));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(returns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(cov(returns), eigen_max);
    // weights = calc_inv(cov(returns), eigen_max)*mean_cols;
    // weights = calc_inv(cov(returns), eigen_max)*mean_cols;
    // // Standard deviation by columns
    // arma::vec sd_cols = mean_cols;
    // for (arma::uword it=0; it < returns.n_cols; it++) {
    //   sd_cols(it) = arma::median(arma::abs((returns.col(it) - sd_cols)));
    // }  // end for
    // sd_cols.replace(0, 1);
    // mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weights = conv_to<vec>::from(arma::sort_index(arma::sort_index(mean_cols)));
    // level;
    weights = (weights - arma::mean(weights));
    break;
  }  // end rankrob
  case meth_od::quantile: {
    // Sum of quantiles for columns
    arma::vec level_s = {con_fi, 1-con_fi};
    weights = conv_to<vec>::from(arma::sum(arma::quantile(returns, level_s, 0), 0));
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
    // arma::vec mean_rows = arma::mean(returns, 1);
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
//' @param \code{coeff} A \emph{numeric} multiplier of the weights.  (The
//'   default is \code{1})
//'   
//' @param \code{bid_offer} A \emph{numeric} bid-offer spread (the default is
//'   \code{0})
//'
//' @param \code{method} A \emph{string} specifying the objective function for
//'   calculating the weights (see Details) (the default is \code{method =
//'   "rank_sharpe"})
//'   
//' @param \code{eigen_thresh} A \emph{numeric} threshold level for discarding
//'   small singular values in order to regularize the inverse of the
//'   \code{returns} matrix (the default is \code{0.001}).
//'   
//' @param \code{eigen_max} An \emph{integer} equal to the number of singular
//'   values used for calculating the regularized inverse of the \code{returns}
//'   matrix (the default is \code{0} - equivalent to \code{eigen_max} equal to
//'   the number of columns of \code{returns}).
//'   
//' @param \code{con_fi} The confidence level for calculating the
//'   quantiles (the default is \code{con_fi = 0.75}).
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
//' @return A column \emph{vector} of strategy returns, with the same length as
//'   the number of rows of \code{returns}.
//'
//' @details 
//'   The function \code{back_test()} performs a backtest simulation of a
//'   rolling portfolio optimization strategy over a \emph{vector} of
//'   \code{endp}.
//'   
//'   It performs a loop over the end points \code{endp}, and subsets the
//'   \emph{matrix} of excess returns \code{excess} along its rows, between the
//'   corresponding end point and the start point. It passes the subset matrix
//'   of excess returns into the function \code{calc_weights()}, which
//'   calculates the optimal portfolio weights. The arguments \code{eigen_max},
//'   \code{alpha}, \code{method}, and \code{scale} are also passed to the
//'   function \code{calc_weights()}.
//'   
//'   The function \code{back_test()} multiplies the weights by the coefficient
//'   \code{coeff} (with default equal to \code{1}), which allows reverting a
//'   strategy if \code{co_eff = -1}.
//'   
//'   The function \code{back_test()} then multiplies the weights times the
//'   future portfolio returns, to calculate the out-of-sample strategy returns.
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
//' re_turns <- na.omit(rutils::etf_env$re_turns[, 1:16])
//' # risk_free is the daily risk-free rate
//' risk_free <- 0.03/260
//' ex_cess <- re_turns - risk_free
//' # Define monthly end points without initial warmpup period
//' end_p <- rutils::calc_endpoints(re_turns, inter_val="months")
//' end_p <- end_p[end_p > 0]
//' len_gth <- NROW(end_p)
//' # Define 12-month look-back interval and start points over sliding window
//' look_back <- 12
//' start_p <- c(rep_len(1, look_back-1), end_p[1:(len_gth-look_back+1)])
//' # Define shrinkage and regularization intensities
//' al_pha <- 0.5
//' eigen_max <- 3
//' # Simulate a monthly rolling portfolio optimization strategy
//' pnl_s <- HighFreq::back_test(ex_cess, re_turns, 
//'                             start_p-1, end_p-1, 
//'                             eigen_max = eigen_max, 
//'                             alpha = al_pha)
//' pnl_s <- xts::xts(pnl_s, index(re_turns))
//' colnames(pnl_s) <- "strat_rets"
//' # Plot dygraph of strategy
//' dygraphs::dygraph(cumsum(pnl_s), 
//'   main="Cumulative Returns of Max Sharpe Portfolio Strategy")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat back_test(const arma::mat& excess, // Portfolio excess returns
                    const arma::mat& returns, // Portfolio returns
                    arma::uvec startp, 
                    arma::uvec endp, 
                    std::string method = "rank_sharpe",
                    double eigen_thresh = 0.001,
                    arma::uword eigen_max = 0,
                    double con_fi = 0.1,
                    double alpha = 0.0,
                    bool scale = true,
                    double vol_target = 0.01,
                    double coeff = 1.0,
                    double bid_offer = 0.0) {
  
  arma::vec weights(returns.n_cols, fill::zeros);
  arma::vec weights_past = zeros(returns.n_cols);
  arma::mat pnl_s = zeros(returns.n_rows, 1);
  
  // Perform loop over the end points
  for (arma::uword it = 1; it < endp.size(); it++) {
    // cout << "it: " << it << endl;
    // Calculate portfolio weights
    weights = coeff*calc_weights(excess.rows(startp(it-1), endp(it-1)), method, eigen_thresh, eigen_max, con_fi, alpha, scale, vol_target);
    // Calculate out-of-sample returns
    pnl_s.rows(endp(it-1)+1, endp(it)) = returns.rows(endp(it-1)+1, endp(it))*weights;
    // Add transaction costs
    pnl_s.row(endp(it-1)+1) -= bid_offer*sum(abs(weights - weights_past))/2;
    weights_past = weights;
  }  // end for
  
  // Return the strategy pnl_s
  return pnl_s;
  
}  // end back_test

