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
//' @param \code{t_series} A single-column \emph{time series} or a
//'   \emph{vector}.
//'
//' @param \code{lagg} An \emph{integer} equal to the number of periods to lag
//'   (the default is \code{lagg = 1}).
//'
//' @param \code{pad_zeros} \emph{Boolean} argument: Should the output be padded
//'   with zeros? (The default is \code{pad_zeros = TRUE}.)
//'
//' @return A column \emph{vector} with the same number of elements as the input
//'   time series.
//'
//' @details The function \code{lag_vec()} applies a lag to the input \emph{time
//'   series} \code{t_series} by shifting its elements by the number equal to
//'   the argument \code{lagg}.  For positive \code{lagg} values, the elements
//'   are shifted forward in time (down), and for negative \code{lagg} values
//'   they are shifted backward (up).
//'   
//'   The output \emph{vector} is padded with either zeros (the default), or
//'   with data from \code{t_series}, so that it has the same number of element
//'   as \code{t_series}.
//'   If the \code{lagg} is positive, then the first element is copied and added
//'   upfront.
//'   If the \code{lagg} is negative, then the last element is copied and added
//'   to the end.
//'   
//'   As a rule, if \code{t_series} contains returns data, then the output
//'   \emph{matrix} should be padded with zeros, to avoid data snooping.
//'   If \code{t_series} contains prices, then the output \emph{matrix} should
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
arma::vec lag_vec(arma::vec& t_series, 
                  arma::sword lagg = 1, 
                  bool pad_zeros = true) {
  
  arma::uword num_rows = (t_series.n_elem-1);
  
  if (lagg > 0) {
    if (pad_zeros) {
      // Pad front with zeros
      return arma::join_cols(arma::zeros(lagg), 
                             t_series.subvec(0, num_rows-lagg));
    } else {
      // Pad front with first element of t_series
      return arma::join_cols(arma::repelem(t_series.subvec(0, 0), lagg, 1), 
                             t_series.subvec(0, num_rows-lagg));
    }  // end if
  } else {
    if (pad_zeros) {
      // Pad back with zeros
      return arma::join_cols(t_series.subvec(-lagg, num_rows), 
                             arma::zeros(-lagg));
    } else {
      // Pad back with last element of t_series
      return arma::join_cols(t_series.subvec(-lagg, num_rows), 
                             arma::repelem(t_series.subvec(num_rows, num_rows), -lagg, 1));
    }  // end if
  }  // end if
  
}  // end lag_vec




////////////////////////////////////////////////////////////
//' Apply a lag to the rows of a \emph{time series} or a \emph{matrix} using
//' \code{RcppArmadillo}.
//' 
//' @param \code{t_series} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lagg} An \emph{integer} equal to the number of periods to lag
//'   (the default is \code{lagg = 1}).
//'
//' @param \code{pad_zeros} \emph{Boolean} argument: Should the output be padded
//'   with zeros? (The default is \code{pad_zeros = TRUE}.)
//'   
//' @return A \emph{matrix} with the same dimensions as the input
//'   argument \code{t_series}.
//'
//' @details The function \code{lag_it()} applies a lag to the input
//'   \emph{matrix} by shifting its rows by the number equal to the argument
//'   \code{lagg}. For positive \code{lagg} values, the rows are shifted forward
//'   (down), and for negative \code{lagg} values they are shifted backward
//'   (up). 
//'   
//'   The output \emph{matrix} is padded with either zeros (the default), or
//'   with rows of data from \code{t_series}, so that it has the same dimensions
//'   as \code{t_series}.
//'   If the \code{lagg} is positive, then the first row is copied and added
//'   upfront.
//'   If the \code{lagg} is negative, then the last row is copied and added
//'   to the end.
//'   
//'   As a rule, if \code{t_series} contains returns data, then the output
//'   \emph{matrix} should be padded with zeros, to avoid data snooping.
//'   If \code{t_series} contains prices, then the output \emph{matrix} should
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
arma::mat lag_it(arma::mat& t_series, 
                 arma::sword lagg = 1, 
                 bool pad_zeros = true) {
  
  arma::uword num_rows = (t_series.n_rows-1);
  arma::uword num_cols = t_series.n_cols;
  
  if (lagg > 0) {
    // Positive lag
    if (pad_zeros) {
      // Pad front with zeros
      return arma::join_cols(arma::zeros(lagg, num_cols), 
                             t_series.rows(0, num_rows-lagg));
    } else {
      // Pad front with first element of t_series
      return arma::join_cols(arma::repmat(t_series.rows(0, 0), lagg, 1), 
                             t_series.rows(0, num_rows-lagg));
    }  // end if
  } else {
    // Negative lag
    if (pad_zeros) {
      // Pad back with zeros
      return arma::join_cols(t_series.rows(-lagg, num_rows), 
                             arma::zeros(-lagg, num_cols));
    } else {
      // Pad back with last element of t_series
      return arma::join_cols(t_series.rows(-lagg, num_rows), 
                             arma::repmat(t_series.rows(num_rows, num_rows), -lagg, 1));
    }  // end if
  }  // end if
  
  // Old code below
  // if (lagg > 0)
  //   // Positive lag
  //   return arma::join_cols(arma::repelem(t_series.row(0), lagg, 1), 
  //                          t_series.rows(0, num_rows-lagg));
  // else
  //   // Negative lag
  //   return arma::join_cols(t_series.rows(-lagg, num_rows), 
  //                          arma::repelem(t_series.row(num_rows), -lagg, 1));
  
}  // end lag_it




////////////////////////////////////////////////////////////
//' Calculate the differences between the neighboring elements of a
//' single-column \emph{time series} or a \emph{vector}.
//' 
//' @param \code{t_series} A single-column \emph{time series} or a \emph{vector}.
//' 
//' @param \code{lagg} An \emph{integer} equal to the number of time periods to
//'   lag when calculating the differences (the default is \code{lagg = 1}).
//'   
//' @param \code{padd} \emph{Boolean} argument: Should the output \emph{vector}
//'   be padded (extended) with zeros, in order to return a \emph{vector} of the
//'   same length as the input? (the default is \code{padd = TRUE})
//'
//' @return A column \emph{vector} containing the differences between the
//'   elements of the input vector.
//'
//' @details The function \code{diff_vec()} calculates the differences between
//'   the input \emph{time series} or \emph{vector} and its lagged version. 
//'   
//'   The argument \code{lagg} specifies the number of lags.  For example, if
//'   \code{lagg=3} then the differences will be taken between each element
//'   minus the element three time periods before it (in the past).  The default
//'   is \code{lagg = 1}.
//' 
//'   The argument \code{padd} specifies whether the output \emph{vector} should
//'   be padded (extended) with zeros at the beginning, in order to return a
//'   \emph{vector} of the same length as the input.  The default is
//'   \code{padd = TRUE}. The padding operation can be time-consuming, because it
//'   requires the copying of data.
//'   
//'   The function \code{diff_vec()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, which makes it several times faster than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random returns
//' re_turns <- rnorm(1e6)
//' # Compare diff_vec() with rutils::diff_it()
//' all.equal(drop(HighFreq::diff_vec(re_turns, lagg=3, padd=TRUE)),
//'   rutils::diff_it(re_turns, lagg=3))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::diff_vec(re_turns, lagg=3, padd=TRUE),
//'   Rcode=rutils::diff_it(re_turns, lagg=3),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec diff_vec(arma::vec& t_series, int lagg = 1, bool padd = true) {
  
  int len_gth = (t_series.n_elem-1);
  
  if (padd)
    // Pad the output with zeros at the beginning
    return (t_series - arma::join_cols(t_series.subvec(0, lagg - 1), 
                                      t_series.subvec(0, len_gth - lagg)));
  else
    // Don't pad the output
    return (t_series.subvec(lagg, len_gth) - t_series.subvec(0, len_gth - lagg));
  
}  // end diff_vec




////////////////////////////////////////////////////////////
//' Calculate the row differences of a a \emph{time series} or a \emph{matrix}
//' using \emph{RcppArmadillo}.
//' 
//' @param \code{t_series} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{lagg} An \emph{integer} equal to the number of rows (time
//'   periods) to lag when calculating the differences (the default is
//'   \code{lagg = 1}).
//'   
//' @param \code{padd} \emph{Boolean} argument: Should the output \emph{matrix}
//'   be padded (extended) with zeros, in order to return a \emph{matrix} with
//'   the same number of rows as the input? (the default is \code{padd = TRUE})
//'
//' @return A \emph{matrix} containing the differences between the rows of the
//'   input \emph{matrix}.
//'
//' @details The function \code{diff_it()} calculates the differences between
//'   the rows of the input \emph{time series} or \emph{matrix} and its lagged
//'   version. The lagged version has its rows shifted down by the number equal
//'   to \code{lagg} rows.
//'   
//'   The argument \code{lagg} specifies the number of lags applied to the rows
//'   of the lagged version. For example, if \code{lagg=3} then the lagged
//'   version will have its rows shifted down by \code{3} rows, and the
//'   differences will be taken between each row minus the row three time
//'   periods before it (in the past). The default is \code{lagg = 1}.
//' 
//'   The argument \code{padd} specifies whether the output \emph{matrix} should
//'   be padded (extended) with the rows of the initial (warmup) period at the
//'   beginning, in order to return a \emph{matrix} with the same number of rows
//'   as the input.  The default is \code{padd = TRUE}. The padding operation
//'   can be time-consuming, because it requires the copying of data.
//'   
//'   The function \code{diff_it()} is implemented in \code{RcppArmadillo}
//'   \code{C++} code, which makes it several times faster than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Create a matrix of random returns
//' re_turns <- matrix(rnorm(5e6), nc=5)
//' # Compare diff_it() with rutils::diff_it()
//' all.equal(HighFreq::diff_it(re_turns, lagg=3, padd=TRUE), 
//'   zoo::coredata(rutils::diff_it(re_turns, lagg=3)), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::diff_it(re_turns, lagg=3, padd=TRUE),
//'   Rcode=rutils::diff_it(re_turns, lagg=3),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat diff_it(arma::mat& t_series, 
                  arma::uword lagg = 1, 
                  bool padd = true) {
  
  arma::uword num_rows = (t_series.n_rows-1);
  // Matrix difference without padding
  arma::mat diff_mat = (t_series.rows(lagg, num_rows) - t_series.rows(0, num_rows - lagg));
  
  if (padd)
    // Pad diff_mat with warmup period at the beginning
    return arma::join_cols(t_series.rows(0, lagg - 1), diff_mat);
  else
    // Don't pad the output
    return diff_mat;
  
}  // end diff_it




////////////////////////////////////////////////////////////
//' Calculate a vector of end points that divides a vector into equal intervals.
//'
//' @param \code{len_gth} An \emph{integer} equal to the length of the vector to
//'   be divide into equal intervals.
//'   
//' @param \code{ste_p} The number of elements in each interval between
//'   neighboring end points.
//' 
//' @param \code{front} \emph{Boolean} argument: if \code{TRUE} then add a stub
//'   interval at the beginning, else add a stub interval at the end.  (default
//'   is \code{TRUE})
//'
//' @return A vector of equally spaced index values representing the end points
//'   (a vector of unsigned integers).
//'
//' @details The end points are a vector of unsigned integers which divide a
//'   vector of length equal to \code{len_gth} into equally spaced intervals. If
//'   a whole number of intervals doesn't fit over the vector, then
//'   \code{calc_endpoints()} adds a stub interval either at the beginning (the
//'   default) or at the end.
//'   The end points are shifted by \code{-1} because indexing starts at
//'   \code{0} in \code{C++} code.
//'
//'   The function \code{calc_endpoints()} is similar to the function
//'   \code{rutils::calc_endpoints()} from package
//'   \href{https://github.com/algoquant/rutils}{rutils}.
//'   
//'   The end points produced by \code{calc_endpoints()} don't include the first
//'   placeholder end point, which is usually equal to zero.
//'   For example, consider the end points for a vector of length \code{20}
//'   divided into intervals of length \code{5}: \code{0, 5, 10, 15, 20}.
//'   In order for all the differences between neighboring end points to be
//'   equal to \code{5}, the first end point must be equal to \code{0}.
//'   The first end point is a placeholder and doesn't correspond to any vector
//'   element.
//'   
//'   This works in \code{R} code because the vector element corresponding to
//'   index \code{0} is empty.  For example, the \code{R} code: \code{(4:1)[c(0,
//'   1)]} produces \code{4}.  So in \code{R} we can select vector elements
//'   using the end points starting at zero.
//'   
//'   In \code{C++} the end points must be shifted by \code{-1} because indexing
//'   starts at \code{0}: \code{-1, 4, 9, 14, 19}.  But there is no vector
//'   element corresponding to index \code{-1}. So in \code{C++} we cannot
//'   select vector elements using the end points starting at \code{-1}. The
//'   solution is to drop the first placeholder end point.
//'   
//' @examples
//' # Calculate end points without a stub interval
//' HighFreq::calc_endpoints(25, 5)
//' # Calculate end points with initial stub interval
//' HighFreq::calc_endpoints(23, 5)
//' # Calculate end points with a stub interval at the end
//' HighFreq::calc_endpoints(23, 5, FALSE)
//'
//' @export
// [[Rcpp::export]]
arma::uvec calc_endpoints(arma::uword len_gth, arma::uword ste_p = 1, bool front = true) {
  
  // Calculate number of intervals that fit over len_gth
  arma::uword num_points = len_gth/ste_p;
  arma::uvec end_p;
  
  if (len_gth == ste_p*num_points) {
    // No stub interval
    // Include the first placeholder end point - legacy code
    // end_p = arma::cumsum(arma::ones<uvec>(num_points));
    // end_p = arma::regspace<uvec>(0, ste_p, len_gth);
    end_p = arma::regspace<uvec>(ste_p, ste_p, len_gth);
  } else {
    // Need to add stub interval
    // Include the first placeholder end point - legacy code
    // end_p = arma::cumsum(arma::ones<uvec>(num_points + 1));
    // end_p = arma::regspace<uvec>(0, ste_p, len_gth + ste_p);
    end_p = arma::regspace<uvec>(ste_p, ste_p, len_gth + ste_p);
    if (front) {
      // Stub interval at beginning
      end_p = end_p - ste_p + len_gth % ste_p;
    } else {
      // Stub interval at end
      // The last end point must be equal to len_gth
      end_p(num_points) = len_gth;
    }  // end if
  }  // end if
  
  // Set the first end point to zero - it's a placeholder
  // end_p(0) = 0;
  // Subtract 1 from end_p because indexing starts at 0
  end_p = end_p - 1;
  return end_p;
  
}  // end calc_endpoints




////////////////////////////////////////////////////////////
//' Calculate a vector of start points equal to the lag of a vector of end
//' points.
//'
//' @param \code{end_points} An \emph{unsigned integer} vector of end
//' points.
//'   
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   lag applied to the end points.
//'   
//' @return An \emph{integer} vector of start points (vector of unsigned
//'   integers), associated with the vector \code{end_points}.
//'
//' @details The start points are equal to the values of the vector
//'   \code{end_points} lagged by an amount equal to \code{look_back}. In
//'   addition, an extra \code{1} is added to them, to avoid data overlaps.  The
//'   lag operation requires adding a beginning warmup interval containing
//'   zeros, so that the vector of start points has the same length as the
//'   \code{end_points}.
//'   
//'   For example, consider the end points for a vector of length \code{25}
//'   divided into equal intervals of length \code{5}: \code{4, 9, 14, 19, 24}.
//'   (In \code{C++} the vector indexing is shifted by \code{-1} and starts at
//'   \code{0} not \code{1}.)
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
arma::uvec calc_startpoints(arma::uvec end_points, arma::uword look_back) {
  
  arma::uword num_points = end_points.n_elem;
  arma::uvec start_p = arma::join_cols(arma::zeros<uvec>(look_back), 
                                       end_points.subvec(0, num_points - look_back - 1) + 1);
  
  return start_p;
  
}  // end calc_startpoints


////////////////////////////////////////////////////////////
//' Multiply in place (without copying) the columns or rows of a \emph{matrix}
//' times a \emph{vector}, element-wise.
//' 
//' @param \code{vec_tor} A \emph{vector}.
//' 
//' @param \code{mat_rix} A \emph{matrix}.
//' 
//' @param \code{by_col} A \emph{Boolean} argument: if \code{TRUE} then multiply
//'   the columns, otherwise multiply the rows. (The default is
//'   \code{by_col = TRUE}.)
//' 
//' @return A single \emph{integer} value, equal to either the number of
//'   \emph{matrix} columns or the number of rows.
//' 
//' @details The function \code{mult_vec_mat()} multiplies the columns or rows
//'   of a \emph{matrix} times a \emph{vector}, element-wise.
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
arma::uword mult_vec_mat(const arma::vec& vec_tor,
                   arma::mat& mat_rix,
                   const bool& by_col = true) {
  arma::uword num_elem = vec_tor.n_elem;
  arma::uword num_rows = mat_rix.n_rows;
  arma::uword num_cols = mat_rix.n_cols;
  
  if ((num_cols == num_rows) && (num_elem == num_rows)) {
    if (by_col) {
      // Multiply each column of mat_rix by vec_tor
      mat_rix.each_col() %= vec_tor;
      return num_rows;
    } else {
      // Multiply each row of mat_rix by vec_tor
      mat_rix.each_row() %= conv_to< rowvec >::from(vec_tor);
      return num_cols;
    }
  } else if (num_elem == num_rows) {
    // Multiply each column of mat_rix by vec_tor
    mat_rix.each_col() %= vec_tor;
    return num_rows;
  } else if (num_elem == num_cols) {
    // Multiply each row of mat_rix by vec_tor
    mat_rix.each_row() %= conv_to< rowvec >::from(vec_tor);
    return num_cols;
  } else 
    stop("Error: Vector length is neither equal to the number of columns nor rows of the matrix!");
    // Return NA_INTEGER;
    
}  // end mult_vec_mat




////////////////////////////////////////////////////////////
//' Calculate the eigen decomposition of the covariance \emph{matrix} of returns
//' using \code{RcppArmadillo}.
//' 
//' @param \code{re_turns} A \emph{time series} or \emph{matrix} of returns
//'   data.
//'
//' @return A list with two elements: a \emph{vector} of eigenvalues 
//'   (named "values"), and a \emph{matrix} of eigenvectors (named
//'   "vectors").
//'
//' @details The function \code{calc_eigen()} first calculates the covariance 
//'   \emph{matrix} of the returns, and then calculates its eigen decomposition.
//'
//' @examples
//' \dontrun{
//' # Create matrix of random returns
//' re_turns <- matrix(rnorm(5e6), nc=5)
//' # Calculate eigen decomposition
//' ei_gen <- HighFreq::calc_eigen(scale(re_turns, scale=FALSE))
//' # Calculate PCA
//' pc_a <- prcomp(re_turns)
//' # Compare PCA with eigen decomposition
//' all.equal(pc_a$sdev^2, drop(ei_gen$values))
//' all.equal(abs(unname(pc_a$rotation)), abs(ei_gen$vectors))
//' # Compare the speed of Rcpp with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_eigen(re_turns),
//'   Rcode=prcomp(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
List calc_eigen(const arma::mat& re_turns) {
  arma::mat eigen_vec;
  arma::vec eigen_val;
  arma::eig_sym(eigen_val, eigen_vec, cov(re_turns));
  // Reverse the order of elements from largest eigenvalue to smallest, similar to R
  return List::create(Named("values") = arma::flipud(eigen_val),
                      Named("vectors") = arma::fliplr(eigen_vec));
}  // end calc_eigen



////////////////////////////////////////////////////////////
//' Calculate the regularized inverse of the covariance \emph{matrix} of returns
//' using \code{RcppArmadillo}.
//' 
//' @param \code{re_turns} A \emph{time series} or \emph{matrix} of returns data.
//' 
//' @param \code{to_l} A \emph{numeric} tolerance level for discarding small
//'   eigenvalues in order to regularize the matrix inverse.  (The default is
//'   \code{0.001})
//'   
//' @param \code{max_eigen} An \emph{integer} equal to the regularization
//'   intensity (the number of eigenvalues and eigenvectors used for calculating
//'   the regularized inverse).
//'
//' @return A \emph{matrix} equal to the regularized inverse. 
//'
//' @details The function calc_inv() calculates the regularized inverse of the
//'   \emph{covariance matrix}, by discarding eigenvectors with small
//'   eigenvalues less than the tolerance level \code{to_l}.
//'   The function \code{calc_inv()} first calculates the covariance
//'   \emph{matrix} of the \code{re_turns}, and then it calculates its
//'   regularized inverse.
//'   If \code{max_eigen} is not specified then it calculates the
//'   regularized inverse using the function \code{arma::pinv()}.
//'   If \code{max_eigen} is specified then it calculates the regularized
//'   inverse using eigen decomposition with only the largest \code{max_eigen}
//'   eigenvalues and their corresponding eigenvectors.
//'
//' @examples
//' \dontrun{
//' # Create random matrix
//' re_turns <- matrix(rnorm(500), nc=5)
//' max_eigen <- 3
//' # Calculate regularized inverse using RcppArmadillo
//' in_verse <- HighFreq::calc_inv(re_turns, max_eigen)
//' # Calculate regularized inverse from eigen decomposition in R
//' ei_gen <- eigen(cov(re_turns))
//' inverse_r <-  ei_gen$vectors[, 1:max_eigen] %*% (t(ei_gen$vectors[, 1:max_eigen]) / ei_gen$values[1:max_eigen])
//' # Compare RcppArmadillo with R
//' all.equal(in_verse, inverse_r)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& re_turns,
                   double to_l = 0.001, 
                   int max_eigen = 0) {
  
  arma::mat cov_mat = cov(re_turns);
  
  if (max_eigen == 0) {
    // Calculate the inverse using arma::pinv()
    return arma::pinv(cov_mat, to_l);
  } else {
    // Calculate the inverse using eigen decomposition
    arma::mat eigen_vec;
    arma::vec eigen_val;
    arma::eig_sym(eigen_val, eigen_vec, cov_mat);
    eigen_vec = eigen_vec.cols(eigen_vec.n_cols-max_eigen, eigen_vec.n_cols-1);
    eigen_val = 1/eigen_val.subvec(eigen_val.n_elem-max_eigen, eigen_val.n_elem-1);
    return eigen_vec*diagmat(eigen_val)*eigen_vec.t();
  }  // end if
  
}  // end calc_inv




////////////////////////////////////////////////////////////
//' Scale (standardize) the columns of a \emph{matrix} of data using
//' \code{RcppArmadillo}.
//' 
//' @param \code{mat_rix} A \emph{matrix} of data.
//' @param \code{use_median} A \emph{Boolean} argument: if \code{TRUE} then the 
//'   centrality (central tendency) is calculated as the \emph{median} and the 
//'   dispersion is calculated as the \emph{median absolute deviation}
//'   (\emph{MAD}).
//'   If \code{use_median} is \code{FALSE} then the centrality is calculated as 
//'   the \emph{mean} and the dispersion is calculated as the \emph{standard
//'   deviation}. (The default is \code{FALSE})
//'
//' @return A \emph{matrix} with the same dimensions as the input
//'   argument \code{mat_rix}.
//'
//' @details The function \code{calc_scaled()} scales (standardizes) the columns
//'   of the \code{mat_rix} argument using \code{RcppArmadillo}.
//'   If the argument \code{use_median} is \code{FALSE} (the default), then it
//'   performs the same calculation as the standard \code{R} function
//'   \code{scale()}, and it calculates the centrality (central tendency) as the
//'   \emph{mean} and the dispersion as the \emph{standard deviation}.
//'   If the argument \code{use_median} is \code{TRUE}, then it calculates the
//'   centrality as the \emph{median} and the dispersion as the \emph{median
//'   absolute deviation} (\emph{MAD}).
//'   
//'   The function \code{calc_scaled()} uses \code{RcppArmadillo} \code{C++}
//'   code and is about \emph{5} times faster than function \code{scale()}, for
//'   a \emph{matrix} with \emph{1,000} rows and \emph{20} columns.
//'   
//' @examples
//' \dontrun{
//' # Create a matrix of random data
//' mat_rix <- matrix(rnorm(20000), nc=20)
//' scale_d <- calc_scaled(mat_rix=mat_rix, use_median=FALSE)
//' scale_d2 <- scale(mat_rix)
//' all.equal(scale_d, scale_d2, check.attributes=FALSE)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=calc_scaled(mat_rix=mat_rix, use_median=FALSE),
//'   Rcode=scale(mat_rix),
//'   times=100))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_scaled(const arma::mat& mat_rix, 
                      const bool use_median=false) {
  arma::mat scaled_mat(mat_rix.n_rows, mat_rix.n_cols);
  arma::vec scale_d(mat_rix.n_rows);
  double cen_ter;
  
  // Perform a loop over the columns
  for (arma::uword it=0; it < mat_rix.n_cols; it++) {
    if (use_median) {
      cen_ter = arma::median(mat_rix.col(it));
      scale_d = (mat_rix.col(it) - cen_ter);
      scale_d = scale_d/arma::median(arma::abs(scale_d));
      scaled_mat.col(it) = scale_d;
    } else {
      cen_ter = arma::mean(mat_rix.col(it));
      scale_d = (mat_rix.col(it) - cen_ter);
      scale_d = scale_d/arma::stddev(scale_d);
      scaled_mat.col(it) = scale_d;
    }  // end if
  }  // end for
  
  return scaled_mat;
  
}  // end calc_scaled




////////////////////////////////////////////////////////////
// Functions for statistics
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
//' Calculate the variance of a a single-column \emph{time series} or a
//' \emph{vector} using \code{RcppArmadillo}.
//' 
//' @param \code{t_series} A single-column \emph{time series} or a \emph{vector}.
//'
//' @return A \emph{numeric} value equal to the variance of the \emph{vector}.
//'
//' @details The function \code{calc_var_vec()} calculates the variance of a
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
double calc_var_vec(arma::vec& t_series) {
  return arma::var(t_series);
}  // end calc_var_vec




////////////////////////////////////////////////////////////
//' Calculate the variance of the columns of a \emph{time series} or a
//' \emph{matrix} using \code{RcppArmadillo}.
//' 
//' @param \code{t_series} A \emph{time series} or a \emph{matrix} of data.
//'
//' @return A row vector equal to the variance of the columns of \code{t_series}
//'   matrix.
//'
//' @details The function \code{calc_var()} calculates the variance of the
//'   columns of a \emph{time series} or a \emph{matrix} of data using
//'   \code{RcppArmadillo} \code{C++} code.
//'   
//'   The function \code{calc_var()} performs the same calculation as the
//'   function \code{colVars()} from package
//'   \href{https://cran.r-project.org/web/packages/matrixStats/index.html}{matrixStats},
//'   but it's much faster because it uses \code{RcppArmadillo} \code{C++} code.
//'
//' @examples
//' \dontrun{
//' # Create a matrix of random returns
//' re_turns <- matrix(rnorm(5e6), nc=5)
//' # Compare calc_var() with standard var()
//' all.equal(drop(HighFreq::calc_var(re_turns)), 
//'   apply(re_turns, 2, var))
//' # Compare calc_var() with matrixStats
//' all.equal(drop(HighFreq::calc_var(re_turns)), 
//'   matrixStats::colVars(re_turns))
//' # Compare the speed of RcppArmadillo with matrixStats and with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_var(re_turns),
//'   matrixStats=matrixStats::colVars(re_turns),
//'   Rcode=apply(re_turns, 2, var),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::rowvec calc_var(arma::mat& t_series) {
  
  return arma::var(t_series);
  
}  // end calc_var




////////////////////////////////////////////////////////////
//' Calculate the variance of an \emph{OHLC time series}, using different range
//' estimators and \code{RcppArmadillo}.
//'
//' @param \code{oh_lc} An \emph{OHLC time series} or a \emph{numeric matrix} of
//'   prices.
//'   
//' @param \code{calc_method} A \emph{character} string representing the range
//'   estimator for calculating the variance.  The estimators include:
//'   \itemize{
//'     \item "close" close-to-close estimator,
//'     \item "rogers_satchell" Rogers-Satchell estimator,
//'     \item "garman_klass" Garman-Klass estimator,
//'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
//'     \item "yang_zhang" Yang-Zhang estimator,
//'    }
//'    (The default is the \code{calc_method = "yang_zhang"}.)
//'    
//' @param \code{lag_close} A \emph{vector} with the lagged \emph{close} prices
//'   of the \emph{OHLC time series}.  This is an optional argument. (The
//'   default is \code{lag_close = 0}.)
//'   
//' @param \code{scal_e} \emph{Boolean} argument: Should the returns be divided
//'   by the time index, the number of seconds in each period? (The default is
//'   \code{scal_e = TRUE}.)
//'
//' @param \code{in_dex} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument. (The default is \code{in_dex = 0}.)
//'   
//' @return A single \emph{numeric} value equal to the variance of the
//'   \emph{OHLC time series}.
//'
//' @details The function \code{calc_var_ohlc()} calculates the variance
//'   from all the different intra-day and day-over-day returns (defined as the
//'   differences of \emph{OHLC} prices), using several different variance
//'   estimation methods.
//'
//'   The default \code{calc_method} is \emph{"yang_zhang"}, which theoretically
//'   has the lowest standard error among unbiased estimators.
//'   The methods \emph{"close"}, \emph{"garman_klass_yz"}, and
//'   \emph{"yang_zhang"} do account for \emph{close-to-open} price jumps, while
//'   the methods \emph{"garman_klass"} and \emph{"rogers_satchell"} do not
//'   account for \emph{close-to-open} price jumps.
//'
//'   If \code{scal_e} is \code{TRUE} (the default), then the returns are
//'   divided by the differences of the time index (which scales the variance to
//'   the units of variance per second squared.) This is useful when calculating
//'   the variance from minutely bar data, because dividing returns by the
//'   number of seconds decreases the effect of overnight price jumps. If the
//'   time index is in days, then the variance is equal to the variance per day
//'   squared.
//'   
//'   The optional argument \code{in_dex} is the time index of the \emph{time
//'   series} \code{oh_lc}. If the time index is in seconds, then the
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
//' # Extract time index of SPY returns
//' in_dex <- c(1, diff(xts::.index(HighFreq::SPY)))
//' # Calculate the variance of SPY returns, with scaling of the returns
//' HighFreq::calc_var_ohlc(HighFreq::SPY, 
//'  calc_method="yang_zhang", scal_e=TRUE, in_dex=in_dex)
//' # Calculate variance without accounting for overnight jumps
//' HighFreq::calc_var_ohlc(HighFreq::SPY, 
//'  calc_method="rogers_satchell", scal_e=TRUE, in_dex=in_dex)
//' # Calculate the variance without scaling the returns
//' HighFreq::calc_var_ohlc(HighFreq::SPY, scal_e=FALSE)
//' # Calculate the variance by passing in the lagged close prices
//' lag_close <- HighFreq::lag_it(HighFreq::SPY[, 4])
//' all.equal(HighFreq::calc_var_ohlc(HighFreq::SPY), 
//'   HighFreq::calc_var_ohlc(HighFreq::SPY, lag_close=lag_close))
//' # Compare with HighFreq::calc_var_ohlc_r()
//' all.equal(HighFreq::calc_var_ohlc(HighFreq::SPY, in_dex=in_dex), 
//'   HighFreq::calc_var_ohlc_r(HighFreq::SPY))
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_var_ohlc(HighFreq::SPY),
//'   Rcode=HighFreq::calc_var_ohlc_r(HighFreq::SPY),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
double calc_var_ohlc(arma::mat& oh_lc, 
                     const std::string& calc_method="yang_zhang", 
                     arma::colvec lag_close=0, 
                     const bool& scal_e = true, 
                     arma::colvec in_dex=0) {
  
  int num_rows = oh_lc.n_rows;
  double co_eff = 0.34/(1.34 + (num_rows+1)/(num_rows-1));
  
  if (!scal_e || (in_dex.n_rows == 1)) {
    in_dex = arma::ones(num_rows);
    // cout << "oh_lc.n_rows = " << num_rows << endl;
    // cout << "in_dex.n_rows = " << in_dex.n_rows << endl;
  }  // end if
  
  // Calculate all the different intra-day and day-over-day returns 
  // (differences of OHLC prices)
  arma::colvec clo_se = oh_lc.col(3);
  arma::colvec open_close(clo_se.n_rows);
  if (lag_close.n_rows == 1) {
    open_close = arma::join_cols(clo_se.subvec(0, 0), clo_se.subvec(0, clo_se.n_elem-2));
    open_close = (oh_lc.col(0) - open_close)/in_dex;
  } else {
    open_close = (oh_lc.col(0) - lag_close)/in_dex;
  }  // end if
  arma::colvec close_open = (clo_se - oh_lc.col(0))/in_dex;
  arma::colvec close_high = (clo_se - oh_lc.col(1))/in_dex;
  arma::colvec close_low = (clo_se - oh_lc.col(2))/in_dex;
  arma::colvec high_low = (oh_lc.col(1) - oh_lc.col(2))/in_dex;
  arma::colvec high_open = (oh_lc.col(1) - oh_lc.col(0))/in_dex;
  arma::colvec low_open = (oh_lc.col(2) - oh_lc.col(0))/in_dex;
  
  if (calc_method == "close") {
    // cout << "Calc method is Close" << endl;
    return arma::var(arma::diff(clo_se));
  } else if (calc_method == "rogers_satchell") {
    // cout << "Calc method is Rogers-Satchell" << endl;
    return -(arma::dot(close_high, high_open) +
             arma::dot(close_low, low_open))/num_rows;
  } else if (calc_method == "garman_klass") {
    // cout << "Calc method is Garman-Klass" << endl;
    return (0.5*arma::dot(high_low, high_low) -
            (2*log(2)-1)*arma::dot(close_open, close_open))/num_rows;
  } else if (calc_method == "garman_klass_yz") {
    // cout << "Calc method is Garman-Klass-YZ" << endl;
    return (0.5*arma::dot(high_low, high_low) -
            (2*log(2)-1)*arma::dot(close_open, close_open))/num_rows + 
            arma::var(open_close);
  } else if (calc_method == "yang_zhang") {
    // cout << "Calc method is Yang-Zhang" << endl;
    return arma::var(open_close) + co_eff*arma::var(close_open) +
      (co_eff-1)*(arma::dot(close_high, high_open) + 
      arma::dot(close_low, low_open))/num_rows;
  } else {
    cout << "Wrong calc method!" << endl;
    return 1;
  }  // end if
  
  // cout << "Calc method is " << calc_method << endl;
  
}  // end calc_var_ohlc



////////////////////////////////////////////////////////////
//' Calculate the ranks of the elements of a single-column \emph{time series} or
//' a \emph{vector} using \code{RcppArmadillo}.
//' 
//' @param \code{vec_tor} A single-column \emph{time series} or a \emph{vector}.
//'
//' @return An \emph{integer vector} with the ranks of the elements of the
//'   \emph{vector}.
//'
//' @details The function \code{calc_ranks()} calculates the ranks of the
//'   elements of a single-column \emph{time series} or a \emph{vector}.
//'   It uses the \code{RcppArmadillo} function \code{arma::sort_index()}.
//'   The function \code{arma::sort_index()} calculates the permutation index to
//'   sort a given vector into ascending order.
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
arma::uvec calc_ranks(const arma::vec& vec_tor) {
  
  return (arma::sort_index(arma::sort_index(vec_tor)) + 1);
  
}  // end calc_ranks




////////////////////////////////////////////////////////////
//' Calculate the Median Absolute Deviations (\emph{MAD}) of the columns of a
//' \emph{time series} or a \emph{matrix} using \code{RcppArmadillo}.
//'
//' @param \code{t_series} A \emph{time series} or a \emph{matrix} of data.
//'
//' @return A single-row matrix with the Median Absolute Deviations \emph{MAD}
//'   of the columns of \code{t_series}.
//'
//' @details The function \code{calc_mad()} calculates the Median Absolute
//'   Deviations \emph{MAD} of the columns of a \emph{time series} or a
//'   \emph{matrix} of data using \code{RcppArmadillo} \code{C++} code.
//'
//'   The function \code{calc_mad()} performs the same calculation as the
//'   function \code{stats::mad()}, but it's much faster because it uses
//'   \code{RcppArmadillo} \code{C++} code.
//'
//' @examples
//' \dontrun{
//' # Calculate VTI returns
//' re_turns <- na.omit(rutils::etf_env$re_turns[ ,"VTI", drop=FALSE])
//' # Compare calc_mad() with stats::mad()
//' all.equal(drop(HighFreq::calc_mad(re_turns)), 
//'   mad(re_turns)/1.4826)
//' # Compare the speed of RcppArmadillo with stats::mad()
//' library(microbenchmark)
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_mad(re_turns),
//'   Rcode=mad(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_mad(arma::mat& t_series) {
  
  // De-median the columns of t_series
  t_series.each_row() -= arma::median(t_series);
  
  return arma::median(arma::abs(t_series));
  
}  // end calc_mad



////////////////////////////////////////////////////////////
// Switch statement in calc_skew() uses C++ enum type.
// This is needed because Rcpp can't map C++ enum type to R variable SEXP.
enum skew_type {Pearson, Quantile, Nonparametric};
// Map string to C++ enum type for switch statement.
skew_type calc_skew_type(const std::string& typ_e) {
  if (typ_e == "Pearson" || typ_e == "pearson" || typ_e == "p") 
    return skew_type::Pearson;
  else if (typ_e == "Quantile" || typ_e == "quantile" || typ_e == "q")
    return skew_type::Quantile;
  else if (typ_e == "Nonparametric" || typ_e == "nonparametric" || typ_e == "n")
    return skew_type::Nonparametric;
  else 
    return skew_type::Pearson;
}  // end calc_skew_type



////////////////////////////////////////////////////////////
//' Calculate the skewness of the columns of a \emph{time series} or a
//' \emph{matrix} using \code{RcppArmadillo}.
//'
//' @param \code{t_series} A \emph{time series} or a \emph{matrix} of data.
//'
//' @param \code{typ_e} A \emph{string} specifying the type of skewness (see
//'   Details). (The default is the \code{typ_e = "pearson"}.)
//'
//' @param \code{al_pha} The confidence level for calculating the quantiles.
//'   (the default is \code{al_pha = 0.25}).
//'
//' @return A single-row matrix with the skewness of the columns of
//'   \code{t_series}.
//'
//' @details The function \code{calc_skew()} calculates the skewness of the
//'   columns of a \emph{time series} or a \emph{matrix} of data using
//'   \code{RcppArmadillo} \code{C++} code.
//'
//'   If \code{typ_e = "pearson"} (the default) then \code{calc_skew()}
//'   calculates the Pearson skewness using the third moment of the data.
//'
//'   If \code{typ_e = "quantile"} then it calculates the skewness using the
//'   differences between the quantiles of the data.
//'
//'   If \code{typ_e = "nonparametric"} then it calculates the skewness as the
//'   difference between the mean of the data minus its median, divided by the
//'   standard deviation.
//'   
//'   The code examples below compare the function \code{calc_skew()} with the
//'   skewness calculated using \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Calculate VTI returns
//' re_turns <- na.omit(rutils::etf_env$re_turns[ ,"VTI", drop=FALSE])
//' # Calculate the Pearson skewness
//' HighFreq::calc_skew(re_turns)
//' # Compare HighFreq::calc_skew() with Pearson skewness
//' calc_skewr <- function(x) {
//'   x <- (x-mean(x)); nr <- NROW(x);
//'   nr*sum(x^3)/(var(x))^1.5/(nr-1)/(nr-2)
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
//' HighFreq::calc_skew(re_turns, typ_e = "quantile", al_pha = 0.1)
//' # Compare HighFreq::calc_skew() with quantile skewness
//' calc_skewq <- function(x) {
//'   	quantile_s <- quantile(x, c(0.25, 0.5, 0.75), type=5)
//'   	(quantile_s[3] + quantile_s[1] - 2*quantile_s[2])/(quantile_s[3] - quantile_s[1])
//' }  # end calc_skewq
//' all.equal(drop(HighFreq::calc_skew(re_turns, typ_e = "quantile")), 
//'   calc_skewq(re_turns), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(re_turns, typ_e = "quantile"),
//'   Rcode=calc_skewq(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' # Calculate the nonparametric skewness
//' HighFreq::calc_skew(re_turns, typ_e = "nonparametric")
//' # Compare HighFreq::calc_skew() with R nonparametric skewness
//' all.equal(drop(HighFreq::calc_skew(re_turns, typ_e = "nonparametric")), 
//'   (mean(re_turns)-median(re_turns))/sd(re_turns), 
//'   check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with R code
//' summary(microbenchmark(
//'   Rcpp=HighFreq::calc_skew(re_turns, typ_e = "nonparametric"),
//'   Rcode=(mean(re_turns)-median(re_turns))/sd(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_skew(arma::mat t_series,
                    const std::string& typ_e = "pearson", 
                    double al_pha = 0.25) {
  
  // switch statement for all the different types of skew
  switch(calc_skew_type(typ_e)) {
  case skew_type::Pearson: {  // Pearson
    double num_rows = t_series.n_rows;
    arma::mat mean_s = arma::mean(t_series);
    arma::mat var_s = arma::var(t_series);
    // De-mean the columns of t_series
    t_series.each_row() -= mean_s;
    return (num_rows/(num_rows-1)/(num_rows-2))*arma::sum(arma::pow(t_series, 3))/arma::pow(var_s, 1.5);
  }  // end pearson
  case skew_type::Quantile: {  // Quantile
    arma::vec prob_s = {al_pha, 0.5, 1.0 - al_pha};
    arma::mat quantile_s = quantile(t_series, prob_s);
    return (quantile_s.row(2) + quantile_s.row(0) - 2*quantile_s.row(1))/(quantile_s.row(2) - quantile_s.row(0));
  }  // end quantile
  case skew_type::Nonparametric: {  // Nonparametric
    return (arma::mean(t_series) - arma::median(t_series))/arma::stddev(t_series);
  }  // end nonparametric
  default : {
    cout << "Invalid typ_e" << endl;
    return 0;
  }  // end default
  }  // end switch
  
}  // end calc_skew



// wippp
// These below are special cases of calc_skew()
// These below need documentation
//' @export
// [[Rcpp::export]]
arma::mat calc_skew_pearson(arma::mat& t_series) {
  
  double num_rows = t_series.n_rows;
  
  arma::mat mean_s = arma::mean(t_series);
  arma::mat var_s = arma::var(t_series);
  
  // De-mean the columns of t_series
  t_series.each_row() -= mean_s;
  
  return (num_rows/(num_rows-1)/(num_rows-2))*arma::sum(arma::pow(t_series, 3))/arma::pow(var_s, 1.5);
  
}  // end calc_skew_pearson



//' @export
// [[Rcpp::export]]
arma::mat calc_skew_quant(arma::mat& t_series, double al_pha = 0.25) {
  
  arma::vec prob_s = {al_pha, 0.5, 1.0 - al_pha};
  arma::mat quantile_s = quantile(t_series, prob_s);
  
  // arma::mat skew_s = (quantile_s.row(2) + quantile_s.row(0) - 2*quantile_s.row(1))/(quantile_s.row(2) - quantile_s.row(0));
  // skew_s /= calc_mad(t_series);
  // skew_s /= (quantile_s.row(2) - quantile_s.row(0));
  
  return (quantile_s.row(2) + quantile_s.row(0) - 2*quantile_s.row(1))/(quantile_s.row(2) - quantile_s.row(0));
  
}  // end calc_skew_quant




//' @export
// [[Rcpp::export]]
arma::mat calc_skew_nonp(arma::mat& t_series) {
  
  return (arma::mean(t_series) - arma::median(t_series))/arma::stddev(t_series);
  
}  // end calc_skew_nonp


// wippp end



////////////////////////////////////////////////////////////
//' Perform multivariate linear regression using \emph{Rcpp}.
//' 
//' @param \code{res_ponse} A \emph{vector} of response data.
//' @param \code{de_sign} A \emph{matrix} of design (predictor i.e.
//'   explanatory) data.
//' 
//' @return A named list with three elements: a \emph{matrix} of coefficients
//'   (named \emph{"coefficients"}), the \emph{z-score} of the last residual
//'   (named \emph{"z_score"}), and a \emph{vector} with the R-squared and
//'   F-statistic (named \emph{"stats"}). The numeric \emph{matrix} of
//'   coefficients named \emph{"coefficients"} containes the alpha and beta
//'   coefficients, and their \emph{t-values} and \emph{p-values}.
//'
//' @details The function \code{calc_lm()} performs the same calculations as the
//'   function \code{lm()} from package \emph{stats}. It uses
//'   \code{RcppArmadillo} \code{C++} code and is about \emph{10} times faster
//'   than \code{lm()}. The code was inspired by this article (but it's not
//'   identical to it):
//'   http://gallery.rcpp.org/articles/fast-linear-model-with-armadillo/
//'
//' @examples
//' \dontrun{
//' # Define design matrix with explanatory variables
//' len_gth <- 100; n_var <- 5
//' de_sign <- matrix(rnorm(n_var*len_gth), nc=n_var)
//' # Response equals linear form plus error terms
//' weight_s <- rnorm(n_var)
//' res_ponse <- -3 + de_sign %*% weight_s + rnorm(len_gth, sd=0.5)
//' # Perform multivariate regression using lm()
//' reg_model <- lm(res_ponse ~ de_sign)
//' sum_mary <- summary(reg_model)
//' # Perform multivariate regression using calc_lm()
//' reg_model_arma <- calc_lm(res_ponse=res_ponse, de_sign=de_sign)
//' reg_model_arma$coefficients
//' # Compare the outputs of both functions
//' all.equal(reg_model_arma$coefficients[, "coeff"], unname(coef(reg_model)))
//' all.equal(unname(reg_model_arma$coefficients), unname(sum_mary$coefficients))
//' all.equal(drop(reg_model_arma$residuals), unname(reg_model$residuals))
//' all.equal(unname(reg_model_arma$stats), c(sum_mary$r.squared, unname(sum_mary$fstatistic[1])))
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List calc_lm(const arma::vec& res_ponse, const arma::mat& de_sign) {
  // Add column for intercept to explanatory matrix
  arma::mat design_p = join_rows(ones(de_sign.n_rows), de_sign);
  int num_rows = de_sign.n_rows, num_cols = design_p.n_cols;
  int deg_free = (num_rows - num_cols);
  
  // fit the model res_ponse ~ de_sign, and calculate alpha and beta coefficients
  arma::colvec co_eff = arma::solve(design_p, res_ponse);
  // Calculate residuals
  arma::colvec resid_uals = res_ponse - design_p*co_eff;
  
  // Calculate TSS, RSS, and ESS
  double tot_sumsq = (num_rows-1)*arma::var(res_ponse);
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
  arma::colvec beta_tvals = co_eff/std_err;
  arma::colvec beta_pvals = 2*Rcpp::pt(-abs(wrap(beta_tvals)), deg_free);
  NumericMatrix coeff_mat = wrap(join_rows(join_rows(join_rows(co_eff, std_err), beta_tvals), beta_pvals));
  Rcpp::colnames(coeff_mat) = Rcpp::CharacterVector::create("coeff", "std_err", "tvals", "pvals");
  
  return Rcpp::List::create(Named("coefficients") = coeff_mat,
                            // Named("residuals") = resid_uals,
                            Named("z_score") = resid_uals(num_rows-1)/arma::stddev(resid_uals),
                            Named("stats") = stat_s);
}  // end calc_lm




////////////////////////////
// Functions for rolling statistics
////////////////////////////


////////////////////////////////////////////////////////////
//' Aggregate a time series of data into a single bar of \emph{OHLC} data.
//'
//' @export
//' @param \code{t_series} A \emph{time series} or a \emph{matrix} with multiple
//'   columns of data.
//'
//' @return A \emph{matrix} containing a single row, with the \emph{open},
//'   \emph{high}, \emph{low}, and \emph{close} values, and also the total
//'   \emph{volume} (if provided as either the second or fifth column of
//'   \code{t_series}).
//'
//' @details The function \code{agg_ohlc()} aggregates a time series of data
//'   into a single bar of \emph{OHLC} data.
//'   It can accept either a single column of data or four columns of
//'   \emph{OHLC} data.
//'   It can also accept an additional column containing the trading volume.
//'   
//' The function \code{agg_ohlc()} calculates the \emph{open} value as equal to
//' the \emph{open} value of the first row of \code{t_series}.
//'   The \emph{high} value as the maximum of the \emph{high} column of
//'   \code{t_series}.
//'   The \emph{low} value as the minimum of the \emph{low} column of
//'   \code{t_series}.
//'   The \emph{close} value as the \emph{close} of the last row of
//'   \code{t_series}.
//'   The \emph{volume} value as the sum of the \emph{volume} column of
//'   \code{t_series}.
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
arma::mat agg_ohlc(const arma::mat& t_series) {
  
  int num_rows = t_series.n_rows;
  int num_cols = t_series.n_cols;
  
  // Number of output columns
  int num_ohlc = num_cols;
  if (num_cols < 4)
    // Add volume column for non-OHLC data
    num_ohlc = 4 + num_cols - 1;
  // Allocate output matrix
  arma::mat oh_lc(1, num_ohlc);
  
  if (num_cols < 4) {
    // Aggregate time series into a single bar of OHLC data.
    oh_lc(0, 0) = t_series(0, 0);
    oh_lc(0, 1) = arma::max(t_series.col(0));
    oh_lc(0, 2) = arma::min(t_series.col(0));
    oh_lc(0, 3) = t_series(num_rows-1, 0);
    if (num_cols == 2) {
      // Aggregate volume data.
      oh_lc(0, 4) = arma::sum(t_series.col(1));
    }  // end if
  } else {
    // Aggregate OHLC time series into a single bar of OHLC data.
    oh_lc(0, 0) = t_series(0, 0);
    oh_lc(0, 1) = arma::max(t_series.col(1));
    oh_lc(0, 2) = arma::min(t_series.col(2));
    oh_lc(0, 3) = t_series(num_rows-1, 3);
    if (num_cols == 5) {
      // Aggregate volume data.
      oh_lc(0, 4) = arma::sum(t_series.col(4));
    }  // end if
  }  // end if
  
  return oh_lc;
  
}  // end agg_ohlc




////////////////////////////////////////////////////////////
//' Count the number of consecutive \code{TRUE} elements in a Boolean vector,
//' and reset the count to zero after every \code{FALSE} element.
//' 
//' @param vec_tor A \emph{Boolean vector} of data.
//'
//' @return An \emph{integer vector} of the same length as the argument
//'   \code{vec_tor}.
//'
//' @details The function \code{roll_count()} calculates the number of
//'   consecutive \code{TRUE} elements in a Boolean vector, and it resets the
//'   count to zero after every \code{FALSE} element.  
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
arma::uvec roll_count(arma::uvec& vec_tor) {
  
  arma::uword len_gth = vec_tor.n_elem;
  arma::uvec count_true(len_gth);
  
  // Initialize count
  count_true[0] = vec_tor[0];
  // Loop over vec_tor
  for (arma::uword it = 1; it < len_gth; it++) {
    if (vec_tor[it])
      // Add count number
      count_true[it] = count_true[it-1] + 1;
    else
      // Reset count to zero
      count_true[it] = vec_tor[it];
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
//' @param \code{t_series} A \emph{time series} or a \emph{matrix} with multiple
//'   columns of data.
//'   
//' @param \emph{end_points} An \emph{integer vector} of end points.
//'
//' @return A \emph{matrix} with \emph{OHLC} data, with the number of rows equal
//'   to the number of \emph{end_points} minus one.
//'   
//' @details The function \code{roll_ohlc()} performs a loop over the
//'   \emph{end_points}, along the rows of the \code{t_series} data. At each
//'   \emph{end_point}, it selects the past rows of \code{t_series} data,
//'   starting at the first bar after the previous \emph{end_point}, and then
//'   calls the function \code{agg_ohlc()} on the selected \code{t_series} data
//'   to calculate the aggregations.
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
//' end_points <- rutils::calc_endpoints(oh_lc, inter_val=25)
//' # Aggregate over end_points:
//' ohlc_agg <- HighFreq::roll_ohlc(t_series=oh_lc, end_points=(end_points-1))
//' # Compare with xts::to.period()
//' ohlc_agg_xts <- .Call("toPeriod", oh_lc, as.integer(end_points), TRUE, NCOL(oh_lc), FALSE, FALSE, colnames(oh_lc), PACKAGE="xts")
//' all.equal(ohlc_agg, coredata(ohlc_agg_xts), check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_ohlc(arma::mat& t_series, arma::uvec& end_points) {
  
  // int num_rows = t_series.n_rows;
  int num_cols = t_series.n_cols;
  arma::uword num_points = end_points.size();
  arma::mat ohlc_agg(num_points-1, num_cols);
  
  // Perform loop over the end points
  for (arma::uword it = 1; it < num_points; it++) {
    // cout << "it: " << it << endl;
    // Aggregate the OHLC
    ohlc_agg.row(it-1) = agg_ohlc(t_series.rows(end_points(it-1)+1, end_points(it)));
  }  // end for
  
  // Return the aggregations
  return ohlc_agg;
  
}  // end roll_ohlc



////////////////////////////////////////////////////////////
//' Calculate the rolling sum over a single-column \emph{time series} or a
//' \emph{vector} using \emph{Rcpp}.
//' 
//' @param \code{t_series} A single-column \emph{time series} or a \emph{vector}.
//' 
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of elements of data used for calculating the sum.
//'
//' @return A column \emph{vector} of the same length as the argument
//'   \code{t_series}.
//'
//' @details The function \code{roll_vec()} calculates a \emph{vector} of
//'   rolling sums, over a \emph{vector} of data, using fast \emph{Rcpp}
//'   \code{C++} code.  The function \code{roll_vec()} is several times faster
//'   than \code{rutils::roll_sum()} which uses vectorized \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random returns
//' re_turns <- rnorm(1e6)
//' # Calculate rolling sums over 11-period lookback intervals
//' sum_rolling <- HighFreq::roll_vec(re_turns, look_back=11)
//' # Compare HighFreq::roll_vec() with rutils::roll_sum()
//' all.equal(HighFreq::roll_vec(re_turns, look_back=11), 
//'          rutils::roll_sum(re_turns, look_back=11))
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
NumericVector roll_vec(NumericVector t_series, int look_back) {
  
  int len_gth = t_series.size();
  NumericVector rolling_sum(len_gth);

  // Warmup period
  rolling_sum[0] = t_series[0];
  for (int it = 1; it < look_back; it++) {
    rolling_sum[it] = rolling_sum[it-1] + t_series[it];
  }  // end for
  
  // Remaining period
  for (int it = look_back; it < len_gth; it++) {
    rolling_sum[it] = rolling_sum[it-1] + t_series[it] - t_series[it-look_back];
  }  // end for
  
  return rolling_sum;
  
}  // end roll_vec




////////////////////////////////////////////////////////////
//' Calculate the rolling weighted sum over a single-column \emph{time series}
//' or a \emph{vector} using \code{RcppArmadillo}.
//' 
//' @param \code{t_series} A single-column \emph{time series} or a \emph{vector}.
//' 
//' @param \code{weight_s} A \emph{vector} of weights.
//'
//' @return A column \emph{vector} of the same length as the argument
//'   \code{t_series}.
//'
//' @details The function \code{roll_vecw()} calculates the rolling weighted sum
//'   of a \emph{vector} over its past values (a convolution with the
//'   \emph{vector} of weights), using \code{RcppArmadillo}. It performs a
//'   similar calculation as the standard \code{R} function
//'   \code{stats::filter(x=t_series, filter=weight_s, method="convolution",
//'   sides=1)}, but it's over \code{6} times faster, and it doesn't produce any
//'   \code{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Create vector from historical prices
//' re_turns <- as.numeric(rutils::etf_env$VTI[, 6])
//' # Create simple weights
//' weight_s <- c(1, rep(0, 10))
//' # Calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_vecw(t_series=re_turns, weight_s=weight_s)
//' # Compare with original
//' all.equal(re_turns, as.numeric(weight_ed))
//' # Second example
//' # Create exponentially decaying weights
//' weight_s <- exp(-0.2*1:11)
//' weight_s <- weight_s/sum(weight_s)
//' # Calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_vecw(t_series=re_turns, weight_s=weight_s)
//' # Calculate rolling weighted sum using filter()
//' filter_ed <- stats::filter(x=re_turns, filter=weight_s, method="convolution", sides=1)
//' # Compare both methods
//' all.equal(filter_ed[-(1:11)], weight_ed[-(1:11)], check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec roll_vecw(const arma::vec& t_series, const arma::vec& weight_s) {
  
  arma::uword len_gth = t_series.n_elem;
  arma::uword look_back = weight_s.n_elem;
  arma::vec rolling_sum(len_gth);
  arma::vec rev_weights = arma::reverse(weight_s);
  // arma::vec rev_weights = weight_s;
  
  // Warmup period
  rolling_sum.subvec(0, look_back-2) = t_series.subvec(0, look_back-2);
  
  // Remaining periods
  for (arma::uword it = look_back-1; it < len_gth; it++) {
    rolling_sum(it) = arma::dot(rev_weights, t_series.subvec(it-look_back+1, it));
  }  // end for
  
  return rolling_sum;
  
}  // end roll_vecw





////////////////////////////////////////////////////////////
//' Calculate the convolutions of the \emph{matrix} columns with a \emph{vector}
//' of weights.
//' 
//' @param \code{mat_rix} A \emph{matrix} of data.
//' @param \code{weight_s} A column \emph{vector} of weights.
//'
//' @return A \emph{matrix} with the same dimensions as the input
//'   argument \code{mat_rix}.
//'
//' @details The function \code{roll_conv()} calculates the convolutions of the
//'   \emph{matrix} columns with a \emph{vector} of weights.  It performs a loop
//'   down over the \emph{matrix} rows and multiplies the past (higher) values
//'   by the weights.  It calculates the rolling weighted sum of the past
//'   values.
//'   
//'   The function \code{roll_conv()} uses the \code{RcppArmadillo} function
//'   \code{arma::conv2()}. It performs a similar calculation to the standard
//'   \code{R} function \code{filter(x=mat_rix, filter=weight_s,
//'   method="convolution", sides=1)}, but it's over \code{6} times faster, and
//'   it doesn't produce any leading \code{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Create matrix from historical prices
//' mat_rix <- na.omit(rutils::etf_env$re_turns[, 1:2])
//' # Create simple weights equal to a 1 value plus zeros
//' weight_s <- matrix(c(1, rep(0, 10)), nc=1)
//' # Calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_conv(mat_rix, weight_s)
//' # Compare with original
//' all.equal(coredata(mat_rix), weight_ed, check.attributes=FALSE)
//' # Second example
//' # Create exponentially decaying weights
//' weight_s <- exp(-0.2*(1:11))
//' weight_s <- matrix(weight_s/sum(weight_s), nc=1)
//' # Calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_conv(mat_rix, weight_s)
//' # Calculate rolling weighted sum using filter()
//' filter_ed <- filter(x=mat_rix, filter=weight_s, method="convolution", sides=1)
//' # Compare both methods
//' all.equal(filter_ed[-(1:11), ], weight_ed[-(1:11), ], check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_conv(arma::mat& mat_rix, arma::mat& weight_s) {
  
  arma::uword look_back = weight_s.n_rows-2;
  arma::uword num_rows = mat_rix.n_rows-1;
  
  // Calculate the convolutions
  arma::mat roll_conv = arma::conv2(mat_rix, weight_s, "full");
  
  // Copy the warmup period
  roll_conv.rows(0, look_back) = mat_rix.rows(0, look_back);
  
  return roll_conv.rows(0, num_rows);
  
}  // end roll_conv




////////////////////////////////////////////////////////////
//' Calculate the convolutions of the \emph{matrix} columns with a \emph{vector}
//' of weights.
//' 
//' @param \code{mat_rix} A \emph{matrix} of data.
//' @param \code{weight_s} A column \emph{vector} of weights.
//'
//' @return A \emph{matrix} with the same dimensions as the input
//'   argument \code{mat_rix}.
//'
//' @details The function \code{roll_conv_ref()} calculates the convolutions of
//'   the \emph{matrix} columns with a \emph{vector} of weights.  It performs a
//'   loop down over the \emph{matrix} rows and multiplies the past (higher)
//'   values by the weights.  It calculates the rolling weighted sum of the past
//'   values.
//'   
//'   The function \code{roll_conv_ref()} accepts a \emph{pointer} to the argument
//'   \code{mat_rix}, and replaces the old \emph{matrix} values with the
//'   weighted sums.
//'   It performs the calculation in place, without copying the \emph{matrix} in
//'   memory (which greatly increases the computation speed).
//'   
//'   The function \code{roll_conv_ref()} uses the \code{RcppArmadillo} function
//'   \code{arma::conv2()}. It performs a similar calculation to the standard
//'   \code{R} function \code{filter(x=mat_rix, filter=weight_s,
//'   method="convolution", sides=1)}, but it's over \code{6} times faster, and
//'   it doesn't produce any leading \code{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Create matrix from historical prices
//' mat_rix <- na.omit(rutils::etf_env$re_turns[, 1:2])
//' # Create simple weights equal to a 1 value plus zeros
//' weight_s <- matrix(c(1, rep(0, 10)), nc=1)
//' # Calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_conv_ref(mat_rix, weight_s)
//' # Compare with original
//' all.equal(coredata(mat_rix), weight_ed, check.attributes=FALSE)
//' # Second example
//' # Create exponentially decaying weights
//' weight_s <- exp(-0.2*(1:11))
//' weight_s <- matrix(weight_s/sum(weight_s), nc=1)
//' # Calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_conv_ref(mat_rix, weight_s)
//' # Calculate rolling weighted sum using filter()
//' filter_ed <- filter(x=mat_rix, filter=weight_s, method="convolution", sides=1)
//' # Compare both methods
//' all.equal(filter_ed[-(1:11), ], weight_ed[-(1:11), ], check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_conv_ref(arma::mat& mat_rix, arma::mat& weight_s) {
  arma::uword look_back = weight_s.n_rows-1;
  arma::uword num_rows = mat_rix.n_rows-1;
  
  // Calculate the convolutions
  arma::mat roll_conv_ref = arma::conv2(mat_rix, weight_s, "full");
  
  // Copy the convolutions
  mat_rix.rows(look_back, num_rows) = roll_conv_ref.rows(look_back, num_rows);
  
  return mat_rix;
}  // end roll_conv_ref




////////////////////////////////////////////////////////////
//' Calculate the rolling weighted sum over a \emph{time series} or a
//' \emph{matrix} using \emph{Rcpp}.
//' 
//' @param \code{t_series} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of data points included in calculating the rolling sum (the default
//'   is \code{look_back = 1}).
//'   
//' @param \code{stu_b} An \emph{integer} value equal to the first stub interval
//'   for calculating the end points.
//' 
//' @param \code{end_points} An \emph{unsigned integer} vector of end
//' points.
//'   
//' @param \code{weight_s} A column \emph{vector} of weights.
//'
//' @return A \emph{matrix} with the same dimensions as the input
//'   argument \code{t_series}.
//'
//' @details The function \code{roll_sum()} calculates the rolling sums over the
//'   columns of the \code{t_series} data.  
//'   The sums are calculated over a number of data points equal to
//'   \code{look_back}.
//'   
//'   The function \code{roll_sum()} returns a \emph{matrix} with the same
//'   dimensions as the input argument \code{t_series}.
//' 
//'   The arguments \code{stu_b}, \code{end_points}, and \code{weight_s} are
//'   optional.
//'   
//'   If either the arguments \code{stu_b} or \code{end_points} are supplied,
//'   then the rolling sums are calculated at the end points. 
//'   
//'   If only the argument \code{stu_b} is supplied, then the end points are
//'   calculated from the \code{stu_b} and \code{look_back} arguments. The first
//'   end point is equal to \code{stu_b} and the end points are spaced
//'   \code{look_back} periods apart.
//'   
//'   If the argument \code{weight_s} is supplied, then weighted sums are
//'   calculated.
//'   Then the function \code{roll_sum()} calculates the rolling weighted sums
//'   of the past values.
//'   
//'   The function \code{roll_sum()} calculates the rolling weighted sums as
//'   convolutions of the \code{t_series} columns with the \emph{vector} of
//'   weights using the \code{RcppArmadillo} function \code{arma::conv2()}.
//'   It performs a similar calculation to the standard \code{R} function
//'   \code{stats::filter(x=t_series, filter=weight_s, method="convolution",
//'   sides=1)}, but it's over \code{6} times faster, and it doesn't produce any
//'   leading \code{NA} values. using fast \emph{RcppArmadillo} \code{C++} code.
//'   The function \code{roll_sum()} is several times faster than
//'   \code{rutils::roll_sum()} which uses vectorized \code{R} code.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Create series of historical returns
//' re_turns <- na.omit(rutils::etf_env$re_turns[, c("VTI", "IEF")])
//' # Define parameters
//' look_back <- 22
//' stu_b <- 21
//' # Calculate rolling sums at each point
//' c_sum <- HighFreq::roll_sum(re_turns, look_back=look_back)
//' r_sum <- rutils::roll_sum(re_turns, look_back=look_back)
//' all.equal(c_sum, coredata(r_sum), check.attributes=FALSE)
//' r_sum <- apply(zoo::coredata(re_turns), 2, cumsum)
//' lag_sum <- rbind(matrix(numeric(2*look_back), nc=2), r_sum[1:(NROW(r_sum) - look_back), ])
//' r_sum <- (r_sum - lag_sum)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' 
//' # Calculate rolling sums at end points
//' c_sum <- HighFreq::roll_sum(re_turns, look_back=look_back, stu_b=stu_b)
//' end_p <- (stu_b + look_back*(0:(NROW(re_turns) %/% look_back)))
//' end_p <- end_p[end_p < NROW(re_turns)]
//' r_sum <- apply(zoo::coredata(re_turns), 2, cumsum)
//' r_sum <- r_sum[end_p+1, ]
//' lag_sum <- rbind(numeric(2), r_sum[1:(NROW(r_sum) - 1), ])
//' r_sum <- (r_sum - lag_sum)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' 
//' # Calculate rolling sums at end points - pass in end_points
//' c_sum <- HighFreq::roll_sum(re_turns, end_points=end_p)
//' all.equal(c_sum, r_sum, check.attributes=FALSE)
//' 
//' # Create exponentially decaying weights
//' weight_s <- exp(-0.2*(1:11))
//' weight_s <- matrix(weight_s/sum(weight_s), nc=1)
//' # Calculate rolling weighted sum
//' c_sum <- HighFreq::roll_sum(re_turns, weight_s=weight_s)
//' # Calculate rolling weighted sum using filter()
//' filter_ed <- filter(x=re_turns, filter=weight_s, method="convolution", sides=1)
//' all.equal(c_sum[-(1:11), ], filter_ed[-(1:11), ], check.attributes=FALSE)
//' 
//' # Calculate rolling weighted sums at end points
//' c_sum <- HighFreq::roll_sum(re_turns, end_points=end_p, weight_s=weight_s)
//' all.equal(c_sum, filter_ed[end_p+1, ], check.attributes=FALSE)
//' 
//' # Create simple weights equal to a 1 value plus zeros
//' weight_s <- matrix(c(1, rep(0, 10)), nc=1)
//' # Calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_sum(re_turns, weight_s)
//' # Compare with original
//' all.equal(coredata(re_turns), weight_ed, check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_sum(arma::mat& t_series,
                   arma::uword look_back = 1,
                   Rcpp::Nullable<int> stu_b = R_NilValue, 
                   Rcpp::Nullable<Rcpp::IntegerVector> end_points = R_NilValue, 
                   Rcpp::Nullable<Rcpp::NumericVector> weight_s = R_NilValue) {
  
  arma::uword num_rows = t_series.n_rows;
  arma::mat cum_sum;
  
  if (weight_s.isNotNull()) {
    // Copy weight_s
    arma::vec weights_vec = Rcpp::as<vec>(weight_s);
    arma::uword num_weights = weights_vec.n_elem;
    // Calculate the weighted averages as convolutions
    cum_sum = arma::conv2(t_series, weights_vec, "full");
    // Copy the warmup period
    // cout << "num_weights = " << num_weights << endl;
    cum_sum.rows(0, num_weights-2) = t_series.rows(0, num_weights-2);
    cum_sum = cum_sum.rows(0, num_rows-1);
    // cout << "cum_sum.n_rows = " << cum_sum.n_rows << endl;
  } else {
    // Calculate cumulative returns
    cum_sum = arma::cumsum(t_series, 0);
  }  // end if
  
  
  // Declare empty end points
  arma::uvec end_p;
  // Update end points
  if (end_points.isNotNull()) {
    // Copy end_points
    end_p = Rcpp::as<uvec>(end_points);
  } else if (stu_b.isNotNull()) {
    // Calculate end points with stu_b
    end_p = arma::regspace<uvec>(Rcpp::as<uword>(stu_b), look_back, num_rows + look_back);
    end_p = end_p.elem(find(end_p < num_rows));
  }  // end if
  
  
  // Calculate the rolling sums
  if (end_p.is_empty() && weight_s.isNotNull()) {
    // Do nothing
    // Return the weighted averages (convolutions) at each point
    // return cum_sum;
  } else if (end_p.is_empty() && !weight_s.isNotNull()) {
    // Return rolling sums at each point
    cum_sum = diff_it(cum_sum, look_back, true);
  } else if (!end_p.is_empty() && weight_s.isNotNull()) {
    // Return the weighted averages (convolutions) at end points
    cum_sum = cum_sum.rows(end_p);
  } else if (!end_p.is_empty() && !weight_s.isNotNull()) {
    // Return the rolling sums at end points
    cum_sum = cum_sum.rows(end_p);
    cum_sum = diff_it(cum_sum, 1, true);
  }  // end if
  
  return cum_sum;
  
}  // end roll_sum




////////////////////////////////////////////////////////////
//' Calculate a \emph{vector} of variance estimates over a rolling look-back
//' interval for a single-column \emph{time series} or a \emph{vector}, using
//' \code{RcppArmadillo}.
//'
//' @param \code{t_series} A single-column \emph{time series} or a \emph{vector}.
//' 
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of \emph{vector} elements used for calculating a single variance
//'   estimate.
//'
//' @return A column \emph{vector} with the same number of elements as the input
//'   argument \code{t_series}.
//'
//' @details The function \code{roll_var_vec()} calculates a \emph{vector} of
//'   variance estimates over a rolling look-back interval for a single-column
//'   \emph{time series} or a \emph{vector}, using \code{RcppArmadillo}
//'   \code{C++} code.
//'   
//'   The function \code{roll_var_vec()} uses an expanding look-back interval in
//'   the initial warmup period, to calculate the same number of elements as the
//'   input argument \code{t_series}.
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
//' # Compare the variance estimates over 11-period lookback intervals
//' all.equal(drop(HighFreq::roll_var_vec(re_turns, look_back=11))[-(1:10)], 
//'   RcppRoll::roll_var(re_turns, n=11))
//' # Compare the speed of RcppArmadillo with RcppRoll
//' library(microbenchmark)
//' summary(microbenchmark(
//'   RcppArmadillo=HighFreq::roll_var_vec(re_turns, look_back=11),
//'   RcppRoll=RcppRoll::roll_var(re_turns, n=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::vec roll_var_vec(arma::vec& t_series, arma::uword look_back=11) {
  arma::uword len_gth = t_series.n_elem;
  arma::vec var_vec = arma::zeros(len_gth);
  
  // Warmup period
  for (arma::uword it = 1; it < look_back; it++) {
    var_vec(it) = arma::var(t_series.subvec(0, it));
  }  // end for
  
  // Remaining period
  for (arma::uword it = look_back; it < len_gth; it++) {
    var_vec(it) = arma::var(t_series.subvec(it-look_back+1, it));
  }  // end for
  
  return var_vec;
  
}  // end roll_var_vec



////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of variance estimates over a rolling look-back
//' interval attached at the end points of a \emph{time series} or a
//' \emph{matrix}.
//'
//' @param \code{t_series} A \emph{time series} or a \emph{matrix}.
//' 
//' @param \code{ste_p} The number of time periods between the end points.
//'
//' @param \code{look_back} The number of end points in the look-back interval.
//'
//' @return A \emph{matrix} with the same number of columns as the input time
//'   series \code{t_series}, and the number of rows equal to the number of end
//'   points.
//'
//' @details The function \code{roll_var()} calculates a \emph{matrix} of
//'   variance estimates over rolling look-back intervals attached at the end
//'   points of the \emph{time series} \code{t_series}.
//'   
//'   The end points are calculated along the rows of \code{t_series} using the
//'   function \code{calc_endpoints()}, with the number of time periods between
//'   the end points equal to \code{ste_p}.
//'   
//'   At each end point, the variance is calculated over a look-back interval
//'   equal to \code{look_back} number of end points.
//'   In the initial warmup period, the variance is calculated over an expanding
//'   look-back interval.
//'   
//'   For example, the rolling variance at \code{25} day end points, with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{ste_p = 25} and \code{look_back = 3}.
//'
//'   The function \code{roll_var()} with the parameter \code{ste_p = 1}
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
//' vari_ance <- HighFreq::roll_var(re_turns, ste_p=25, look_back=3)
//' # Compare the variance estimates over 11-period lookback intervals
//' all.equal(HighFreq::roll_var(re_turns, look_back=11)[-(1:10), ], 
//'   drop(RcppRoll::roll_var(re_turns, n=11)), check.attributes=FALSE)
//' # Compare the speed of RcppArmadillo with RcppRoll
//' library(microbenchmark)
//' summary(microbenchmark(
//'   RcppArmadillo=HighFreq::roll_var(re_turns, look_back=11),
//'   RcppRoll=RcppRoll::roll_var(re_turns, n=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_var(arma::mat& t_series, arma::uword ste_p = 1, arma::uword look_back = 1) {
  
  // Calculate end points
  arma::uword num_rows = t_series.n_rows;
  arma::uvec end_p = calc_endpoints(num_rows, ste_p);
  // Start points equal to end points lagged by look_back
  arma::uvec start_p = calc_startpoints(end_p, look_back);
  // Allocate variance matrix
  arma::uword num_points = end_p.n_elem;
  arma::mat vari_ance = arma::zeros(num_points, t_series.n_cols);

  // Perform loop over the end_points
  for (arma::uword ep = 0; ep < num_points; ep++) {
    // Calculate variance
    if (end_p(ep) > start_p(ep)) {
      vari_ance.row(ep) = arma::var(t_series.rows(start_p(ep), end_p(ep)));
    }  // end if
  }  // end for

  // Old code below
  // Warmup period
  // for (arma::uword it = 1; it < look_back; it++) {
  //   vari_ance.row(it) = arma::var(t_series.rows(0, it));
  // }  // end for
  
  // Remaining period
  // for (arma::uword it = look_back; it < num_rows; it++) {
  //   vari_ance.row(it) = arma::var(t_series.rows(it-look_back+1, it));
  // }  // end for
  
  return vari_ance;
  
}  // end roll_var



////////////////////////////////////////////////////////////
//' Calculate a \emph{vector} of variance estimates over a rolling look-back
//' interval attached at the end points of a \emph{time series} or a
//' \emph{matrix} with \emph{OHLC} price data.
//' 
//' @param \code{oh_lc} A \emph{time series} or a \emph{matrix} with \emph{OHLC}
//'   price data.
//'   
//' @param \code{ste_p} The number of time periods between the end points.
//'
//' @param \code{look_back} The number of end points in the look-back interval.
//'   
//' @param \code{calc_method} A \emph{character} string representing the range
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
//' @param \code{scal_e} \emph{Boolean} argument: Should the returns be divided
//'   by the time index, the number of seconds in each period? (The default is
//'   \code{scal_e = TRUE}.)
//'   
//' @param \code{in_dex} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument. (The default is \code{in_dex=0}.)
//'
//' @return A column \emph{vector} of variance estimates, with the number of
//'   rows equal to the number of end points.
//'
//' @details The function \code{roll_var_ohlc()} calculates a \emph{vector} of
//'   variance estimates over a rolling look-back interval attached at the end
//'   points of the \emph{time series} \code{oh_lc}.
//'   
//'   The end points are calculated along the rows of \code{oh_lc} using the
//'   function \code{calc_endpoints()}, with the number of time periods between
//'   the end points equal to \code{ste_p}.
//'   
//'   The function \code{roll_var_ohlc()} performs a loop over the end points,
//'   subsets the previous (past) rows of \code{oh_lc}, and passes them into the
//'   function \code{calc_var_ohlc()}.
//' 
//'   At each end point, the variance is calculated over a look-back interval
//'   equal to \code{look_back} number of end points.
//'   In the initial warmup period, the variance is calculated over an expanding
//'   look-back interval.
//'   
//'   For example, the rolling variance at daily end points with an \code{11}
//'   day look-back, can be calculated using the parameters \code{ste_p = 1} and
//'   \code{look_back = 11} (Assuming the \code{oh_lc} data has daily
//'   frequency.)
//' 
//'   Similarly, the rolling variance at \code{25} day end points with a
//'   \code{75} day look-back, can be calculated using the parameters
//'   \code{ste_p = 25} and \code{look_back = 3} (because \code{3*25 = 75}).
//' 
//'   The function \code{roll_var_ohlc()} calculates the variance from all the
//'   different intra-day and day-over-day returns (defined as the differences
//'   between \emph{OHLC} prices), using several different variance estimation
//'   methods.
//'   
//'   The default \code{calc_method} is \emph{"yang_zhang"}, which theoretically
//'   has the lowest standard error among unbiased estimators.
//'   The methods \emph{"close"}, \emph{"garman_klass_yz"}, and
//'   \emph{"yang_zhang"} do account for \emph{close-to-open} price jumps, while
//'   the methods \emph{"garman_klass"} and \emph{"rogers_satchell"} do not
//'   account for \emph{close-to-open} price jumps.
//'
//'   If \code{scal_e} is \code{TRUE} (the default), then the returns are
//'   divided by the differences of the time index (which scales the variance to
//'   the units of variance per second squared.) This is useful when calculating
//'   the variance from minutely bar data, because dividing returns by the
//'   number of seconds decreases the effect of overnight price jumps. If the
//'   time index is in days, then the variance is equal to the variance per day
//'   squared.
//'   
//'   The optional argument \code{in_dex} is the time index of the \emph{time
//'   series} \code{oh_lc}. If the time index is in seconds, then the
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
//' # Extract time index of SPY returns
//' oh_lc <- HighFreq::SPY
//' in_dex <- c(1, diff(xts::.index(oh_lc)))
//' # Rolling variance at minutely end points, with a 21 minute look-back
//' var_rolling <- HighFreq::roll_var_ohlc(oh_lc, 
//'                               ste_p=1, look_back=21, 
//'                               calc_method="yang_zhang", 
//'                               in_dex=in_dex, scal_e=TRUE)
//' # Daily OHLC prices
//' oh_lc <- rutils::etf_env$VTI
//' in_dex <- c(1, diff(xts::.index(oh_lc)))
//' # Rolling variance at 5 day end points, with a 20 day look-back (20=4*5)
//' var_rolling <- HighFreq::roll_var_ohlc(oh_lc, 
//'                               ste_p=5, look_back=4, 
//'                               calc_method="yang_zhang", 
//'                               in_dex=in_dex, scal_e=TRUE)
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
//'   HighFreq::calc_var_ohlc(sub_ohlc, lag_close=sub_close, scal_e=TRUE, in_dex=sub_index)
//' })  # end sapply
//' var_rollingr <- c(0, var_rollingr)
//' all.equal(drop(var_rolling), var_rollingr)
//' }
//' @export
// [[Rcpp::export]]
arma::vec roll_var_ohlc(arma::mat& oh_lc, 
                        arma::uword ste_p = 1, 
                        arma::uword look_back = 1, 
                        const std::string& calc_method = "yang_zhang", 
                        const bool& scal_e = true, 
                        arma::colvec in_dex = 0) {

  // Calculate end points
  arma::uword num_rows = oh_lc.n_rows;
  arma::uvec end_p = calc_endpoints(num_rows, ste_p);
  // Start points equal to end points lagged by look_back
  arma::uvec start_p = calc_startpoints(end_p, look_back);
  
  // Allocate variance matrix
  arma::uword num_points = end_p.n_elem;
  arma::vec vari_ance = arma::zeros(num_points);
  
  // Extract OHLC close prices
  arma::colvec clo_se = oh_lc.col(3);
  arma::colvec lag_close = lag_it(clo_se);

  // Set the time index to 1 if scal_e = FALSE
  if (!scal_e || (in_dex.n_rows == 1)) {
    in_dex = arma::ones(num_rows);
  }  // end if

  // Define data subsets over look-back intervals
  arma::mat sub_ohlc;
  arma::colvec sub_close;
  arma::colvec sub_index;
  
  // Perform loop over the end_points
  for (arma::uword ep = 0; ep < num_points; ep++) {
    if (end_p(ep) > start_p(ep)) {
      sub_ohlc = oh_lc.rows(start_p(ep), end_p(ep));
      sub_close = lag_close.rows(start_p(ep), end_p(ep));
      sub_index = in_dex.subvec(start_p(ep), end_p(ep));
      // Calculate variance
      vari_ance.row(ep) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, scal_e, sub_index);
    }  // end if
  }  // end for

  // Old code below

  // Warmup period
  // for (arma::uword it = 1; it < look_back; it++) {
  //   arma::mat sub_ohlc = oh_lc.rows(0, it);
  //   arma::colvec sub_close = lag_close.rows(0, it);
  //   arma::colvec sub_index = in_dex.subvec(0, it);
  //   vari_ance(it) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, scal_e, sub_index);
  // }  // end for
  
  // Remaining period
  // for (arma::uword it = look_back; it < num_rows; it++) {
  //   arma::mat sub_ohlc = oh_lc.rows(it-look_back+1, it);
  //   arma::colvec sub_close = lag_close.rows(it-look_back+1, it);
  //   arma::colvec sub_index = in_dex.subvec(it-look_back+1, it);
  //   vari_ance(it) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, scal_e, sub_index);
  // }  // end for
  
  return vari_ance;
  
}  // end roll_var_ohlc



////////////////////////////////////////////////////////////
//' Perform a rolling scaling (standardization) of the columns of a
//' \emph{matrix} of data using \code{RcppArmadillo}.
//' 
//' @param mat_rix A \emph{matrix} of data.
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
//'   deviation}. (The default is \code{use_median = FALSE})
//'
//' @return A \emph{matrix} with the same dimensions as the input argument
//'   \code{mat_rix}.
//'
//' @details The function \code{roll_scale()} performs a rolling scaling
//'   (standardization) of the columns of the \code{mat_rix} argument using
//'   \code{RcppArmadillo}.
//'   The function \code{roll_scale()} performs a loop over the rows of 
//'   \code{mat_rix}, subsets a number of previous (past) rows equal to 
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
//' rolled_scaled2 <- roll_scale(mat_rix=mat_rix, look_back = look_back, use_median=FALSE)
//' all.equal(rolled_scaled[-1, ], rolled_scaled2[-1, ])
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat roll_scale(const arma::mat& mat_rix, 
                     const arma::uword& look_back,
                     const bool use_median=false) {
  
  arma::uword num_rows = mat_rix.n_rows;
  arma::mat scaled_mat(num_rows, mat_rix.n_cols);
  arma::mat sub_mat;
  
  // Warmup period
  scaled_mat.row(0) = mat_rix.row(0);
  for (arma::uword it = 1; it < look_back; it++) {
    sub_mat = mat_rix.rows(0, it);
    sub_mat = calc_scaled(sub_mat, use_median);
    scaled_mat.row(it) = sub_mat.row(sub_mat.n_rows-1);
  }  // end for
  
  // Remaining periods
  for (arma::uword it = look_back; it < num_rows; it++) {
    sub_mat = mat_rix.rows(it-look_back+1, it);
    sub_mat = calc_scaled(sub_mat, use_median);
    scaled_mat.row(it) = sub_mat.row(sub_mat.n_rows-1);
  }  // end for
  
  return scaled_mat;
}  // end roll_scale



////////////////////////////////////////////////////////////
//' Perform rolling regressions over the rows of the design matrix, and
//' calculate a \emph{vector} of z-scores of the residuals.
//' 
//' @param \code{res_ponse} A \emph{vector} of response data.
//' 
//' @param \code{de_sign} A \emph{matrix} of design (predictor i.e.
//'   explanatory) data.
//'   
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of elements of data used for calculating the regressions.
//'
//' @return A column \emph{vector} of the same length as the number of rows of
//'   \code{de_sign}.
//'
//' @details The function \code{roll_zscores()} performs rolling regressions
//'   along the rows of the design \emph{matrix} \code{de_sign}, using the
//'   function \code{calc_lm()}.
//'   
//'   The function \code{roll_zscores()} performs a loop over the rows of
//'   \code{de_sign}, and it subsets \code{de_sign} and \code{res_ponse} over a
//'   number of previous (past) rows equal to \code{look_back}.  It performs a
//'   regression on the subset data, and calculates the \emph{z-score} of the
//'   last residual value for each regression. It returns a numeric
//'   \emph{vector} of the \emph{z-scores}.
//'   
//' @examples
//' \dontrun{
//' # Calculate Z-scores from rolling time series regression using RcppArmadillo
//' look_back <- 11
//' clo_se <- as.numeric(Cl(rutils::etf_env$VTI))
//' date_s <- xts::.index(rutils::etf_env$VTI)
//' z_scores <- HighFreq::roll_zscores(res_ponse=clo_se, 
//'  de_sign=matrix(as.numeric(date_s), nc=1), 
//'  look_back=look_back)
//' # Define design matrix with explanatory variables
//' len_gth <- 100; n_var <- 5
//' de_sign <- matrix(rnorm(n_var*len_gth), nc=n_var)
//' # response equals linear form plus error terms
//' weight_s <- rnorm(n_var)
//' res_ponse <- -3 + de_sign %*% weight_s + rnorm(len_gth, sd=0.5)
//' # Calculate Z-scores from rolling multivariate regression using RcppArmadillo
//' look_back <- 11
//' z_scores <- HighFreq::roll_zscores(res_ponse=res_ponse, de_sign=de_sign, look_back=look_back)
//' # Calculate z-scores in R from rolling multivariate regression using lm()
//' z_scores_r <- sapply(1:NROW(de_sign), function(ro_w) {
//'   if (ro_w == 1) return(0)
//'   start_point <- max(1, ro_w-look_back+1)
//'   sub_response <- res_ponse[start_point:ro_w]
//'   sub_design <- de_sign[start_point:ro_w, ]
//'   reg_model <- lm(sub_response ~ sub_design)
//'   resid_uals <- reg_model$residuals
//'   resid_uals[NROW(resid_uals)]/sd(resid_uals)
//' })  # end sapply
//' # Compare the outputs of both functions
//' all.equal(unname(z_scores[-(1:look_back)]), 
//'   unname(z_scores_r[-(1:look_back)]))
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec roll_zscores(const arma::vec& res_ponse, 
                       const arma::mat& de_sign, 
                       const arma::uword& look_back) {
  
  arma::uword num_rows = de_sign.n_rows;
  arma::vec z_scores(num_rows);
  arma::vec sub_response;
  arma::mat sub_design;

  // Warmup period
  for (arma::uword it = 1; it < look_back; it++) {
    sub_response = res_ponse.subvec(0, it);
    sub_design = de_sign.rows(0, it);
    z_scores(it) = calc_lm(sub_response, sub_design)["z_score"];
  }  // end for
  
  // Remaining periods
  for (arma::uword it = look_back; it < num_rows; it++) {
    sub_response = res_ponse.subvec(it-look_back+1, it);
    sub_design = de_sign.rows(it-look_back+1, it);
    z_scores(it) = calc_lm(sub_response, sub_design)["z_score"];
  }  // end for
  
  return z_scores;
}  // end roll_zscores




////////////////////////////
// Functions for simulation
////////////////////////////



////////////////////////////////////////////////////////////
//' Simulate a \emph{GARCH} process using \emph{Rcpp}.
//' 
//' @param \code{om_ega} Parameter proportional to the long-term average level
//'   of variance.
//' @param \code{al_pha} The weight associated with recent realized variance
//'   updates.
//' @param \code{be_ta} The weight associated with the past variance estimates.
//' @param \code{in_nov} A \emph{vector} of innovations (random numbers).
//' 
//' @return A \emph{matrix} with two columns: the simulated returns and
//'   variance, and with the same number of rows as the length of the argument 
//'   \code{in_nov}.
//'
//' @details The function \code{sim_garch()} simulates a \emph{GARCH} process
//'   using fast \emph{Rcpp} \code{C++} code.
//'
//' @examples
//' \dontrun{
//' # Define the GARCH model parameters
//' om_ega <- 0.01
//' al_pha <- 0.5
//' be_ta <- 0.2
//' # Simulate the GARCH process using Rcpp
//' garch_rcpp <- sim_garch(om_ega=om_ega, al_pha=al_pha, be_ta=be_ta, in_nov=rnorm(10000))
//' }
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix sim_garch(double om_ega, 
                        double al_pha, 
                        double be_ta, 
                        NumericVector in_nov) {
  
  int len_gth = in_nov.size();
  NumericVector vari_ance(len_gth);
  NumericVector re_turns(len_gth);
  vari_ance[0] = om_ega/(1-al_pha-be_ta);
  re_turns[0] = sqrt(vari_ance[0])*in_nov[0];
  
  for (int it = 1; it < len_gth; it++) {
    re_turns[it] = sqrt(vari_ance[it-1])*in_nov[it];
    vari_ance[it] = om_ega + al_pha*pow(re_turns[it], 2) + be_ta*vari_ance[it-1];
  }  // end for
  
  return cbind(re_turns, vari_ance);
  
}  // end sim_garch



////////////////////////////////////////////////////////////
//' Simulate an \emph{Ornstein-Uhlenbeck} process using \emph{Rcpp}.
//' 
//' @param \code{eq_price} The equilibrium price. 
//' @param \code{vol_at} The volatility of returns.
//' @param \code{the_ta} The strength of mean reversion.
//' @param \code{in_nov} A \emph{vector} of innovations (random numbers).
//' 
//' @return A column \emph{vector} representing the \emph{time series} of
//'   prices, with the same length as the argument \code{in_nov}.
//'
//' @details The function \code{sim_ou()} simulates an \emph{Ornstein-Uhlenbeck}
//'   process using fast \emph{Rcpp} \code{C++} code.
//'   It returns a column \emph{vector} representing the \emph{time series} of
//'   prices.
//'
//' @examples
//' \dontrun{
//' # Define the Ornstein-Uhlenbeck model parameters
//' eq_price <- 5.0
//' vol_at <- 0.01
//' the_ta <- 0.01
//' # Simulate Ornstein-Uhlenbeck process using Rcpp
//' price_s <- HighFreq::sim_ou_rcpp(eq_price=eq_price, vol_at=vol_at, the_ta=the_ta, in_nov=rnorm(1000))
//' }
//' 
//' @export
// [[Rcpp::export]]
NumericVector sim_ou(double eq_price, 
                     double vol_at, 
                     double the_ta, 
                     NumericVector in_nov) {
  int len_gth = in_nov.size();
  NumericVector price_s(len_gth);
  NumericVector re_turns(len_gth);
  
  price_s[0] = eq_price;
  for (int it = 1; it < len_gth; it++) {
    re_turns[it] = the_ta*(eq_price - price_s[it-1]) + vol_at*in_nov[it-1];
    price_s[it] = price_s[it-1] * exp(re_turns[it]);
  }  // end for
  return price_s;
}  // end sim_ou



////////////////////////////////////////////////////////////
//' Recursively filter a \emph{vector} of innovations through a \emph{vector} of
//' \emph{ARIMA} coefficients.
//' 
//' @param \code{in_nov} A \emph{vector} of innovations (random numbers).
//' @param \code{co_eff} A \emph{vector} of \emph{ARIMA} coefficients.
//'
//' @return A column \emph{vector} of the same length as the argument
//'   \code{in_nov}.
//'
//' @details The function \code{sim_arima()} recursively filters a \emph{vector}
//'   of innovations through a \emph{vector} of \emph{ARIMA} coefficients, using
//'   \code{RcppArmadillo} \code{C++} code.
//'   It performs the same calculation as the standard \code{R} function
//'   \code{filter(x=in_nov, filter=co_eff, method="recursive")}, but it's over
//'   \code{6} times faster.
//'   
//' @examples
//' \dontrun{
//' # Create vector of innovations
//' in_nov <- rnorm(100)
//' # Create ARIMA coefficients
//' co_eff <- c(-0.8, 0.2)
//' # Calculate recursive filter using filter()
//' filter_ed <- filter(in_nov, filter=co_eff, method="recursive")
//' # Calculate recursive filter using RcppArmadillo
//' ari_ma <- HighFreq::sim_arima(in_nov, rev(co_eff))
//' # Compare the two methods
//' all.equal(as.numeric(ari_ma), as.numeric(filter_ed))
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec sim_arima(const arma::vec& in_nov, const arma::vec& co_eff) {
  arma::uword len_gth = in_nov.n_elem;
  arma::uword look_back = co_eff.n_elem;
  arma::vec ari_ma(len_gth);
  
  // Warmup period
  ari_ma(0) = in_nov(0);
  ari_ma(1) = in_nov(1) + co_eff(look_back-1) * ari_ma(0);
  for (arma::uword it=2; it < look_back-1; it++) {
    ari_ma(it) = in_nov(it) + arma::dot(co_eff.subvec(look_back-it, look_back-1), ari_ma.subvec(0, it-1));
  }  // end for
  
  // Remaining periods
  for (arma::uword it = look_back; it < len_gth; it++) {
    ari_ma(it) = in_nov(it) + arma::dot(co_eff, ari_ma.subvec(it-look_back, it-1));
  }  // end for
  
  return ari_ma;
}  // end sim_arima



////////////////////////////
// Functions for backtests
////////////////////////////


////////////////////////////////////////////////////////////
//' Calculate the optimal portfolio weights for different objective functions.
//' 
//' @param \code{re_turns} A \emph{time series} or a \emph{matrix} of returns
//'   data (the returns in excess of the risk-free rate).
//'   
//' @param \code{typ_e} A \emph{string} specifying the objective for calculating
//'   the weights (see Details).
//'   
//' @param \code{to_l} A \emph{numeric} tolerance level for discarding small
//'   eigenvalues in order to regularize the matrix inverse.  (The default is
//'   \code{0.001})
//'   
//' @param \code{max_eigen} An \emph{integer} equal to the number of
//'   eigenvectors used for calculating the regularized inverse of the
//'   covariance \emph{matrix} (the default is the number of columns of
//'   \code{re_turns}).
//'   
//' @param \code{al_pha} The shrinkage intensity between \code{0} and \code{1}.
//'   (the default is \code{0}).
//' 
//' @param \code{scal_e} A \emph{Boolean} specifying whether the weights should
//'   be scaled (the default is \code{scal_e = TRUE}).
//'
//' @param \code{vo_l} A \emph{numeric} volatility target for scaling the
//'   weights.  (The default is \code{0.001})
//'   
//' @return A column \emph{vector} of the same length as the number of columns
//'   of \code{re_turns}.
//'
//' @details The function \code{calc_weights()} calculates the optimal portfolio
//'   weights for different objective functions, using \code{RcppArmadillo}
//'   \code{C++} code.
//' 
//'   If \code{typ_e = "max_sharpe"} (the default) then \code{calc_weights()}
//'   calculates the weights of the maximum Sharpe portfolio, by multiplying the
//'   inverse of the covariance \emph{matrix} times the mean column returns.
//'   
//'   If \code{typ_e = "min_var"} then it calculates the weights of the minimum
//'   variance portfolio under linear constraints.
//'   
//'   If \code{typ_e = "min_varpca"} then it calculates the weights of the
//'   minimum variance portfolio under quadratic constraints (which is the
//'   highest order principal component).
//' 
//'   If \code{typ_e = "rank"} then it calculates the weights as the ranks
//'   (order index) of the trailing Sharpe ratios of the portfolio assets.
//' 
//'   If \code{scal_e = TRUE} (the default) then the weights are scaled so that
//'   the resulting portfolio has a volatility equal to \code{vo_l}.
//'   
//'   \code{calc_weights()} applies dimensional regularization to calculate the
//'   inverse of the covariance \emph{matrix} of returns from its eigen
//'   decomposition, using the function \code{arma::eig_sym()}.
//'   
//'   In addition, it applies shrinkage to the \emph{vector} of mean column
//'   returns, by shrinking it to its common mean value.
//'   The shrinkage intensity \code{al_pha} determines the amount of shrinkage 
//'   that is applied, with \code{al_pha = 0} representing no shrinkage (with 
//'   the estimator of mean returns equal to the means of the columns of 
//'   \code{re_turns}), and \code{al_pha = 1} representing complete shrinkage 
//'   (with the estimator of mean returns equal to the single mean of all the
//'   columns of \code{re_turns})
//' 
//' @examples
//' \dontrun{
//' # Calculate covariance matrix of ETF returns
//' re_turns <- na.omit(rutils::etf_env$re_turns[, 1:16])
//' ei_gen <- eigen(cov(re_turns))
//' # Calculate regularized inverse of covariance matrix
//' max_eigen <- 3
//' eigen_vec <- ei_gen$vectors[, 1:max_eigen]
//' eigen_val <- ei_gen$values[1:max_eigen]
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
//' weight_s <- drop(HighFreq::calc_weights(re_turns, max_eigen=max_eigen, al_pha=al_pha))
//' all.equal(weight_s, weights_r)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec calc_weights(const arma::mat& re_turns, 
                       const std::string& typ_e = "max_sharpe",
                       double to_l = 0.001,
                       int max_eigen = 0,
                       const double& pro_b = 0.1,
                       const double& al_pha = 0.0,
                       const bool scal_e = true,
                       double vo_l = 0.01) {
  // Initialize
  arma::vec weight_s(re_turns.n_cols);
  if (max_eigen == 0)  max_eigen = re_turns.n_cols;
  
  // Calculate weights depending on typ_e
  if (typ_e == "max_sharpe") {
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(re_turns, 0));
    // Shrink mean_cols to the mean of re_turns
    mean_cols = ((1-al_pha)*mean_cols + al_pha*arma::mean(mean_cols));
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(re_turns, max_eigen);
    // weight_s = calc_inv(re_turns, max_eigen=max_eigen)*mean_cols;
    weight_s = calc_inv(re_turns, to_l=to_l, max_eigen=max_eigen)*mean_cols;
  } else if (typ_e == "max_sharpe_median") {
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::median(re_turns, 0));
    // Shrink mean_cols to the mean of re_turns
    mean_cols = ((1-al_pha)*mean_cols + al_pha*arma::median(mean_cols));
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(re_turns, max_eigen);
    weight_s = calc_inv(re_turns, to_l=to_l, max_eigen=max_eigen)*mean_cols;
  } else if (typ_e == "min_var") {
    // Apply regularized inverse to unit vector
    weight_s = calc_inv(re_turns, to_l=to_l, max_eigen=max_eigen)*arma::ones(re_turns.n_cols);
  } else if (typ_e == "min_varpca") {
    // Calculate highest order principal component
    arma::vec eigen_val;
    arma::mat eigen_vec;
    arma::eig_sym(eigen_val, eigen_vec, cov(re_turns));
    weight_s = eigen_vec.col(0);
  } else if (typ_e == "rank") {
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(re_turns, 0));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(re_turns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(mean_cols)));
    weight_s = (weight_s - arma::mean(weight_s));
  } else if (typ_e == "rankrob") {
    // Median returns by columns
    arma::vec mean_cols = arma::trans(arma::median(re_turns, 0));
    // mean_cols = ((1-al_pha)*mean_cols + al_pha*arma::mean(mean_cols));
    // Standard deviation by columns
    arma::vec sd_cols = arma::trans(arma::stddev(re_turns, 0));
    sd_cols.replace(0, 1);
    mean_cols = mean_cols/sd_cols;
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(re_turns, max_eigen);
    // weight_s = calc_inv(re_turns, max_eigen)*mean_cols;
    // weight_s = calc_inv(re_turns, max_eigen)*mean_cols;
    // // Standard deviation by columns
    // arma::vec sd_cols = mean_cols;
    // for (arma::uword it=0; it < re_turns.n_cols; it++) {
    //   sd_cols(it) = arma::median(arma::abs((re_turns.col(it) - sd_cols)));
    // }  // end for
    // sd_cols.replace(0, 1);
    // mean_cols = mean_cols/sd_cols;
    // Weights equal to ranks of Sharpe
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(mean_cols)));
    // pro_b;
    weight_s = (weight_s - arma::mean(weight_s));
  } else if (typ_e == "quan_tile") {
    // Sum of quantiles for columns
    arma::vec prob_s = {pro_b, 1-pro_b};
    weight_s = conv_to< vec >::from(arma::sum(arma::quantile(re_turns, prob_s, 0), 0));
    // Weights equal to ranks
    weight_s = conv_to< vec >::from(arma::sort_index(arma::sort_index(weight_s)));
    weight_s = (weight_s - arma::mean(weight_s));
  } else {
    cout << "Warning: Incorrect typ_e argument: " << typ_e << endl;
    return arma::ones(re_turns.n_cols);
  }  // end if

  if (scal_e == TRUE) {
    // return weight_s/sqrt(sum(square(weight_s)));
    // return weight_s/sum(weight_s);
    // Returns of equally weighted portfolio
    // arma::vec mean_rows = arma::mean(re_turns, 1);
    // Returns of weighted portfolio
    // arma::vec returns_portf = re_turns*weight_s;
    // Scale weight_s to equally weighted portfolio and return them
    // return weight_s*arma::stddev(arma::mean(re_turns, 1))/arma::stddev(re_turns*weight_s);
    // Scale weight_s so the resulting portfolio has a volatility equal to vo_l
    return weight_s*vo_l/arma::stddev(re_turns*weight_s);
  }  // end if
  
  return weight_s;
}  // end calc_weights



////////////////////////////////////////////////////////////
//' Simulate (backtest) a rolling portfolio optimization strategy, using
//' \code{RcppArmadillo}.
//' 
//' @param \code{ex_cess} A \emph{time series} or a \emph{matrix} of excess
//'   returns data (the returns in excess of the risk-free rate).
//'   
//' @param \code{re_turns} A \emph{time series} or a \emph{matrix} of returns
//'   data (the returns in excess of the risk-free rate).
//'   
//' @param \code{start_points} An \emph{integer vector} of start points.
//' 
//' @param \code{end_points} An \emph{integer vector} of end points.
//' 
//' @param \code{typ_e} A \emph{string} specifying the objective for calculating
//'   the weights (see Details).
//'   
//' @param \code{to_l} A \emph{numeric} tolerance level for discarding small
//'   eigenvalues in order to regularize the matrix inverse.  (The default is
//'   \code{0.001})
//'   
//' @param \code{max_eigen} An \emph{integer} equal to the number of
//'   eigenvectors used for calculating the regularized inverse of the
//'   covariance \emph{matrix} (the default is the number of columns of
//'   \code{re_turns}).
//'   
//' @param \code{al_pha} The shrinkage intensity between \code{0} and \code{1}.
//'   (the default is \code{0}).
//'   
//' @param \code{scal_e} A \emph{Boolean} specifying whether the weights should
//'   be scaled (the default is \code{scal_e = TRUE}).
//'   
//' @param \code{vo_l} A \emph{numeric} volatility target for scaling the
//'   weights.  (The default is \code{0.001})
//'   
//' @param \code{co_eff} A \emph{numeric} multiplier of the weights.  (The
//'   default is \code{1})
//'   
//' @param \code{bid_offer} A \emph{numeric} bid-offer spread.  (The default is
//'   \code{0})
//'
//' @return A column \emph{vector} of strategy returns, with the same length as
//'   the number of rows of \code{re_turns}.
//'
//' @details The function \code{back_test()} performs a backtest simulation of a
//'   rolling portfolio optimization strategy over a \emph{vector} of
//'   \code{end_points}.
//'   
//'   It performs a loop over the \code{end_points}, and subsets the
//'   \emph{matrix} of excess returns \code{ex_cess} along its rows, between the
//'   corresponding end point and the start point. It passes the subset matrix
//'   of excess returns into the function \code{calc_weights()}, which
//'   calculates the optimal portfolio weights. The arguments \code{max_eigen},
//'   \code{al_pha}, \code{typ_e}, and \code{scal_e} are also passed to the
//'   function \code{calc_weights()}.
//'   
//'   The function \code{back_test()} multiplies the weights by the coefficient
//'   \code{co_eff} (with default equal to \code{1}), which allows reverting a
//'   strategy if \code{co_eff = -1}.
//'   
//'   The function \code{back_test()} then multiplies the weights times the
//'   future portfolio returns, to calculate the out-of-sample strategy returns.
//'   
//'   The function \code{back_test()} calculates the transaction costs by
//'   multiplying the bid-offer spread \code{bid_offer} times the absolute
//'   difference between the current weights minus the weights from the previous
//'   period. It then subtracts the transaction costs from the out-of-sample
//'   strategy returns.
//'   
//'   The function \code{back_test()} returns a \emph{time series} (column
//'   \emph{vector}) of strategy returns, of the same length as the number of
//'   rows of \code{re_turns}.
//'
//' @examples
//' \dontrun{
//' # Calculate the ETF daily excess returns
//' re_turns <- na.omit(rutils::etf_env$re_turns[, 1:16])
//' # risk_free is the daily risk-free rate
//' risk_free <- 0.03/260
//' ex_cess <- re_turns - risk_free
//' # Define monthly end_points without initial warmpup period
//' end_points <- rutils::calc_endpoints(re_turns, inter_val="months")
//' end_points <- end_points[end_points>50]
//' len_gth <- NROW(end_points)
//' # Define 12-month look_back interval and start_points over sliding window
//' look_back <- 12
//' start_points <- c(rep_len(1, look_back-1), end_points[1:(len_gth-look_back+1)])
//' # Define shrinkage and regularization intensities
//' al_pha <- 0.5
//' max_eigen <- 3
//' # Simulate a monthly rolling portfolio optimization strategy
//' pnl_s <- HighFreq::back_test(ex_cess, re_turns, 
//'                             start_points-1, end_points-1, 
//'                             max_eigen = max_eigen, 
//'                             al_pha = al_pha)
//' pnl_s <- xts::xts(pnl_s, index(re_turns))
//' colnames(pnl_s) <- "strat_rets"
//' # Plot dygraph of strategy
//' dygraphs::dygraph(cumsum(pnl_s), 
//'   main="Cumulative Returns of Max Sharpe Portfolio Strategy")
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat back_test(const arma::mat& ex_cess, // Portfolio excess returns
                    const arma::mat& re_turns, // Portfolio returns
                    const arma::uvec& start_points, 
                    const arma::uvec& end_points, 
                    const std::string& typ_e = "max_sharpe",
                    double to_l = 0.001,
                    int max_eigen = 0,
                    const double& pro_b = 0.1,
                    const double& al_pha = 0,
                    const bool& scal_e = true,
                    double vo_l = 0.01,
                    const double& co_eff = 1.0,
                    const double& bid_offer = 0.0) {
  
  arma::vec weight_s(re_turns.n_cols);
  arma::vec weights_past = zeros(re_turns.n_cols);
  arma::vec pnl_s = zeros(re_turns.n_rows);
  
  // Perform loop over the end_points
  for (arma::uword it = 1; it < end_points.size(); it++) {
    // cout << "it: " << it << endl;
    // Calculate portfolio weights
    weight_s = co_eff*calc_weights(ex_cess.rows(start_points(it-1), end_points(it-1)), typ_e, to_l, max_eigen, pro_b, al_pha, scal_e, vo_l);
    // Calculate out-of-sample returns
    pnl_s.subvec(end_points(it-1)+1, end_points(it)) = re_turns.rows(end_points(it-1)+1, end_points(it))*weight_s;
    // Add transaction costs
    pnl_s.row(end_points(it-1)+1) -= bid_offer*sum(abs(weight_s - weights_past))/2;
    weights_past = weight_s;
  }  // end for
  
  // Return the strategy pnl_s
  return pnl_s;
  
}  // end back_test
