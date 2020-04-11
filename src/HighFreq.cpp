// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
using namespace std;
using namespace Rcpp;
using namespace arma;

////////////////////////////
// Rcpp and RcppArmadillo functions for package HighFreq
////////////////////////////


////////////////////////////
// Functions for matrix algebra
////////////////////////////


////////////////////////////////////////////////////////////
//' Apply a lag to a \emph{vector} or a single-column \emph{time series}
//' using \code{RcppArmadillo}.
//' 
//' @param \code{vec_tor} A \emph{vector} or a single-column \emph{time series}.
//' @param \code{lagg} An \emph{integer} equal to the number of periods to lag
//'   (the default is \code{lagg=1}).
//'
//' @return A column \emph{vector} with the same number of elements as the input
//'   vector.
//'
//' @details The function \code{lag_vec()} applies a lag to the input
//'   \emph{vector} by shifting its elements by the number equal to the argument
//'   \code{lagg}. For positive \code{lagg} values, the elements are shifted
//'   forward, and for negative \code{lagg} values they are shifted backward.
//'   The output \emph{vector} is padded with either the first or the last
//'   element, to maintain its original length.
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
//'   rcpp=HighFreq::lag_vec(re_turns),
//'   rcode=rutils::lag_it(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec lag_vec(arma::vec& vec_tor, int lagg=1) {
  
  int len_gth = (vec_tor.n_elem-1);
  
  if (lagg > 0)
    return arma::join_cols(arma::repelem(vec_tor.subvec(0, 0), lagg, 1), 
                           vec_tor.subvec(0, len_gth-lagg));
  else
    return arma::join_cols(vec_tor.subvec(-lagg, len_gth), 
                           arma::repelem(vec_tor.subvec(len_gth, len_gth), -lagg, 1));
  
}  // end lag_vec




////////////////////////////////////////////////////////////
//' Apply a lag to a \emph{matrix} or \emph{time series} using
//' \code{RcppArmadillo}.
//' 
//' @param \code{mat_rix} A \emph{matrix} or \emph{time series}.
//' @param \code{lagg} An \emph{integer} equal to the number of periods to lag
//'   (the default is \code{lagg=1}).
//'
//' @return A \emph{matrix} with the same dimensions as the input
//'   argument \code{mat_rix}.
//'
//' @details The function \code{lag_it()} applies a lag to the input
//'   \emph{matrix} by shifting its rows by the number equal to the argument
//'   \code{lagg}. For positive \code{lagg} values, the rows are shifted forward
//'   (down), and for negative \code{lagg} values they are shifted backward
//'   (up). The output \emph{matrix} is padded with either the first or the last
//'   row, to maintain it original dimensions. The function \code{lag_it()} can
//'   be applied to vectors in the form of single-column matrices.
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
//'   rcpp=HighFreq::lag_it(re_turns),
//'   rcode=rutils::lag_it(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat lag_it(arma::mat& mat_rix, int lagg=1) {
  int num_rows = (mat_rix.n_rows-1);
  
  if (lagg > 0)
    return arma::join_cols(arma::repelem(mat_rix.row(0), lagg, 1), 
                           mat_rix.rows(0, num_rows-lagg));
  else
    return arma::join_cols(mat_rix.rows(-lagg, num_rows), 
                           arma::repelem(mat_rix.row(num_rows), -lagg, 1));
  
}  // end lag_it




////////////////////////////////////////////////////////////
//' Calculate the differences of a \emph{vector} or a single-column
//' \emph{time series} using \emph{RcppArmadillo}.
//' 
//' @param \code{vec_tor} A \emph{vector} or single-column \emph{time series}.
//' @param \code{lagg} An \emph{integer} equal to the number of time periods to
//'   lag when calculating the differences (the default is \code{lagg=1}).
//' @param \code{padd} \emph{Boolean} argument: Should the output \emph{vector}
//'   be padded (extended) with zeros, in order to return a \emph{vector} of the
//'   same length as the input? (the default is \code{padd=FALSE})
//'
//' @return A column \emph{vector} containing the differences of the input
//'   vector.
//'
//' @details The function \code{diff_vec()} calculates the differences between
//'   the input \emph{vector} or \emph{time series} and its lagged version. 
//'   
//'   The argument \code{lagg} specifies the number of lags.  For example, if
//'   \code{lagg=3} then the differences will be taken between each element
//'   minus the element three time periods before it (in the past).  The default
//'   is \code{lagg=1}.
//' 
//'   The argument \code{padd} specifies whether the output \emph{vector} should
//'   be padded (extended) with zeros at the beginning, in order to return a
//'   \emph{vector} of the same length as the input.  The default is
//'   \code{padd=FALSE}. The padding operation is time-consuming, so that
//'   \code{padd=FALSE} can be twice as fast as \code{padd=TRUE}.
//'   
//'   The function \code{diff_vec()} is implemented in \code{RcppArmadillo}
//'   code, which makes it slightly faster than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random returns
//' re_turns <- rnorm(1e6)
//' # Compare diff_vec() with rutils::diff_it()
//' all.equal(drop(HighFreq::diff_vec(re_turns, padd=TRUE)),
//'   rutils::diff_it(re_turns))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   rcpp=HighFreq::diff_vec(re_turns, padd=TRUE),
//'   rcode=rutils::diff_it(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec diff_vec(arma::vec& vec_tor, int lagg=1, const bool& padd=false) {
  
  if (padd)
    return arma::join_cols(arma::zeros(lagg), 
                           arma::diff(vec_tor, lagg));
  else
    return arma::diff(vec_tor, lagg);
  
}  // end diff_vec




////////////////////////////////////////////////////////////
//' Calculate the row differences of a \emph{matrix} or a \emph{time
//' series} using \emph{RcppArmadillo}.
//' 
//' @param \code{mat_rix} A \emph{matrix} or \emph{time series}.
//' @param \code{lagg} An \emph{integer} equal to the number of rows (time
//'   periods) to lag when calculating the differences (the default is
//'   \code{lagg=1}).
//' @param \code{padd} \emph{Boolean} argument: Should the output \emph{matrix}
//'   be padded (extended) with zeros, in order to return a \emph{matrix} with
//'   the same number of rows as the input? (the default is \code{padd=FALSE})
//'
//' @return A \emph{matrix} containing the differences of the input
//'   \emph{matrix}.
//'
//' @details The function \code{diff_it()} calculates the differences between
//'   the rows of the input \emph{matrix} or \emph{time series} and its lagged
//'   version. The lagged version has its rows shifted down by the number equal
//'   to \code{lagg} rows.
//'   
//'   The argument \code{lagg} specifies the number of lags applied to the rows
//'   of the lagged version. For example, if \code{lagg=3} then the lagged
//'   version will have its rows shifted down by \code{3} rows, and the
//'   differences will be taken between each row minus the row three time
//'   periods before it (in the past). The default is \code{lagg=1}.
//' 
//'   The argument \code{padd} specifies whether the output \emph{matrix} should
//'   be padded (extended) with rows of zeros at the beginning, in order to
//'   return a \emph{matrix} with the same number of rows as the input.  The
//'   default is \code{padd=FALSE}. The padding operation is time-consuming, so
//'   that \code{padd=FALSE} can be twice as fast as \code{padd=TRUE}.
//'   
//'   The function \code{diff_it()} is implemented in \code{RcppArmadillo}
//'   code, which makes it slightly faster than \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Create a matrix of random returns
//' re_turns <- matrix(rnorm(5e6), nc=5)
//' # Compare diff_it() with rutils::diff_it()
//' all.equal(HighFreq::diff_it(re_turns, padd=TRUE),
//'   rutils::diff_it(re_turns))
//' # Compare the speed of RcppArmadillo with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   rcpp=HighFreq::diff_it(re_turns, padd=TRUE),
//'   rcode=rutils::diff_it(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat diff_it(arma::mat& mat_rix, int lagg=1, const bool& padd=false) {
  
  if (padd)
    return arma::join_cols(arma::zeros(lagg, mat_rix.n_cols), 
                           arma::diff(mat_rix, lagg));
  else
    return arma::diff(mat_rix, lagg);
  
}  // end diff_it




////////////////////////////////////////////////////////////
//' Multiply the columns or rows of a \emph{matrix} times a \emph{vector},
//' element-wise.
//' 
//' @param \code{vec_tor} A \emph{vector}.
//' @param \code{mat_rix} A \emph{matrix}.
//' @param \code{by_col} A \emph{Boolean} argument: if \code{TRUE} then multiply
//'   the columns, otherwise multiply the rows. (The default is
//'   \code{by_col=TRUE}.)
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
//'   The function \code{mult_vec_mat()} uses \code{RcppArmadillo}, so when
//'   multiplying large \emph{matrix} columns it's several times faster than
//'   vectorized \code{R} code, and it's even much faster compared to \code{R}
//'   when multiplying the \emph{matrix} rows.
//'   
//' @examples
//' \dontrun{
//' # Multiply matrix columns using R
//' mat_rix <- matrix(round(runif(25e4), 2), nc=5e2)
//' vec_tor <- round(runif(5e2), 2)
//' prod_uct <- vec_tor*mat_rix
//' # Multiply the matrix in place
//' mult_vec_mat(vec_tor, mat_rix)
//' all.equal(prod_uct, mat_rix)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     rcpp=mult_vec_mat(vec_tor, mat_rix),
//'     rcode=vec_tor*mat_rix,
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' 
//' # Multiply matrix rows using R
//' mat_rix <- matrix(round(runif(25e4), 2), nc=5e2)
//' vec_tor <- round(runif(5e2), 2)
//' prod_uct <- t(vec_tor*t(mat_rix))
//' # Multiply the matrix in place
//' mult_vec_mat(vec_tor, mat_rix, by_col=FALSE)
//' all.equal(prod_uct, mat_rix)
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'     rcpp=mult_vec_mat(vec_tor, mat_rix, by_col=FALSE),
//'     rcode=t(vec_tor*t(mat_rix)),
//'     times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::uword mult_vec_mat(const arma::vec& vec_tor,
                   arma::mat& mat_rix,
                   const bool& by_col=true) {
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
//' @param \code{mat_rix} A numeric \emph{matrix} or \emph{time series} of
//'   returns data.
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
//'   rcpp=HighFreq::calc_eigen(re_turns),
//'   rcode=prcomp(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
List calc_eigen(const arma::mat& mat_rix) {
  arma::mat eigen_vec;
  arma::vec eigen_val;
  arma::eig_sym(eigen_val, eigen_vec, cov(mat_rix));
  // Reverse the order of elements from largest eigenvalue to smallest, similar to R
  return List::create(Named("values") = arma::flipud(eigen_val),
                      Named("vectors") = arma::fliplr(eigen_vec));
}  // end calc_eigen



////////////////////////////////////////////////////////////
//' Calculate the regularized inverse of the covariance \emph{matrix} of returns
//' using \code{RcppArmadillo}.
//' 
//' @param \code{mat_rix} A \emph{matrix} of returns data.
//' @param \code{max_eigen} An \emph{integer} equal to the regularization
//'   intensity (the number of eigenvalues and eigenvectors used for calculating
//'   the regularized inverse).
//'
//' @return A \emph{matrix} equal to the regularized inverse. 
//'
//' @details The function calc_inv() calculates the regularized inverse of the
//'   \emph{covariance matrix}, by truncating the number of eigenvectors to
//'   \code{max_eigen}. The function \code{calc_inv()} first calculates the
//'   covariance \emph{matrix} of the \code{mat_rix}, and then it calculates the
//'   regularized inverse from the truncated eigen decomposition. It uses only
//'   the largest \code{max_eigen} eigenvalues and their corresponding
//'   eigenvectors.
//'
//' @examples
//' \dontrun{
//' # Create random matrix
//' mat_rix <- matrix(rnorm(500), nc=5)
//' max_eigen <- 3
//' # Calculate regularized inverse using RcppArmadillo
//' in_verse <- HighFreq::calc_inv(mat_rix, max_eigen)
//' # Calculate regularized inverse from eigen decomposition in R
//' ei_gen <- eigen(cov(mat_rix))
//' inverse_r <-  ei_gen$vectors[, 1:max_eigen] %*% (t(ei_gen$vectors[, 1:max_eigen]) / ei_gen$values[1:max_eigen])
//' # Compare RcppArmadillo with R
//' all.equal(in_verse, inverse_r)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::mat calc_inv(const arma::mat& mat_rix, const arma::uword& max_eigen) {
  arma::mat eigen_vec;
  arma::vec eigen_val;
  
  arma::eig_sym(eigen_val, eigen_vec, cov(mat_rix));
  eigen_vec = eigen_vec.cols(eigen_vec.n_cols-max_eigen, eigen_vec.n_cols-1);
  eigen_val = 1/eigen_val.subvec(eigen_val.n_elem-max_eigen, eigen_val.n_elem-1);
  // arma::mat eigen_valmat = diagmat(eigen_val);
  
  return eigen_vec*diagmat(eigen_val)*eigen_vec.t();
  
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
//'   The function \code{calc_scaled()} uses \code{RcppArmadillo} and is about
//'   \emph{5} times faster than function \code{scale()}, for a \emph{matrix}
//'   with \emph{1,000} rows and \emph{20} columns.
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
//'   rcpp=calc_scaled(mat_rix=mat_rix, use_median=FALSE),
//'   rcode=scale(mat_rix),
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




////////////////////////////
// Functions for statistics
////////////////////////////


////////////////////////////////////////////////////////////
//' Calculate the variance of a \emph{vector} or a single-column \emph{time
//' series} using \code{RcppArmadillo}.
//' 
//' @param \code{vec_tor} A \emph{vector} or a single-column \emph{time series}.
//'
//' @return A \emph{numeric} value equal to the variance of the \emph{vector}.
//'
//' @details The function \code{calc_var_vec()} calculates the variance of a
//'   \emph{vector} using \code{RcppArmadillo}, so it's significantly faster
//'   than the \code{R} function \code{var()}.
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
//'   rcpp=HighFreq::calc_var_vec(re_turns),
//'   rcode=var(re_turns),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
double calc_var_vec(arma::vec& vec_tor) {
  return arma::var(vec_tor);
}  // end calc_var_vec




////////////////////////////////////////////////////////////
//' Calculate the variance of the columns of a \emph{matrix} or \emph{time
//' series} using \code{RcppArmadillo}.
//' 
//' @param \code{mat_rix} A \emph{matrix} or a \emph{time series}.
//'
//' @return A row vector equal to the variance of the \emph{matrix} columns.
//'
//' @details The function \code{calc_var()} calculates the variance of the
//'   columns of a \emph{matrix} using \code{RcppArmadillo}. 
//'   
//'   The function \code{calc_var()} performs the same calculation as the
//'   function \code{colVars()} from package
//'   \href{https://cran.r-project.org/web/packages/matrixStats/index.html}{matrixStats},
//'   but it's much faster because it uses \code{RcppArmadillo}.
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
//'   rcpp=HighFreq::calc_var(re_turns),
//'   matrixStats=matrixStats::colVars(re_turns),
//'   rcode=apply(re_turns, 2, var),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::rowvec calc_var(arma::mat& mat_rix) {
  return arma::var(mat_rix);
}  // end calc_var




////////////////////////////////////////////////////////////
//' Calculate the variance of an \emph{OHLC time series}, using different range
//' estimators and \code{RcppArmadillo}.
//'
//' @param \code{oh_lc} An \emph{OHLC time series} or a \emph{numeric matrix} of
//'   prices.
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
//' @param \code{lag_close} A \emph{vector} with the lagged \emph{close} prices
//'   of the \emph{OHLC time series}.  This is an optional argument. (The
//'   default is \code{lag_close=0}.)
//' @param \code{in_dex} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument. (The default is \code{in_dex=0}.)
//' @param \code{scal_e} \emph{Boolean} argument: Should the returns be divided
//'   by the number of seconds in each period? (The default is
//'   \code{scal_e=TRUE}.)
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
//'   The optional argument \code{in_dex} is the time index of the \emph{time
//'   series} \code{oh_lc}. If the time index is in seconds, then the
//'   differences of the index are equal to the number of seconds in each time
//'   period.  If the time index is in days, then the differences are equal to
//'   the number of days in each time period.
//'   
//'   If \code{scal_e} is \code{TRUE} (the default), then the returns are
//'   divided by the differences of the time index (which scales the variance to
//'   the units of variance per second squared.) This is useful when calculating
//'   the variance from minutely bar data, because dividing returns by the
//'   number of seconds decreases the effect of overnight price jumps. If the
//'   time index is in days, then the variance is equal to the variance per day
//'   squared.
//'   
//'   The optional argument \code{lag_close} are the lagged \emph{close} prices
//'   of the \emph{OHLC time series}.  Passing in the lagged \emph{close} prices
//'   speeds up the calculation, so it's useful for rolling calculations.
//'   
//'   The function \code{calc_var_ohlc()} is implemented in \code{RcppArmadillo}
//'   code, and it's over \code{10} times faster than \code{calc_var_ohlc_r()},
//'   which is implemented in \code{R} code.
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
//'   rcpp=HighFreq::calc_var_ohlc(HighFreq::SPY),
//'   rcode=HighFreq::calc_var_ohlc_r(HighFreq::SPY),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
double calc_var_ohlc(arma::mat& oh_lc, 
                     const std::string& calc_method="yang_zhang", 
                     arma::colvec lag_close=0, 
                     arma::colvec in_dex=0, 
                     const bool& scal_e=true) {
  
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
//' Calculate the ranks of the elements of a \emph{vector} or a single-column
//' \emph{time series} using \code{RcppArmadillo}.
//' 
//' @param \code{vec_tor} A \emph{vector} or a single-column \emph{time series}.
//'
//' @return An \emph{integer vector} with the ranks of the elements of the
//'   \emph{vector}.
//'
//' @details The function \code{calc_ranks()} calculates the ranks of the
//'   elements of a \emph{vector} or a single-column \emph{time series}.
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
//'   rcpp=calc_ranks(da_ta),
//'   rcode=rank(da_ta),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::uvec calc_ranks(const arma::vec& vec_tor) {
  return (arma::sort_index(arma::sort_index(vec_tor)) + 1);
}  // end calc_ranks




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
//'   \code{RcppArmadillo} and is about \emph{10} times faster than \code{lm()}.
//'   The code was inspired by this article (but it's not identical to it):
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
//' Calculate the rolling sum over a \emph{vector} or a single-column \emph{time
//' series} using \emph{Rcpp}.
//' 
//' @param \code{vec_tor} A \emph{vector} or a single-column \emph{time series}.
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of elements of data used for calculating the sum.
//'
//' @return A column \emph{vector} of the same length as the argument
//'   \code{vec_tor}.
//'
//' @details The function \code{roll_sum()} calculates a \emph{vector} of 
//'   rolling sums, over a \emph{vector} of data, using \emph{Rcpp}.  The
//'   function \code{roll_sum()} is several times faster than
//'   \code{rutils::roll_sum()} which uses vectorized \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Create a vector of random returns
//' re_turns <- rnorm(1e6)
//' # Calculate rolling sums over 11-period lookback intervals
//' sum_rolling <- HighFreq::roll_sum(re_turns, look_back=11)
//' # Compare HighFreq::roll_sum() with rutils::roll_sum()
//' all.equal(HighFreq::roll_sum(re_turns, look_back=11), 
//'          rutils::roll_sum(re_turns, look_back=11))
//' # Compare the speed of Rcpp with R code
//' library(microbenchmark)
//' summary(microbenchmark(
//'   rcpp=HighFreq::roll_sum(re_turns, look_back=11),
//'   rcode=rutils::roll_sum(re_turns, look_back=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' 
//' @export
// [[Rcpp::export]]
NumericVector roll_sum(NumericVector vec_tor, int look_back) {
  int len_gth = vec_tor.size();
  NumericVector rolling_sum(len_gth);

  // Warmup period
  rolling_sum[0] = vec_tor[0];
  for (int it=1; it < look_back; it++) {
    rolling_sum[it] = rolling_sum[it-1] + vec_tor[it];
  }  // end for
  
  // Remaining period
  for (int it=look_back; it < len_gth; it++) {
    rolling_sum[it] = rolling_sum[it-1] + vec_tor[it] - vec_tor[it-look_back];
  }  // end for
  
  return rolling_sum;
}  // end roll_sum




////////////////////////////////////////////////////////////
//' Calculate the rolling weighted sum over a \emph{vector} or a single-column
//' \emph{time series} using \code{RcppArmadillo}.
//' 
//' @param \code{vec_tor} A \emph{vector} or a single-column \emph{time series}.
//' @param \code{weight_s} A \emph{vector} of weights.
//'
//' @return A column \emph{vector} of the same length as the argument
//'   \code{vec_tor}.
//'
//' @details The function \code{roll_wsum()} calculates the rolling weighted sum
//'   of a \emph{vector} over its past values (a convolution with the
//'   \emph{vector} of weights), using \code{RcppArmadillo}. It performs a
//'   similar calculation as the standard \code{R} function
//'   \code{stats::filter(x=vec_tor, filter=weight_s, method="convolution",
//'   sides=1)}, but it's over \code{6} times faster, and it doesn't produce any
//'   \code{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # Create vector from historical prices
//' vec_tor <- as.numeric(rutils::etf_env$VTI[, 6])
//' # Create simple weights
//' weight_s <- c(1, rep(0, 10))
//' # Calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_wsum(vec_tor=vec_tor, weight_s=weight_s)
//' # Compare with original
//' all.equal(vec_tor, as.numeric(weight_ed))
//' # Second example
//' # Create exponentially decaying weights
//' weight_s <- exp(-0.2*1:11)
//' weight_s <- weight_s/sum(weight_s)
//' # Calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_wsum(vec_tor=vec_tor, weight_s=weight_s)
//' # Calculate rolling weighted sum using filter()
//' filter_ed <- stats::filter(x=vec_tor, filter=weight_s, method="convolution", sides=1)
//' # Compare both methods
//' all.equal(filter_ed[-(1:11)], weight_ed[-(1:11)], check.attributes=FALSE)
//' }
//' 
//' @export
// [[Rcpp::export]]
arma::vec roll_wsum(const arma::vec& vec_tor, const arma::vec& weight_s) {
  arma::uword len_gth = vec_tor.n_elem;
  arma::uword look_back = weight_s.n_elem;
  arma::vec rolling_sum(len_gth);
  arma::vec rev_weights = arma::reverse(weight_s);
  // arma::vec rev_weights = weight_s;
  
  // Warmup period
  rolling_sum.subvec(0, look_back-2) = vec_tor.subvec(0, look_back-2);
  
  // Remaining periods
  for (arma::uword it=look_back-1; it < len_gth; it++) {
    rolling_sum(it) = arma::dot(rev_weights, vec_tor.subvec(it-look_back+1, it));
  }  // end for
  
  return rolling_sum;
}  // end roll_wsum





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
//' Calculate a \emph{vector} of variance estimates over a rolling look-back
//' interval for a \emph{vector} or a single-column \emph{time series}, using
//' \code{RcppArmadillo}.
//'
//' @param \code{vec_tor} A \emph{vector} or a single-column \emph{time series}.
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of \emph{vector} elements used for calculating a single variance
//'   estimate.
//'
//' @return A column \emph{vector} with the same number of elements as the input
//'   argument \code{vec_tor}.
//'
//' @details The function \code{roll_var_vec()} calculates a \code{vec_tor} of
//'   variance estimates over a rolling look-back interval for a \emph{vector}
//'   or a single-column \emph{time series}, using \code{RcppArmadillo}.
//'   
//'   The function \code{roll_var_vec()} uses an expanding look-back interval in
//'   the initial warmup period, to calculate the same number of elements as the
//'   input argument \code{vec_tor}.
//'
//'   The function \code{roll_var_vec()} performs the same calculation as the
//'   function \code{roll_var()} from package
//'   \href{https://cran.r-project.org/web/packages/RcppRoll/index.html}{RcppRoll},
//'   but it's several times faster because it uses \code{RcppArmadillo}.
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
arma::vec roll_var_vec(arma::vec& vec_tor, arma::uword look_back=11) {
  arma::uword len_gth = vec_tor.n_elem;
  arma::vec var_vec = arma::zeros(len_gth);
  
  // Warmup period
  for (arma::uword it=1; it < look_back; it++) {
    var_vec(it) = arma::var(vec_tor.subvec(0, it));
  }  // end for
  
  // Remaining period
  for (arma::uword it=look_back; it < len_gth; it++) {
    var_vec(it) = arma::var(vec_tor.subvec(it-look_back+1, it));
  }  // end for
  
  return var_vec;
  
}  // end roll_var_vec




////////////////////////////////////////////////////////////
//' Calculate a \emph{matrix} of variance estimates over a rolling look-back
//' interval for a \emph{time series} or a \emph{matrix}, using
//' \code{RcppArmadillo}.
//'
//' @param \code{mat_rix} A \emph{matrix} or a \emph{time series}.
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of time periods (\emph{matrix} rows) used for calculating a single
//'   variance estimate.
//'
//' @return A \emph{matrix} with the same number of rows and columns as the
//'   input argument \code{mat_rix}.
//'
//' @details The function \code{roll_var()} calculates a \code{mat_rix} of
//'   variance estimates over a rolling look-back interval for a \emph{time
//'   series} or a \emph{matrix}, using \code{RcppArmadillo}.
//'
//'   The function \code{roll_var()} uses an expanding look-back interval in the
//'   initial warmup period, to calculate the same number of rows as the input
//'   argument \code{mat_rix}.
//'
//'   The function \code{roll_var()} performs the same calculation as the
//'   function \code{roll_var()} from package
//'   \href{https://cran.r-project.org/web/packages/RcppRoll/index.html}{RcppRoll},
//'   but it's several times faster because it uses \code{RcppArmadillo}.
//'
//' @examples
//' \dontrun{
//' # Create a matrix of random returns
//' re_turns <- matrix(rnorm(5e3), nc=5)
//' # Compare the variance estimates over 11-period lookback intervals
//' all.equal(HighFreq::roll_var(re_turns, look_back=11)[-(1:10), ], 
//'   RcppRoll::roll_var(re_turns, n=11))
//' # Compare the speed of RcppArmadillo with RcppRoll
//' library(microbenchmark)
//' summary(microbenchmark(
//'   RcppArmadillo=HighFreq::roll_var(re_turns, look_back=11),
//'   RcppRoll=RcppRoll::roll_var(re_turns, n=11),
//'   times=10))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_var(arma::mat& mat_rix, arma::uword look_back=11) {
  arma::uword num_rows = mat_rix.n_rows;
  arma::mat var_mat = arma::zeros(num_rows, mat_rix.n_cols);
  
  // Warmup period
  for (arma::uword it=1; it < look_back; it++) {
    var_mat.row(it) = arma::var(mat_rix.rows(0, it));
  }  // end for
  
  // Remaining period
  for (arma::uword it=look_back; it < num_rows; it++) {
    var_mat.row(it) = arma::var(mat_rix.rows(it-look_back+1, it));
  }  // end for
  
  return var_mat;
  
}  // end roll_var


////////////////////////////////////////////////////////////
//' Calculate a \emph{vector} of variance estimates over a rolling look-back
//' interval for an \emph{OHLC time series}, using different range estimators
//' and \code{RcppArmadillo}.
//'
//' @param \code{oh_lc} An \emph{OHLC time series} or a \emph{numeric matrix} of
//'   prices.
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
//' @param \code{in_dex} A \emph{vector} with the time index of the \emph{time
//'   series}.  This is an optional argument. (The default is \code{in_dex=0}.)
//' @param \code{scal_e} \emph{Boolean} argument: Should the returns be divided
//'   by the number of seconds in each period? (The default is
//'   \code{scal_e=TRUE}.)
//' @param \code{look_back} The length of the look-back interval, equal to the
//'   number of time periods (\code{oh_lc} rows) used for calculating a single
//'   variance estimate.
//'
//' @return A column \emph{vector} of the same length as the number of rows of
//'   \code{oh_lc}.
//'
//' @details The function \code{roll_var_ohlc()} performs a loop over the rows
//'   of \code{oh_lc}, subsets a number of previous (past) rows equal to
//'   \code{look_back}, and passes them into the function
//'   \code{calc_var_ohlc()}. It uses an expanding look-back interval in the
//'   initial warmup period, to calculate the same number of elements as the
//'   number of rows in the input argument \code{oh_lc}.
//' 
//'   The function \code{roll_var_ohlc()} calculates the variance from all the
//'   different intra-day and day-over-day returns (defined as the differences
//'   of \emph{OHLC} prices), using several different variance estimation
//'   methods.
//'   
//'   The default \code{calc_method} is \emph{"yang_zhang"}, which theoretically
//'   has the lowest standard error among unbiased estimators.
//'   The methods \emph{"close"}, \emph{"garman_klass_yz"}, and
//'   \emph{"yang_zhang"} do account for \emph{close-to-open} price jumps, while
//'   the methods \emph{"garman_klass"} and \emph{"rogers_satchell"} do not
//'   account for \emph{close-to-open} price jumps.
//'
//'   The optional argument \code{in_dex} is the time index of the \emph{time
//'   series} \code{oh_lc}. If the time index is in seconds, then the
//'   differences of the index are equal to the number of seconds in each time
//'   period.  If the time index is in days, then the differences are equal to
//'   the number of days in each time period.
//'   
//'   If \code{scal_e} is \code{TRUE} (the default), then the returns are
//'   divided by the differences of the time index (which scales the variance to
//'   the units of variance per second squared.) This is useful when calculating
//'   the variance from minutely bar data, because dividing returns by the
//'   number of seconds decreases the effect of overnight price jumps. If the
//'   time index is in days, then the variance is equal to the variance per day
//'   squared.
//'   
//'   The function \code{roll_var_ohlc()} is implemented in \code{RcppArmadillo}
//'   code, so it's many times faster than the equivalent \code{R} code.
//'
//' @examples
//' \dontrun{
//' # Extract time index of SPY returns
//' oh_lc <- HighFreq::SPY
//' in_dex <- c(1, diff(xts::.index(HighFreq::SPY)))
//' # Calculate the rolling variance of SPY returns, with scaling of the returns
//' var_rolling <- roll_var_ohlc(oh_lc, 
//'                               calc_method="yang_zhang", 
//'                               in_dex=in_dex,
//'                               scal_e=TRUE, 
//'                               look_back=21)
//' }
//' @export
// [[Rcpp::export]]
arma::vec roll_var_ohlc(arma::mat& oh_lc, 
                        const std::string& calc_method="yang_zhang", 
                        arma::colvec in_dex=0, 
                        const bool& scal_e=true, 
                        arma::uword look_back=11) {

  arma::uword num_rows = oh_lc.n_rows;
  arma::vec var_vec = arma::zeros(num_rows);
  arma::colvec clo_se = oh_lc.col(3);
  arma::colvec lag_close = lag_it(clo_se);

  if (!scal_e || (in_dex.n_rows == 1)) {
    in_dex = arma::ones(num_rows);
  }  // end if
  
  // Warmup period
  for (arma::uword it=1; it < look_back; it++) {
    arma::mat sub_ohlc = oh_lc.rows(0, it);
    arma::colvec sub_close = lag_close.rows(0, it);
    arma::colvec sub_index = in_dex.subvec(0, it);
    var_vec(it) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, sub_index, scal_e);
  }  // end for
  
  // Remaining period
  for (arma::uword it=look_back; it < num_rows; it++) {
    arma::mat sub_ohlc = oh_lc.rows(it-look_back+1, it);
    arma::colvec sub_close = lag_close.rows(it-look_back+1, it);
    arma::colvec sub_index = in_dex.subvec(it-look_back+1, it);
    var_vec(it) = calc_var_ohlc(sub_ohlc, calc_method, sub_close, sub_index, scal_e);
  }  // end for
  
  return var_vec;
  
}  // end roll_var_ohlc




////////////////////////////////////////////////////////////
//' Perform a rolling scaling (standardization) of the columns of a
//' \emph{matrix} of data using \code{RcppArmadillo}.
//' 
//' @param mat_rix A \emph{matrix} of data.
//' @param look_back The length of the look-back interval, equal to the number 
//'   of rows of data used in the scaling.
//' @param use_median A \emph{Boolean} argument: if \code{TRUE} then the 
//'   centrality (central tendency) is calculated as the \emph{median} and the 
//'   dispersion is calculated as the \emph{median absolute deviation}
//'   (\emph{MAD}).
//'   If \code{use_median} is \code{FALSE} then the centrality is calculated as 
//'   the \emph{mean} and the dispersion is calculated as the \emph{standard
//'   deviation}. (The default is \code{use_median=FALSE})
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
//' rolled_scaled <- roll::roll_scale(data=mat_rix, width=look_back, min_obs=1)
//' rolled_scaled2 <- roll_scale(mat_rix=mat_rix, look_back=look_back, use_median=FALSE)
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
  for (arma::uword it=1; it < look_back; it++) {
    sub_mat = mat_rix.rows(0, it);
    sub_mat = calc_scaled(sub_mat, use_median);
    scaled_mat.row(it) = sub_mat.row(sub_mat.n_rows-1);
  }  // end for
  // Remaining periods
  for (arma::uword it=look_back; it < num_rows; it++) {
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
//' @param res_ponse A \emph{vector} of response data.
//' @param de_sign A \emph{matrix} of design (predictor i.e.
//'   explanatory) data.
//' @param look_back The length of the look-back interval, equal to the number
//'   of elements of data used for calculating the regressions.
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
//'   if (ro_w==1) return(0)
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
  for (arma::uword it=1; it < look_back; it++) {
    sub_response = res_ponse.subvec(0, it);
    sub_design = de_sign.rows(0, it);
    z_scores(it) = calc_lm(sub_response, sub_design)["z_score"];
  }  // end for
  
  // Remaining periods
  for (arma::uword it=look_back; it < num_rows; it++) {
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
//' @param om_ega Parameter proportional to the long-term average level of
//'   variance.
//' @param al_pha The weight associated with recent realized variance updates.
//' @param be_ta The weight associated with the past variance estimates.
//' @param in_nov A \emph{vector} of innovations (random numbers).
//' 
//' @return A \emph{matrix} with two columns: the simulated returns and
//'   variance, and with the same number of rows as the length of the argument 
//'   \code{in_nov}.
//'
//' @details The function \code{sim_garch()} simulates a \emph{GARCH} process
//'   using \emph{Rcpp}.
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
  
  for (int it=1; it < len_gth; it++) {
    re_turns[it] = sqrt(vari_ance[it-1])*in_nov[it];
    vari_ance[it] = om_ega + al_pha*pow(re_turns[it], 2) + be_ta*vari_ance[it-1];
  }  // end for
  return cbind(re_turns, vari_ance);
}  // end sim_garch



////////////////////////////////////////////////////////////
//' Simulate an \emph{Ornstein-Uhlenbeck} process using \emph{Rcpp}.
//' 
//' @param eq_price The equilibrium price. 
//' @param vol_at The volatility of returns.
//' @param the_ta The strength of mean reversion.
//' @param in_nov A \emph{vector} of innovations (random numbers).
//' 
//' @return A column \emph{vector} representing the \emph{time series} of
//'   prices, with the same length as the argument \code{in_nov}.
//'
//' @details The function \code{sim_ou()} simulates an \emph{Ornstein-Uhlenbeck}
//'   process using \emph{Rcpp}, and returns A column \emph{vector} representing the 
//'   \emph{time series} of prices.
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
  for (int it=1; it < len_gth; it++) {
    re_turns[it] = the_ta*(eq_price - price_s[it-1]) + vol_at*in_nov[it-1];
    price_s[it] = price_s[it-1] * exp(re_turns[it]);
  }  // end for
  return price_s;
}  // end sim_ou



////////////////////////////////////////////////////////////
//' Recursively filter a \emph{vector} of innovations through a \emph{vector} of
//' \emph{ARIMA} coefficients.
//' 
//' @param in_nov A \emph{vector} of innovations (random numbers).
//' @param co_eff A \emph{vector} of \emph{ARIMA} coefficients.
//'
//' @return A column \emph{vector} of the same length as the argument
//'   \code{in_nov}.
//'
//' @details The function \code{sim_arima()} recursively filters a \emph{vector}
//'   of innovations through a \emph{vector} of \emph{ARIMA} coefficients, using
//'   \code{RcppArmadillo}.
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
  for (arma::uword it=look_back; it < len_gth; it++) {
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
//' @param re_turns A \emph{matrix} of excess returns data (the returns
//'   in excess of the risk-free rate).
//' @param typ_e A \emph{string} specifying the objective for calculating the
//'   weights (see Details).
//' @param max_eigen An \emph{integer} equal to the number of eigenvectors used
//'   for calculating the regularized inverse of the covariance \emph{matrix}
//'   (the default is the number of columns of \code{re_turns}).
//' @param al_pha The shrinkage intensity (the default is \code{0}).
//' @param scal_e A \emph{Boolean} specifying whether the weights should be
//'   scaled (the default is \code{scal_e=TRUE}).
//'
//' @return A column \emph{vector} of the same length as the number of columns
//'   of \code{re_turns}.
//'
//' @details The function \code{calc_weights()} calculates the optimal portfolio
//'   weights for different objective functions, using \code{RcppArmadillo}.
//' 
//'   If \code{typ_e == "max_sharpe"} (the default) then \code{calc_weights()}
//'   calculates the weights of the maximum Sharpe portfolio, by multiplying the
//'   inverse of the covariance \emph{matrix} times the mean column returns.
//'   
//'   If \code{typ_e == "min_var"} then it calculates the weights of the minimum
//'   variance portfolio under linear constraints.
//'   
//'   If \code{typ_e == "min_varpca"} then it calculates the weights of the
//'   minimum variance portfolio under quadratic constraints (which is the
//'   highest order principal component).
//' 
//'   If \code{typ_e == "rank"} then it calculates the weights as the ranks
//'   (order index) of the trailing Sharpe ratios of the portfolio assets.
//' 
//'   If \code{scal_e == TRUE} (the default) then \code{calc_weights()} scales
//'   the weights so that the resulting portfolio has the same volatility as the
//'   equally weighted portfolio.
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
                       int max_eigen = 1,
                       const double& pro_b = 0.1,
                       const double& al_pha = 0.0,
                       const bool scal_e = true) {
  // Initialize
  arma::vec weight_s(re_turns.n_cols);
  if (max_eigen == 1)  max_eigen = re_turns.n_cols;
  
  // Calculate weights depending on typ_e
  if (typ_e == "max_sharpe") {
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::mean(re_turns, 0));
    // Shrink mean_cols to the mean of re_turns
    mean_cols = ((1-al_pha)*mean_cols + al_pha*arma::mean(mean_cols));
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(re_turns, max_eigen);
    weight_s = calc_inv(re_turns, max_eigen)*mean_cols;
  } else if (typ_e == "max_sharpe_median") {
    // Mean returns by columns
    arma::vec mean_cols = arma::trans(arma::median(re_turns, 0));
    // Shrink mean_cols to the mean of re_turns
    mean_cols = ((1-al_pha)*mean_cols + al_pha*arma::median(mean_cols));
    // Apply regularized inverse
    // arma::mat in_verse = calc_inv(re_turns, max_eigen);
    weight_s = calc_inv(re_turns, max_eigen)*mean_cols;
  } else if (typ_e == "min_var") {
    // Apply regularized inverse to unit vector
    weight_s = calc_inv(re_turns, max_eigen)*arma::ones(re_turns.n_cols);
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
    // Returns of equally weighted portfolio
    // arma::vec mean_rows = arma::mean(re_turns, 1);
    // Returns of weighted portfolio
    // arma::vec returns_portf = re_turns*weight_s;
    // Scale weight_s to equally weighted portfolio and return them
    // Return weight_s/sqrt(sum(square(weight_s)));
    // Return weight_s/sum(weight_s);
    return weight_s*arma::stddev(arma::mean(re_turns, 1))/arma::stddev(re_turns*weight_s);
  }  // end if
  
  return weight_s;
}  // end calc_weights



////////////////////////////////////////////////////////////
//' Simulate (backtest) a rolling portfolio optimization strategy, using
//' \code{RcppArmadillo}.
//' 
//' @param ex_cess A \emph{matrix} of excess returns data (the returns
//'   in excess of the risk-free rate).
//' @param re_turns A \emph{matrix} of excess returns data (the returns
//'   in excess of the risk-free rate).
//' @param start_points An \emph{integer vector} of start points.
//' @param end_points An \emph{integer vector} of end points.
//' @param typ_e A \emph{string} specifying the objective for calculating the
//'   weights (see Details).
//' @param max_eigen An \emph{integer} equal to the number of eigenvectors used
//'   for calculating the regularized inverse of the covariance \emph{matrix}
//'   (the default is the number of columns of \code{re_turns}).
//' @param al_pha A numeric shrinkage intensity.  (The default is \code{0})
//' @param scal_e A \emph{Boolean} specifying whether the weights should be
//'   scaled (the default is \code{scal_e=TRUE}).
//' @param co_eff A numeric multiplier of the weights.  (The default is
//'   \code{1})
//' @param bid_offer A numeric bid-offer spread.  (The default is \code{0})
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
                    const arma::uword& max_eigen = 1,
                    const double& pro_b = 0.1,
                    const double& al_pha = 0,
                    const bool& scal_e = true,
                    const double& co_eff = 1.0,
                    const double& bid_offer = 0.0) {
  arma::vec pnl_s = zeros(re_turns.n_rows);
  arma::vec weights_past = zeros(re_turns.n_cols);
  arma::vec weight_s(re_turns.n_cols);
  
  // Perform loop over the end_points
  for (arma::uword it=1; it < end_points.size(); it++) {
    // cout << "it: " << it << endl;
    // Calculate portfolio weights
    weight_s = co_eff*calc_weights(ex_cess.rows(start_points(it-1), end_points(it-1)), typ_e, max_eigen, pro_b, al_pha, scal_e);
    // Calculate out-of-sample returns
    pnl_s.subvec(end_points(it-1)+1, end_points(it)) = re_turns.rows(end_points(it-1)+1, end_points(it))*weight_s;
    // Add transaction costs
    pnl_s.row(end_points(it-1)+1) -= bid_offer*sum(abs(weight_s - weights_past))/2;
    weights_past = weight_s;
  }  // end for
  // Return the strategy returns
  return pnl_s;
}  // end back_test
