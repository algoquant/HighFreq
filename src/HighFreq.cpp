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


//' Calculate the eigen decomposition of the covariance matrix of returns using
//' \emph{RcppArmadillo}.
//' 
//' @param mat_rix A numeric \emph{matrix} of returns data.
//'
//' @return A list with two elements: a numeric \emph{vector} of eigenvalues 
//'   (named "values"), and a numeric \emph{matrix} of eigenvectors (named
//'   "vectors").
//'
//' @details The function \code{calc_eigen()} first calculates the covariance 
//'   matrix of the \code{mat_rix}, and then calculates its eigen decomposition.
//'
//' @examples
//' \dontrun{
//' # Create random matrix
//' mat_rix <- matrix(rnorm(500), nc=5)
//' # Calculate eigen decomposition
//' ei_gen <- HighFreq::calc_eigen(scale(mat_rix, scale=FALSE))
//' # Calculate PCA
//' pc_a <- prcomp(mat_rix)
//' # Compare PCA with eigen decomposition
//' all.equal(pc_a$sdev^2, drop(ei_gen$values))
//' all.equal(abs(unname(pc_a$rotation)), abs(ei_gen$vectors))
//' }
//' @export
// [[Rcpp::export]]
List calc_eigen(const arma::mat& mat_rix) {
  arma::mat eigen_vec;
  arma::vec eigen_val;
  arma::eig_sym(eigen_val, eigen_vec, cov(mat_rix));
  // reverse the order of elements from largest eigenvalue to smallest, similar to R
  return List::create(Named("values") = arma::flipud(eigen_val),
                      Named("vectors") = arma::fliplr(eigen_vec));
}  // end calc_eigen



// The function calc_inv() calculates the regularized inverse 
// of the covariance matrix, by truncating the number of 
// eigenvectors to max_eigen.

//' Calculate the regularized inverse of the covariance matrix of returns using 
//' \emph{RcppArmadillo}.
//' 
//' @param mat_rix A numeric \emph{matrix} of returns data.
//' @param max_eigen An \emph{integer} equal to the regularization intensity
//'   (the number of eigenvalues and eigenvectors used for calculating the 
//'   regularized inverse).
//'
//' @return A numeric \emph{matrix} equal to the regularized inverse. 
//'
//' @details The function \code{calc_inv()} first calculates the covariance 
//'   matrix of the \code{mat_rix}, and then it calculates the regularized
//'   inverse from the truncated eigen decomposition.
//'   It uses only the largest \code{max_eigen} eigenvalues and their
//'   corresponding eigenvectors.
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



//' Scale (standardize) the columns of a \emph{matrix} of data using
//' \emph{RcppArmadillo}.
//' 
//' @param mat_rix A numeric \emph{matrix} of data.
//' @param use_median A \emph{Boolean} argument: if \code{TRUE} then the 
//'   centrality (central tendency) is calculated as the \emph{median} and the 
//'   dispersion is calculated as the \emph{median absolute deviation}
//'   (\emph{MAD}).
//'   If \code{use_median} is \code{FALSE} then the centrality is calculated as 
//'   the \emph{mean} and the dispersion is calculated as the \emph{standard
//'   deviation}. (The default is \code{FALSE})
//'
//' @return A numeric \emph{matrix} with the same dimensions as the input
//'   argument \code{mat_rix}.
//'
//' @details The function \code{calc_scaled()} scales (standardizes) the columns
//'   of the \code{mat_rix} argument using \emph{RcppArmadillo}.
//'   If the argument \code{use_median} is \code{FALSE} (the default), then it
//'   performs a similar calculation as the standard \emph{R} function
//'   \code{scale()}, and it calculates the centrality (central tendency) as the
//'   \emph{mean} and the dispersion as the \emph{standard deviation}.
//'   If the argument \code{use_median} is \code{TRUE}, then it calculates the
//'   centrality as the \emph{median} and the dispersion as the \emph{median
//'   absolute deviation} (\emph{MAD}).
//'   
//'   The function \code{calc_scaled()} uses \emph{RcppArmadillo} and is about
//'   \emph{5} times faster than function \code{scale()}, for a matrix with
//'   \emph{1,000} rows and \emph{20} columns.
//'   
//' @examples
//' \dontrun{
//' mat_rix <- matrix(rnorm(20000), nc=20)
//' scale_d <- calc_scaled(mat_rix=mat_rix, use_median=FALSE)
//' scale_d2 <- scale(mat_rix)
//' all.equal(scale_d, scale_d2, check.attributes=FALSE)
//' library(microbenchmark)
//' summary(microbenchmark(
//'   pure_r=scale(mat_rix),
//'   rcpp=calc_scaled(mat_rix=mat_rix, use_median=FALSE),
//'   times=100))[, c(1, 4, 5)]  # end microbenchmark summary
//' }
//' @export
// [[Rcpp::export]]
arma::mat calc_scaled(const arma::mat& mat_rix, 
                      const bool use_median=false) {
  arma::mat scaled_mat(mat_rix.n_rows, mat_rix.n_cols);
  arma::vec scale_d(mat_rix.n_rows);
  double cen_ter;
  
  // perform a loop over the columns
  for (arma::uword it = 0; it < mat_rix.n_cols; it++) {
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


//' Calculate the variance of a vector using \emph{Rcpp}.
//' 
//' @param vec_tor A numeric \emph{vector} of data.
//'
//' @return A single \emph{numeric} value.
//'
//' @details The function \code{vari_ance()} calculates the variance of a vector
//'   using \emph{Rcpp}. The function \code{vari_ance()} is slightly faster than
//'   the \emph{R} function \code{var()}.
//' 
//' @examples
//' \dontrun{
//' # calculate variance 
//' HighFreq::vari_ance(rnorm(1000))
//' }
//' @export
// [[Rcpp::export]]
double vari_ance(NumericVector vec_tor) {
  return sum(pow(vec_tor - sum(vec_tor)/vec_tor.size(), 2))/(vec_tor.size()-1);
}  // end vari_ance



//' Perform multivariate linear regression using RcppArmadillo.
//' 
//' @param res_ponse A numeric \emph{vector} of response data.
//' @param de_sign A numeric \emph{matrix} of design (predictor i.e.
//'   explanatory) data.
//' 
//' @return A named list with three elements: a numeric \emph{matrix} of 
//'   coefficients (named \emph{"coefficients"}), the \emph{z-score} of the last
//'   residual (named \emph{"z_score"}), and a numeric \emph{vector} with the 
//'   R-squared and F-statistic (named \emph{"stats"}). The numeric 
//'   \emph{matrix} of coefficients named \emph{"coefficients"} containes the 
//'   alpha and beta coefficients, and their \emph{t-values} and 
//'   \emph{p-values}.
//'
//' @details The function \code{calc_lm()} performs the same calculations as the
//'   function \code{lm()} from package \emph{stats}.
//'   It uses \emph{RcppArmadillo} and is about \emph{10} times faster than
//'   \code{lm()}.
//'   The code was inspired by this article (but it's not identical to it):
//'   http://gallery.rcpp.org/articles/fast-linear-model-with-armadillo/
//'
//' @examples
//' \dontrun{
//' # Define design matrix with explanatory variables
//' len_gth <- 100; n_var <- 5
//' de_sign <- matrix(rnorm(n_var*len_gth), nc=n_var)
//' # response equals linear form plus error terms
//' weight_s <- rnorm(n_var)
//' res_ponse <- -3 + de_sign %*% weight_s + rnorm(len_gth, sd=0.5)
//' # perform multivariate regression using lm()
//' reg_model <- lm(res_ponse ~ de_sign)
//' sum_mary <- summary(reg_model)
//' # perform multivariate regression using calc_lm()
//' reg_model_arma <- calc_lm(res_ponse=res_ponse, de_sign=de_sign)
//' reg_model_arma$coefficients
//' # compare the outputs of both functions
//' all.equal(reg_model_arma$coefficients[, "coeff"], unname(coef(reg_model)))
//' all.equal(unname(reg_model_arma$coefficients), unname(sum_mary$coefficients))
//' all.equal(drop(reg_model_arma$residuals), unname(reg_model$residuals))
//' all.equal(unname(reg_model_arma$stats), c(sum_mary$r.squared, unname(sum_mary$fstatistic[1])))
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List calc_lm(const arma::vec& res_ponse, const arma::mat& de_sign) {
  // add column for intercept to explanatory matrix
  arma::mat design_p = join_rows(ones(de_sign.n_rows), de_sign);
  int num_rows = de_sign.n_rows, num_cols = design_p.n_cols;
  int deg_free = (num_rows - num_cols);
  
  // fit the model res_ponse ~ de_sign, and calculate alpha and beta coefficients
  arma::colvec co_eff = arma::solve(design_p, res_ponse);
  // calculate residuals
  arma::colvec resid_uals = res_ponse - design_p*co_eff;
  
  // calculate TSS, RSS, and ESS
  double tot_sumsq = (num_rows-1)*arma::var(res_ponse);
  double res_sumsq = arma::dot(resid_uals, resid_uals);
  double exp_sumsq = tot_sumsq - res_sumsq;
  
  // calculate R-squared and F-statistic
  double r_squared = exp_sumsq/tot_sumsq;
  double f_stat = (exp_sumsq*deg_free)/(res_sumsq*(num_cols-1));
  // arma::rowvec stat_s=join_horiz(r_squared, f_stat);
  Rcpp::NumericVector stat_s(2);
  stat_s(0) = r_squared;
  stat_s(1) = f_stat;
  stat_s.attr("names") = Rcpp::CharacterVector::create("R-squared", "F-statistic");
  
  // calculate standard errors of beta coefficients
  arma::colvec std_err = arma::sqrt(res_sumsq/deg_free*arma::diagvec(arma::pinv(arma::trans(design_p)*design_p)));
  // calculate t-values and p-values of beta coefficients
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

//' Calculate the rolling sum over a vector using \emph{Rcpp}.
//' 
//' @param vec_tor A numeric \emph{vector} of data.
//' @param look_back The length of the look-back interval, equal to the number 
//'   of elements of data used for calculating the sum.
//'
//' @return A numeric \emph{vector} of the same length as the argument
//'   \code{vec_tor}.
//'
//' @details The function \code{roll_sum()} calculates a \emph{vector} of 
//'   rolling sums, over a \emph{vector} of data, using \emph{Rcpp}.  The
//'   function \code{roll_sum()} is over \emph{6} times faster than
//'   \code{rutils::roll_sum()} which uses vectorized \emph{R}.
//'
//' @examples
//' \dontrun{
//' # calculate rolling sums over 11-period intervals
//' sum_rolling <- HighFreq::roll_sum(rnorm(1000), look_back=11)
//' }
//' @export
// [[Rcpp::export]]
NumericVector roll_sum(NumericVector vec_tor, int look_back) {
  int len_gth = vec_tor.size();
  NumericVector rolling_sum(len_gth);

  // startup period
  rolling_sum[0] = vec_tor[0];
  for (int it = 1; it < look_back; it++) {
    rolling_sum[it] = rolling_sum[it-1] + vec_tor[it];
  }  // end for
  
  for (int it = look_back; it < len_gth; it++) {
    rolling_sum[it] = rolling_sum[it-1] + vec_tor[it] - vec_tor[it-look_back];
  }  // end for
  
  return rolling_sum;
}  // end roll_sum



//' Calculate the rolling weighted sum over a vector of data using
//' \emph{RcppArmadillo}.
//' 
//' @param vec_tor A numeric \emph{vector} of data.
//' @param wei_ghts A numeric \emph{vector} of weights.
//'
//' @return A numeric \emph{vector} of the same length as the argument
//'   \code{vec_tor}.
//'
//' @details The function \code{roll_wsum()} calculates the rolling weighted sum
//'   of a vector over its past values (a convolution with the \emph{vector} of 
//'   weights), using \emph{RcppArmadillo}. It performs a similar calculation as
//'   the standard \emph{R} function \code{filter(x=vec_tor, filter=wei_ghts, 
//'   method="convolution", sides=1)}, but it's about \emph{6} times faster, and it 
//'   doesn't produce any \emph{NA} values.
//'   
//' @examples
//' \dontrun{
//' # First example
//' # create vector from historical prices
//' vec_tor <- as.numeric(rutils::env_etf$VTI[, 6])
//' # create simple weights
//' wei_ghts <- c(1, rep(0, 10))
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_wsum(vec_tor=vec_tor, wei_ghts=rev(wei_ghts))
//' # compare with original
//' all.equal(vec_tor, as.numeric(weight_ed))
//' # Second example
//' # create exponentially decaying weights
//' wei_ghts <- exp(-0.2*1:11)
//' wei_ghts <- wei_ghts/sum(wei_ghts)
//' # calculate rolling weighted sum
//' weight_ed <- HighFreq::roll_wsum(vec_tor=vec_tor, wei_ghts=rev(wei_ghts))
//' # calculate rolling weighted sum using filter()
//' filter_ed <- filter(x=vec_tor, filter=wei_ghts, method="convolution", sides=1)
//' # compare both methods
//' all.equal(as.numeric(filter_ed[-(1:11)]), as.numeric(weight_ed[-(1:11)]))
//' }
//' @export
// [[Rcpp::export]]
arma::vec roll_wsum(const arma::vec& vec_tor, const arma::vec& wei_ghts) {
  uword len_gth = vec_tor.n_elem;
  uword look_back = wei_ghts.n_elem;
  arma::vec rolling_sum(len_gth);
  // arma::vec rev_weights = arma::reverse(wei_ghts);
  arma::vec rev_weights = wei_ghts;
  
  // startup period
  rolling_sum.subvec(0, look_back-2) = vec_tor.subvec(0, look_back-2);
  
  // remaining periods
  for (uword it = look_back-1; it < len_gth; it++) {
    rolling_sum(it) = arma::dot(rev_weights, vec_tor.subvec(it-look_back+1, it));
  }  // end for
  
  return rolling_sum;
}  // end roll_wsum



//' Calculate a time series of variance estimates over a rolling look-back
//' interval for an \emph{OHLC} time series of prices, using different range
//' estimators for variance.
//' 
//' Currently only works for vectors
//'
//' @param oh_lc An \emph{OHLC} time series of prices in \emph{xts} format.
//' @param calc_method \emph{character} string representing method for estimating
//'   variance.  The methods include:
//'   \itemize{
//'     \item "close" close to close,
//'     \item "garman_klass" Garman-Klass,
//'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
//'     \item "rogers_satchell" Rogers-Satchell,
//'     \item "yang_zhang" Yang-Zhang,
//'    }
//'    (The default is \code{"yang_zhang"})
//' @param look_back The size of the look-back interval, equal to the number of rows
//'   of data used for calculating the variance.
//' @param sca_le \emph{Boolean} argument: should the returns be divided by the
//'   number of seconds in each period? (The default is \code{TRUE})
//'
//' @return An \emph{xts} time series with a single column and the same number of
//'   rows as the argument \code{oh_lc}.
//'
//' @details The function \code{roll_var()} calculates a time series of rolling 
//'   variance estimates of percentage returns, from over a
//'   \emph{vector} of returns, using several different variance estimation
//'   methods based on the range of \emph{OHLC} prices.
//'
//'   If \code{sca_le} is \code{TRUE} (the default), then the variance is divided
//'   by the squared differences of the time index (which scales the variance to
//'   units of variance per second squared.) This is useful for example, when
//'   calculating intra-day variance from minutely bar data, because dividing
//'   returns by the number of seconds decreases the effect of overnight price
//'   jumps.
//'
//'   If \code{sca_le} is \code{TRUE} (the default), then the variance is
//'   expressed in the scale of the time index of the \emph{OHLC} time series.
//'   For example, if the time index is in seconds, then the variance is given in
//'   units of variance per second squared.  If the time index is in days, then
//'   the variance is equal to the variance per day squared.
//'
//'   The time index of the \code{oh_lc} time series is assumed to be in
//'   \emph{POSIXct} format, so that its internal value is equal to the number of
//'   seconds that have elapsed since the \emph{epoch}.
//'
//'   The methods \code{"close"}, \code{"garman_klass_yz"}, and
//'   \code{"yang_zhang"} do account for close-to-open price jumps, while the
//'   methods \code{"garman_klass"} and \code{"rogers_satchell"} do not account
//'   for close-to-open price jumps.
//'
//'   The default method is \code{"yang_zhang"}, which theoretically has the
//'   lowest standard error among unbiased estimators.
//'
//'   The function \code{roll_var()} performs the same calculations as the
//'   function \code{volatility()} from package
//'   \href{https://cran.r-project.org/web/packages/TTR/index.html}{TTR}, but
//'   it's a little faster because it uses function RcppRoll::roll_sd(), and it
//'   performs less data validation.
//'
//' @examples
//' \dontrun{
//' # create minutely OHLC time series of random prices
//' oh_lc <- HighFreq::random_ohlc()
//' # calculate variance estimates for oh_lc over a 21 period interval
//' var_rolling <- HighFreq::roll_var(oh_lc, look_back=21)
//' # calculate variance estimates for SPY
//' var_rolling <- HighFreq::roll_var(HighFreq::SPY, calc_method="yang_zhang")
//' # calculate SPY variance without accounting for overnight jumps
//' var_rolling <- HighFreq::roll_var(HighFreq::SPY, calc_method="rogers_satchell")
//' }
//' @export
// [[Rcpp::export]]

NumericVector roll_var(NumericVector vec_tor, int look_back) {
  int len_gth = vec_tor.size();
  NumericVector var_vec(len_gth);
  NumericVector roll_mean(len_gth);
  
  var_vec[0] = 0;
  for (int it = 1; it < vec_tor.size(); it++) {
    var_vec[it] = vari_ance(vec_tor[Range(std::max(it-look_back+1, 0), it)]);
  }  // end for
  
  return var_vec;
}  // end roll_var



//' Perform a rolling scaling (standardization) of the columns of a
//' \emph{matrix} of data using \emph{RcppArmadillo}.
//' 
//' @param mat_rix A numeric \emph{matrix} of data.
//' @param look_back The length of the look-back interval, equal to the number 
//'   of rows of data used in the scaling.
//' @param use_median A \emph{Boolean} argument: if \code{TRUE} then the 
//'   centrality (central tendency) is calculated as the \emph{median} and the 
//'   dispersion is calculated as the \emph{median absolute deviation}
//'   (\emph{MAD}).
//'   If \code{use_median} is \code{FALSE} then the centrality is calculated as 
//'   the \emph{mean} and the dispersion is calculated as the \emph{standard
//'   deviation}. (The default is \code{FALSE})
//'
//' @return A numeric \emph{matrix} with the same dimensions as the input
//'   argument \code{mat_rix}.
//'
//' @details The function \code{roll_scale()} performs a rolling scaling
//'   (standardization) of the columns of the \code{mat_rix} argument using
//'   \emph{RcppArmadillo}.
//'   The function \code{roll_scale()} performs a loop over the rows of 
//'   \code{mat_rix}, then subsets a number of previous (past) rows equal to 
//'   \code{look_back}, and scales the subset matrix.  It assigns the last row
//'   of the scaled subset matrix to the return matrix.
//'   
//'   If the argument \code{use_median} is \code{FALSE} (the default), then it
//'   performs a similar calculation as the function \code{roll::roll_scale()}.
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
//' @export
// [[Rcpp::export]]
arma::mat roll_scale(const arma::mat& mat_rix, 
                     const arma::uword& look_back,
                     const bool use_median=false) {
  arma::uword num_rows = mat_rix.n_rows;
  arma::mat scaled_mat(num_rows, mat_rix.n_cols);
  arma::mat sub_mat;
  
  // startup period
  scaled_mat.row(0) = mat_rix.row(0);
  for (uword it = 1; it < look_back; it++) {
    sub_mat = mat_rix.rows(0, it);
    sub_mat = calc_scaled(sub_mat, use_median);
    scaled_mat.row(it) = sub_mat.row(sub_mat.n_rows-1);
  }  // end for
  // remaining periods
  for (uword it = look_back; it < num_rows; it++) {
    sub_mat = mat_rix.rows(it-look_back+1, it);
    sub_mat = calc_scaled(sub_mat, use_median);
    scaled_mat.row(it) = sub_mat.row(sub_mat.n_rows-1);
  }  // end for
  
  return scaled_mat;
}  // end roll_scale



//' Perform rolling regressions over the rows of the design matrix, and
//' calculate a \emph{vector} of z-scores of the residuals.
//' 
//' @param res_ponse A numeric \emph{vector} of response data.
//' @param de_sign A numeric \emph{matrix} of design (predictor i.e.
//'   explanatory) data.
//' @param look_back The length of the look-back interval, equal to the number 
//'   of elements of data used for calculating the regressions.
//'
//' @return A numeric \emph{vector} of the same length as the number of rows of
//'   \code{de_sign}.
//'
//' @details The function \code{roll_zscores()} performs rolling regressions 
//'   along the rows of the design matrix \code{de_sign}, using function
//'   \code{calc_lm()}. 
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
//' # calculate Z-scores from rolling time series regression using RcppArmadillo
//' look_back <- 11
//' clo_se <- as.numeric(Cl(rutils::env_etf$VTI))
//' date_s <- xts::.index(rutils::env_etf$VTI)
//' z_scores <- HighFreq::roll_zscores(res_ponse=clo_se, 
//'  de_sign=matrix(as.numeric(date_s), nc=1), 
//'  look_back=look_back)
//' # Define design matrix with explanatory variables
//' len_gth <- 100; n_var <- 5
//' de_sign <- matrix(rnorm(n_var*len_gth), nc=n_var)
//' # response equals linear form plus error terms
//' weight_s <- rnorm(n_var)
//' res_ponse <- -3 + de_sign %*% weight_s + rnorm(len_gth, sd=0.5)
//' # calculate Z-scores from rolling multivariate regression using RcppArmadillo
//' look_back <- 11
//' z_scores <- HighFreq::roll_zscores(res_ponse=res_ponse, de_sign=de_sign, look_back=look_back)
//' # calculate z-scores in R from rolling multivariate regression using lm()
//' z_scores_r <- sapply(1:NROW(de_sign), function(ro_w) {
//'   if (ro_w==1) return(0)
//'   start_point <- max(1, ro_w-look_back+1)
//'   sub_response <- res_ponse[start_point:ro_w]
//'   sub_design <- de_sign[start_point:ro_w, ]
//'   reg_model <- lm(sub_response ~ sub_design)
//'   resid_uals <- reg_model$residuals
//'   resid_uals[NROW(resid_uals)]/sd(resid_uals)
//' })  # end sapply
//' # compare the outputs of both functions
//' all.equal(unname(z_scores[-(1:look_back)]), 
//'   unname(z_scores_r[-(1:look_back)]))
//' }
//' @export
// [[Rcpp::export]]
arma::vec roll_zscores(const arma::vec& res_ponse, 
                       const arma::mat& de_sign, 
                       const arma::uword& look_back) {
  arma::uword num_rows = de_sign.n_rows;
  arma::vec z_scores(num_rows);
  arma::vec sub_response;
  arma::mat sub_design;
  Rcpp::List lm_list;
  
  // startup period
  for (uword it = 1; it < look_back; it++) {
    sub_response = res_ponse.subvec(0, it);
    sub_design = de_sign.rows(0, it);
    lm_list = calc_lm(sub_response, sub_design);
    z_scores[it] = lm_list["z_score"];
  }  // end for
  
  // remaining periods
  for (uword it = look_back; it < num_rows; it++) {
    sub_response = res_ponse.subvec(it-look_back+1, it);
    sub_design = de_sign.rows(it-look_back+1, it);
    lm_list = calc_lm(sub_response, sub_design);
    z_scores[it] = lm_list["z_score"];
  }  // end for
  
  return z_scores;
}  // end roll_zscores




////////////////////////////
// Functions for simulation
////////////////////////////



//' Simulate a \emph{GARCH} process using \emph{Rcpp}.
//' 
//' @param om_ega Parameter proportional to the long-term average level of variance.
//' @param al_pha The weight associated with recent realized variance updates.
//' @param be_ta The weight associated with the past variance estimates.
//' @param in_nov A numeric \emph{vector} of innovations (random numbers).
//' 
//' @return A numeric \emph{matrix} with two columns: the simulated returns and
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



//' Simulate an \emph{Ornstein-Uhlenbeck} process using \emph{Rcpp}.
//' 
//' @param eq_price The equilibrium price. 
//' @param vol_at The volatility of returns.
//' @param the_ta The strength of mean reversion.
//' @param in_nov A numeric \emph{vector} of innovations (random numbers).
//' 
//' @return A numeric \emph{vector} representing the time series of prices, with
//'   the same length as the argument \code{in_nov}.
//'
//' @details The function \code{sim_ou()} simulates an \emph{Ornstein-Uhlenbeck}
//'   process using \emph{Rcpp}, and returns a \emph{vector} representing the 
//'   time series of prices.
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



//' Recursively filter a vector of innovations through a vector of \emph{ARIMA} 
//' coefficients.
//' 
//' @param in_nov A numeric \emph{vector} of innovations (random numbers).
//' @param co_eff A numeric \emph{vector} of \emph{ARIMA} coefficients.
//'
//' @return A numeric \emph{vector} of the same length as the argument
//'   \code{in_nov}.
//'
//' @details The function \code{sim_arima()} recursively filters a vector of
//'   innovations through a vector of \emph{ARIMA} coefficients, using 
//'   \emph{RcppArmadillo}.
//'   It performs the same calculation as the standard \emph{R} function 
//'   \code{filter(x=in_nov, filter=co_eff, method="recursive")}, but it's about
//'   \emph{6} times faster.
//'   
//' @examples
//' \dontrun{
//' # create vector of innovations
//' in_nov <- rnorm(100)
//' # create ARIMA coefficients
//' co_eff <- c(-0.8, 0.2)
//' # calculate recursive filter using filter()
//' filter_ed <- filter(in_nov, filter=co_eff, method="recursive")
//' # calculate recursive filter using RcppArmadillo
//' ari_ma <- HighFreq::sim_arima(in_nov, rev(co_eff))
//' # compare the two methods
//' all.equal(as.numeric(ari_ma), as.numeric(filter_ed))
//' }
//' @export
// [[Rcpp::export]]
arma::vec sim_arima(const arma::vec& in_nov, const arma::vec& co_eff) {
  uword len_gth = in_nov.n_elem;
  uword look_back = co_eff.n_elem;
  arma::vec ari_ma(len_gth);
  
  // startup period
  ari_ma(0) = in_nov(0);
  ari_ma(1) = in_nov(1) + co_eff(look_back-1) * ari_ma(0);
  for (uword it = 2; it < look_back-1; it++) {
    ari_ma(it) = in_nov(it) + arma::dot(co_eff.subvec(look_back-it, look_back-1), ari_ma.subvec(0, it-1));
  }  // end for
  
  // remaining periods
  for (uword it = look_back; it < len_gth; it++) {
    ari_ma(it) = in_nov(it) + arma::dot(co_eff, ari_ma.subvec(it-look_back, it-1));
  }  // end for
  
  return ari_ma;
}  // end sim_arima



////////////////////////////
// Functions for backtests
////////////////////////////


//' Calculate the portfolio weights with the maximum Sharpe ratio.
//' 
//' @param re_turns A numeric \emph{matrix} of excess returns data (the returns
//'   in excess of the risk-free rate).
//' @param max_eigen An \emph{integer} equal to the regularization intensity
//'   (the number of eigenvalues and eigenvectors used for calculating the 
//'   regularized inverse).
//' @param al_pha The shrinkage intensity.  (The default is \code{0})
//'
//' @return A numeric \emph{vector} of the same length as the number of columns
//'   of \code{re_turns}.
//'
//' @details The function \code{calc_weights()} calculates the scaled weights of
//'   the portfolio with the maximum Sharpe ratio, using \emph{RcppArmadillo}.
//'   
//'   It first calculates the regularized inverse of the covariance matrix of
//'   returns using function \code{HighFreq::calc_inv()}.
//'   It then estimates the vector of mean returns and applies shrinkage to it,
//'   by shrinking it to its common mean value.
//'   The shrinkage intensity \code{al_pha} determines the amount of shrinkage 
//'   that is applied, with \code{al_pha = 0} representing no shrinkage (with 
//'   the estimator of mean returns equal to the means of the columns of 
//'   \code{re_turns}), and \code{al_pha = 1} representing complete shrinkage 
//'   (with the estimator of mean returns equal to the single mean of all the
//'   columns of \code{re_turns})
//'   
//'   The function \code{calc_weights()} calculates the weights by multiplying 
//'   the inverse of the covariance matrix times the estimator of the mean
//'   returns. It finally scales the weights by their sum, so that the sum of
//'   the weights is equal to \code{1}.
//' 
//' @examples
//' \dontrun{
//' # Calculate ETF prices
//' sym_bols <- colnames(rutils::env_etf$price_s)
//' sym_bols <- sym_bols[!(sym_bols=="VXX")]
//' price_s <- rutils::env_etf$price_s[, sym_bols]
//' # Carry forward non-NA prices
//' price_s <- zoo::na.locf(price_s)
//' price_s <- na.omit(price_s)
//' # Calculate simple ETF returns
//' re_turns <- rutils::diff_it(price_s)
//' # Calculate covariance matrix
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
//' weights_r <- drop(weight_s/sum(weight_s))
//' # Calculate weights using RcppArmadillo
//' weight_s <- drop(HighFreq::calc_weights(re_turns, max_eigen, al_pha=al_pha))
//' all.equal(weight_s, weights_r)
//' }
//' @export
// [[Rcpp::export]]
arma::vec calc_weights(const arma::mat& re_turns, 
                       const arma::uword& max_eigen,
                       const double& al_pha = 0
) {
  arma::mat in_verse = calc_inv(re_turns, max_eigen);
  arma::vec weight_s = arma::trans(arma::mean(re_turns, 0));
  
  // shrink weight_s to the mean of weight_s
  weight_s = ((1-al_pha)*weight_s + al_pha*arma::mean(weight_s));
  // apply regularized inverse
  weight_s = in_verse*weight_s;
  // scale weight_s and return them
  // return weight_s/sqrt(sum(square(weight_s)));
  return weight_s/sum(weight_s);
}  // end calc_weights



//' Simulate (backtest) a rolling portfolio optimization strategy, using
//' \emph{RcppArmadillo}.
//' 
//' @param ex_cess A numeric \emph{matrix} of excess returns data (the returns
//'   in excess of the risk-free rate).
//' @param re_turns A numeric \emph{matrix} of excess returns data (the returns
//'   in excess of the risk-free rate).
//' @param start_points An integer \emph{vector} of start points.
//' @param end_points An integer \emph{vector} of end points.
//' @param max_eigen An \emph{integer} equal to the regularization intensity
//'   (the number of eigenvalues and eigenvectors used for calculating the 
//'   regularized inverse).
//' @param al_pha The shrinkage intensity.  (The default is \code{0})
//'
//' @return A numeric \emph{vector} of strategy returns, with the same length as
//'   the number of rows of \code{re_turns}.
//'
//' @details The function \code{roll_portf()} calculates a \emph{vector} of 
//'   returns for a rolling portfolio optimization strategy. The function 
//'   \code{roll_portf()} performs a loop over the \code{end_points}, and 
//'   subsets the \emph{matrix} of excess returns \code{ex_cess} along its rows,
//'   between the corresponding end point and the start point. It passes the
//'   subset matrix of excess returns into the function \code{calc_weights()},
//'   which calculates the optimal portfolio weights. It then multiplies the
//'   weights times the future portfolio returns, to calculate the out-of-sample
//'   strategy returns.
//'
//' @examples
//' \dontrun{
//' # Calculate ETF prices
//' sym_bols <- colnames(rutils::env_etf$price_s)
//' sym_bols <- sym_bols[!(sym_bols=="VXX")]
//' price_s <- rutils::env_etf$price_s[, sym_bols]
//' # Carry forward non-NA prices
//' price_s <- zoo::na.locf(price_s)
//' price_s <- na.omit(price_s)
//' # Calculate simple ETF returns
//' re_turns <- rutils::diff_it(price_s)
//' # Calculate the daily excess returns
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
//' strat_rets <- HighFreq::roll_portf(ex_cess, re_turns, 
//'                                   start_points-1, end_points-1, 
//'                                   max_eigen, al_pha)
//' strat_rets <- xts(strat_rets, index(re_turns))
//' colnames(strat_rets) <- "strat_rets"
//' # Plot dygraph of strategy
//' dygraphs::dygraph(cumsum(strat_rets), 
//'   main="Cumulative Returns of Max Sharpe Portfolio Strategy")
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_portf(const arma::mat& ex_cess, // portfolio excess returns
                     const arma::mat& re_turns, // portfolio returns
                     const arma::uvec& start_points, 
                     const arma::uvec& end_points, 
                     const arma::uword& max_eigen,
                     const double& al_pha = 0
) {
  arma::vec sre_turns = zeros(re_turns.n_rows);
  arma::vec weight_s(re_turns.n_cols);
  
  // initial period
  // sre_turns.subvec(0, end_points[1]) = re_turns.rows(0, end_points[1]);
  // perform a loop over the end_points
  for (arma::uword it = 1; it < end_points.size(); it++) {
    // cout << "it: " << it << endl;
    // subset the returns
    // arma::mat sub_returns = ex_cess.rows(start_points[it-1], end_points[it-1]);
    // calculate portfolio weights
    weight_s = calc_weights(ex_cess.rows(start_points[it-1], end_points[it-1]), max_eigen, al_pha);
    // sub_returns = re_turns.rows(end_points[it-1]+1, end_points[it]);
    sre_turns.subvec(end_points[it-1]+1, end_points[it]) = re_turns.rows(end_points[it-1]+1, end_points[it])*weight_s;
    // arma::mat foo = re_turns.rows(end_points[it-1]+1, end_points[it])*weight_s;
  }  // end for
  // return the strategy returns
  return sre_turns;
}  // end roll_portf
