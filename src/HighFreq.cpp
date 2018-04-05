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

// double vari_ance(NumericVector vec_tor);



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
//'   function \code{roll_sum()} is over six times faster than
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

  // warmup period
  rolling_sum[0] = vec_tor[0];
  for (int it = 1; it < look_back; it++) {
    rolling_sum[it] = rolling_sum[it-1] + vec_tor[it];
  }  // end for
  
  for (int it = look_back; it < len_gth; it++) {
    rolling_sum[it] = rolling_sum[it-1] + vec_tor[it] - vec_tor[it-look_back];
  }  // end for
  
  return rolling_sum;
}  // end roll_sum



//' Calculate the rolling weighted sum over a vector using \emph{RcppArmadillo}.
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
//'   method="convolution", sides=1)}, but it's about six times faster, and it 
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
  
  // warmup period
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
//'    (default is \code{"yang_zhang"})
//' @param look_back The size of the look-back interval, equal to the number of rows
//'   of data used for calculating the variance.
//' @param sca_le \emph{Boolean} argument: should the returns be divided by the
//'   number of seconds in each period? (default is \code{TRUE})
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



////////////////////////////
// Functions for matrix algebra
////////////////////////////


// The function inv_reg() calculates the regularized inverse 
// of the covariance matrix, by truncating the number of 
// eigen-vectors to max_eigen.
//' Calculate a time series of variance estimates over a rolling look-back interval
//' for an \emph{OHLC} time series of prices, using different range estimators
//' for variance.
//' 
//' @param vec_tor A numeric \emph{vector} of data.
//' @param look_back The length of the look-back interval, equal to the number of rows
//'   of data used for calculating the variance.
//'
//' @return A numeric \emph{vector} of the same length as the argument \code{vec_tor}.
//'
//' @details The function \code{inv_reg()} calculates a \emph{vector} of rolling 
//'   variance estimates, from over a \emph{vector} of returns, using \emph{Rcpp}.
//'
//' @examples
//' \dontrun{
//' # create minutely OHLC time series of random prices
//' vec_tor <- HighFreq::random_ohlc()
//' # calculate variance estimates for vec_tor over a 21 period interval
//' var_rolling <- HighFreq::inv_reg(vec_tor, look_back=21)
//' # calculate variance estimates for SPY
//' var_rolling <- HighFreq::inv_reg(HighFreq::SPY, calc_method="yang_zhang")
//' # calculate SPY variance without accounting for overnight jumps
//' var_rolling <- HighFreq::inv_reg(HighFreq::SPY, calc_method="rogers_satchell")
//' }
//' @export
// [[Rcpp::export]]
arma::mat inv_reg(const arma::mat& re_turns, const arma::uword& max_eigen) {
  arma::mat eigen_vec;
  arma::vec eigen_val;
  
  arma::eig_sym(eigen_val, eigen_vec, cov(re_turns));
  eigen_vec = eigen_vec.cols(eigen_vec.n_cols-max_eigen, eigen_vec.n_cols-1);
  eigen_val = 1/eigen_val.subvec(eigen_val.n_elem-max_eigen, eigen_val.n_elem-1);
  // arma::mat eigen_valmat = diagmat(eigen_val);
  
  return eigen_vec*diagmat(eigen_val)*eigen_vec.t();
  
}  // end inv_reg


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
//'   six times faster.
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
  
  // warmup period
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


// The function sharpe_weights_reg() calculates the maximum 
// Sharpe ratio portfolio weights for the matrix re_turns.
// It uses the regularized inverse of the covariance matrix.
//' Calculate a time series of variance estimates over a rolling look-back interval
//' for an \emph{OHLC} time series of prices, using different range estimators
//' for variance.
//' 
//' @param vec_tor A numeric \emph{vector} of data.
//' @param look_back The length of the look-back interval, equal to the number of rows
//'   of data used for calculating the variance.
//'
//' @return A numeric \emph{vector} of the same length as the argument \code{vec_tor}.
//'
//' @details The function \code{sharpe_weights_reg()} calculates a \emph{vector} of rolling 
//'   variance estimates, from over a \emph{vector} of returns, using \emph{Rcpp}.
//'
//' @examples
//' \dontrun{
//' # create minutely OHLC time series of random prices
//' vec_tor <- HighFreq::random_ohlc()
//' # calculate variance estimates for vec_tor over a 21 period interval
//' var_rolling <- HighFreq::sharpe_weights_reg(vec_tor, look_back=21)
//' # calculate variance estimates for SPY
//' var_rolling <- HighFreq::sharpe_weights_reg(HighFreq::SPY, calc_method="yang_zhang")
//' # calculate SPY variance without accounting for overnight jumps
//' var_rolling <- HighFreq::sharpe_weights_reg(HighFreq::SPY, calc_method="rogers_satchell")
//' }
//' @export
// [[Rcpp::export]]
arma::vec sharpe_weights_reg(const arma::mat& re_turns, 
                             const arma::vec alpha_s, 
                             const arma::vec alphas_b, 
                             const arma::uword& max_eigen) {
  arma::mat in_verse = inv_reg(re_turns, max_eigen);
  arma::vec weight_s = arma::trans(arma::mean(re_turns, 0));
  arma::vec mean_s(weight_s.n_elem);
  mean_s.fill(arma::mean(weight_s));
  
  // shrink weight_s to the mean of weight_s
  weight_s = (alphas_b % weight_s + alpha_s % mean_s);
  // apply regularized inverse
  weight_s = in_verse*weight_s;
  return weight_s/sqrt(sum(square(weight_s)));
}  // end sharpe_weights_reg



// The function roll_portf() performs a loop over the 
// end_points, subsets the re_turns matrix, and calculates 
// the PCA variances using eigen decomposition.
//' Calculate a time series of variance estimates over a rolling look-back interval
//' for an \emph{OHLC} time series of prices, using different range estimators
//' for variance.
//' 
//' @param vec_tor A numeric \emph{vector} of data.
//' @param look_back The length of the look-back interval, equal to the number of rows
//'   of data used for calculating the variance.
//'
//' @return A numeric \emph{vector} of the same length as the argument \code{vec_tor}.
//'
//' @details The function \code{roll_portf()} calculates a \emph{vector} of rolling 
//'   variance estimates, from over a \emph{vector} of returns, using \emph{Rcpp}.
//'
//' @examples
//' \dontrun{
//' # create minutely OHLC time series of random prices
//' vec_tor <- HighFreq::random_ohlc()
//' # calculate variance estimates for vec_tor over a 21 period interval
//' var_rolling <- HighFreq::roll_portf(vec_tor, look_back=21)
//' # calculate variance estimates for SPY
//' var_rolling <- HighFreq::roll_portf(HighFreq::SPY, calc_method="yang_zhang")
//' # calculate SPY variance without accounting for overnight jumps
//' var_rolling <- HighFreq::roll_portf(HighFreq::SPY, calc_method="rogers_satchell")
//' }
//' @export
// [[Rcpp::export]]
arma::mat roll_portf(const arma::mat& ex_cess, // portfolio returns
                     const arma::mat& re_turns, // portfolio returns
                     const arma::uvec& start_points, 
                     const arma::uvec& end_points, 
                     const double& al_pha, 
                     const arma::uword& max_eigen) {
  arma::vec sre_turns = zeros(re_turns.n_rows);
  arma::vec weight_s(re_turns.n_cols);
  arma::vec alpha_s(re_turns.n_cols);
  alpha_s.fill(al_pha);
  arma::vec alphas_b(re_turns.n_cols);
  alphas_b.fill(1-al_pha);
  
  // perform a loop over the end_points
  for (arma::uword i = 1; i < end_points.size(); i++) {
    // subset the returns
    arma::mat sub_returns = ex_cess.rows(start_points[i-1], end_points[i-1]);
    // calculate portfolio weights
    weight_s = sharpe_weights_reg(sub_returns, alpha_s, alphas_b, max_eigen);
    // sub_returns = re_turns.rows(end_points[i-1]+1, end_points[i]);
    sre_turns.subvec(end_points[i-1]+1, end_points[i]) = re_turns.rows(end_points[i-1]+1, end_points[i])*weight_s;
    // arma::mat foo = re_turns.rows(end_points[i-1]+1, end_points[i])*weight_s;
  }  // end for
  // return the strategy returns
  return sre_turns;
}  // end roll_portf

