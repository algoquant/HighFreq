##########################################################################
#' Calculate a random \emph{TAQ} time series of prices and trading volumes, in
#' \emph{xts} format.
#'
#' Calculate a \emph{TAQ} time series of random prices following geometric
#' Brownian motion, combined with random trading volumes.
#'
#' @export
#' @param \code{vol_at} The volatility per period of the \code{in_dex} time index
#'   (default is \code{6.5e-05} per second, or about \code{0.01=1.0\%} per day).
#' @param \code{dri_ft} The drift per period of the \code{in_dex} time index (default
#'   is 0.0).
#' @param \code{in_dex} The time index for the \emph{TAQ} time series.
#' @param bid_offer The bid-offer spread expressed as a fraction of the prices
#'   (default is 0.001=10bps).
#'
#' @return An \emph{xts} time series, with time index equal to the input
#'   \code{in_dex} time index, and with four columns containing the bid, ask,
#'   and trade prices, and the trade volume.
#'
#' @details The function \code{random_taq()} calculates an \emph{xts} time
#'   series with four columns containing random prices following geometric
#'   Brownian motion: the bid, ask, and trade prices, combined with random trade
#'   volume data.
#'   If \code{in_dex} isn't supplied as an argument, then by default it's
#'   equal to the secondly index over the two previous calendar days.
#'
#' @examples
#' # Create secondly TAQ time series of random prices
#' ta_q <- HighFreq::random_taq()
#' # Create random TAQ time series from SPY index
#' ta_q <- HighFreq::random_taq(in_dex=index(HighFreq::SPY["2012-02-13/2012-02-15"]))

random_taq <- function(vol_at=6.5e-5, dri_ft=0.0,
  in_dex=seq(from=as.POSIXct(paste(Sys.Date()-3, "09:30:00")),
             to=as.POSIXct(paste(Sys.Date()-1, "16:00:00")), by="1 sec"),
  bid_offer=0.001, ...) {
  len_gth <- NROW(in_dex)
  # Create xts of random prices following geometric Brownian motion
  ta_q <- xts(exp(cumsum(vol_at*rnorm(len_gth) + dri_ft - vol_at^2/2)),
              order.by=in_dex)
  # Create vector of random bid-offer spreads
  bid_offer <- bid_offer*(1 + runif(len_gth))/2
  # Create TAQ data from bid and offer prices
  ta_q <- merge(ta_q*(1-bid_offer), ta_q*(1+bid_offer))
  # Add traded price to TAQ data
  r_unif <- runif(len_gth)
  ta_q <- merge(ta_q, r_unif*ta_q[, 1] + (1-r_unif)*ta_q[, 2])
  # Add trade volume column
  ta_q <- merge(ta_q, sample(x=10*(2:18), size=len_gth, replace=TRUE))
  colnames(ta_q) <- c("Bid.Price", "Ask.Price", "Trade.Price", "Volume")
  ta_q
}  # end random_taq




##########################################################################
#' Calculate a random \emph{OHLC} time series of prices and trading volumes, in
#' \emph{xts} format.
#'
#' Calculate a random \emph{OHLC} time series either by simulating random prices
#' following geometric Brownian motion, or by randomly sampling from an input
#' time series.
#'
#' @export
#' @param \code{oh_lc} An \emph{OHLC} time series of prices and trading volumes, in
#'   \emph{xts} format (default is \emph{NULL}).
#' @param \code{vol_at} The volatility per period of the \code{in_dex} time index
#'   (default is \code{6.5e-05} per second, or about \code{0.01=1.0\%} per day).
#' @param \code{dri_ft} The drift per period of the \code{in_dex} time index (default
#'   is 0.0).
#' @param \code{in_dex} The time index for the \emph{OHLC} time series.
#' @param \code{re_duce} \emph{Boolean} argument: should \code{oh_lc} time series be
#'   transformed to reduced form? (default is \code{TRUE})
#'
#' @return An \emph{xts} time series with the same dimensions and the same time
#'   index as the input \code{oh_lc} time series.
#'
#' @details If the input \code{oh_lc} time series is \emph{NULL} (the default),
#'   then the function \code{random_ohlc()} simulates a minutely \emph{OHLC}
#'   time series of random prices following geometric Brownian motion, over the
#'   two previous calendar days.
#'
#'   If the input \code{oh_lc} time series is not \emph{NULL}, then the rows of
#'   \code{oh_lc} are randomly sampled, to produce a random time series.
#'
#'   If \code{re_duce} is \code{TRUE} (the default), then the \code{oh_lc} time
#'   series is first transformed to reduced form, then randomly sampled, and
#'   finally converted to standard form.
#'
#'   Note: randomly sampling from an intraday time series over multiple days
#'   will cause the overnight price jumps to be re-arranged into intraday price
#'   jumps.  This will cause moment estimates to become inflated compared to the
#'   original time series.
#'
#' @examples
#' # Create minutely synthetic OHLC time series of random prices
#' oh_lc <- HighFreq::random_ohlc()
#' # Create random time series from SPY by randomly sampling it
#' oh_lc <- HighFreq::random_ohlc(oh_lc=HighFreq::SPY["2012-02-13/2012-02-15"])

random_ohlc <- function(oh_lc=NULL, re_duce=TRUE, vol_at=6.5e-5, dri_ft=0.0,
    in_dex=seq(from=as.POSIXct(paste(Sys.Date()-3, "09:30:00")),
      to=as.POSIXct(paste(Sys.Date()-1, "16:00:00")), by="1 sec"), ...) {
  if (is.null(oh_lc)) {
    len_gth <- NROW(in_dex)
    # Create xts of random prices following geometric Brownian motion
    x_ts <- xts(exp(cumsum(vol_at*rnorm(len_gth) + dri_ft - vol_at^2/2)), order.by=in_dex)
    # Add trade volume column
    x_ts <- merge(x_ts, volume=sample(x=10*(2:18), size=len_gth, replace=TRUE))
    # Aggregate to minutes OHLC data
    to.period(x=x_ts, period="minutes")
  } else {
    oh_lc <- log(oh_lc)  # transform to normal
    if (re_duce)  # Calculate reduced form of oh_lc
      oh_lc <- rutils::diff_ohlc(oh_lc)
    # randomly sample from the rows of oh_lc
    oh_lc <- xts(coredata(oh_lc)[c(1, sample(x=2:NROW(oh_lc), replace=TRUE)), ], order.by=index(oh_lc))
    # Return standard form of randomized oh_lc
    exp(rutils::diff_ohlc(oh_lc, re_duce=FALSE))
  }
}  # end random_ohlc




##########################################################################
#' Remove overnight close-to-open price jumps from an \emph{OHLC} time series,
#' by adding adjustment terms to its prices.
#'
#' @export
#' @param \code{oh_lc} An \emph{OHLC} time series of prices and trading volumes, in
#'   \emph{xts} format.
#'
#' @return An \emph{OHLC} time series with the same dimensions and the same time
#'   index as the input \code{oh_lc} time series.
#'
#' @details The function \code{remove_jumps()} removes the overnight
#'   close-to-open price jumps from an \emph{OHLC} time series, by adjusting its
#'   prices so that the first \emph{Open} price of the day is equal to the last
#'   \emph{Close} price of the previous day.
#'
#'   The function \code{remove_jumps()} adds adjustment terms to all the
#'   \emph{OHLC} prices, so that intra-day returns and volatilities are not
#'   affected.
#'
#'   The function \code{remove_jumps()} identifies overnight periods as those
#'   that are greater than 60 seconds. This assumes that intra-day periods
#'   between neighboring rows of data are 60 seconds or less.
#'
#'   The time index of the \code{oh_lc} time series is assumed to be in
#'   \emph{POSIXct} format, so that its internal value is equal to the number of
#'   seconds that have elapsed since the \emph{epoch}.
#'
#' @examples
#' # Remove overnight close-to-open price jumps from SPY data
#' oh_lc <- remove_jumps(HighFreq::SPY)

remove_jumps <- function(oh_lc) {
  # find time index of the periods greater than 60 seconds
  which_periods <- which(c(1, diff(xts::.index(oh_lc))) > 60)
  # Calculate cumulative sum of overnight price jumps
  jump_s <- numeric(NROW(oh_lc))
  jump_s[which_periods] <- as.numeric(oh_lc[which_periods, 1]) - as.numeric(oh_lc[which_periods-1, 4])
  jump_s <- cumsum(jump_s)
  # subtract overnight price jumps from OHLC
  oh_lc[, 1:4] <- coredata(oh_lc[, 1:4]) - jump_s
  oh_lc
}  # end remove_jumps




##########################################################################
#' Calculate single period percentage returns from either \emph{TAQ} or
#' \emph{OHLC} prices.
#'
#' @export
#' @param \code{x_ts} An \emph{xts} time series of either \emph{TAQ} or \emph{OHLC} data.
#' @param \code{lagg} An integer equal to the number of time periods of lag. (default
#'   is 1)
#' @param \code{col_umn} The column number to extract from the \emph{OHLC} data.
#'   (default is \code{4}, or the \emph{Close} prices column)
#' @param \code{scal_e} \emph{Boolean} argument: should the returns be divided by the
#'   number of seconds in each period? (default is \code{TRUE})
#'
#' @return A single-column \emph{xts} time series of returns.
#'
#' @details The function \code{ohlc_returns()} calculates the percentage returns
#'   for either \emph{TAQ} or \emph{OHLC} data, defined as the difference of log
#'   prices.  Multi-period returns can be calculated by setting the \code{lag}
#'   parameter to values greater than \code{1} (the default).
#'
#'   If \code{scal_e} is \code{TRUE} (the default), then the returns are divided
#'   by the differences of the time index (which scales the returns to units of
#'   returns per second.)
#'
#'   The time index of the \code{x_ts} time series is assumed to be in
#'   \emph{POSIXct} format, so that its internal value is equal to the number of
#'   seconds that have elapsed since the \emph{epoch}.
#'
#'   If \code{scal_e} is \code{TRUE} (the default), then the returns are
#'   expressed in the scale of the time index of the \code{x_ts} time series.
#'   For example, if the time index is in seconds, then the returns are given in
#'   units of returns per second.  If the time index is in days, then the
#'   returns are equal to the returns per day.
#'
#'   The function \code{ohlc_returns()} identifies the \code{x_ts} time series as
#'   \emph{TAQ} data when it has six columns, otherwise assumes it's \emph{OHLC}
#'   data. By default, for \emph{OHLC} data, it differences the \emph{Close}
#'   prices, but can also difference other prices depending on the value of
#'   \code{col_umn}.
#'
#' @examples
#' # Calculate secondly returns from TAQ data
#' re_turns <- HighFreq::ohlc_returns(x_ts=HighFreq::SPY_TAQ)
#' # Calculate close to close returns
#' re_turns <- HighFreq::ohlc_returns(x_ts=HighFreq::SPY)
#' # Calculate open to open returns
#' re_turns <- HighFreq::ohlc_returns(x_ts=HighFreq::SPY, col_umn=1)

ohlc_returns <- function(x_ts, lagg=1, col_umn=4, scal_e=TRUE) {
  # Return NULL if no data
  if (is.null(x_ts))  return(NULL)
  # Calculate mid prices
  if (NCOL(x_ts)==6)  # TAQ data has 6 columns
    re_turns <- 0.5 * (x_ts[, "Bid.Price"] + x_ts[, "Ask.Price"])
  else
    re_turns <- x_ts[, col_umn]  # OHLC data
  # Calculate returns
  re_turns <- rutils::diff_it(log(re_turns), lagg=lagg)
  if (scal_e)
    re_turns <- re_turns / c(rep(1, lagg), diff(xts::.index(re_turns), lagg=lagg))
  re_turns[1:lagg, ] <- 0
  # Colnames(re_turns) <- paste0(rutils::get_name(colnames(x_ts)[1]), ".returns")
  re_turns
}  # end ohlc_returns




##########################################################################
#' Calculate a \emph{Boolean} vector that identifies extreme tail values in a
#' single-column \emph{xts} time series or vector, over a rolling look-back
#' interval.
#'
#' @export
#' @param \code{x_ts} A single-column \emph{xts} time series, or a \emph{numeric} or
#'   \emph{Boolean} vector.
#' @param \code{look_back} The number of data points in rolling look-back interval for 
#'   estimating rolling quantile.
#' @param \code{vol_mult} The quantile multiplier.
#'
#' @return A \emph{Boolean} vector with the same number of rows as the input
#'   time series or vector.
#'
#' @details The function \code{which_extreme()} calculates a \emph{Boolean}
#'   vector, with \code{TRUE} for values that belong to the extreme tails
#'   of the distribution of values.
#'
#'   The function \code{which_extreme()} applies a version of the Hampel median
#'   filter to identify extreme values, but instead of using the median absolute
#'   deviation (MAD), it uses the \code{0.9} quantile values calculated over a
#'   rolling look-back interval.
#'
#'   Extreme values are defined as those that exceed the product of the
#'   multiplier times the rolling quantile. Extreme values belong to the fat
#'   tails of the recent (trailing) distribution of values, so they are present
#'   only when the trailing distribution of values has fat tails.  If the
#'   trailing distribution of values is closer to normal (without fat tails),
#'   then there are no extreme values.
#'
#'   The quantile multiplier \code{vol_mult} controls the threshold at which
#'   values are identified as extreme. Smaller quantile multiplier values will
#'   cause more values to be identified as extreme.
#'
#' @examples
#' # Create local copy of SPY TAQ data
#' ta_q <- HighFreq::SPY_TAQ
#' # scrub quotes with suspect bid-offer spreads
#' bid_offer <- ta_q[, "Ask.Price"] - ta_q[, "Bid.Price"]
#' sus_pect <- which_extreme(bid_offer, look_back=51, vol_mult=3)
#' # Remove suspect values
#' ta_q <- ta_q[!sus_pect]

which_extreme <- function(x_ts, look_back=51, vol_mult=2) {
# Calculate volatility as rolling quantile
  quan_tile <- caTools::runquantile(x=abs(as.numeric(x_ts)), k=look_back,
                        probs=0.9, endrule="constant", align="center")
#  quan_tile <- xts(quan_tile, order.by=index(x_ts))
#  colnames(quan_tile) <- "volat"
# Carry forward non-zero volatility values
  quan_tile[quan_tile==0] <- NA
  quan_tile[1] <- 1
  quan_tile <- rutils::na_locf(quan_tile)
#  quan_tile <- rutils::na_locf(quan_tile, fromLast=TRUE)

# extreme value if x_ts greater than scaled volatility
  ex_treme <- (abs(x_ts) > 2*vol_mult*quan_tile)
  ex_treme[1] <- FALSE
#  colnames(ex_treme) <- "suspect"

  cat("date:", format(as.Date(index(first(x_ts)))), "\tfound", sum(ex_treme), "extreme values\n")
  ex_treme
}  # end which_extreme




##########################################################################
#' Calculate a \emph{Boolean} vector that identifies isolated jumps (spikes) in
#' a single-column \emph{xts} time series or vector, over a rolling interval.
#'
#' @export
#' @inheritParams which_extreme
#'
#' @return A \emph{Boolean} vector with the same number of rows as the input
#'   time series or vector.
#'
#' @details The function \code{which_jumps()} calculates a \emph{Boolean}
#'   vector, with \code{TRUE} for values that are isolated jumps (spikes).
#'
#'   The function \code{which_jumps()} applies a version of the Hampel median
#'   filter to identify jumps, but instead of using the median absolute
#'   deviation (MAD), it uses the \code{0.9} quantile of returns calculated over
#'   a rolling interval.
#'   This is in contrast to function \code{which_extreme()}, which applies a
#'   Hampel filter to the values themselves, instead of the returns.
#'   Returns are defined as simple differences between neighboring values.
#'
#'   Jumps (or spikes), are defined as isolated values that are very different
#'   from the neighboring values, either before or after.  Jumps create
#'   pairs of large neighboring returns of opposite sign.
#'
#'   Jumps (spikes) must satisfy two conditions:
#'   \enumerate{
#'     \item Neighboring returns both exceed a multiple of the rolling quantile,
#'     \item The sum of neighboring returns doesn't exceed that multiple.
#'   }
#'
#'   The quantile multiplier \code{vol_mult} controls the threshold at which
#'   values are identified as jumps. Smaller quantile multiplier values will
#'   cause more values to be identified as jumps.
#'
#' @examples
#' # Create local copy of SPY TAQ data
#' ta_q <- SPY_TAQ
#' # Calculate mid prices
#' mid_prices <- 0.5 * (ta_q[, "Bid.Price"] + ta_q[, "Ask.Price"])
#' # Replace whole rows containing suspect price jumps with NA, and perform locf()
#' ta_q[which_jumps(mid_prices, look_back=31, vol_mult=1.0), ] <- NA
#' ta_q <- xts:::na.locf.xts(ta_q)

which_jumps <- function(x_ts, look_back=51, vol_mult=2) {
# Calculate simple returns
  re_turns <- rutils::diff_it(as.numeric(x_ts))
#  re_turns[1] <- 0
#  colnames(re_turns) <- "diffs"
  rets_advanced <- rutils::lag_it(re_turns, -1)
#  rets_advanced[NROW(rets_advanced)] <- 0
#  colnames(rets_advanced) <- "rets_advanced"

# Calculate volatility as the rolling quantile of returns
  quan_tile <- caTools::runquantile(x=abs(re_turns), k=look_back,
                        probs=0.9, endrule="constant", align="center")
#  quan_tile <- xts(quan_tile, order.by=index(re_turns))
#  colnames(quan_tile) <- "volat"
# Carry forward non-zero quan_tile values
  quan_tile[quan_tile==0] <- NA
  quan_tile[1] <- 1
  quan_tile <- rutils::na_locf(quan_tile)
#  quan_tile <- rutils::na_locf(quan_tile, fromLast=TRUE)

# value is suspect if abs re_turns greater than quan_tile,
# and if abs sum of re_turns less than quan_tile
  sus_pect <- ((abs(re_turns) > vol_mult*quan_tile) &
      (abs(rets_advanced) > vol_mult*quan_tile) &
      (abs(re_turns+rets_advanced) < 2*vol_mult*quan_tile))
  sus_pect[1] <- FALSE
#  colnames(sus_pect) <- "suspect"
# Cat("Parsing", deparse(substitute(ta_q)), "\n")
# Cat("Parsing", strsplit(deparse(substitute(ta_q)), split="[.]")[[1]][4], "on date:", format(to_day), "\tfound", sum(sus_pect), "suspect prices\n")
  cat("date:", format(as.Date(index(first(x_ts)))), "\tfound", sum(sus_pect), "jump prices\n")
  sus_pect
}  # end which_jumps




##########################################################################
#' Scrub a single day of \emph{TAQ} data in \emph{xts} format, without
#' aggregation.
#'
#' @export
#' @inheritParams which_extreme
#' 
#' @param \code{ta_q} \emph{TAQ} A time series in \emph{xts} format.
#' 
#' @param \code{tzone} The timezone to convert.
#'
#' @return A \emph{TAQ} time series in \emph{xts} format.
#'
#' @details The function \code{scrub_taq()} performs the same scrubbing
#'   operations as \code{scrub_agg}, except it doesn't aggregate, and returns
#'   the \emph{TAQ} data in \emph{xts} format.
#'
#' @examples
# scrub a single day of TAQ data without aggregating it
#' ta_q <- HighFreq::scrub_taq(ta_q=HighFreq::SPY_TAQ, look_back=11, vol_mult=1)
#' # Create random TAQ prices and scrub them
#' ta_q <- HighFreq::random_taq()
#' ta_q <- HighFreq::scrub_taq(ta_q=ta_q)
#' ta_q <- HighFreq::scrub_taq(ta_q=ta_q, look_back=11, vol_mult=1)

scrub_taq <- function(ta_q, look_back=51, vol_mult=2, tzone="America/New_York") {
# Convert timezone of index to New_York
  index(ta_q) <- lubridate::with_tz(time=index(ta_q), tzone=tzone)
# subset data to NYSE trading hours
  ta_q <- ta_q["T09:30:00/T16:00:00", ]
# Return NULL if no data
  if (NROW(ta_q)==0)  return(NULL)
#  to_day <- as.Date(index(first(ta_q)))

# Remove duplicate time stamps using duplicated()
  ta_q <- ta_q[!duplicated(index(ta_q)), ]

# scrub quotes with suspect bid-offer spreads
  bid_offer <- ta_q[, "Ask.Price"] - ta_q[, "Bid.Price"]
#  bid_offer <- na.omit(bid_offer)
  sus_pect <- which_extreme(bid_offer, look_back=look_back, vol_mult=vol_mult)
# Remove suspect values
  ta_q <- ta_q[!sus_pect]
# Replace suspect values
# ta_q[sus_pect, "Bid.Price"] <- ta_q[sus_pect, "Trade.Price"]
# ta_q[sus_pect, "Ask.Price"] <- ta_q[sus_pect, "Trade.Price"]

# scrub quotes with suspect price jumps
# Calculate mid prices
  mid_prices <- 0.5 * (ta_q[, "Bid.Price"] + ta_q[, "Ask.Price"])
#  mid_prices <- na.omit(mid_prices)
#  colnames(mid_prices) <- "Mid.Price"
# Replace NA volumes with zero
  ta_q[is.na(ta_q[, "Volume"]), "Volume"] <- 0
# Replace whole rows containing suspect price jumps with NA, and perform locf()
  ta_q[which_jumps(mid_prices, look_back=look_back, vol_mult=vol_mult), ] <- NA
  rutils::na_locf(ta_q)
}  # end scrub_taq




##########################################################################
#' Scrub a single day of \emph{TAQ} data, aggregate it, and convert to
#' \emph{OHLC} format.
#'
#' @export
#' @inheritParams scrub_taq
#' @param \code{period} The aggregation period.
#'
#' @return A \emph{OHLC} time series in \emph{xts} format.
#'
#' @details The function \code{scrub_agg()} performs:
#' \itemize{
#'   \item index timezone conversion,
#'   \item data subset to trading hours,
#'   \item removal of duplicate time stamps,
#'   \item scrubbing of quotes with suspect bid-offer spreads,
#'   \item scrubbing of quotes with suspect price jumps,
#'   \item cbinding of mid prices with volume data,
#'   \item aggregation to OHLC using function \code{to.period()} from package \emph{xts},
#' }
#' Valid 'period' character strings include: "minutes", "3 min", "5 min", "10
#' min", "15 min", "30 min", and "hours". The time index of the output time
#' series is rounded up to the next integer multiple of 'period'.
#'
#' @examples
#' # Create random TAQ prices
#' ta_q <- HighFreq::random_taq()
#' # Aggregate to ten minutes OHLC data
#' oh_lc <- HighFreq::scrub_agg(ta_q, period="10 min")
#' chart_Series(oh_lc, name="random prices")
#' # scrub and aggregate a single day of SPY TAQ data to OHLC
#' oh_lc <- HighFreq::scrub_agg(ta_q=HighFreq::SPY_TAQ)
#' chart_Series(oh_lc, name=sym_bol)

scrub_agg <- function(ta_q, look_back=51, vol_mult=2,
                      period="minutes", tzone="America/New_York") {
# Convert timezone of index to New_York
  index(ta_q) <- lubridate::with_tz(time=index(ta_q), tzone=tzone)
# subset data to NYSE trading hours
  ta_q <- ta_q["T09:30:00/T16:00:00", ]
# Return NULL if no data
  if (NROW(ta_q)==0)  return(NULL)
#  to_day <- as.Date(index(first(ta_q)))

# Remove duplicate time stamps using duplicated()
  ta_q <- ta_q[!duplicated(index(ta_q)), ]

# scrub quotes with suspect bid-offer spreads
  bid_offer <- ta_q[, "Ask.Price"] - ta_q[, "Bid.Price"]
#  bid_offer <- na.omit(bid_offer)
  sus_pect <- which_extreme(bid_offer, look_back=look_back, vol_mult=vol_mult)
# Remove suspect values
  ta_q <- ta_q[!sus_pect]
# Replace suspect values
# ta_q[sus_pect, "Bid.Price"] <- ta_q[sus_pect, "Trade.Price"]
# ta_q[sus_pect, "Ask.Price"] <- ta_q[sus_pect, "Trade.Price"]

# scrub quotes with suspect price jumps
# Calculate mid prices
  mid_prices <- 0.5 * (ta_q[, "Bid.Price"] + ta_q[, "Ask.Price"])
#  mid_prices <- na.omit(mid_prices)
  colnames(mid_prices) <- "Mid.Price"
# Replace whole rows containing suspect price jumps with NA, and perform locf()
  mid_prices[which_jumps(mid_prices, look_back=look_back, vol_mult=vol_mult)] <- NA
  mid_prices <- rutils::na_locf(mid_prices)
#  mid_prices <- rutils::na_locf(mid_prices, fromLast=TRUE)
# Cbind mid_prices with volume data, and replace NA volumes with zero
  mid_prices <- cbind(mid_prices, ta_q[index(mid_prices), "Volume"])
  mid_prices[is.na(mid_prices[, "Volume"]), "Volume"] <- 0

# Aggregate to OHLC and cumulative volume data
  mid_prices <- switch(period,
                       "minutes"={sec_onds <- 60; to.period(x=mid_prices, period=period)},
                       "3 min"={sec_onds <- 3*60; to.minutes3(x=mid_prices)},
                       "5 min"={sec_onds <- 5*60; to.minutes5(x=mid_prices)},
                       "10 min"={sec_onds <- 10*60; to.minutes10(x=mid_prices)},
                       "15 min"={sec_onds <- 15*60; to.minutes15(x=mid_prices)},
                       "30 min"={sec_onds <- 30*60; to.minutes30(x=mid_prices)},
                       "hours"={sec_onds <- 60*60; to.period(x=mid_prices, period=period)}
                       )  # end switch
# round up times to next period
  index(mid_prices) <- align.time(x=index(mid_prices), n=sec_onds)
  mid_prices
}  # end scrub_agg




##########################################################################
#' Load, scrub, aggregate, and rbind multiple days of \emph{TAQ} data for a
#' single symbol, and save the \emph{OHLC} time series to a single
#' \sQuote{\code{*.RData}} file.
#'
#' @export
#' @param \code{sym_bol} A \emph{character} string representing symbol or ticker.
#' @param \code{data_dir} A \emph{character} string representing directory containing
#'   input \sQuote{\code{*.RData}} files.
#' @param \code{output_dir} A \emph{character} string representing directory containing
#'   output \sQuote{\code{*.RData}} files.
#' @inheritParams scrub_agg
#'
#' @return An \emph{OHLC} time series in \emph{xts} format.
#'
#' @details The function \code{save_scrub_agg()} loads multiple days of
#'   \emph{TAQ} data, then scrubs, aggregates, and rbinds them into a
#'   \emph{OHLC} time series, and finally saves it to a single
#'   \sQuote{\code{*.RData}} file. The \emph{OHLC} time series is stored in a
#'   variable named \sQuote{\code{symbol}}, and then it's saved to a file named
#'   \sQuote{\code{symbol.RData}} in the \sQuote{\code{output_dir}} directory.
#'   The \emph{TAQ} data files are assumed to be stored in separate directories
#'   for each \sQuote{\code{symbol}}. Each \sQuote{\code{symbol}} has its own
#'   directory (named \sQuote{\code{symbol}}) in the \sQuote{\code{data_dir}}
#'   directory. Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \emph{TAQ}
#'   data.
#' @examples
#' \dontrun{
#' # set data directories
#' data_dir <- "C:/Develop/data/hfreq/src/"
#' output_dir <- "C:/Develop/data/hfreq/scrub/"
#' sym_bol <- "SPY"
#' # Aggregate SPY TAQ data to 15-min OHLC bar data, and save the data to a file
#' save_scrub_agg(sym_bol=sym_bol, data_dir=data_dir, output_dir=output_dir, period="15 min")
#' }

save_scrub_agg <- function(sym_bol,
                      data_dir="E:/mktdata/sec/",
                      output_dir="E:/output/data/",
                      look_back=51,
                      vol_mult=2,
                      period="minutes",
                      tzone="America/New_York") {
# Create path to directory containing *.RData files
  file_dir <- file.path(data_dir, sym_bol)
# get list of *.RData files
  file_list <- list.files(file_dir)
# Create paths to *.RData files
  file_names <- file.path(file_dir, file_list)

# load TAQ data one by one, scrub and aggregate it, return list of xts
da_ta <- lapply(file_names, function(file_name) {
  cat("loading", sym_bol, "from file: ", file_name, "\n")
  sym_bol <- load(file_name)
  scrub_agg(get(sym_bol),
            look_back=look_back,
            vol_mult=vol_mult,
            period=period, tzone=tzone)
})  # end sapply

# Recursively "rbind" the list into a single xts
  da_ta <- rutils::do_call_rbind(da_ta)
# assign column names, i.e. "symbol.High"
  colnames(da_ta) <- sapply(strsplit(colnames(da_ta), split="[.]"),
                           function(strng) paste(sym_bol, strng[-1], sep="."))

# Copy the xts data to a variable with the name 'sym_bol'
  assign(sym_bol, da_ta)

# save the xts data to a file in the output_dir
  save(list=sym_bol, file=file.path(output_dir, paste0(sym_bol, ".RData")))
  invisible(sym_bol)

}  # end save_scrub_agg




##########################################################################
#' Load and scrub multiple days of \emph{TAQ} data for a single symbol, and save
#' it to multiple \sQuote{\code{*.RData}} files.
#'
#' @export
#' @inheritParams save_scrub_agg
#'
#' @return a \emph{TAQ} time series in \emph{xts} format.
#'
#' @details The function \code{save_taq()} loads multiple days of \emph{TAQ}
#'   data, scrubs it, and saves the scrubbed TAQ data to individual
#'   \sQuote{\code{*.RData}} files. It uses the same file names for output as
#'   the input file names. The \emph{TAQ} data files are assumed to be stored in
#'   separate directories for each \sQuote{\code{symbol}}. Each
#'   \sQuote{\code{symbol}} has its own directory (named \sQuote{\code{symbol}})
#'   in the \sQuote{\code{data_dir}} directory.
#'   Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \emph{TAQ}
#'   data.
#'
#' @examples
#' \dontrun{
#' save_taq("SPY")
#' }

save_taq <- function(sym_bol,
                      data_dir="E:/mktdata/sec/",
                      output_dir="E:/output/data/",
                      look_back=51,
                      vol_mult=2,
                      tzone="America/New_York") {
# Create path to directory containing *.RData files
  data_dir <- file.path(data_dir, sym_bol)
# get list of *.RData files
  file_names <- list.files(data_dir)
# Create path to directory for writing *.RData files
  output_dir <- file.path(output_dir, sym_bol)

# load TAQ data one-by-one, scrub, and save
  dummy_data <- sapply(file_names, function(file_name) {
    cat("loading", sym_bol, "from file: ", file_name, "\n")
    file_name_in <- file.path(data_dir, file_name)
    sym_bol <- load(file_name_in)
    file_name_out <- file.path(output_dir, file_name)
# save the xts data to a file in the output_dir
    ta_q <- scrub_taq(get(sym_bol), look_back=look_back, vol_mult=vol_mult, tzone=tzone)
    if (!is.null(ta_q)) {
      assign(sym_bol, ta_q)
      save(list=sym_bol, file=file_name_out)
      cat("finished saving", sym_bol, "to file: ", file_name, "\n")
    }
    file_name
  })  # end sapply

  invisible(sym_bol)

}  # end save_taq




##########################################################################
#' Load, scrub, aggregate, and rbind multiple days of \emph{TAQ} data for a
#' single symbol. Calculate returns and save them to a single
#' \sQuote{\code{*.RData}} file.
#'
#' @export
#' @inheritParams save_scrub_agg
#'
#' @return A time series of returns and volume in \emph{xts} format.
#'
#' @details The function \code{save_rets} loads multiple days of \emph{TAQ}
#'   data, then scrubs, aggregates, and rbinds them into a \emph{OHLC} time
#'   series.  It then calculates returns using function \code{ohlc_returns()}, and
#'   stores them in a variable named \sQuote{\code{symbol.rets}}, and saves them
#'   to a file called \sQuote{\code{symbol.rets.RData}}.
#'   The \emph{TAQ} data files are assumed to be stored in separate directories
#'   for each \sQuote{\code{symbol}}. Each \sQuote{\code{symbol}} has its own
#'   directory (named \sQuote{\code{symbol}}) in the \sQuote{\code{data_dir}}
#'   directory. Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \emph{TAQ}
#'   data.
#'
#' @examples
#' \dontrun{
#' save_rets("SPY")
#' }

save_rets <- function(sym_bol,
                      data_dir="E:/mktdata/sec/",
                      output_dir="E:/output/data/",
                      look_back=51,
                      vol_mult=2,
                      period="minutes",
                      tzone="America/New_York") {
# Create path to directory containing *.RData files
  file_dir <- file.path(data_dir, sym_bol)
# get list of *.RData files
  file_list <- list.files(file_dir)
# Create paths to *.RData files
  file_names <- file.path(file_dir, file_list)

# load TAQ data into list
  ta_q <- sapply(file_names, function(file_name) {
    cat("loading", sym_bol, "from file: ", file_name, "\n")
    sym_bol <- load(file_name)
    get(sym_bol)
  })

# scrub and aggregate the TAQ data
  oh_lc <- lapply(ta_q, scrub_agg,
                      look_back=look_back,
                      vol_mult=vol_mult,
                      period=period,
                      tzone=tzone)

# Calculate returns
  oh_lc <- lapply(oh_lc, ohlc_returns)

# Recursively "rbind" the list into a single xts
  oh_lc <- rutils::do_call_rbind(oh_lc)
# assign column names, i.e. "symbol.rets"
  colnames(oh_lc) <-
    c(paste(sym_bol, "rets", sep="."), paste(sym_bol, "vol", sep="."))

# Copy the xts data to a variable with the name 'sym_bol'
  sym_bol_rets <- paste(sym_bol, "rets", sep=".")
  assign(sym_bol_rets, oh_lc)

# save the xts data to a file in the output_dir
  save(list=eval(sym_bol_rets),
       file=file.path(output_dir, paste0(sym_bol_rets, ".RData")))
  invisible(sym_bol_rets)

}  # end save_rets




##########################################################################
#' Load \emph{OHLC} time series data for a single symbol, calculate its returns,
#' and save them to a single \sQuote{\code{*.RData}} file, without aggregation.
#'
#' @export
#' @inheritParams save_scrub_agg
#'
#' @return A time series of returns and volume in \emph{xts} format.
#'
#' @details The function \code{save_rets_ohlc()} loads \emph{OHLC} time series
#'   data from a single file.  It then calculates returns using function
#'   \code{ohlc_returns()}, and stores them in a variable named
#'   \sQuote{\code{symbol.rets}}, and saves them to a file called
#'   \sQuote{\code{symbol.rets.RData}}.
#'
#' @examples
#' \dontrun{
#' save_rets_ohlc("SPY")
#' }

save_rets_ohlc <- function(sym_bol,
                      data_dir="E:/output/data/",
                      output_dir="E:/output/data/") {
# Create path to directory containing sym_bol.RData file
  file_name <- file.path(data_dir, paste0(sym_bol, ".RData"))
# load OHLC data
  cat("loading", sym_bol, "from file: ", file_name, "\n")
  sym_bol <- load(file_name)

# Calculate returns
  da_ta <- ohlc_returns(get(sym_bol))

# Copy the xts data to a variable with the name 'sym_bol'
  sym_bol_rets <- paste(sym_bol, "rets", sep=".")
  assign(sym_bol_rets, da_ta)

# save the xts data to a file in the output_dir
  cat("saving", sym_bol, "to file: ", paste0(sym_bol_rets, ".RData"), "\n")
  save(list=eval(sym_bol_rets), file=file.path(output_dir, paste0(sym_bol_rets, ".RData")))
  invisible(sym_bol_rets)

}  # end save_rets_ohlc




##########################################################################
#' Calculate the Value at Risk (\emph{VaR}) or the Conditional Value at Risk
#' (\emph{CVaR}) of an \emph{xts} \emph{time series} of returns, using \code{R}
#' code.
#' 
#' @param \code{tseries} An \emph{xts} \emph{time series} of returns with
#'   multiple columns.
#'   
#' @param \code{method} A \emph{string} specifying the type of risk measure
#'   (the default is \code{method = "var"} - see Details).
#'    
#' @param \code{con_fi} The confidence level for calculating the
#'   quantile (the default is \code{con_fi = pnorm(-2) = 0.02275}).
#'
#' @return A vector with the risk measures of the columns of the input
#'   \emph{time series} \code{tseries}.
#'
#' @details 
#'   The function \code{calc_cvar()} calculates the Value at Risk (\emph{VaR})
#'   or the Conditional Value at Risk (\emph{CVaR}) of an \emph{xts} \emph{time
#'   series} of returns, using \code{R}
#'   
#'   The Value at Risk (\emph{VaR}) and the Conditional Value at Risk
#'   (\emph{CVaR}) are measures of the tail risk of returns.
#'
#'   If \code{method = "var"} then \code{calc_cvar()} calculates the Value at
#'   Risk (\emph{VaR}) as the quantile of the returns as follows:
#'   \deqn{
#'     \alpha = \int_{-\infty}^{\mathrm{VaR}(\alpha)} \mathrm{f}(r) \, \mathrm{d}r
#'   }
#'   Where \eqn{\alpha} is the confidence level for calculating the quantile,
#'   and \eqn{\mathrm{f}(r)} is the probability density (distribution) of
#'   returns.
#'   
#'   If \code{method = "cvar"} then \code{calc_cvar()} calculates the Value at
#'   Risk (\emph{VaR}) as the Expected Tail Loss (\emph{ETL}) of the returns as
#'   follows:
#'   \deqn{
#'     \mathrm{CVaR} = \frac{1}{\alpha} \int_{0}^\alpha \mathrm{VaR}(p) \, \mathrm{d}p
#'   }
#'   Where \eqn{\alpha} is the confidence level for calculating the quantile.
#'   
#'   
#' @examples
#' \dontrun{
#' # Calculate VTI and XLF returns
#' re_turns <- na.omit(rutils::etf_env$re_turns[, c("VTI", "XLF")])
#' # Calculate VaR
#' all.equal(HighFreq::calc_cvar(re_turns), 
#'   sapply(re_turns, quantile, probs=pnorm(-2)), check.attributes=FALSE)
#' # Calculate CVaR
#' all.equal(HighFreq::calc_cvar(re_turns, method="cvar", con_fi=0.02), 
#'   sapply(re_turns, function(x) mean(x[x < quantile(x, 0.02)])), 
#'   check.attributes=FALSE)
#' }
#' 
#' @export
calc_cvar <- function(tseries, method = "var", con_fi = pnorm(-2)) {

  # Switch for the different risk methods
  risk <- switch(method,
                 "var"={sapply(tseries, quantile, probs=con_fi)},
                 # Calculate CVaR as expected loss
                 "cvar"={sapply(tseries, function(x) mean(x[x < quantile(x, con_fi)]))}
  )  # end switch
  
  risk
  
}  # end calc_cvar




##########################################################################
#' Calculate a time series of point estimates of variance for an \emph{OHLC}
#' time series, using different range estimators for variance.
#'
#' Calculates the point variance estimates from individual rows of \emph{OHLC}
#' prices (rows of data), using the squared differences of \emph{OHLC} prices at
#' each point in time, without averaging them over time.
#'
#' @export
#' @param \code{oh_lc} An \emph{OHLC} time series of prices in \emph{xts} format.
#' 
#' @param \code{method} A \emph{character} string representing the method for
#'   estimating variance.  The methods include:
#'   \itemize{
#'     \item "close" close to close,
#'     \item "garman_klass" Garman-Klass,
#'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
#'     \item "rogers_satchell" Rogers-Satchell,
#'     \item "yang_zhang" Yang-Zhang,
#'    }
#'    (default is \code{"yang_zhang"})
#' @param \code{scal_e} \emph{Boolean} argument: should the returns be divided by the
#'   number of seconds in each period? (default is \code{TRUE})
#'
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#'
#' @details The function \code{ohlc_variance()} calculates a time series of point
#'   variance estimates of percentage returns, from \emph{OHLC} prices, without
#'   averaging them over time. For example, the method \code{"close"} simply
#'   calculates the squares of the differences of the log \emph{Close} prices.
#'
#'   The other methods calculate the squares of other possible differences of
#'   the log \emph{OHLC} prices.  This way the point variance estimates only
#'   depend on the price differences within individual rows of data (and
#'   possibly from the neighboring rows.)
#'   All the methods are implemented assuming zero drift, since the calculations
#'   are performed only for a single row of data, at a single point in time.
#'
#'   The user can choose from several different variance estimation methods. The
#'   methods \code{"close"}, \code{"garman_klass_yz"}, and \code{"yang_zhang"}
#'   do account for close-to-open price jumps, while the methods
#'   \code{"garman_klass"} and \code{"rogers_satchell"} do not account for
#'   close-to-open price jumps. The default method is \code{"yang_zhang"}, which
#'   theoretically has the lowest standard error among unbiased estimators.
#'
#'   The point variance estimates can be passed into function \code{roll_vwap()}
#'   to perform averaging, to calculate rolling variance estimates.  This is
#'   appropriate only for the methods \code{"garman_klass"} and
#'   \code{"rogers_satchell"}, since they don't require subtracting the rolling
#'   mean from the point variance estimates.
#'
#'   The point variance estimates can also be considered to be technical
#'   indicators, and can be used as inputs into trading models.
#'
#'   If \code{scal_e} is \code{TRUE} (the default), then the variance is divided
#'   by the squared differences of the time index (which scales the variance to
#'   units of variance per second squared.) This is useful for example, when
#'   calculating intra-day variance from minutely bar data, because dividing
#'   returns by the number of seconds decreases the effect of overnight price
#'   jumps.
#'
#'   If \code{scal_e} is \code{TRUE} (the default), then the variance is
#'   expressed in the scale of the time index of the \emph{OHLC} time series.
#'   For example, if the time index is in seconds, then the variance is given in
#'   units of variance per second squared.  If the time index is in days, then
#'   the variance is equal to the variance per day squared.
#'
#'   The time index of the \code{oh_lc} time series is assumed to be in
#'   \emph{POSIXct} format, so that its internal value is equal to the number of
#'   seconds that have elapsed since the \emph{epoch}.
#'
#'   The function \code{ohlc_variance()} performs similar calculations to the
#'   function \code{volatility()} from package
#'   \href{https://cran.r-project.org/web/packages/TTR/index.html}{TTR}, but it
#'   assumes zero drift, and doesn't calculate a running sum using
#'   \code{runSum()}.  It's also a little faster because it performs less data
#'   validation.
#'
#' @examples
#' # Create minutely OHLC time series of random prices
#' oh_lc <- HighFreq::random_ohlc()
#' # Calculate variance estimates for oh_lc
#' var_running <- HighFreq::ohlc_variance(oh_lc)
#' # Calculate variance estimates for SPY
#' var_running <- HighFreq::ohlc_variance(HighFreq::SPY, method="yang_zhang")
#' # Calculate SPY variance without overnight jumps
#' var_running <- HighFreq::ohlc_variance(HighFreq::SPY, method="rogers_satchell")

ohlc_variance <- function(oh_lc, method="yang_zhang", scal_e=TRUE) {
  sym_bol <- rutils::get_name(colnames(oh_lc)[1])
  # oh_lc <- log(oh_lc[, 1:4])
  vari_ance <- switch(method,
         "close"={rutils::diff_it(oh_lc[, 4])^2},
         "garman_klass"={0.5*(oh_lc[, 2]-oh_lc[, 3])^2 -
                         (2*log(2)-1)*(oh_lc[, 4]-oh_lc[, 1])^2},
         "rogers_satchell"={(oh_lc[, 2]-oh_lc[, 4])*(oh_lc[, 2]-oh_lc[, 1]) +
                            (oh_lc[, 3]-oh_lc[, 4])*(oh_lc[, 3]-oh_lc[, 1])},
         "garman_klass_yz"={(oh_lc[, 1]-rutils::lag_it(oh_lc[, 4]))^2 +
                            0.5*(oh_lc[, 2]-oh_lc[, 3])^2 -
                            (2*log(2)-1)*(oh_lc[, 4]-oh_lc[, 1])^2},
         "yang_zhang"={co_eff <- 0.34/2.34
                            (oh_lc[, 1]-rutils::lag_it(oh_lc[, 4]))^2 +
                            co_eff*(oh_lc[, 1]-oh_lc[, 4])^2 +
                            (1-co_eff)*((oh_lc[, 2]-oh_lc[, 4])*(oh_lc[, 2]-oh_lc[, 1]) +
                               (oh_lc[, 3]-oh_lc[, 4])*(oh_lc[, 3]-oh_lc[, 1]))}
  )  # end switch
  if (scal_e)
    vari_ance <- vari_ance/c(1, diff(xts::.index(oh_lc)))^2
  vari_ance[1, ] <- 0
  vari_ance <- rutils::na_locf(vari_ance)
  # Colnames(vari_ance) <- paste0(sym_bol, ".Variance")
  vari_ance
}  # end ohlc_variance




##########################################################################
#' Calculate time series of point skew estimates from a \emph{OHLC} time series,
#' assuming zero drift.
#'
#' @export
#' @param \code{oh_lc} An \emph{OHLC} time series of prices in \emph{xts} format.
#' 
#' @param \code{method} A \emph{character} string representing method for
#'   estimating skew.
#'
#' @return A time series of point skew estimates.
#'
#' @details The function \code{ohlc_skew()} calculates a time series of skew
#'   estimates from \emph{OHLC} prices, one for each row of \emph{OHLC} data.
#'   The skew estimates are expressed in the time scale of the index of the
#'   \emph{OHLC} time series.
#'   For example, if the time index is in seconds, then the skew is given in
#'   units of skew per second.  If the time index is in days, then the skew is
#'   equal to the skew per day.
#'
#'   Currently only the \code{"close"} skew estimation method is correct
#'   (assuming zero drift), while the \code{"rogers_satchell"} method produces a
#'   skew-like indicator, proportional to the skew. The default method is
#'   \code{"rogers_satchell"}.
#'
#' @examples
#' # Calculate time series of skew estimates for SPY
#' sk_ew <- HighFreq::ohlc_skew(HighFreq::SPY)

ohlc_skew <- function(oh_lc, method="rogers_satchell") {
  sym_bol <- rutils::get_name(colnames(oh_lc)[1])
  # oh_lc <- log(oh_lc[, 1:4])
  sk_ew <- switch(method,
                  "close"={rutils::diff_it(oh_lc[, 4])^3},
                  "garman_klass"={0.5*(oh_lc[, 2]-oh_lc[, 3])^3 -
                      (2*log(2)-1)*(oh_lc[, 4]-oh_lc[, 1])^3},
                  "rogers_satchell"={
                    (oh_lc[, 2]-oh_lc[, 4])*(oh_lc[, 2]-oh_lc[, 1])*(oh_lc[, 2]-0.5*(oh_lc[, 4] + oh_lc[, 1])) +
                      (oh_lc[, 3]-oh_lc[, 4])*(oh_lc[, 3]-oh_lc[, 1])*(oh_lc[, 3]-0.5*(oh_lc[, 4] + oh_lc[, 1]))},
                  "garman_klass_yz"={(oh_lc[, 1]-rutils::lag_it(oh_lc[, 4]))^3 +
                      0.5*(oh_lc[, 2]-oh_lc[, 3])^3 -
                      (2*log(2)-1)*(oh_lc[, 4]-oh_lc[, 1])^3},
                  "yang_zhang"={c_o <- oh_lc[, 1]-rutils::lag_it(oh_lc[, 4]);
                  o_c <- oh_lc[, 1]-oh_lc[, 4];
                  (c_o-sum(c_o)/NROW(c_o))^3 +
                    0.67*(o_c-sum(o_c)/NROW(o_c))^3 +
                    0.33*((oh_lc[, 2]-oh_lc[, 4])*(oh_lc[, 2]-oh_lc[, 1]) +
                            (oh_lc[, 3]-oh_lc[, 4])*(oh_lc[, 3]-oh_lc[, 1]))}
  )  # end switch
  sk_ew <- sk_ew/c(1, diff(xts::.index(oh_lc)))^3
  sk_ew[1, ] <- 0
  sk_ew <- rutils::na_locf(sk_ew)
  # Colnames(sk_ew) <- paste0(sym_bol, ".Skew")
  sk_ew
}  # end ohlc_skew




##########################################################################
#' Calculate time series of point Sharpe-like statistics for each row of a
#' \emph{OHLC} time series.
#'
#' @export
#' @param \code{oh_lc} An \emph{OHLC} time series of prices in \emph{xts} format.
#' 
#' @param \code{method} A \emph{character} string representing method for
#'   estimating the Sharpe-like exponent.
#'
#' @return An \emph{xts} time series with the same number of rows as the
#'   argument \code{oh_lc}.
#'
#' @details The function \code{ohlc_sharpe()} calculates Sharpe-like statistics
#'   for each row of a \emph{OHLC} time series.
#'   The Sharpe-like statistic is defined as the ratio of the difference between
#'   \emph{Close} minus \emph{Open} prices divided by the difference between
#'   \emph{High} minus \emph{Low} prices.
#'   This statistic may also be interpreted as something like a \emph{Hurst
#'   exponent} for a single row of data.
#'   The motivation for the Sharpe-like statistic is the notion that if prices
#'   are trending in the same direction inside a given time bar of data, then
#'   this statistic is close to either 1 or -1.
#'
#' @examples
#' # Calculate time series of running Sharpe ratios for SPY
#' sharpe_running <- ohlc_sharpe(HighFreq::SPY)

ohlc_sharpe <- function(oh_lc, method="close") {
  sharpe_ratio <- switch(method,
                   "close"={(oh_lc[, 4]-oh_lc[, 1])/(oh_lc[, 2]-oh_lc[, 3])},
                   "method2"={(oh_lc[, 4]-oh_lc[, 1])/(oh_lc[, 2]-oh_lc[, 3])}
  )  # end switch
  sharpe_ratio <- ifelse(oh_lc[, 2]==oh_lc[, 3], 0, sharpe_ratio)
  # Colnames(sharpe_ratio) <- paste0(rutils::get_name(colnames(oh_lc)[1]), ".Sharpe")
  sharpe_ratio
}  # end ohlc_sharpe




##########################################################################
#' Calculate the aggregation (weighted average) of a statistical estimator over
#' a \emph{OHLC} time series using \code{R} code.
#'
#' @export
#' @param \code{oh_lc} An \emph{OHLC} time series of prices and trading volumes,
#'   in \emph{xts} format.
#'
#' @param \code{calc_bars} A \emph{character} string representing a function
#'   for calculating statistics for individual \emph{OHLC} bars.
#'
#' @param \code{weight_ed} \emph{Boolean} argument: should estimate be weighted
#'   by the trading volume? (default is \code{TRUE})
#'
#' @param ... additional parameters to the function \code{calc_bars}.
#'
#' @return A single \emph{numeric} value equal to the volume weighted average of
#'   an estimator over the time series.
#'
#' @details The function \code{agg_stats_r()} calculates a single number
#'   representing the volume weighted average of statistics of individual
#'   \emph{OHLC} bars.
#'   It first calls the function \code{calc_bars} to calculate a vector of
#'   statistics for the \emph{OHLC} bars.
#'   For example, the statistic may simply be the difference between the
#'   \emph{High} minus \emph{Low} prices.  In this case the function
#'   \code{calc_bars} would calculate a vector of \emph{High} minus
#'   \emph{Low} prices.
#'   The function \code{agg_stats_r()} then calculates a trade volume
#'   weighted average of the vector of statistics.
#'   
#'   The function \code{agg_stats_r()} is implemented in \code{R} code.
#'
#' @examples
#' # Calculate weighted average variance for SPY (single number)
#' vari_ance <- agg_stats_r(oh_lc=HighFreq::SPY, calc_bars="ohlc_variance")
#' # Calculate time series of daily skew estimates for SPY
#' skew_daily <- apply.daily(x=HighFreq::SPY, FUN=agg_stats_r, calc_bars="ohlc_skew")

agg_stats_r <- function(oh_lc, calc_bars="ohlc_variance", weight_ed=TRUE, ...) {
  
# Match "calc_bars" with moment function
  calc_bars <- match.fun(calc_bars)
  agg_regations <- calc_bars(oh_lc, ...)
  
# Weight the estimates by volume
  if (weight_ed) {
    agg_regations <- oh_lc[, 5]*agg_regations
    agg_regations <- sum(agg_regations)/sum(oh_lc[, 5])
  } else
    agg_regations <- sum(agg_regations)/NROW(agg_regations)
  
  agg_regations
  
}  # end agg_stats_r




##########################################################################
#' Calculate the volume-weighted average price of an \emph{OHLC} time series
#' over a rolling look-back interval.
#'
#' Performs the same operation as function \code{VWAP()} from package
#' \href{https://cran.r-project.org/web/packages/TTR/index.html}{TTR},
#' but using vectorized functions, so it's a little faster.
#'
#' @export
#' @param \code{oh_lc} An \emph{OHLC} time series of prices in \emph{xts} format.
#' 
#' @param \code{close} A time series of close prices.
#' 
#' @param \code{look_back} The size of the look-back interval, equal to the number of 
#'   rows of data used for calculating the average price.
#'
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#'
#' @details The function \code{roll_vwap()} calculates the volume-weighted
#'   average closing price, defined as the sum of the prices multiplied by
#'   trading volumes in the look-back interval, divided by the sum of trading
#'   volumes in the interval. If the argument \code{close} is passed in explicitly,
#'   then its volume-weighted average value over time is calculated.
#'
#' @examples
#' # Calculate and plot rolling volume-weighted average closing prices (VWAP)
#' prices_rolling <- roll_vwap(oh_lc=HighFreq::SPY["2013-11"], look_back=11)
#' chart_Series(HighFreq::SPY["2013-11-12"], name="SPY prices")
#' add_TA(prices_rolling["2013-11-12"], on=1, col="red", lwd=2)
#' legend("top", legend=c("SPY prices", "VWAP prices"),
#' bg="white", lty=c(1, 1), lwd=c(2, 2),
#' col=c("black", "red"), bty="n")
#' # Calculate running returns
#' returns_running <- ohlc_returns(x_ts=HighFreq::SPY)
#' # Calculate the rolling volume-weighted average returns
#' roll_vwap(oh_lc=HighFreq::SPY, close=returns_running, look_back=11)

roll_vwap <- function(oh_lc, close=oh_lc[, 4, drop=FALSE], look_back) {
  roll_vwap <- rutils::roll_sum(x_ts=close*oh_lc[, 5, drop=FALSE], look_back=look_back)
  volume_rolling <- rutils::roll_sum(x_ts=oh_lc[, 5, drop=FALSE], look_back=look_back)
  # roll_vwap <- HighFreq::roll_sum(tseries=close*oh_lc[, 5, drop=FALSE], look_back=look_back)
  # volume_rolling <- HighFreq::roll_sum(tseries=oh_lc[, 5, drop=FALSE], look_back=look_back)
  roll_vwap <- ifelse(volume_rolling > 0, roll_vwap/volume_rolling, 0)
  # roll_vwap[is.na(roll_vwap)] <- 0
  # Colnames(roll_vwap) <- paste0(rutils::get_name(colnames(oh_lc)[1]), ".VWAP")
  # Colnames(roll_vwap) <- colnames(oh_lc)
  roll_vwap
}  # end roll_vwap




##########################################################################
#' Calculate a vector of statistics over an \emph{OHLC} time series, and
#' calculate a rolling mean over the statistics.
#'
#' @export
#' @param \code{oh_lc} An \emph{OHLC} time series of prices and trading volumes, in
#'   \emph{xts} format.
#'   
#' @param \code{calc_stats} The name of the function for estimating statistics of a single
#'   row of \emph{OHLC} data, such as volatility, skew, and higher moments.
#'   
#' @param \code{look_back} The size of the look-back interval, equal to the number of
#'   rows of data used for calculating the rolling mean.
#'   
#' @param \code{weight_ed} \emph{Boolean} argument: should statistic be weighted by
#'   trade volume? (default \code{TRUE})
#'   
#' @param ... additional parameters to the function \code{calc_stats}.
#'
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#'
#' @details The function \code{roll_stats()} calculates a vector of statistics
#'   over an \emph{OHLC} time series, such as volatility, skew, and higher
#'   moments.  The statistics could also be any other aggregation of a single
#'   row of \emph{OHLC} data, for example the \emph{High} price minus the
#'   \emph{Low} price squared.  The length of the vector of statistics is equal
#'   to the number of rows of the argument \code{oh_lc}. Then it calculates a
#'   trade volume weighted rolling mean over the vector of statistics over and
#'   calculate statistics.
#'
#' @examples
#' # Calculate time series of rolling variance and skew estimates
#' var_rolling <- roll_stats(oh_lc=HighFreq::SPY, look_back=21)
#' skew_rolling <- roll_stats(oh_lc=HighFreq::SPY, calc_stats="ohlc_skew", look_back=21)
#' skew_rolling <- skew_rolling/(var_rolling)^(1.5)
#' skew_rolling[1, ] <- 0
#' skew_rolling <- rutils::na_locf(skew_rolling)

roll_stats <- function(oh_lc, calc_stats="ohlc_variance", look_back=11, weight_ed=TRUE, ...) {
  
# Match "calc_stats" with moment function
  calc_stats <- match.fun(calc_stats)
  agg_regations <- calc_stats(oh_lc, ...)
  
# Weight by volume
  if (weight_ed) {
    agg_regations <- oh_lc[, 5]*agg_regations
    volume_rolling <- rutils::roll_sum(oh_lc[, 5], look_back=look_back)
    agg_regations <- rutils::roll_sum(agg_regations, look_back=look_back)/volume_rolling
    agg_regations[is.na(agg_regations)] <- 0
  } else
    agg_regations <- rutils::roll_sum(agg_regations, look_back=look_back)/look_back
  # Colnames(agg_regations) <- paste(rutils::get_name(colnames(oh_lc)[1]), "Vol", sep=".")
  
  agg_regations
  
}  # end roll_stats




##########################################################################
#' Calculate the variance of an \emph{OHLC} time series, using different range
#' estimators for variance.
#'
#' @export
#' @param \code{oh_lc} An \emph{OHLC} time series of prices in \emph{xts} format.
#' 
#' @param \code{method} A \emph{character} string representing the method for
#'   estimating variance.  The methods include:
#'   \itemize{
#'     \item "close" close to close,
#'     \item "garman_klass" Garman-Klass,
#'     \item "garman_klass_yz" Garman-Klass with account for close-to-open price jumps,
#'     \item "rogers_satchell" Rogers-Satchell,
#'     \item "yang_zhang" Yang-Zhang,
#'    }
#'    (default is \code{"yang_zhang"})
#'    
#' @param \code{scal_e} \emph{Boolean} argument: should the returns be divided by the
#'   number of seconds in each period? (default is \code{TRUE})
#'
#' @return A single \emph{numeric} value equal to the variance.
#'
#' @details The function \code{calc_var_ohlc_r()} calculates the variance
#'   from all the different intra-day and day-over-day returns (defined as the
#'   differences of \emph{OHLC} prices), using several different variance
#'   estimation methods.
#'
#'   The default method is \code{"yang_zhang"}, which theoretically has the
#'   lowest standard error among unbiased estimators.
#'   The methods \code{"close"}, \code{"garman_klass_yz"}, and
#'   \code{"yang_zhang"} do account for close-to-open price jumps, while the
#'   methods \code{"garman_klass"} and \code{"rogers_satchell"} do not account
#'   for close-to-open price jumps.
#'
#'   If \code{scal_e} is \code{TRUE} (the default), then the returns are divided
#'   by the differences of the time index (which scales the variance to the
#'   units of variance per second squared.) This is useful when calculating the
#'   variance from minutely bar data, because dividing returns by the number of
#'   seconds decreases the effect of overnight price jumps. If the time index is
#'   in days, then the variance is equal to the variance per day squared.
#'   
#'   The function \code{calc_var_ohlc_r()} is implemented in \code{R} code.
#'
#' @examples
#' # Calculate the variance of SPY returns
#' HighFreq::calc_var_ohlc_r(HighFreq::SPY, method="yang_zhang")
#' # Calculate variance without accounting for overnight jumps
#' HighFreq::calc_var_ohlc_r(HighFreq::SPY, method="rogers_satchell")
#' # Calculate the variance without scaling the returns
#' HighFreq::calc_var_ohlc_r(HighFreq::SPY, scal_e=FALSE)

calc_var_ohlc_r <- function(oh_lc, method="yang_zhang", scal_e=TRUE) {
  
  # oh_lc <- log(oh_lc[, 1:4])
  num_rows <- NROW(oh_lc)
  
  # Define the time index for scaling the returns
  if (scal_e)
    in_dex <- c(1, diff(xts::.index(oh_lc)))
  else
    in_dex <- rep(1, num_rows)
  
  # Coerce oh_lc to matrix
  if (is.xts(oh_lc))
    oh_lc <- coredata(oh_lc)
  
  # Calculate all the different intra-day and day-over-day returns 
  # (differences of OHLC prices)
  close_close <- rutils::diff_it(oh_lc[, 4])/in_dex
  open_close <- (oh_lc[, 1]-rutils::lag_it(oh_lc[, 4]))/in_dex
  close_open <- (oh_lc[, 4]-oh_lc[, 1])/in_dex
  close_high <- (oh_lc[, 4]-oh_lc[, 2])/in_dex
  close_low <- (oh_lc[, 4]-oh_lc[, 3])/in_dex
  high_low <- (oh_lc[, 2]-oh_lc[, 3])/in_dex
  high_open <- (oh_lc[, 2]-oh_lc[, 1])/in_dex
  low_open <- (oh_lc[, 3]-oh_lc[, 1])/in_dex
  
  switch(method,
         "close"={var(close_close)},
         "rogers_satchell"={-sum(close_high*high_open + close_low*low_open)/num_rows},
         "garman_klass"={sum(0.5*high_low^2 - (2*log(2)-1)*close_open^2)/num_rows},
         "garman_klass_yz"={sum(0.5*high_low^2 - (2*log(2)-1)*close_open^2)/num_rows + 
             var(open_close)},
         "yang_zhang"={co_eff <- 0.34/(1.34 + (num_rows+1)/(num_rows-1))
         var(open_close) + co_eff*var(close_open) +
           (co_eff-1)*sum(close_high*high_open + close_low*low_open)/num_rows}
  )  # end switch
}  # end calc_var_ohlc_r




##########################################################################
#' Calculate a time series of Sharpe ratios over a rolling look-back interval
#' for an \emph{OHLC} time series.
#'
#' @export
#' @param \code{oh_lc} An \emph{OHLC} time series of prices in \emph{xts} format.
#' 
#' @param \code{look_back} The size of the look-back interval, equal to the number of
#'   rows of data used for aggregating the \emph{OHLC} prices.
#'
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#'
#' @details The function \code{roll_sharpe()} calculates the rolling Sharpe 
#'   ratio defined as the ratio of percentage returns over the look-back
#'   interval, divided by the average volatility of percentage returns.
#'
#' @examples
#' # Calculate rolling Sharpe ratio over SPY
#' sharpe_rolling <- roll_sharpe(oh_lc=HighFreq::SPY, look_back=11)

roll_sharpe <- function(oh_lc, look_back=11) {
  re_turns <- ohlc_returns(oh_lc, lag=look_back, scal_e=FALSE)
  var_rolling <- sqrt(HighFreq::roll_var_ohlc(oh_lc, look_back=look_back, scal_e=FALSE))
  sharpe_rolling <- ifelse(var_rolling==0,
                           1.0,
                           re_turns/var_rolling)
  # Colnames(sharpe_rolling) <- paste0(rutils::get_name(colnames(oh_lc)[1]), ".Sharpe")
  rutils::na_locf(sharpe_rolling)
}  # end roll_sharpe




##########################################################################
#' Calculate a time series of \emph{Hurst} exponents over a rolling look-back
#' interval.
#'
#' @export
#' @param \code{oh_lc} An \emph{OHLC} time series of prices in \emph{xts} format.
#' 
#' @param \code{look_back} The size of the look-back interval, equal to the number of 
#'   rows of data used for aggregating the \emph{OHLC} prices.
#'
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{oh_lc}.
#'
#' @details The function \code{roll_hurst()} calculates a time series of
#'   \emph{Hurst} exponents from \emph{OHLC} prices, over a rolling look-back
#'   interval.
#'
#'   The \emph{Hurst} exponent is defined as the logarithm of the ratio of the
#'   price range, divided by the standard deviation of returns, and divided by
#'   the logarithm of the interval length.
#'
#'   The function \code{roll_hurst()} doesn't use the same definition as the
#'   rescaled range definition of the \emph{Hurst} exponent.
#'   First, because the price range is calculated using \emph{High} and
#'   \emph{Low} prices, which produces bigger range values, and higher
#'   \emph{Hurst} exponent estimates.
#'   Second, because the \emph{Hurst} exponent is estimated using a single
#'   aggregation interval, instead of multiple intervals in the rescaled range
#'   definition.
#'
#'   The rationale for using a different definition of the \emph{Hurst} exponent
#'   is that it's designed to be a technical indicator for use as input into
#'   trading models, rather than an estimator for statistical analysis.
#'
#' @examples
#' # Calculate rolling Hurst for SPY in March 2009
#' hurst_rolling <- roll_hurst(oh_lc=HighFreq::SPY["2009-03"], look_back=11)
#' chart_Series(hurst_rolling["2009-03-10/2009-03-12"], name="SPY hurst_rolling")

roll_hurst <- function(oh_lc, look_back=11) {
  ran_ge <- c(rep(0, look_back-1), (RcppRoll::roll_max(x=oh_lc[, 2], n=look_back) +
               RcppRoll::roll_max(x=-oh_lc[, 3], n=look_back)))
  var_rolling <- sqrt(HighFreq::roll_var_ohlc(oh_lc, look_back=look_back, scal_e=FALSE))
  hurst_rolling <- ifelse((var_rolling==0) | (ran_ge==0),
                          0.5,
                          log(ran_ge/var_rolling)/log(look_back))
  # Colnames(hurst_rolling) <- paste0(rutils::get_name(colnames(oh_lc)[1]), ".Hurst")
  rutils::na_locf(hurst_rolling)
}  # end roll_hurst




##########################################################################
#' Apply an aggregation function over a rolling look-back interval and the end
#' points of an \emph{OHLC} time series, using \code{R} code.
#'
#' @export
#' @param \code{x_ts} An \emph{OHLC} time series of prices and trading volumes, in
#'   \emph{xts} format.
#'   
#' @param \code{agg_fun} The name of the aggregation function to be applied over a
#'   rolling look-back interval.
#'   
#' @param \code{look_back} The number of end points in the look-back interval used for
#'   applying the aggregation function (including the current row).
#'   
#' @param \code{by_columns} \emph{Boolean} argument: should the function
#'   \code{agg_fun()} be applied column-wise (individually), or should it be
#'   applied to all the columns combined? (default is \code{FALSE})
#'   
#' @param \code{out_xts} \emph{Boolean} argument: should the output be coerced into an
#'   \emph{xts} series? (default is \code{TRUE})
#'   
#' @param \code{end_points} An integer vector of end points.
#' 
#' @param ... additional parameters to the function \code{agg_fun}.
#'
#' @return Either an \emph{xts} time series with the number of rows equal to the
#'   length of argument \code{end_points}, or a list the length of argument
#'   \code{end_points}.
#'
#' @details The function \code{roll_apply()} applies an aggregation function 
#'   over a rolling look-back interval attached at the end points of an
#'   \emph{OHLC} time series.
#'
#'   The function \code{roll_apply()} is implemented in \code{R} code.
#'
#'
#'   \code{HighFreq::roll_apply()} performs similar operations to the functions
#'   \code{rollapply()} and \code{period.apply()} from package 
#'   \href{https://cran.r-project.org/web/packages/xts/index.html}{xts}, and 
#'   also the function \code{apply.rolling()} from package 
#'   \href{https://cran.r-project.org/web/packages/PerformanceAnalytics/index.html}{PerformanceAnalytics}.
#'   (The function \code{rollapply()} isn't exported from the package 
#'   \href{https://cran.r-project.org/web/packages/xts/index.html}{xts}.)
#'
#'   But \code{HighFreq::roll_apply()} is faster because it performs less 
#'   type-checking and skips other overhead. Unlike the other functions, 
#'   \code{roll_apply()} doesn't produce any leading \emph{NA} values.
#'
#'   The function \code{roll_apply()} can be called in two different ways,
#'   depending on the argument \code{end_points}.
#'   If the argument \code{end_points} isn't explicitly passed to 
#'   \code{roll_apply()}, then the default value is used, and 
#'   \code{roll_apply()} performs aggregations over overlapping intervals at
#'   each point in time.
#'   
#'   If the argument \code{end_points} is explicitly passed to 
#'   \code{roll_apply()}, then \code{roll_apply()} performs aggregations over 
#'   intervals attached at the end_points.  If look_back=2 then the aggregations
#'   are performed over non-overlapping intervals, otherwise they are performed
#'   over overlapping intervals.
#'
#'   If the argument \code{out_xts} is \code{TRUE} (the default) then the output
#'   is coerced into an \emph{xts} series, with the number of rows equal to the 
#'   length of argument \code{end_points}.  Otherwise a list is returned, with
#'   the length equal to the length of argument \code{end_points}.
#'   
#'   If \code{out_xts} is \code{TRUE} and the aggregation function 
#'   \code{agg_fun()} returns a single value, then \code{roll_apply()} returns 
#'   an \emph{xts} time series with a single column. If \code{out_xts} is
#'   \code{TRUE} and if \code{agg_fun()} returns a vector of values, then
#'   \code{roll_apply()} returns an \emph{xts} time series with multiple
#'   columns, equal to the length of the vector returned by the aggregation
#'   function \code{agg_fun()}.
#'
#' @examples
#' # extract a single day of SPY data
#' oh_lc <- HighFreq::SPY["2012-02-13"]
#' inter_val <- 11  # number of data points between end points
#' look_back <- 4  # number of end points in look-back interval
#' # Calculate the rolling sums of oh_lc columns over a rolling look-back interval
#' agg_regations <- roll_apply(oh_lc, agg_fun=sum, look_back=look_back, by_columns=TRUE)
#' # Apply a vector-valued aggregation function over a rolling look-back interval
#' agg_function <- function(oh_lc)  c(max(oh_lc[, 2]), min(oh_lc[, 3]))
#' agg_regations <- roll_apply(oh_lc, agg_fun=agg_function, look_back=look_back)
#' # Define end points at 11-minute intervals (HighFreq::SPY is minutely bars)
#' end_points <- rutils::end_points(oh_lc, inter_val=inter_val)
#' # Calculate the sums of oh_lc columns over end_points using non-overlapping intervals
#' agg_regations <- roll_apply(oh_lc, agg_fun=sum, end_points=end_points, by_columns=TRUE)
#' # Apply a vector-valued aggregation function over the end_points of oh_lc
#' # using overlapping intervals
#' agg_regations <- roll_apply(oh_lc, agg_fun=agg_function,
#'                             look_back=5, end_points=end_points)

roll_apply <- function(x_ts, agg_fun, look_back=2, end_points=seq_along(x_ts), 
                       by_columns=FALSE, out_xts=TRUE, ...) {
  # Match "agg_fun" with some aggregation function
  agg_fun <- match.fun(agg_fun)
  len_gth <- NROW(end_points)
  # Define start_points as lag of end_points
  start_points <- c(rep_len(1, look_back-1), end_points[1:(len_gth-look_back+1)])
  # Define list of look-back intervals for aggregations over past
  look_backs <- lapply(seq_along(end_points), 
                       function(in_dex) {
                         start_points[in_dex]:end_points[in_dex]
                       })  # end lapply
  # Perform aggregations over length of end_points
  if (by_columns) {
    # Perform individual aggregations by columns
    agg_regations <- lapply(x_ts, function(col_umn)
      lapply(look_backs, function(look_back)
        agg_fun(x_ts[look_back], ...)
      ))  # end lapply
  } else {  # not by_columns
    agg_regations <- lapply(look_backs, function(look_back)
      agg_fun(x_ts[look_back], ...)
    )  # end lapply
  }  # end if
  
  if (out_xts) {
    # Coerce agg_regations into matrix and transpose it
    if (is.null(dim(agg_regations)))
      agg_regations <- t(agg_regations)
    agg_regations <- t(agg_regations)
    # Coerce agg_regations into xts series
    xts(agg_regations, order.by=index(x_ts[end_points]))
  } else
    agg_regations
  
}  # end roll_apply




##########################################################################
#' Perform a backtest simulation of a trading strategy (model) over a vector of
#' end points along a time series of prices.
#'
#' @export
#' @param \code{x_ts} A time series of prices, asset returns, trading volumes, and
#'   other data, in \emph{xts} format.
#'   
#' @param \code{train_func} The name of the function for training (calibrating) a
#'   forecasting model, to be applied over a rolling look-back interval.
#'   
#' @param \code{trade_func} The name of the trading model function, to be applied over
#'   a rolling look-forward interval.
#'   
#' @param \code{look_back} The size of the look-back interval, equal to the number of
#'   rows of data used for training the forecasting model.
#'   
#' @param \code{look_forward} The size of the look-forward interval, equal to the number
#'   of rows of data used for trading the strategy.
#'   
#' @param \code{end_points} A vector of end points along the rows of the \code{x_ts}
#'   time series, given as either integers or dates.
#'   
#' @param ... additional parameters to the functions \code{train_func()} and
#'   \code{trade_func()}.
#'
#' @return An \emph{xts} time series with the number of rows equal to the number
#'   of end points minus two.
#'
#' @details The function \code{roll_backtest()} performs a rolling backtest 
#'   simulation of a trading strategy over a vector of end points. At each end 
#'   point, it trains (calibrates) a forecasting model using past data taken 
#'   from the \code{x_ts} time series over the look-back interval, and applies
#'   the forecasts to the \code{trade_func()} trading model, using out-of-sample
#'   future data from the look-forward interval.
#'   
#'   The function \code{trade_func()} should simulate the trading model, and it 
#'   should return a named list with at least two elements: a named vector of 
#'   performance statistics, and an \emph{xts} time series of out-of-sample 
#'   returns.  The list returned by \code{trade_func()} can also have additional
#'   elements, like the in-sample calibrated model statistics, etc.
#'
#'   The function \code{roll_backtest()} returns a named list containing the
#'   lists returned by function \code{trade_func()}.  The list names are equal
#'   to the \emph{end_points} dates. 
#'   The number of list elements is equal to the number of \emph{end_points}
#'   minus two (because the first and last end points can't be included in the
#'   backtest).
#'
#' @examples
#' \dontrun{
#' # Combine two time series of prices
#' price_s <- cbind(rutils::etf_env$XLU, rutils::etf_env$XLP)
#' look_back <- 252
#' look_forward <- 22
#' # Define end points
#' end_points <- rutils::calc_endpoints(price_s, look_forward)
#' # Perform back-test
#' back_test <- roll_backtest(end_points=end_points,
#'     look_forward=look_forward,
#'     look_back=look_back,
#'     train_func = train_model,
#'     trade_func = trade_model,
#'     model_params = model_params,
#'     trading_params = trading_params,
#'     x_ts=price_s)
#' }

roll_backtest <- function(x_ts,
                          train_func, trade_func,
                          look_back=look_forward,
                          look_forward,
                          end_points=rutils::calc_endpoints(x_ts, look_forward),
                          ...) {
  # Match train and trade function names to functions
  train_func <- match.fun(train_func)
  trade_func <- match.fun(trade_func)

  # Convert end_points dates to integer
  if (inherits(end_points, c("Date", "POSIXt"))) {
    end_dates <- end_points
    end_points <- match(end_points, index(x_ts))
  } else {
    end_dates <- index(x_ts[end_points])
  }  # end if

  # Define integer back_points and fwd_points from integer end_points
  back_points <- end_points - look_back + 1
  back_points[back_points < 1] <- 1
  
  fwd_points <- end_points + look_forward
  fwd_points[fwd_points > NROW(x_ts)] <- NROW(x_ts)

  # Perform backtest over length of end_points
  backtest_range <- 2:(NROW(end_points)-1)
  back_test <- lapply(backtest_range, function(in_dex) {
    trained_model <- 
      train_func(x_ts[back_points[in_dex]:end_points[in_dex]], ...)
    trade_func(x_ts[(end_points[in_dex]+1):fwd_points[in_dex]],
               trained_model, ...)
  })  # end lapply

  names(back_test) <- end_dates[backtest_range]
  back_test
  # Coerce back_test into matrix and transpose it
  # if (is.null(dim(back_test)))
  #   back_test <- t(back_test)
  # back_test <- t(back_test)
  # Coerce back_test into xts series
  # xts(back_test, order.by=index(x_ts[end_points[backtest_range]]))
}  # end roll_backtest




##########################################################################
#' Perform seasonality aggregations over a single-column \emph{xts} time series.
#'
#' @export
#' @param \code{x_ts} A single-column \emph{xts} time series.
#' 
#' @param \code{in_dex} A vector of \emph{character} strings representing points in
#'   time, of the same length as the argument \code{x_ts}.
#'
#' @return An \emph{xts} time series with mean aggregations over the seasonality
#'   interval.
#'
#' @details The function \code{season_ality()} calculates the mean of values
#'   observed at the same points in time specified by the argument
#'   \code{in_dex}. An example of a daily seasonality aggregation is the average
#'   price of a stock between 9:30AM and 10:00AM every day, over many days. The
#'   argument \code{in_dex} is passed into function \code{tapply()}, and must be
#'   the same length as the argument \code{x_ts}.
#'
#' @examples
#' # Calculate running variance of each minutely OHLC bar of data
#' x_ts <- ohlc_variance(HighFreq::SPY)
#' # Remove overnight variance spikes at "09:31"
#' in_dex <- format(index(x_ts), "%H:%M")
#' x_ts <- x_ts[!in_dex=="09:31", ]
#' # Calculate daily seasonality of variance
#' var_seasonal <- season_ality(x_ts=x_ts)
#' chart_Series(x=var_seasonal, name=paste(colnames(var_seasonal),
#'   "daily seasonality of variance"))

season_ality <- function(x_ts, in_dex=format(zoo::index(x_ts), "%H:%M")) {
# Aggregate the mean
  agg_regation <- tapply(X=x_ts, INDEX=in_dex, FUN=mean)
# Coerce from array to named vector
  agg_regation <- structure(as.vector(agg_regation), names=names(agg_regation))
# Coerce to xts
  agg_regation <- xts(x=agg_regation,
      order.by=as.POSIXct(paste(Sys.Date(), names(agg_regation))))
  colnames(agg_regation) <- colnames(x_ts)
  agg_regation
}  # end season_ality


