##########################################################################
#' Simulate a random \emph{TAQ} time series of prices and trading volumes, in
#' \emph{xts} format.
#'
#' Simulate a \emph{TAQ} time series of random prices following geometric
#' Brownian motion, combined with random trading volumes.
#'
#' @param \code{volat} The volatility per period of the \code{datev} time index
#'   (default is \code{6.5e-05} per second, or about \code{0.01=1.0\%} per day).
#' @param \code{drift} The drift per period of the \code{datev} time index (default
#'   is 0.0).
#' @param \code{datev} The time index for the \emph{TAQ} time series.
#' @param bidask The bid-ask spread expressed as a fraction of the prices
#'   (default is 0.001=10bps).
#'
#' @return An \emph{xts} time series, with time index equal to the input
#'   \code{datev} time index, and with four columns containing the bid, ask,
#'   and trade prices, and the trade volume.
#'
#' @details The function \code{random_taq()} calculates an \emph{xts} time
#'   series with four columns containing random prices following geometric
#'   Brownian motion: the bid, ask, and trade prices, combined with random trade
#'   volume data.
#'   If \code{datev} isn't supplied as an argument, then by default it's
#'   equal to the secondly index over the two previous calendar days.
#'
#' @examples
#' # Create secondly TAQ time series of random prices
#' taq <- HighFreq::random_taq()
#' # Create random TAQ time series from SPY index
#' taq <- HighFreq::random_taq(datev=index(HighFreq::SPY["2012-02-13/2012-02-15"]))
#'
#' @export

random_taq <- function(volat=6.5e-5, drift=0.0,
  datev=seq(from=as.POSIXct(paste(Sys.Date()-3, "09:30:00")),
             to=as.POSIXct(paste(Sys.Date()-1, "16:00:00")), by="1 sec"),
  bidask=0.001, ...) {
  
  nrows <- NROW(datev)
  # Create xts of random prices following geometric Brownian motion
  taq <- xts(exp(cumsum(volat*rnorm(nrows) + drift - volat^2/2)),
              order.by=datev)
  # Create vector of random bid-ask spreads
  bidask <- bidask*(1 + runif(nrows))/2
  # Create TAQ data from bid and offer prices
  taq <- merge(taq*(1-bidask), taq*(1+bidask))
  # Add traded price to TAQ data
  pricer <- runif(nrows)
  taq <- merge(taq, pricer*taq[, 1] + (1-pricer)*taq[, 2])
  # Add trade volume column
  taq <- merge(taq, sample(x=10*(2:18), size=nrows, replace=TRUE))
  colnames(taq) <- c("Bid.Price", "Ask.Price", "Trade.Price", "Volume")
  return(taq)
  
}  # end random_taq




##########################################################################
#' Simulate a random \emph{OHLC} time series of prices and trading volumes, in
#' \emph{xts} format.
#'
#' Simulate a random \emph{OHLC} time series either by simulating random prices
#' following geometric Brownian motion, or by randomly sampling from an input
#' time series.
#'
#' @param \code{ohlc} An \emph{OHLC} time series of prices and trading volumes,
#'   in \emph{xts} format (default is \emph{NULL}).
#'   
#' @param \code{volat} The volatility per period of the \code{datev} time index
#'   (default is \code{6.5e-05} per second, or about \code{0.01=1.0\%} per day).
#'   
#' @param \code{drift} The drift per period of the \code{datev} time index
#'   (default is 0.0).
#'   
#' @param \code{datev} The time index for the \emph{OHLC} time series.
#' 
#' @param \code{reducit} \emph{Boolean} argument: should \code{ohlc} time series
#'   be transformed to reduced form? (default is \code{TRUE})
#'
#' @return An \emph{xts} time series with the same dimensions and the same time
#'   index as the input \code{ohlc} time series.
#'
#' @details If the input \code{ohlc} time series is \emph{NULL} (the default),
#'   then the function \code{random_ohlc()} simulates a minutely \emph{OHLC}
#'   time series of random prices following geometric Brownian motion, over the
#'   two previous calendar days.
#'
#'   If the input \code{ohlc} time series is not \emph{NULL}, then the rows of
#'   \code{ohlc} are randomly sampled, to produce a random time series.
#'
#'   If \code{reducit} is \code{TRUE} (the default), then the \code{ohlc} time
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
#' ohlc <- HighFreq::random_ohlc()
#' # Create random time series from SPY by randomly sampling it
#' ohlc <- HighFreq::random_ohlc(ohlc=HighFreq::SPY["2012-02-13/2012-02-15"])
#'
#' @export

random_ohlc <- function(ohlc=NULL, reducit=TRUE, volat=6.5e-5, drift=0.0,
    datev=seq(from=as.POSIXct(paste(Sys.Date()-3, "09:30:00")),
      to=as.POSIXct(paste(Sys.Date()-1, "16:00:00")), by="1 sec"), ...) {
  
  if (is.null(ohlc)) {
    nrows <- NROW(datev)
    # Create xts of random prices following geometric Brownian motion
    xtsv <- xts(exp(cumsum(volat*rnorm(nrows) + drift - volat^2/2)), order.by=datev)
    # Add trade volume column
    xtsv <- merge(xtsv, volume=sample(x=10*(2:18), size=nrows, replace=TRUE))
    # Aggregate to minutes OHLC data
    return(to.period(x=xtsv, period="minutes"))
  } else {
    ohlc <- log(ohlc)  # transform to normal
    if (reducit)  # Calculate reduced form of ohlc
      ohlc <- rutils::diffohlc(ohlc)
    # randomly sample from the rows of ohlc
    ohlc <- xts(coredata(ohlc)[c(1, sample(x=2:NROW(ohlc), replace=TRUE)), ], order.by=index(ohlc))
    # Return standard form of randomized ohlc
    return(exp(rutils::diffohlc(ohlc, reducit=FALSE)))
  } # end if
  
}  # end random_ohlc




##########################################################################
#' Remove overnight close-to-open price jumps from an \emph{OHLC} time series,
#' by adding adjustment terms to its prices.
#'
#' @param \code{ohlc} An \emph{OHLC} time series of prices and trading volumes, in
#'   \emph{xts} format.
#'
#' @return An \emph{OHLC} time series with the same dimensions and the same time
#'   index as the input \code{ohlc} time series.
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
#'   The time index of the \code{ohlc} time series is assumed to be in
#'   \emph{POSIXct} format, so that its internal value is equal to the number of
#'   seconds that have elapsed since the \emph{epoch}.
#'
#' @examples
#' # Remove overnight close-to-open price jumps from SPY data
#' ohlc <- remove_jumps(HighFreq::SPY)
#'
#' @export

remove_jumps <- function(ohlc) {
  
  # Find time index of the periods greater than 60 seconds
  timev <- which(c(1, diff(xts::.index(ohlc))) > 60)
  # Calculate cumulative sum of overnight price jumps
  jumpv <- numeric(NROW(ohlc))
  jumpv[timev] <- as.numeric(ohlc[timev, 1]) - as.numeric(ohlc[timev-1, 4])
  jumpv <- cumsum(jumpv)
  # subtract overnight price jumps from OHLC
  ohlc[, 1:4] <- coredata(ohlc[, 1:4]) - jumpv
  return(ohlc)
  
}  # end remove_jumps




##########################################################################
#' Calculate single period percentage returns from either \emph{TAQ} or
#' \emph{OHLC} prices.
#'
#' @param \code{xtsv} An \emph{xts} time series of either \emph{TAQ} or
#'   \emph{OHLC} data.
#' 
#' @param \code{lagg} An integer equal to the number of time periods of lag
#'   (default is 1).
#'   
#' @param \code{colnum} The column number to extract from the \emph{OHLC} data
#'   (default is \code{4}, or the \emph{Close} prices column).
#'   
#' @param \code{scalit} \emph{Boolean} argument: should the returns be divided
#'   by the number of seconds in each period? (default is \code{TRUE})
#'
#' @return A single-column \emph{xts} time series of returns.
#'
#' @details The function \code{ohlc_returns()} calculates the percentage returns
#'   for either \emph{TAQ} or \emph{OHLC} data, defined as the difference of log
#'   prices.  Multi-period returns can be calculated by setting the \code{lag}
#'   parameter to values greater than \code{1} (the default).
#'
#'   If \code{scalit} is \code{TRUE} (the default), then the returns are divided
#'   by the differences of the time index (which scales the returns to units of
#'   returns per second.)
#'
#'   The time index of the \code{xtsv} time series is assumed to be in
#'   \emph{POSIXct} format, so that its internal value is equal to the number of
#'   seconds that have elapsed since the \emph{epoch}.
#'
#'   If \code{scalit} is \code{TRUE} (the default), then the returns are
#'   expressed in the scale of the time index of the \code{xtsv} time series.
#'   For example, if the time index is in seconds, then the returns are given in
#'   units of returns per second.  If the time index is in days, then the
#'   returns are equal to the returns per day.
#'
#'   The function \code{ohlc_returns()} identifies the \code{xtsv} time series as
#'   \emph{TAQ} data when it has six columns, otherwise assumes it's \emph{OHLC}
#'   data. By default, for \emph{OHLC} data, it differences the \emph{Close}
#'   prices, but can also difference other prices depending on the value of
#'   \code{colnum}.
#'
#' @examples
#' # Calculate secondly returns from TAQ data
#' returns <- HighFreq::ohlc_returns(xtsv=HighFreq::SPY_TAQ)
#' # Calculate close to close returns
#' returns <- HighFreq::ohlc_returns(xtsv=HighFreq::SPY)
#' # Calculate open to open returns
#' returns <- HighFreq::ohlc_returns(xtsv=HighFreq::SPY, colnum=1)
#' 
#' @export

ohlc_returns <- function(xtsv, lagg=1, colnum=4, scalit=TRUE) {
  
  # Return NULL if no data
  if (is.null(xtsv))  return(NULL)
  # Calculate mid prices
  if (NCOL(xtsv)==6)  # TAQ data has 6 columns
    retp <- 0.5 * (xtsv[, "Bid.Price"] + xtsv[, "Ask.Price"])
  else
    retp <- xtsv[, colnum]  # OHLC data
  # Calculate log returns
  retp <- rutils::diffit(log(retp), lagg=lagg)
  if (scalit)
    retp <- retp / c(rep(1, lagg), diff(xts::.index(retp), lagg=lagg))
  retp[1:lagg, ] <- 0
  # Colnames(retp) <- paste0(rutils::get_name(colnames(xtsv)[1]), ".retp")
  return(retp)
  
}  # end ohlc_returns




##########################################################################
#' Calculate a \emph{Boolean} vector that identifies extreme tail values in a
#' single-column \emph{xts} time series or vector, over a rolling look-back
#' interval.
#'
#' @param \code{xtsv} A single-column \emph{xts} time series, or a
#'   \emph{numeric} or \emph{Boolean} vector.
#'   
#' @param \code{lookb} The number of data points in rolling look-back interval
#'   for estimating rolling quantile.
#'   
#' @param \code{volm} The quantile multiplier.
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
#'   The quantile multiplier \code{volm} controls the threshold at which
#'   values are identified as extreme. Smaller quantile multiplier values will
#'   cause more values to be identified as extreme.
#'
#' @examples
#' # Create local copy of SPY TAQ data
#' taq <- HighFreq::SPY_TAQ
#' # scrub quotes with suspect bid-ask spreads
#' bidask <- taq[, "Ask.Price"] - taq[, "Bid.Price"]
#' jumpv <- which_extreme(bidask, lookb=51, volm=3)
#' # Remove suspect values
#' taq <- taq[!jumpv]
#' 
#' @export

which_extreme <- function(xtsv, lookb=51, volm=2) {
  
# Calculate volatility as rolling quantile
  quantv <- caTools::runquantile(x=abs(as.numeric(xtsv)), k=lookb,
                        probs=0.9, endrule="constant", align="center")
#  quantv <- xts(quantv, order.by=index(xtsv))
#  colnames(quantv) <- "volat"
# Carry forward non-zero volatility values
  quantv[quantv==0] <- NA
  quantv[1] <- 1
  quantv <- rutils::na_locf(quantv)
#  quantv <- rutils::na_locf(quantv, fromLast=TRUE)

# extreme value if xtsv greater than scaled volatility
  exv <- (abs(xtsv) > 2*volm*quantv)
  exv[1] <- FALSE
#  colnames(exv) <- "suspect"

  cat("date:", format(as.Date(index(first(xtsv)))), "\tfound", sum(exv), "extreme values\n")
  return(exv)
  
}  # end which_extreme




##########################################################################
#' Calculate a \emph{Boolean} vector that identifies isolated jumps (spikes) in
#' a single-column \emph{xts} time series or vector, over a rolling interval.
#'
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
#'   The quantile multiplier \code{volm} controls the threshold at which
#'   values are identified as jumps. Smaller quantile multiplier values will
#'   cause more values to be identified as jumps.
#'
#' @examples
#' # Create local copy of SPY TAQ data
#' taq <- SPY_TAQ
#' # Calculate mid prices
#' pricem <- 0.5 * (taq[, "Bid.Price"] + taq[, "Ask.Price"])
#' # Replace whole rows containing suspect price jumps with NA, and perform locf()
#' taq[which_jumps(pricem, lookb=31, volm=1.0), ] <- NA
#' taq <- xts:::na.locf.xts(taq)
#' 
#' @export

which_jumps <- function(xtsv, lookb=51, volm=2) {
  
# Calculate simple returns
  rets <- rutils::diffit(as.numeric(xtsv))
#  rets[1] <- 0
#  colnames(rets) <- "diffs"
  retadv <- rutils::lagit(rets, -1)
#  retadv[NROW(retadv)] <- 0
#  colnames(retadv) <- "retadv"

# Calculate volatility as the rolling quantile of returns
  quantv <- caTools::runquantile(x=abs(rets), k=lookb,
                        probs=0.9, endrule="constant", align="center")
#  quantv <- xts(quantv, order.by=index(rets))
#  colnames(quantv) <- "volat"
# Carry forward non-zero quantv values
  quantv[quantv==0] <- NA
  quantv[1] <- 1
  quantv <- rutils::na_locf(quantv)
#  quantv <- rutils::na_locf(quantv, fromLast=TRUE)

# value is suspect if abs returns greater than quantv,
# and if abs sum of returns less than quantv
  jumpv <- ((abs(rets) > volm*quantv) &
      (abs(retadv) > volm*quantv) &
      (abs(rets+retadv) < 2*volm*quantv))
  jumpv[1] <- FALSE
#  colnames(jumpv) <- "suspect"
# Cat("Parsing", deparse(substitute(taq)), "\n")
# Cat("Parsing", strsplit(deparse(substitute(taq)), split="[.]")[[1]][4], "on date:", format(todayd), "\tfound", sum(jumpv), "suspect prices\n")
  cat("date:", format(as.Date(index(first(xtsv)))), "\tfound", sum(jumpv), "jump prices\n")
  return(jumpv)
  
}  # end which_jumps




##########################################################################
#' Scrub a single day of \emph{TAQ} data in \emph{xts} format, without
#' aggregation.
#'
#' @inheritParams which_extreme
#' 
#' @param \code{taq} \emph{TAQ} A time series in \emph{xts} format.
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
#' taq <- HighFreq::scrub_taq(taq=HighFreq::SPY_TAQ, lookb=11, volm=1)
#' # Create random TAQ prices and scrub them
#' taq <- HighFreq::random_taq()
#' taq <- HighFreq::scrub_taq(taq=taq)
#' taq <- HighFreq::scrub_taq(taq=taq, lookb=11, volm=1)
#' 
#' @export

scrub_taq <- function(taq, lookb=51, volm=2, tzone="America/New_York") {
  
# Convert timezone of index to New_York
  index(taq) <- lubridate::with_tz(time=index(taq), tzone=tzone)
# subset data to NYSE trading hours
  taq <- taq["T09:30:00/T16:00:00", ]
# Return NULL if no data
  if (NROW(taq)==0)  return(NULL)
#  todayd <- as.Date(index(first(taq)))

# Remove duplicate time stamps using duplicated()
  taq <- taq[!duplicated(index(taq)), ]

# scrub quotes with suspect bid-ask spreads
  bidask <- taq[, "Ask.Price"] - taq[, "Bid.Price"]
#  bidask <- na.omit(bidask)
  jumpv <- which_extreme(bidask, lookb=lookb, volm=volm)
# Remove suspect values
  taq <- taq[!jumpv]
# Replace suspect values
# taq[jumpv, "Bid.Price"] <- taq[jumpv, "Trade.Price"]
# taq[jumpv, "Ask.Price"] <- taq[jumpv, "Trade.Price"]

# scrub quotes with suspect price jumps
# Calculate mid prices
  pricem <- 0.5 * (taq[, "Bid.Price"] + taq[, "Ask.Price"])
#  pricem <- na.omit(pricem)
#  colnames(pricem) <- "Mid.Price"
# Replace NA volumes with zero
  taq[is.na(taq[, "Volume"]), "Volume"] <- 0
# Replace whole rows containing suspect price jumps with NA, and perform locf()
  taq[which_jumps(pricem, lookb=lookb, volm=volm), ] <- NA
  return(rutils::na_locf(taq))
  
}  # end scrub_taq




##########################################################################
#' Scrub a single day of \emph{TAQ} data, aggregate it, and convert to
#' \emph{OHLC} format.
#'
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
#'   \item scrubbing of quotes with suspect bid-ask spreads,
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
#' taq <- HighFreq::random_taq()
#' # Aggregate to ten minutes OHLC data
#' ohlc <- HighFreq::scrub_agg(taq, period="10 min")
#' chart_Series(ohlc, name="random prices")
#' # scrub and aggregate a single day of SPY TAQ data to OHLC
#' ohlc <- HighFreq::scrub_agg(taq=HighFreq::SPY_TAQ)
#' chart_Series(ohlc, name=symbol)
#' 
#' @export

scrub_agg <- function(taq, lookb=51, volm=2,
                      period="minutes", tzone="America/New_York") {
# Convert timezone of index to New_York
  index(taq) <- lubridate::with_tz(time=index(taq), tzone=tzone)
# subset data to NYSE trading hours
  taq <- taq["T09:30:00/T16:00:00", ]
# Return NULL if no data
  if (NROW(taq)==0)  return(NULL)
#  todayd <- as.Date(index(first(taq)))

# Remove duplicate time stamps using duplicated()
  taq <- taq[!duplicated(index(taq)), ]

# scrub quotes with suspect bid-ask spreads
  bidask <- taq[, "Ask.Price"] - taq[, "Bid.Price"]
#  bidask <- na.omit(bidask)
  jumpv <- which_extreme(bidask, lookb=lookb, volm=volm)
# Remove suspect values
  taq <- taq[!jumpv]
# Replace suspect values
# taq[jumpv, "Bid.Price"] <- taq[jumpv, "Trade.Price"]
# taq[jumpv, "Ask.Price"] <- taq[jumpv, "Trade.Price"]

# scrub quotes with suspect price jumps
# Calculate mid prices
  pricem <- 0.5 * (taq[, "Bid.Price"] + taq[, "Ask.Price"])
#  pricem <- na.omit(pricem)
  colnames(pricem) <- "Mid.Price"
# Replace whole rows containing suspect price jumps with NA, and perform locf()
  pricem[which_jumps(pricem, lookb=lookb, volm=volm)] <- NA
  pricem <- rutils::na_locf(pricem)
#  pricem <- rutils::na_locf(pricem, fromLast=TRUE)
# Cbind pricem with volume data, and replace NA volumes with zero
  pricem <- cbind(pricem, taq[index(pricem), "Volume"])
  pricem[is.na(pricem[, "Volume"]), "Volume"] <- 0

# Aggregate to OHLC and cumulative volume data
  pricem <- switch(period,
                   "minutes"={numsec <- 60; to.period(x=pricem, period=period)},
                   "3 min"={numsec <- 3*60; to.minutes3(x=pricem)},
                   "5 min"={numsec <- 5*60; to.minutes5(x=pricem)},
                   "10 min"={numsec <- 10*60; to.minutes10(x=pricem)},
                   "15 min"={numsec <- 15*60; to.minutes15(x=pricem)},
                   "30 min"={numsec <- 30*60; to.minutes30(x=pricem)},
                   "hours"={numsec <- 60*60; to.period(x=pricem, period=period)}
  )  # end switch
  # round up times to next period
  index(pricem) <- align.time(x=index(pricem), n=numsec)
  return(pricem)
  
}  # end scrub_agg




##########################################################################
#' Load, scrub, aggregate, and rbind multiple days of \emph{TAQ} data for a
#' single symbol, and save the \emph{OHLC} time series to a single
#' \sQuote{\code{*.RData}} file.
#'
#' @param \code{symbol} A \emph{character} string representing symbol or ticker.
#' 
#' @param \code{dirin} A \emph{character} string representing directory
#'   containing input \sQuote{\code{*.RData}} files.
#'   
#' @param \code{dirout} A \emph{character} string representing directory
#'   containing output \sQuote{\code{*.RData}} files.
#'   
#' @inheritParams scrub_agg
#'
#' @return An \emph{OHLC} time series in \emph{xts} format.
#'
#' @details The function \code{save_scrub_agg()} loads multiple days of
#'   \emph{TAQ} data, then scrubs, aggregates, and rbinds them into a
#'   \emph{OHLC} time series, and finally saves it to a single
#'   \sQuote{\code{*.RData}} file. The \emph{OHLC} time series is stored in a
#'   variable named \sQuote{\code{symbol}}, and then it's saved to a file named
#'   \sQuote{\code{symbol.RData}} in the \sQuote{\code{dirout}} directory.
#'   The \emph{TAQ} data files are assumed to be stored in separate directories
#'   for each \sQuote{\code{symbol}}. Each \sQuote{\code{symbol}} has its own
#'   directory (named \sQuote{\code{symbol}}) in the \sQuote{\code{dirin}}
#'   directory. Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \emph{TAQ}
#'   data.
#' @examples
#' \dontrun{
#' # set data directories
#' dirin <- "C:/Develop/data/hfreq/src/"
#' dirout <- "C:/Develop/data/hfreq/scrub/"
#' symbol <- "SPY"
#' # Aggregate SPY TAQ data to 15-min OHLC bar data, and save the data to a file
#' save_scrub_agg(symbol=symbol, dirin=dirin, dirout=dirout, period="15 min")
#' }
#' 
#' @export

save_scrub_agg <- function(symbol,
                      dirin="E:/mktdata/sec/",
                      dirout="E:/output/data/",
                      lookb=51,
                      volm=2,
                      period="minutes",
                      tzone="America/New_York") {
  
# Create path to directory containing *.RData files
  pathin <- file.path(dirin, symbol)
# Get list of *.RData files
  filelist <- list.files(pathin)
# Create paths to *.RData files
  filell <- file.path(pathin, filelist)

# load TAQ data one by one, scrub and aggregate it, return list of xts
datav <- lapply(filell, function(filen) {
  cat("loading", symbol, "from file: ", filen, "\n")
  symbol <- load(filen)
  scrub_agg(get(symbol),
            lookb=lookb,
            volm=volm,
            period=period, 
            tzone=tzone)
})  # end sapply

# Recursively "rbind" the list into a single xts
  datav <- rutils::do_call_rbind(datav)
# assign column names, i.e. "symbol.High"
  colnames(datav) <- sapply(strsplit(colnames(datav), split="[.]"),
                           function(strng) paste(symbol, strng[-1], sep="."))

# Copy the xts data to a variable with the name 'symbol'
  assign(symbol, datav)

# save the xts data to a file in the dirout
  save(list=symbol, file=file.path(dirout, paste0(symbol, ".RData")))
  invisible(symbol)

}  # end save_scrub_agg




##########################################################################
#' Load and scrub multiple days of \emph{TAQ} data for a single symbol, and save
#' it to multiple \sQuote{\code{*.RData}} files.
#'
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
#'   in the \sQuote{\code{dirin}} directory.
#'   Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \emph{TAQ}
#'   data.
#'
#' @examples
#' \dontrun{
#' save_taq("SPY")
#' }
#' 
#' @export

save_taq <- function(symbol,
                     dirin="E:/mktdata/sec/",
                     dirout="E:/output/data/",
                     lookb=51,
                     volm=2,
                     tzone="America/New_York") {
  
# Create path to directory containing *.RData files
  dirin <- file.path(dirin, symbol)
  # Create path to directory for writing *.RData files
  dirout <- file.path(dirout, symbol)
  # Get list of *.RData files
  filell <- list.files(dirin)

# load TAQ data one-by-one, scrub, and save
  outputs <- sapply(filell, function(filen) {
    cat("loading", symbol, "from file: ", filen, "\n")
    filin <- file.path(dirin, filen)
    symbol <- load(filin)
    filout <- file.path(dirout, filen)
# save the xts data to a file in the dirout
    taq <- scrub_taq(get(symbol), lookb=lookb, volm=volm, tzone=tzone)
    if (!is.null(taq)) {
      assign(symbol, taq)
      save(list=symbol, file=filout)
      cat("finished saving", symbol, "to file: ", filen, "\n")
    } # end if
    filen
  })  # end sapply

  invisible(symbol)

}  # end save_taq




##########################################################################
#' Load, scrub, aggregate, and rbind multiple days of \emph{TAQ} data for a
#' single symbol. Calculate the returns and save them to a single
#' \sQuote{\code{*.RData}} file.
#'
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
#'   directory (named \sQuote{\code{symbol}}) in the \sQuote{\code{dirin}}
#'   directory. Each \sQuote{\code{symbol}} directory contains multiple daily
#'   \sQuote{\code{*.RData}} files, each file containing one day of \emph{TAQ}
#'   data.
#'
#' @examples
#' \dontrun{
#' save_rets("SPY")
#' }
#' 
#' @export

save_rets <- function(symbol,
                      dirin="E:/mktdata/sec/",
                      dirout="E:/output/data/",
                      lookb=51,
                      volm=2,
                      period="minutes",
                      tzone="America/New_York") {
  
# Create path to directory containing *.RData files
  pathin <- file.path(dirin, symbol)
# Get list of *.RData files
  filelist <- list.files(pathin)
# Create paths to *.RData files
  filell <- file.path(pathin, filelist)

# load TAQ data into list
  taq <- sapply(filell, function(filen) {
    cat("loading", symbol, "from file: ", filen, "\n")
    symbol <- load(filen)
    get(symbol)
  })

# scrub and aggregate the TAQ data
  ohlc <- lapply(taq, scrub_agg,
                 lookb=lookb,
                 volm=volm,
                 period=period,
                 tzone=tzone)

# Calculate returns
  ohlc <- lapply(ohlc, ohlc_returns)

# Recursively "rbind" the list into a single xts
  ohlc <- rutils::do_call_rbind(ohlc)
# assign column names, i.e. "symbol.rets"
  colnames(ohlc) <- c(paste(symbol, "rets", sep="."), paste(symbol, "vol", sep="."))

# Copy the xts data to a variable with the name 'symbol'
  symbolr <- paste(symbol, "rets", sep=".")
  assign(symbolr, ohlc)

# save the xts data to a file in the dirout
  save(list=eval(symbolr),
       file=file.path(dirout, paste0(symbolr, ".RData")))
  
  invisible(symbolr)

}  # end save_rets




##########################################################################
#' Load \emph{OHLC} time series data for a single symbol, calculate its returns,
#' and save them to a single \sQuote{\code{*.RData}} file, without aggregation.
#'
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
#' 
#' @export

save_rets_ohlc <- function(symbol,
                      dirin="E:/output/data/",
                      dirout="E:/output/data/") {
  
# Create path to directory containing symbol.RData file
  filen <- file.path(dirin, paste0(symbol, ".RData"))
# load OHLC data
  cat("loading", symbol, "from file: ", filen, "\n")
  symbol <- load(filen)

# Calculate returns
  datav <- ohlc_returns(get(symbol))

# Copy the xts data to a variable with the name 'symbol'
  symbolr <- paste(symbol, "rets", sep=".")
  assign(symbolr, datav)

# save the xts data to a file in the dirout
  cat("saving", symbol, "to file: ", paste0(symbolr, ".RData"), "\n")
  save(list=eval(symbolr), file=file.path(dirout, paste0(symbolr, ".RData")))
  invisible(symbolr)

}  # end save_rets_ohlc




##########################################################################
#' Calculate the Value at Risk (\emph{VaR}) or the Conditional Value at Risk
#' (\emph{CVaR}) of an \emph{xts} \emph{time series} of returns, using \code{R}
#' code.
#' 
#' @param \code{timeser} An \emph{xts} \emph{time series} of returns with
#'   multiple columns.
#'   
#' @param \code{method} A \emph{string} specifying the type of risk measure
#'   (the default is \code{method = "var"} - see Details).
#'    
#' @param \code{confi} The confidence level for calculating the
#'   quantile (the default is \code{confi = pnorm(-2) = 0.02275}).
#'
#' @return A vector with the risk measures of the columns of the input
#'   \emph{time series} \code{timeser}.
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
#' returns <- na.omit(rutils::etfenv$returns[, c("VTI", "XLF")])
#' # Calculate VaR
#' all.equal(HighFreq::calc_cvar(returns), 
#'   sapply(returns, quantile, probs=pnorm(-2)), check.attributes=FALSE)
#' # Calculate CVaR
#' all.equal(HighFreq::calc_cvar(returns, method="cvar", confi=0.02), 
#'   sapply(returns, function(x) mean(x[x < quantile(x, 0.02)])), 
#'   check.attributes=FALSE)
#' }
#' 
#' @export

calc_cvar <- function(timeser, method = "var", confi = pnorm(-2)) {

  # Switch for the different risk methods
  riskv <- switch(method, "var"={sapply(timeser, quantile, probs=confi)},
                 # Calculate CVaR as expected loss
                 "cvar"={sapply(timeser, function(x) mean(x[x < quantile(x, confi)]))}
  )  # end switch
  
  return(riskv)
  
}  # end calc_cvar




##########################################################################
#' Calculate a time series of point estimates of variance for an \emph{OHLC}
#' time series, using different range estimators for variance.
#'
#' Calculates the point variance estimates from individual rows of \emph{OHLC}
#' prices (rows of data), using the squared differences of \emph{OHLC} prices at
#' each point in time, without averaging them over time.
#'
#' @param \code{ohlc} An \emph{OHLC} time series of prices in \emph{xts} format.
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
#' @param \code{scalit} \emph{Boolean} argument: should the returns be divided by the
#'   number of seconds in each period? (default is \code{TRUE})
#'
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{ohlc}.
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
#'   If \code{scalit} is \code{TRUE} (the default), then the variance is divided
#'   by the squared differences of the time index (which scales the variance to
#'   units of variance per second squared.) This is useful for example, when
#'   calculating intra-day variance from minutely bar data, because dividing
#'   returns by the number of seconds decreases the effect of overnight price
#'   jumps.
#'
#'   If \code{scalit} is \code{TRUE} (the default), then the variance is
#'   expressed in the scale of the time index of the \emph{OHLC} time series.
#'   For example, if the time index is in seconds, then the variance is given in
#'   units of variance per second squared.  If the time index is in days, then
#'   the variance is equal to the variance per day squared.
#'
#'   The time index of the \code{ohlc} time series is assumed to be in
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
#' ohlc <- HighFreq::random_ohlc()
#' # Calculate variance estimates for ohlc
#' var_running <- HighFreq::ohlc_variance(ohlc)
#' # Calculate variance estimates for SPY
#' var_running <- HighFreq::ohlc_variance(HighFreq::SPY, method="yang_zhang")
#' # Calculate SPY variance without overnight jumps
#' var_running <- HighFreq::ohlc_variance(HighFreq::SPY, method="rogers_satchell")
#' 
#' @export

ohlc_variance <- function(ohlc, method="yang_zhang", scalit=TRUE) {
  
  symbol <- rutils::get_name(colnames(ohlc)[1])
  # ohlc <- log(ohlc[, 1:4])
  varv <- switch(method,
     "close"={rutils::diffit(ohlc[, 4])^2},
     "garman_klass"={0.5*(ohlc[, 2]-ohlc[, 3])^2 -
         (2*log(2)-1)*(ohlc[, 4]-ohlc[, 1])^2},
     "rogers_satchell"={(ohlc[, 2]-ohlc[, 4])*(ohlc[, 2]-ohlc[, 1]) +
         (ohlc[, 3]-ohlc[, 4])*(ohlc[, 3]-ohlc[, 1])},
     "garman_klass_yz"={(ohlc[, 1]-rutils::lagit(ohlc[, 4]))^2 +
         0.5*(ohlc[, 2]-ohlc[, 3])^2 -
         (2*log(2)-1)*(ohlc[, 4]-ohlc[, 1])^2},
     "yang_zhang"={coeff <- 0.34/2.34
     (ohlc[, 1]-rutils::lagit(ohlc[, 4]))^2 +
       coeff*(ohlc[, 1]-ohlc[, 4])^2 +
       (1-coeff)*((ohlc[, 2]-ohlc[, 4])*(ohlc[, 2]-ohlc[, 1]) +
                    (ohlc[, 3]-ohlc[, 4])*(ohlc[, 3]-ohlc[, 1]))}
  )  # end switch
  if (scalit)
    varv <- varv/c(1, diff(xts::.index(ohlc)))^2
  varv[1, ] <- 0
  varv <- rutils::na_locf(varv)
  # Colnames(varv) <- paste0(symbol, ".varv")
  return(varv)
  
}  # end ohlc_variance




##########################################################################
#' Calculate time series of point skew estimates from a \emph{OHLC} time series,
#' assuming zero drift.
#'
#' @param \code{ohlc} An \emph{OHLC} time series of prices in \emph{xts} format.
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
#' skew <- HighFreq::ohlc_skew(HighFreq::SPY)
#' 
#' @export

ohlc_skew <- function(ohlc, method="rogers_satchell") {
  
  symbol <- rutils::get_name(colnames(ohlc)[1])
  # ohlc <- log(ohlc[, 1:4])
  skewv <- switch(method,
    "close"={rutils::diffit(ohlc[, 4])^3},
    "garman_klass"={0.5*(ohlc[, 2]-ohlc[, 3])^3 -
        (2*log(2)-1)*(ohlc[, 4]-ohlc[, 1])^3},
    "rogers_satchell"={
      (ohlc[, 2]-ohlc[, 4])*(ohlc[, 2]-ohlc[, 1])*(ohlc[, 2]-0.5*(ohlc[, 4] + ohlc[, 1])) +
        (ohlc[, 3]-ohlc[, 4])*(ohlc[, 3]-ohlc[, 1])*(ohlc[, 3]-0.5*(ohlc[, 4] + ohlc[, 1]))},
    "garman_klass_yz"={(ohlc[, 1]-rutils::lagit(ohlc[, 4]))^3 +
        0.5*(ohlc[, 2]-ohlc[, 3])^3 -
        (2*log(2)-1)*(ohlc[, 4]-ohlc[, 1])^3},
    "yang_zhang"={c_o <- ohlc[, 1]-rutils::lagit(ohlc[, 4]);
    o_c <- ohlc[, 1]-ohlc[, 4];
    (c_o-sum(c_o)/NROW(c_o))^3 +
      0.67*(o_c-sum(o_c)/NROW(o_c))^3 +
      0.33*((ohlc[, 2]-ohlc[, 4])*(ohlc[, 2]-ohlc[, 1]) +
              (ohlc[, 3]-ohlc[, 4])*(ohlc[, 3]-ohlc[, 1]))}
  )  # end switch
  
  skewv <- skewv/c(1, diff(xts::.index(ohlc)))^3
  skewv[1, ] <- 0
  skewv <- rutils::na_locf(skewv)
  # Colnames(skewv) <- paste0(symbol, ".Skew")
  return(skewv)
  
}  # end ohlc_skew




##########################################################################
#' Calculate time series of point Sharpe-like statistics for each row of a
#' \emph{OHLC} time series.
#'
#' @param \code{ohlc} An \emph{OHLC} time series of prices in \emph{xts} format.
#' 
#' @param \code{method} A \emph{character} string representing method for
#'   estimating the Sharpe-like exponent.
#'
#' @return An \emph{xts} time series with the same number of rows as the
#'   argument \code{ohlc}.
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
#' 
#' @export

ohlc_sharpe <- function(ohlc, method="close") {
  
  sharper <- switch(method,
   "close"={(ohlc[, 4]-ohlc[, 1])/(ohlc[, 2]-ohlc[, 3])},
   "method2"={(ohlc[, 4]-ohlc[, 1])/(ohlc[, 2]-ohlc[, 3])}
  )  # end switch
  sharper <- ifelse(ohlc[, 2]==ohlc[, 3], 0, sharper)
  # Colnames(sharper) <- paste0(rutils::get_name(colnames(ohlc)[1]), ".Sharpe")
  return(sharper)
  
}  # end ohlc_sharpe




##########################################################################
#' Calculate the aggregation (weighted average) of a statistical estimator over
#' a \emph{OHLC} time series using \code{R} code.
#'
#' @param \code{ohlc} An \emph{OHLC} time series of prices and trading volumes,
#'   in \emph{xts} format.
#'
#' @param \code{calc_bars} A \emph{character} string representing a function
#'   for calculating statistics for individual \emph{OHLC} bars.
#'
#' @param \code{weighted} \emph{Boolean} argument: should estimate be weighted
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
#' varv <- agg_stats_r(ohlc=HighFreq::SPY, calc_bars="ohlc_variance")
#' # Calculate time series of daily skew estimates for SPY
#' skew_daily <- apply.daily(x=HighFreq::SPY, FUN=agg_stats_r, calc_bars="ohlc_skew")
#' 
#' @export

agg_stats_r <- function(ohlc, calc_bars="ohlc_variance", weighted=TRUE, ...) {
  
# Match "calc_bars" with moment function
  calc_bars <- match.fun(calc_bars)
  aggvar <- calc_bars(ohlc, ...)
  
# Weight the estimates by volume
  if (weighted) {
    aggvar <- ohlc[, 5]*aggvar
    aggvar <- sum(aggvar)/sum(ohlc[, 5])
  } else
    aggvar <- sum(aggvar)/NROW(aggvar)
  
  return(aggvar)
  
}  # end agg_stats_r




##########################################################################
#' Calculate the volume-weighted average price of an \emph{OHLC} time series
#' over a rolling look-back interval.
#'
#' Performs the same operation as function \code{VWAP()} from package
#' \href{https://cran.r-project.org/web/packages/TTR/index.html}{TTR},
#' but using vectorized functions, so it's a little faster.
#'
#' @param \code{ohlc} An \emph{OHLC} time series of prices in \emph{xts} format.
#' 
#' @param \code{close} A time series of close prices.
#' 
#' @param \code{lookb} The size of the look-back interval, equal to the number of 
#'   rows of data used for calculating the average price.
#'
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{ohlc}.
#'
#' @details The function \code{roll_vwap()} calculates the volume-weighted
#'   average closing price, defined as the sum of the prices multiplied by
#'   trading volumes in the look-back interval, divided by the sum of trading
#'   volumes in the interval. If the argument \code{close} is passed in explicitly,
#'   then its volume-weighted average value over time is calculated.
#'
#' @examples
#' # Calculate and plot rolling volume-weighted average closing prices (VWAP)
#' priceroll <- roll_vwap(ohlc=HighFreq::SPY["2013-11"], lookb=11)
#' chart_Series(HighFreq::SPY["2013-11-12"], name="SPY prices")
#' add_TA(priceroll["2013-11-12"], on=1, col="red", lwd=2)
#' legend("top", legend=c("SPY prices", "VWAP prices"),
#' bg="white", lty=c(1, 1), lwd=c(2, 2),
#' col=c("black", "red"), bty="n")
#' # Calculate running returns
#' returns_running <- ohlc_returns(xtsv=HighFreq::SPY)
#' # Calculate the rolling volume-weighted average returns
#' roll_vwap(ohlc=HighFreq::SPY, close=returns_running, lookb=11)
#' 
#' @export

roll_vwap <- function(ohlc, close=ohlc[, 4, drop=FALSE], lookb) {
  
  vwapr <- rutils::roll_sum(xtsv=close*ohlc[, 5, drop=FALSE], lookb=lookb)
  volume_rolling <- rutils::roll_sum(xtsv=ohlc[, 5, drop=FALSE], lookb=lookb)
  # vwapr <- HighFreq::roll_sum(timeser=close*ohlc[, 5, drop=FALSE], lookb=lookb)
  # volume_rolling <- HighFreq::roll_sum(timeser=ohlc[, 5, drop=FALSE], lookb=lookb)
  vwapr <- ifelse(volume_rolling > 0, vwapr/volume_rolling, 0)
  # vwapr[is.na(vwapr)] <- 0
  # Colnames(vwapr) <- paste0(rutils::get_name(colnames(ohlc)[1]), ".VWAP")
  # Colnames(vwapr) <- colnames(ohlc)
  return(vwapr)
  
}  # end roll_vwap




##########################################################################
#' Calculate a vector of statistics over an \emph{OHLC} time series, and
#' calculate a rolling mean over the statistics.
#'
#' @param \code{ohlc} An \emph{OHLC} time series of prices and trading volumes, in
#'   \emph{xts} format.
#'   
#' @param \code{calc_stats} The name of the function for estimating statistics of a single
#'   row of \emph{OHLC} data, such as volatility, skew, and higher moments.
#'   
#' @param \code{lookb} The size of the look-back interval, equal to the number of
#'   rows of data used for calculating the rolling mean.
#'   
#' @param \code{weighted} \emph{Boolean} argument: should statistic be weighted by
#'   trade volume? (default \code{TRUE})
#'   
#' @param ... additional parameters to the function \code{calc_stats}.
#'
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{ohlc}.
#'
#' @details The function \code{roll_stats()} calculates a vector of statistics
#'   over an \emph{OHLC} time series, such as volatility, skew, and higher
#'   moments.  The statistics could also be any other aggregation of a single
#'   row of \emph{OHLC} data, for example the \emph{High} price minus the
#'   \emph{Low} price squared.  The length of the vector of statistics is equal
#'   to the number of rows of the argument \code{ohlc}. Then it calculates a
#'   trade volume weighted rolling mean over the vector of statistics over and
#'   calculate statistics.
#'
#' @examples
#' # Calculate time series of rolling variance and skew estimates
#' varoll <- roll_stats(ohlc=HighFreq::SPY, lookb=21)
#' skewroll <- roll_stats(ohlc=HighFreq::SPY, calc_stats="ohlc_skew", lookb=21)
#' skewroll <- skewroll/(varoll)^(1.5)
#' skewroll[1, ] <- 0
#' skewroll <- rutils::na_locf(skewroll)
#' 
#' @export

roll_stats <- function(ohlc, calc_stats="ohlc_variance", lookb=11, weighted=TRUE, ...) {
  
# Match "calc_stats" with moment function
  calc_stats <- match.fun(calc_stats)
  aggvar <- calc_stats(ohlc, ...)
  
# Weight by volume
  if (weighted) {
    aggvar <- ohlc[, 5]*aggvar
    volume_rolling <- rutils::roll_sum(ohlc[, 5], lookb=lookb)
    aggvar <- rutils::roll_sum(aggvar, lookb=lookb)/volume_rolling
    aggvar[is.na(aggvar)] <- 0
  } else
    aggvar <- rutils::roll_sum(aggvar, lookb=lookb)/lookb
  # Colnames(aggvar) <- paste(rutils::get_name(colnames(ohlc)[1]), "Vol", sep=".")
  
  return(aggvar)
  
}  # end roll_stats




##########################################################################
#' Calculate the variance of an \emph{OHLC} time series, using different range
#' estimators for variance.
#'
#' @param \code{ohlc} An \emph{OHLC} time series of prices in \emph{xts} format.
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
#' @param \code{scalit} \emph{Boolean} argument: should the returns be divided by the
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
#'   If \code{scalit} is \code{TRUE} (the default), then the returns are divided
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
#' HighFreq::calc_var_ohlc_r(HighFreq::SPY, scalit=FALSE)
#' 
#' @export

calc_var_ohlc_r <- function(ohlc, method="yang_zhang", scalit=TRUE) {
  
  # ohlc <- log(ohlc[, 1:4])
  nrows <- NROW(ohlc)
  
  # Define the time index for scaling the returns
  if (scalit)
    timed <- c(1, diff(xts::.index(ohlc)))
  else
    timed <- rep(1, nrows)
  
  # Coerce ohlc to matrix
  if (is.xts(ohlc))
    ohlc <- coredata(ohlc)
  
  # Calculate all the different intra-day and day-over-day returns 
  # (differences of OHLC prices)
  retc <- rutils::diffit(ohlc[, 4])/timed
  reton <- (ohlc[, 1]-rutils::lagit(ohlc[, 4]))/timed
  retd <- (ohlc[, 4]-ohlc[, 1])/timed
  closhi <- (ohlc[, 4]-ohlc[, 2])/timed
  closlo <- (ohlc[, 4]-ohlc[, 3])/timed
  hilo <- (ohlc[, 2]-ohlc[, 3])/timed
  hiop <- (ohlc[, 2]-ohlc[, 1])/timed
  lowop <- (ohlc[, 3]-ohlc[, 1])/timed
  
  switch(method,
         "close"={var(retc)},
         "rogers_satchell"={-sum(closhi*hiop + closlo*lowop)/nrows},
         "garman_klass"={sum(0.5*hilo^2 - (2*log(2)-1)*retd^2)/nrows},
         "garman_klass_yz"={sum(0.5*hilo^2 - (2*log(2)-1)*retd^2)/nrows + 
             var(reton)},
         "yang_zhang"={coeff <- 0.34/(1.34 + (nrows+1)/(nrows-1))
         var(reton) + coeff*var(retd) +
           (coeff-1)*sum(closhi*hiop + closlo*lowop)/nrows}
  )  # end switch
  
}  # end calc_var_ohlc_r




##########################################################################
#' Calculate a time series of Sharpe ratios over a rolling look-back interval
#' for an \emph{OHLC} time series.
#'
#' @param \code{ohlc} An \emph{OHLC} time series of prices in \emph{xts} format.
#' 
#' @param \code{lookb} The size of the look-back interval, equal to the number of
#'   rows of data used for aggregating the \emph{OHLC} prices.
#'
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{ohlc}.
#'
#' @details The function \code{roll_sharpe()} calculates the rolling Sharpe 
#'   ratio defined as the ratio of percentage returns over the look-back
#'   interval, divided by the average volatility of percentage returns.
#'
#' @examples
#' # Calculate rolling Sharpe ratio over SPY
#' sharper <- roll_sharpe(ohlc=HighFreq::SPY, lookb=11)
#' 
#' @export

roll_sharpe <- function(ohlc, lookb=11) {
  
  retp <- ohlc_returns(ohlc, lag=lookb, scalit=FALSE)
  varoll <- sqrt(HighFreq::roll_var_ohlc(ohlc, lookb=lookb, scalit=FALSE))
  sharper <- ifelse(varoll == 0, 1.0, retp/varoll)
  # colnames(sharper) <- paste0(rutils::get_name(colnames(ohlc)[1]), ".Sharpe")
  return(rutils::na_locf(sharper))
  
}  # end roll_sharpe




##########################################################################
#' Calculate a time series of \emph{Hurst} exponents over a rolling look-back
#' interval.
#'
#' @param \code{ohlc} An \emph{OHLC} time series of prices in \emph{xts} format.
#' 
#' @param \code{lookb} The size of the look-back interval, equal to the number
#'   of rows of data used for aggregating the \emph{OHLC} prices.
#'
#' @return An \emph{xts} time series with a single column and the same number of
#'   rows as the argument \code{ohlc}.
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
#' hurstr <- roll_hurst(ohlc=HighFreq::SPY["2009-03"], lookb=11)
#' chart_Series(hurstr["2009-03-10/2009-03-12"], name="SPY hurstr")
#' 
#' @export

roll_hurst <- function(ohlc, lookb=11) {
  
  rangev <- c(rep(0, lookb-1), (RcppRoll::roll_max(x=ohlc[, 2], n=lookb) +
               RcppRoll::roll_max(x=-ohlc[, 3], n=lookb)))
  varoll <- sqrt(HighFreq::roll_var_ohlc(ohlc, lookb=lookb, scalit=FALSE))
  hurstr <- ifelse((varoll==0) | (rangev==0), 0.5, log(rangev/varoll)/log(lookb))
  # Colnames(hurstr) <- paste0(rutils::get_name(colnames(ohlc)[1]), ".Hurst")
  return(rutils::na_locf(hurstr))
  
}  # end roll_hurst




##########################################################################
#' Apply an aggregation function over a rolling look-back interval and the end
#' points of an \emph{OHLC} time series, using \code{R} code.
#'
#' @param \code{xtsv} An \emph{OHLC} time series of prices and trading volumes,
#'   in \emph{xts} format.
#'   
#' @param \code{agg_fun} The name of the aggregation function to be applied over
#'   a rolling look-back interval.
#'   
#' @param \code{lookb} The number of end points in the look-back interval used
#'   for applying the aggregation function (including the current row).
#'   
#' @param \code{by_columns} \emph{Boolean} argument: should the function
#'   \code{agg_fun()} be applied column-wise (individually), or should it be
#'   applied to all the columns combined? (default is \code{FALSE})
#'   
#' @param \code{out_xts} \emph{Boolean} argument: should the output be coerced
#'   into an \emph{xts} series? (default is \code{TRUE})
#'   
#' @param \code{endpoints} An integer vector of end points.
#' 
#' @param ... additional parameters to the function \code{agg_fun}.
#'
#' @return Either an \emph{xts} time series with the number of rows equal to the
#'   length of argument \code{endpoints}, or a list the length of argument
#'   \code{endpoints}.
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
#'   depending on the argument \code{endpoints}.
#'   If the argument \code{endpoints} isn't explicitly passed to 
#'   \code{roll_apply()}, then the default value is used, and 
#'   \code{roll_apply()} performs aggregations over overlapping intervals at
#'   each point in time.
#'   
#'   If the argument \code{endpoints} is explicitly passed to 
#'   \code{roll_apply()}, then \code{roll_apply()} performs aggregations over 
#'   intervals attached at the endpoints.  If lookb=2 then the aggregations
#'   are performed over non-overlapping intervals, otherwise they are performed
#'   over overlapping intervals.
#'
#'   If the argument \code{out_xts} is \code{TRUE} (the default) then the output
#'   is coerced into an \emph{xts} series, with the number of rows equal to the 
#'   length of argument \code{endpoints}.  Otherwise a list is returned, with
#'   the length equal to the length of argument \code{endpoints}.
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
#' ohlc <- HighFreq::SPY["2012-02-13"]
#' interval <- 11  # number of data points between end points
#' lookb <- 4  # number of end points in look-back interval
#' # Calculate the rolling sums of ohlc columns over a rolling look-back interval
#' aggvar <- roll_apply(ohlc, agg_fun=sum, lookb=lookb, by_columns=TRUE)
#' # Apply a vector-valued aggregation function over a rolling look-back interval
#' agg_function <- function(ohlc)  c(max(ohlc[, 2]), min(ohlc[, 3]))
#' aggvar <- roll_apply(ohlc, agg_fun=agg_function, lookb=lookb)
#' # Define end points at 11-minute intervals (HighFreq::SPY is minutely bars)
#' endpoints <- rutils::endpoints(ohlc, interval=interval)
#' # Calculate the sums of ohlc columns over endpoints using non-overlapping intervals
#' aggvar <- roll_apply(ohlc, agg_fun=sum, endpoints=endpoints, by_columns=TRUE)
#' # Apply a vector-valued aggregation function over the endpoints of ohlc
#' # using overlapping intervals
#' aggvar <- roll_apply(ohlc, agg_fun=agg_function,
#'                             lookb=5, endpoints=endpoints)
#'                             
#' @export

roll_apply <- function(xtsv, agg_fun, lookb=2, endpoints=seq_along(xtsv), 
                       by_columns=FALSE, out_xts=TRUE, ...) {
  
  # Match "agg_fun" with some aggregation function
  agg_fun <- match.fun(agg_fun)
  nrows <- NROW(endpoints)
  # Define startp as lag of endpoints
  startp <- c(rep_len(1, lookb-1), endpoints[1:(nrows-lookb+1)])
  # Define list of look-back intervals for aggregations over past
  lookbs <- lapply(seq_along(endpoints), function(endp) {
    startp[endp]:endpoints[endp]
  })  # end lapply
  # Perform aggregations over length of endpoints
  if (by_columns) {
    # Perform individual aggregations by columns
    aggvar <- lapply(xtsv, function(colnum)
      lapply(lookbs, function(lookb)
        agg_fun(xtsv[lookb], ...)
      ))  # end lapply
  } else {  # not by_columns
    aggvar <- lapply(lookbs, function(lookb)
      agg_fun(xtsv[lookb], ...)
    )  # end lapply
  }  # end if
  
  if (out_xts) {
    # Coerce aggvar into matrix and transpose it
    if (is.null(dim(aggvar)))
      aggvar <- t(aggvar)
    aggvar <- t(aggvar)
    # Coerce aggvar into xts series
    xts(aggvar, order.by=index(xtsv[endpoints]))
  } else
    aggvar
  
}  # end roll_apply




##########################################################################
#' Perform a backtest simulation of a trading strategy (model) over a vector of
#' end points along a time series of prices.
#'
#' @param \code{xtsv} A time series of prices, asset returns, trading volumes,
#'   and other data, in \emph{xts} format.
#'   
#' @param \code{train_func} The name of the function for training (calibrating)
#'   a forecasting model, to be applied over a rolling look-back interval.
#'   
#' @param \code{trade_func} The name of the trading model function, to be
#'   applied over a rolling look-forward interval.
#'   
#' @param \code{lookb} The size of the look-back interval, equal to the number
#'   of rows of data used for training the forecasting model.
#'   
#' @param \code{look_forward} The size of the look-forward interval, equal to
#'   the number of rows of data used for trading the strategy.
#'   
#' @param \code{endpoints} A vector of end points along the rows of the
#'   \code{xtsv} time series, given as either integers or dates.
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
#'   from the \code{xtsv} time series over the look-back interval, and applies
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
#'   listv returned by function \code{trade_func()}.  The list names are equal
#'   to the \emph{endpoints} dates. 
#'   The number of list elements is equal to the number of \emph{endpoints}
#'   minus two (because the first and last end points can't be included in the
#'   backtest).
#'
#' @examples
#' \dontrun{
#' # Combine two time series of prices
#' prices <- cbind(rutils::etfenv$XLU, rutils::etfenv$XLP)
#' lookb <- 252
#' look_forward <- 22
#' # Define end points
#' endpoints <- rutils::calc_endpoints(prices, look_forward)
#' # Perform back-test
#' back_test <- roll_backtest(endpoints=endpoints,
#'     look_forward=look_forward,
#'     lookb=lookb,
#'     train_func = train_model,
#'     trade_func = trade_model,
#'     model_params = model_params,
#'     trading_params = trading_params,
#'     xtsv=prices)
#' }
#' 
#' @export

roll_backtest <- function(xtsv,
                          train_func, trade_func,
                          lookb=look_forward,
                          look_forward,
                          endpoints=rutils::calc_endpoints(xtsv, look_forward),
                          ...) {
  # Match train and trade function names to functions
  train_func <- match.fun(train_func)
  trade_func <- match.fun(trade_func)

  # Convert endpoints dates to integer
  if (inherits(endpoints, c("Date", "POSIXt"))) {
    endds <- endpoints
    endpoints <- match(endpoints, index(xtsv))
  } else {
    endds <- index(xtsv[endpoints])
  }  # end if

  # Define integer backp and fwdp from integer endpoints
  backp <- endpoints - lookb + 1
  backp[backp < 1] <- 1
  
  fwdp <- endpoints + look_forward
  fwdp[fwdp > NROW(xtsv)] <- NROW(xtsv)

  # Perform backtest over length of endpoints
  backtest_range <- 2:(NROW(endpoints)-1)
  back_test <- lapply(backtest_range, function(endp) {
    trained_model <- train_func(xtsv[backp[endp]:endpoints[endp]], ...)
    trade_func(xtsv[(endpoints[endp]+1):fwdp[endp]], trained_model, ...)
  })  # end lapply

  names(back_test) <- endds[backtest_range]
  return(back_test)
  # Coerce back_test into matrix and transpose it
  # if (is.null(dim(back_test)))
  #   back_test <- t(back_test)
  # back_test <- t(back_test)
  # Coerce back_test into xts series
  # xts(back_test, order.by=index(xtsv[endpoints[backtest_range]]))
  
}  # end roll_backtest




##########################################################################
#' Calculate the intraday seasonality aggregations over a single-column
#' \emph{xts} time series.
#'
#' @param \code{xtsv} A single-column \emph{xts} time series.
#' 
#' @param \code{endp} A vector of \emph{character} strings representing points in
#'   time, of the same length as the argument \code{xtsv}.
#'
#' @return An \emph{xts} time series with mean aggregations over the seasonality
#'   interval.
#'
#' @details The function \code{calc_season()} calculates the mean of values
#'   observed at the same points in time specified by the argument
#'   \code{endp}. 
#'   An example of a daily seasonality aggregation is the average price of a
#'   stock between 9:30AM and 10:00AM every day, over many days. The argument
#'   \code{endp} is passed into function \code{tapply()}, and must be the same
#'   length as the argument \code{xtsv}.
#'
#' @examples
#' # Calculate running variance of each minutely OHLC bar of data
#' xtsv <- ohlc_variance(HighFreq::SPY)
#' # Remove overnight variance spikes at "09:31"
#' endp <- format(index(xtsv), "%H:%M")
#' xtsv <- xtsv[!endp=="09:31", ]
#' # Calculate daily seasonality of variance
#' var_seasonal <- calc_season(xtsv=xtsv)
#' chart_Series(x=var_seasonal, name=paste(colnames(var_seasonal),
#'   "daily seasonality of variance"))
#'   
#' @export

calc_season <- function(xtsv, endp=format(zoo::index(xtsv), "%H:%M")) {
  
# Aggregate the mean
  aggv <- tapply(X=xtsv, INDEX=endp, FUN=mean)
# Coerce from array to named vector
  aggv <- structure(as.vector(aggv), names=names(aggv))
# Coerce to xts
  aggv <- xts(x=aggv, order.by=as.POSIXct(paste(Sys.Date(), names(aggv))))
  colnames(aggv) <- colnames(xtsv)
  aggv
  
}  # end calc_season


