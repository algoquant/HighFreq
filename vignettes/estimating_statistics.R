## ----eval=FALSE----------------------------------------------------------
#  # load HighFreq to load SPY data
#  library(HighFreq)
#  # rolling average prices
#  look_back <- 10
#  prices_rolling <- rutils::roll_sum(Cl(HighFreq::SPY), look_back=look_back)/look_back
#  colnames(prices_rolling) <- "SPY.Prices"
#  chart_Series(HighFreq::SPY["2013-11-12"], name="SPY Prices")
#  add_TA(prices_rolling["2013-11-12"], on=1, col="red", lwd=2)
#  legend("top", legend=c("SPY prices", "average prices"),
#         bg="white", lty=c(1, 1), lwd=c(2, 2),
#         col=c("black", "red"), bty="n")

## ----eval=FALSE, echo=(-(1:2))-------------------------------------------
#  # load HighFreq
#  library(HighFreq)
#  # rolling volume-weighted average prices
#  look_back <- 10
#  prices_rolling <- roll_vwap(oh_lc=HighFreq::SPY["2013-11-12"], look_back=look_back)
#  colnames(prices_rolling) <- "SPY.Prices"
#  chart_Series(HighFreq::SPY["2013-11-12"], name="SPY VWAP Prices")
#  add_TA(prices_rolling["2013-11-12"], on=1, col="red", lwd=2)
#  legend("top", legend=c("SPY prices", "VWAP prices"),
#         bg="white", lty=c(1, 1), lwd=c(2, 2),
#         col=c("black", "red"), bty="n")

## ----eval=FALSE, echo=(-(1:2))-------------------------------------------
#  # load HighFreq
#  library(HighFreq)
#  # calculate variance of SPY using method yang_zhang
#  # scale from minutely to daily frequency and also apply factor to compensate for secondly units
#  vari_ance <- (6.5*60*60^2)*HighFreq::calc_variance(HighFreq::SPY, calc_method="yang_zhang")
#  # calculate variance of SPY without accounting for overnight jumps
#  vari_ance <- (6.5*60*60^2)*HighFreq::calc_variance(HighFreq::SPY, calc_method="rogers_satchell")
#  # calcuate daily intraday volatilities
#  var_daily <- (6.5*60*60^2)*period.apply(x=HighFreq::SPY, INDEX=end_days, HighFreq::calc_variance)
#  index(var_daily) <- lubridate::floor_date(index(var_daily), "day")

## ----eval=FALSE, echo=(-(1:2))-------------------------------------------
#  # load HighFreq
#  library(HighFreq)
#  # calculate running variance using method rogers_satchell
#  # scale from minutely to daily frequency and also apply factor to compensate for secondly units
#  var_running <- (6.5*60*60^2)*run_variance(oh_lc=HighFreq::SPY,
#                                            calc_method="rogers_satchell")
#  # calculate rolling volume-weighted average daily variance
#  look_back <- 21
#  var_rolling <- roll_vwap(oh_lc=HighFreq::SPY, x_ts=var_running, look_back=look_back)
#  colnames(var_rolling) <- "SPY.Var"
#  
#  # calculate rolling daily variance using roll_variance()
#  var_rolling <- (6.5*60*60^2)*roll_variance(oh_lc=HighFreq::SPY,
#                                            calc_method="rogers_satchell",
#                                            look_back=look_back)
#  
#  # calculate rolling volume-weighted average skew indicator
#  skew_running <- run_skew(oh_lc=HighFreq::SPY)
#  skew_rolling <- roll_vwap(oh_lc=HighFreq::SPY, x_ts=skew_running, look_back=look_back)
#  skew_rolling <- ifelse(var_rolling>0, skew_rolling/(var_rolling)^(1.5), 0)
#  colnames(skew_rolling) <- "SPY.Skew"
#  chart_Series(skew_rolling["2013-11-12"], name="SPY Rolling Skew")

## ----eval=FALSE, echo=(-(1:2))-------------------------------------------
#  # load HighFreq
#  library(HighFreq)
#  # calculate rolling volume-weighted average variance and skew
#  look_back <- 21
#  var_rolling <- roll_moment(oh_lc=HighFreq::SPY, look_back=look_back)
#  skew_rolling <- roll_moment(oh_lc=HighFreq::SPY, mo_ment="run_skew", look_back=look_back)
#  skew_rolling <- ifelse(var_rolling>0, skew_rolling/(var_rolling)^(1.5), 0)
#  chart_Series(skew_rolling["2013-11-12"], name="SPY Rolling Skew")

## ----eval=FALSE, echo=(-(1:2))-------------------------------------------
#  # load HighFreq
#  library(HighFreq)
#  # calculate daily average open to close variance
#  var_daily <- (6.5*60*60^2)*xts::apply.daily(x=HighFreq::SPY, FUN=agg_regate,
#                                mo_ment="run_variance", calc_method="rogers_satchell")
#  colnames(var_daily) <- "SPY.Var"
#  chart_Series(100*sqrt(var_daily["/2010"]), name="SPY daily standard deviation")
#  
#  # calculate daily average skew
#  skew_daily <- xts::apply.daily(x=HighFreq::SPY, FUN=agg_regate, mo_ment="run_skew")
#  skew_daily <- skew_daily/(var_daily)^(1.5)
#  colnames(skew_daily) <- "SPY.Skew"
#  inter_val <- "2013-10/2013-12"
#  chart_Series(skew_daily[inter_val], name="SPY Skew")

## ----eval=FALSE, echo=(-(1:2))-------------------------------------------
#  # load HighFreq
#  library(HighFreq)
#  volume_seasonal <- season_ality(Vo(HighFreq::SPY))
#  colnames(volume_seasonal) <- "SPY.volume_seasonal"
#  chart_Series(volume_seasonal, name="SPY intraday seasonality of volume")
#  var_seasonal <- season_ality((6.5*60*60^2)*run_variance(oh_lc=HighFreq::SPY))
#  colnames(var_seasonal) <- "SPY.var_seasonal"
#  chart_Series(var_seasonal, name="SPY intraday seasonality of variance")

## ----eval=FALSE, echo=(-(1:2))-------------------------------------------
#  # load HighFreq
#  library(HighFreq)
#  # extract a single day of SPY data
#  x_ts <- SPY["2012-02-13"]
#  look_back <- 11
#  # calculate the rolling sums of the columns of x_ts
#  agg_regations <- roll_apply(x_ts, agg_fun=sum, look_back=look_back, by_columns=TRUE)
#  # define a vector-valued aggregation function
#  agg_function <- function(x_ts)  c(max(x_ts[, 2]), min(x_ts[, 3]))
#  # apply the aggregation function over a rolling window
#  agg_regations <- roll_apply(x_ts, agg_fun=agg_function, look_back=look_back)
#  # define end points at 11-minute intervals (SPY is minutely bars)
#  end_points <- rutils::end_points(x_ts, inter_val=look_back)
#  # calculate the rolling sums of the columns of x_ts over end_points
#  agg_regations <- roll_apply(x_ts, agg_fun=sum, look_back=2, end_points=end_points, by_columns=TRUE)
#  # apply the vector-valued aggregation function over the end_points of x_ts
#  agg_regations <- roll_apply(x_ts, agg_fun=agg_function, look_back=2, end_points=end_points)

