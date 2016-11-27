## ----eval=FALSE----------------------------------------------------------
#  # load HighFreq to load SPY data
#  library(HighFreq)
#  sym_bol <- "SPY"
#  # rolling average prices
#  win_dow <- 10
#  prices_rolling <- rutils::roll_sum(Cl(get(sym_bol)), win_dow=win_dow)/win_dow
#  colnames(prices_rolling) <- paste0(sym_bol, ".Prices")
#  chart_Series(get(sym_bol)["2013-11-12"], name=paste(sym_bol, "Prices"))
#  add_TA(prices_rolling["2013-11-12"], on=1, col="red", lwd=2)
#  legend("top", legend=c("SPY prices", "average prices"),
#         bg="white", lty=c(1, 1), lwd=c(2, 2),
#         col=c("black", "red"), bty="n")

## ----eval=FALSE----------------------------------------------------------
#  # load HighFreq
#  library(HighFreq)
#  sym_bol <- "SPY"
#  # rolling volume-weighted average prices
#  win_dow <- 10
#  prices_rolling <- roll_vwap(oh_lc=get(sym_bol)["2013-11-12"], win_dow=win_dow)
#  colnames(prices_rolling) <- paste0(sym_bol, ".Prices")
#  chart_Series(get(sym_bol)["2013-11-12"], name=paste(sym_bol, "VWAP Prices"))
#  add_TA(prices_rolling["2013-11-12"], on=1, col="red", lwd=2)
#  legend("top", legend=c("SPY prices", "VWAP prices"),
#         bg="white", lty=c(1, 1), lwd=c(2, 2),
#         col=c("black", "red"), bty="n")

## ----eval=FALSE----------------------------------------------------------
#  # load HighFreq
#  library(HighFreq)
#  # calculate running variance using method rogers_satchell
#  # scale to daily frequency and also apply factor for secondly units
#  var_running <- (6.5*60*60^2)*run_variance(oh_lc=SPY,
#                                            calc_method="rogers_satchell")
#  # calculate rolling volume-weighted average daily variance
#  win_dow <- 21
#  var_rolling <- roll_vwap(oh_lc=SPY, x_ts=var_running, win_dow=win_dow)
#  colnames(var_rolling) <- paste0(sym_bol, ".Var")
#  
#  # calculate rolling daily variance using roll_variance()
#  var_rolling <- (6.5*60*60^2)*roll_variance(oh_lc=SPY,
#                                            calc_method="rogers_satchell",
#                                            win_dow=win_dow)
#  
#  # calculate rolling volume-weighted average skew indicator
#  skew_running <- run_skew(oh_lc=SPY)
#  skew_rolling <- roll_vwap(oh_lc=SPY, x_ts=skew_running, win_dow=win_dow)
#  skew_rolling <- ifelse(var_rolling>0, skew_rolling/(var_rolling)^(1.5), 0)
#  colnames(skew_rolling) <- paste0(sym_bol, ".Skew")
#  chart_Series(skew_rolling["2013-11-12"], name=paste(sym_bol, "Rolling Skew"))

## ----eval=FALSE----------------------------------------------------------
#  # load HighFreq
#  library(HighFreq)
#  # calculate rolling volume-weighted average variance and skew
#  win_dow <- 21
#  var_rolling <- roll_moment(oh_lc=SPY, win_dow=win_dow)
#  skew_rolling <- roll_moment(oh_lc=SPY, mo_ment="run_skew", win_dow=win_dow)
#  skew_rolling <- ifelse(var_rolling>0, skew_rolling/(var_rolling)^(1.5), 0)
#  chart_Series(skew_rolling["2013-11-12"], name=paste(sym_bol, "Rolling Skew"))

## ----eval=FALSE----------------------------------------------------------
#  # load HighFreq
#  library(HighFreq)
#  sym_bol <- SPY
#  # calculate daily open to close variance and skew
#  var_daily <- (6.5*60*60)*xts::apply.daily(x=get(sym_bol), FUN=agg_regate,
#                                mo_ment="run_variance", calc_method="rogers_satchell")
#  colnames(var_daily) <- paste0(sym_bol, ".Var")
#  skew_daily <- xts::apply.daily(x=get(sym_bol), FUN=agg_regate, mo_ment="run_skew")
#  skew_daily <- skew_daily/(var_daily)^(1.5)
#  colnames(skew_daily) <- paste0(sym_bol, ".Skew")
#  inter_val <- "2013-10/2013-12"
#  chart_Series(skew_daily[inter_val], name=paste(sym_bol, "Skew"))

## ----eval=FALSE----------------------------------------------------------
#  # load HighFreq
#  library(HighFreq)
#  sym_bol <- "SPY"
#  volume_seasonal <- season_ality(Vo(get(sym_bol)))
#  colnames(volume_seasonal) <- paste0(sym_bol, ".volume_seasonal")
#  chart_Series(volume_seasonal, name=paste(sym_bol, "intraday seasonality of volume"))
#  var_seasonal <- season_ality((6.5*60*60)*run_variance(oh_lc=get(sym_bol)))
#  colnames(var_seasonal) <- paste0(sym_bol, ".var_seasonal")
#  chart_Series(var_seasonal, name=paste(sym_bol, "intraday seasonality of variance"))

## ----eval=FALSE----------------------------------------------------------
#  # load HighFreq
#  library(HighFreq)
#  # extract a single day of SPY data
#  x_ts <- SPY["2012-02-13"]
#  win_dow <- 11
#  # calculate the rolling sums of the columns of x_ts
#  agg_regations <- roll_apply(x_ts, agg_fun=sum, win_dow=win_dow, by_columns=TRUE)
#  # define a vector-valued aggregation function
#  agg_function <- function(x_ts)  c(max(x_ts[, 2]), min(x_ts[, 3]))
#  # apply the aggregation function over a rolling window
#  agg_regations <- roll_apply(x_ts, agg_fun=agg_function, win_dow=win_dow)
#  # define end points at 11-minute intervals (SPY is minutely bars)
#  end_points <- rutils::end_points(x_ts, inter_val=win_dow)
#  # calculate the rolling sums of the columns of x_ts over end_points
#  agg_regations <- roll_apply(x_ts, agg_fun=sum, win_dow=2, end_points=end_points, by_columns=TRUE)
#  # apply the vector-valued aggregation function over the end_points of x_ts
#  agg_regations <- roll_apply(x_ts, agg_fun=agg_function, win_dow=2, end_points=end_points)

