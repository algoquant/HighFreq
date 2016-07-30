################################################
###
### Demos for managing high frequency data using package HighFreq
###
################################################

# Set the time-zone to New_York
Sys.setenv(TZ="America/New_York")
# setwd("C:/Develop/data")
# search()  # get search path
options(digits.secs=6)
options(digits=7)
# suppress spurious timezone warning messages
options(xts_check_TZ=FALSE)
options(stringsAsFactors=FALSE)
options(max.print=80)
rm(list=ls())  # remove all objects

# install HighFreq from github
install.packages("devtools")
library(devtools)
install_github(repo="algoquant/HighFreq")

# install HighFreq using install.packages()
install.packages(pkgs="C:/Develop/R/HighFreq", repos=NULL, type="source")
install.packages(pkgs="C:/Develop/R/HighFreq", repos=NULL, type="source", lib="C:/Users/Jerzy/Downloads")
install.packages(pkgs="C:/Develop/R/HighFreq", repos=NULL, type="source", lib="C:/Users/Jerzy/Documents/R/win-library/3.2")

# load HighFreq and attach the data
library(HighFreq)
# data(hf_data)

# set data directories
data_dir <- "E:/mktdata/sec/"
output_dir <- "E:/output/data/"

# set data directories
data_dir <- "C:/Develop/data/hfreq/src/"
output_dir <- "C:/Develop/data/hfreq/scrub/"


###########
# process single day of data

# load list of symbols from file in cwd
# sym_bols <- read.csv(file="etf_list_hf.csv")
# sym_bols <- sym_bols[[1]]
# the file "hf_data.RData" is part of "HighFreq" package, and contains "sym_bols"
data("hf_data")

# define sym_bol
sym_bol <- "SPY"
inter_val <- "2013-11-11/2013-11-15"

# load a single day of TAQ data
sym_bol <- load(
  file.path(data_dir, 
            paste0(sym_bol, "/2014.05.02.", sym_bol, ".RData")))

### scrub a single day of TAQ data (don't aggregate)
taq_data <- scrub_TAQ(taq_data=get(sym_bol))

# calculate returns
returns_running <- run_returns(x_ts=taq_data)

### scrub and aggregate a single day of TAQ data to OHLC
ohlc_data <- scrub_agg(taq_data=get(sym_bol))
chartSeries(ohlc_data, name=sym_bol, theme=chartTheme("white"))


###########
# process TAQ data using package HighFreq: load TAQ data, aggregate to OHLC, and save to file

# aggregate TAQ data for a single symbol, and save to file
save_scrub_agg(sym_bol, data_dir=data_dir, output_dir=output_dir, period="15 min")

save_TAQ(sym_bol, data_dir=data_dir, output_dir=output_dir)

# calculate returns for a single symbol, and save to file
save_rets(sym_bol, data_dir=data_dir, output_dir=output_dir, period="15 min")

# aggregate data for list of symbols, and save to multiple files
sapply(head(sym_bols), save_scrub_agg, data_dir=data_dir, output_dir=output_dir, period="15 min")

# calculate returns for list of symbols, and save to file
sapply(head(sym_bols), save_rets, data_dir=data_dir, output_dir=output_dir, period="15 min")

# load processed OHLC data for a single symbol
# load(file=paste0(sym_bol, ".RData"))
load(file.path(output_dir, paste0(sym_bol, ".RData")))
# load(file="SPY.RData")
# plot OHLC data
chartSeries(get(sym_bol), name=sym_bol, theme=chartTheme("white"))


###########
### estimating volatility, skewness, and kurtosis using package TTR

# function volatility() from package TTR
vol_at <- volatility(OHLC=get(sym_bol)[inter_val], 
                     calc="yang.zhang", n=10)
vol_at <- volatility(OHLC=get(sym_bol)[inter_val], 
                     calc="rogers.satchell", n=10)
vol_at <- volatility(OHLC=get(sym_bol)[inter_val], 
                     calc="garman.klass", n=10)


# estimating rolling aggregations and moments using package TTR
vol_at <- runMedian(x=get(sym_bol)[inter_val], n=100)
vol_at <- runSD(x=get(sym_bol)[inter_val], n=100)
chart_xts(vol_at)


###########
# estimating rolling moments using package HighFreq
library(HighFreq)

win_dow <- 10

### daily open to close variance and skew
daily_volume <- xts::apply.daily(x=get(sym_bol)[, 5], FUN=sum)
colnames(daily_volume) <- paste0(na_me(get(sym_bol)), ".Volume")
var_daily <- xts::apply.daily(x=get(sym_bol), FUN=agg_regate, mo_ment="run_variance", calc_method="rogers.satchell")
colnames(var_daily) <- paste0(na_me(get(sym_bol)), ".Var")
skew_daily <- xts::apply.daily(x=get(sym_bol), FUN=agg_regate, mo_ment="run_skew")
skew_daily <- skew_daily/(var_daily)^(1.5)
colnames(skew_daily) <- paste0(na_me(get(sym_bol)), ".Skew")

# daily Sharpe
sharpe_daily <- xts::apply.daily(x=get(sym_bol), FUN=agg_regate, mo_ment="run_sharpe")
colnames(sharpe_daily) <- paste0(na_me(get(sym_bol)), ".Sharpe")
chart_Series(sharpe_daily[inter_val], name=paste(sym_bol, "Sharpe"))

# simple autocorrelation
returns_running <- run_returns(x_ts=get(sym_bol))
hurst_daily <- returns_running*lag(returns_running)/returns_running^2
hurst_daily[1, ] <- 0
hurst_daily[is.nan(hurst_daily), ] <- 0
hurst_daily <- xts::apply.daily(x=hurst_daily, FUN=sum)

library(PerformanceAnalytics)
hurst_daily <- xts::apply.daily(x=hurst_daily, FUN=HurstIndex)

hurst_daily <- xts::apply.daily(x=hurst_daily, FUN=sum)

colnames(hurst_daily) <- paste0(sym_bol, ".autocorr")

x11(width=6, height=4)
inter_val <- "2013-06-01/"
chart_Series(var_daily[inter_val], name=paste(sym_bol, "variance"))
chart_Series(rutils::roll_sum(hurst_daily, win_dow=10)[-(1:10)]/10, name=paste(sym_bol, "Hurst"))
abline(h=0.5, col="blue", lwd=2)

# minutely returns
returns_running <- run_returns(x_ts=get(sym_bol)[inter_val])

# rolling over end_points
n_row <- NROW(returns_running)
num_agg <- n_row %/% 10
end_points <- c(0, n_row-10*num_agg + 10*(0:num_agg))

# rolling prices
prices_rolling <- rutils::roll_sum(get(sym_bol)[, 1], win_dow=win_dow)/win_dow
colnames(prices_rolling) <- paste0(sym_bol, ".Rets")
chart_Series(get(sym_bol)["2013-11-12", ], name=paste(sym_bol, "Prices"))
add_TA(prices_rolling["2013-11-12"], on=1, col="red", lwd=2)

### rolling volume-weighted returns
returns_rolling <- roll_vwap(oh_lc=SPY, x_ts=returns_running, win_dow=win_dow)
returns_rolling <- returns_rolling[end_points, ]

### calculate SPY variance without overnight jumps
var_running <- run_variance(oh_lc=SPY, calc_method="rogers.satchell")
# calculate rolling SPY variance
var_rolling <- roll_vwap(oh_lc=SPY, x_ts=var_running, win_dow=win_dow)
var_rolling <- var_rolling[end_points, ]

### calculate running and rolling SPY skew
skew_running <- run_skew(oh_lc=SPY)
skew_rolling <- roll_vwap(oh_lc=SPY, x_ts=skew_running, win_dow=win_dow)

# calculate daily seasonality of skew
skew_seasonal <- season_ality(skew_running)

inter_val <- "2014-03-01/"
chart_Series(var_rolling[inter_val], name=colnames(var_rolling))


### calculate running and rolling Sharpe ratio
sharpe_running <- run_sharpe(oh_lc=SPY)
sharpe_rolling <- roll_vwap(oh_lc=SPY, x_ts=sharpe_running, win_dow=win_dow)


### calculate rolling Sharpe Hurst exponent using ratio of range variance estimators

hurst_rolling <- roll_hurst(oh_lc=SPY, win_dow=win_dow)
# calculate hurst_rolling by applying roll_hurst() over argument off_set
hurst_rolling <- lapply(0:(win_dow-1), roll_hurst, oh_lc=SPY, win_dow=win_dow)
hurst_rolling <- rutils::do_call_rbind(hurst_rolling)
# merge and remove duplicate trailing values
dupli_cated <- hurst_rolling[duplicated(index(hurst_rolling)), ]
dupli_cates <- rutils::do_call_rbind(lapply(unique(index(dupli_cated)), 
                                            function(dupli_cate) 
                                              xts(mean(hurst_rolling[dupli_cate, ]), 
                                                  order.by=dupli_cate)))
dupli_cated <- index(hurst_rolling) == unique(index(dupli_cated))
hurst_rolling <- hurst_rolling[!dupli_cated, ]
hurst_rolling <- rbind(hurst_rolling, dupli_cates)
tail(hurst_rolling, 12)
dim(hurst_rolling)

# remove daily warmup periods
hurst_rolling <- hurst_rolling["T09:41:00/T16:00:00"]

# remove stub periods with low frequency of observations
# calculate frequency of time periods
period_freq <- table(in_dex)
period_freq <- structure(as.vector(period_freq), names=names(period_freq))
# remove stub periods
period_freq <- period_freq[period_freq > mean(period_freq)]
period_freq <- period_freq[-length(period_freq)]
hurst_rolling <- hurst_rolling[in_dex %in% names(period_freq), ]
chart_Series(x=hurst_rolling["2012-02-13"], 
             name=paste(colnames(hurst_rolling), "10-minute aggregations"))

### rolling autocorrelation
roll_autocorr <- returns_running*lag(returns_running)
roll_autocorr[1, ] <- 0
roll_autocorr <- rutils::roll_sum(roll_autocorr[, 1], win_dow=10)/10
roll_autocorr[1:10, ] <- 0
colnames(roll_autocorr) <- paste0(sym_bol, ".autocorr")
roll_autocorr <- roll_autocorr[end_points, ]


### design matrix called SPY_design containing columns of aggregations

win_dow <- 5
returns_running <- run_returns(x_ts=SPY)
returns_rolling <- roll_vwap(oh_lc=SPY, x_ts=returns_running, win_dow=win_dow)
colnames(returns_running) <- "returns"
colnames(returns_rolling) <- "returns.WA5"
SPY_design <- cbind(returns_running, returns_rolling)

var_running <- run_variance(oh_lc=SPY)
var_running <- run_variance(oh_lc=SPY, calc_method="rogers.satchell")
var_rolling <- roll_vwap(oh_lc=SPY, x_ts=var_running, win_dow=win_dow)
colnames(var_running) <- "variance"
colnames(var_rolling) <- "variance.WA5"
SPY_design <- cbind(SPY_design, var_running, var_rolling)

skew_running <- run_skew(oh_lc=SPY)
skew_rolling <- roll_vwap(oh_lc=SPY, x_ts=skew_running, win_dow=win_dow)
colnames(skew_running) <- "skew"
colnames(skew_rolling) <- "skew.WA5"
SPY_design <- cbind(SPY_design, skew_running, skew_rolling)

sharpe_running <- run_sharpe(oh_lc=SPY)
sharpe_rolling <- roll_vwap(oh_lc=SPY, x_ts=sharpe_running, win_dow=win_dow)
colnames(sharpe_running) <- "sharpe_running"
colnames(sharpe_rolling) <- "sharpe_running.WA5"
SPY_design <- cbind(SPY_design, sharpe_running, sharpe_rolling)

sharpe_rolling <- roll_sharpe(oh_lc=SPY, win_dow=win_dow)
colnames(sharpe_rolling) <- "sharpe_rolling"
x11()
chart_Series(sharpe_rolling["2013-11-12", ], name=paste(sym_bol, "sharpe_rolling"))
SPY_design <- cbind(SPY_design, sharpe_rolling)

hurst_rolling <- roll_hurst(oh_lc=SPY, win_dow=win_dow)
colnames(hurst_rolling) <- "hurst_rolling"
chart_Series(hurst_rolling["2013-11-12", ], name=paste(sym_bol, "hurst_rolling"))
SPY_design <- cbind(SPY_design, hurst_rolling)

SPY_design <- cbind(SPY_design, Vo(SPY))

colnames(SPY_design)
# SPY_design <- cbind(SPY_design[, -1], SPY_design[, "SPY.Volume"])

head(SPY_design["2013-11-12", ], 7)
save(SPY_design, file="C:/Develop/data/SPY_design.RData")
load(file="C:/Develop/data/SPY_design.RData")


### rolling lm using package roll

library(roll)

returns_running <- rutils::lag_xts(SPY_design[, "returns"], k=-1)
tail(cbind(returns_running, SPY_design[, "returns"]))

lm_coef <- function(x, y) {
  lm(x ~ y)$coef
}  # end lm_coef

roll_lm <- rollapply(data = returns, width = 252,
                       FUN = lm_coef, by.column = FALSE,
                       align = "right")


### daily seasonality

# calculate daily seasonality of volume
volume_seasonal <- season_ality(Vo(get(sym_bol)[inter_val]))
colnames(volume_seasonal) <- paste0(na_me(get(sym_bol)), ".volume_seasonal")
season_data <- volume_seasonal

# calculate daily seasonality of variance
var_seasonal <- season_ality(run_variance(oh_lc=get(sym_bol)[inter_val, 1:4], calc_method="rogers.satchell"))
colnames(var_seasonal) <- paste0(na_me(get(sym_bol)), ".var_seasonal")
season_data <- var_seasonal

# calculate daily seasonality of variance
# calculate variance of each minutely OHLC bar of data
x_ts <- run_variance(get("SPY"))
# remove overnight variance spikes at "09:31"
x_ts <- x_ts["T09:32:00/T16:00:00"]
# calculate daily seasonality of variance
var_seasonal <- season_ality(x_ts=x_ts)
chart_Series(x=var_seasonal, 
             name=paste(colnames(var_seasonal), "daily seasonality"))

season_illiquid <- 1e6*sqrt(var_seasonal/volume_seasonal)
foo <- 1e6*ifelse(Vo(get(sym_bol))==0, 0, sqrt(run_variance(oh_lc=get(sym_bol))/Vo(get(sym_bol))))
season_illiquid <- season_ality(foo)
season_illiquid <- season_ality(sqrt(run_variance(oh_lc=get(sym_bol))/Vo(get(sym_bol))))
colnames(season_illiquid) <- paste0(na_me(get(sym_bol)), ".season_illiquid")
season_data <- season_illiquid

# calculate daily seasonality of skew
skew_seasonal <- season_ality(run_skew(oh_lc=get(sym_bol)[inter_val, 1:4], calc_method="rogers.satchell"))
# skew_seasonal <- skew_seasonal/(var_seasonal)^(1.5)
colnames(skew_seasonal) <- paste0(na_me(get(sym_bol)), ".skew_seasonal")
season_data <- skew_seasonal

# calculate daily seasonality of Hurst exponent
hurst_seasonal <- season_ality(hurst_rolling)
colnames(hurst_seasonal) <- paste0(colnames(hurst_rolling), ".seasonal")
# plot without daily warmup period
chart_Series(x=hurst_seasonal[-(1:10), ], 
             name=paste(colnames(hurst_seasonal), "daily seasonality"))
# below is for run_sharpe() which isn't really true Hurst
hurst_seasonal <- season_ality(run_sharpe(oh_lc=get(sym_bol)[inter_val, 1:4]))
hurst_seasonal <- hurst_seasonal[-NROW(hurst_seasonal)]
colnames(hurst_seasonal) <- paste0(rutils::na_me(get(sym_bol)), ".seasonal")
season_data <- hurst_seasonal

season_autocorr <- returns_running*lag(returns_running)
season_autocorr[1, ] <- 0
season_autocorr <- returns_running*lag(rutils::roll_sum(returns_running, win_dow=5)/5)
season_autocorr[1:5, ] <- 0
season_autocorr <- season_ality(season_autocorr)
colnames(season_autocorr) <- paste0(na_me(get(sym_bol)), ".season_autocorr")
season_data <- rutils::roll_sum(season_autocorr, win_dow=5)/5

# daily Hurst exponents
hurst_daily <- xts::apply.daily(x=SPY, FUN=agg_regate, mo_ment="run_sharpe")
hurst_daily <- xts::apply.daily(x=SPY, FUN=function(x, ...) abs(agg_regate(oh_lc=x, ...)), mo_ment="run_sharpe")
colnames(hurst_daily) <- paste0(rutils::na_me(get(sym_bol)), ".Hurst.daily")
chart_Series(roll_sum(hurst_daily, 10)[-(1:10)]/10, name=paste(sym_bol, "Hurst"))
abline(h=0.5, col="blue", lwd=2)


## aggregate SPY to 10-minute bars
na_mes <- colnames(get(sym_bol))
SPY <- to.minutes10(get(sym_bol))
colnames(SPY) <- col_names
head(SPY)

## plot
season_data <- season_data[-(NROW(season_data))]
x11(width=6, height=4)
chart_Series(x=season_data, 
             name=paste(colnames(season_data), "daily seasonality"))
# or, adjust y-axis range
rutils::chart_xts(season_data,
          name=paste(colnames(season_data), "daily seasonality"), 
          ylim=c(min(season_data), max(season_data)/2))
# or, adjust y-axis range by hand
plot_theme <- chart_theme()
plot_theme$format.labels <- "%H:%M"
ch_ob <- chart_Series(x=season_data, 
                      name=paste(colnames(season_data), "daily seasonality"), 
                      theme=plot_theme, plot=FALSE)
y_lim <- ch_ob$get_ylim()
y_lim[[2]] <- structure(c(y_lim[[2]][1], y_lim[[2]][2]/2), fixed=TRUE)
ch_ob$set_ylim(y_lim)
plot(ch_ob)


## regressions
da_ta <- cbind(roll_autocorr, var_rolling)
da_ta <- cbind(foo, var_rolling)

da_ta <- cbind(hurst_daily, var_daily)
da_ta <- cbind(hurst_seasonal, var_seasonal)
da_ta <- cbind(hurst_daily, skew_daily)

for_mula <- as.formula(paste(colnames(da_ta)[1], 
                             paste(paste(colnames(da_ta)[-1], 
                                         collapse=" + ")), sep="~"))
for_mula <- as.formula(paste(colnames(da_ta)[1], 
                             paste(paste(colnames(da_ta)[-1], 
                                         collapse=" + "), "- 1"), sep="~"))
l_m <- lm(for_mula, data=da_ta)
summary(l_m)
# scatterplot of da_ta
x11()
plot(for_mula, data=da_ta)
# plot of da_ta
x11()
inter_val <- "2014-04-07/2014-04-10"
ran_ge <- range(da_ta[inter_val, 1])
chart_Series(da_ta[inter_val, 1]/(ran_ge[2]-ran_ge[1]), name=paste(sym_bol, "data"))
ran_ge <- range(da_ta[inter_val, 2])
add_TA(da_ta[inter_val, 2]/(ran_ge[2]-ran_ge[1]), on=1, col="blue", lwd=2)

# rolling volume-weighted skewness
skew_rolling <- roll_moment(oh_lc=get(sym_bol)[inter_val], mo_ment="run_skew", win_dow=10)
skew_rolling <- skew_rolling/(var_rolling)^(1.5)
skew_rolling[1, ] <- 0
skew_rolling <- na.locf(skew_rolling)
tail(skew_rolling, 11)


## Bootstrap of standard errors of all the methods in function run_variance()

set.seed(1121)  # reset random number generator

# create time index of one second intervals for a single day
in_dex <- seq(from=as.POSIXct("2016-01-01 00:00:00"),
              to=as.POSIXct("2016-01-30 00:00:00"), by="1 sec")
# create xts of random prices
returns_running <- rnorm(length(in_dex), sd=0.001)
x_ts <- xts(exp(cumsum(returns_running)), order.by=in_dex)
x_ts <- xts(exp(cumsum(rnorm(length(in_dex), sd=0.001))), order.by=in_dex)
# aggregate to minutes OHLC data
oh_lc <- xts::to.period(x=x_ts, period="minutes")
# calculate variance estimates
vari_ance <- run_variance(oh_lc)
sum(vari_ance)

# calculate variance using all the different estimators in run_variance()
meth_ods <- c("close", "garman.klass", "rogers.satchell", "garman.klass_yz", "yang.zhang")
sapply(meth_ods, function(meth_od) {
  sum(run_variance(oh_lc, calc_method=meth_od))
})

# calculate standard errors of variance estimators using bootstrap
boot_strap <- sapply(1:100, function(x) {
  boot_sample <- sample(returns_running, replace=TRUE)
  x_ts <- xts(exp(cumsum(boot_sample)), order.by=in_dex)
# aggregate to minutes OHLC data
  oh_lc <- xts::to.period(x=x_ts, period="minutes")
# calculate variance estimates
  sapply(meth_ods, function(meth_od) {
    sum(run_variance(oh_lc, calc_method=meth_od))
  })  # end sapply
})  # end sapply

# calculate standard errors of variance estimators using bootstrap - parallel version
boot_strap <- sapply(1:100, function(x) {
  boot_sample <- sample(returns_running, replace=TRUE)
  x_ts <- xts(exp(cumsum(boot_sample)), order.by=in_dex)
  # aggregate to minutes OHLC data
  oh_lc <- xts::to.period(x=x_ts, period="minutes")
  # calculate variance estimates
  sapply(meth_ods, function(meth_od) {
    sum(run_variance(oh_lc, calc_method=meth_od))
  })  # end sapply
})  # end sapply


head(t(boot_strap))
boot_errors <- apply(boot_strap, MARGIN=1, sd)

save(boot_strap, file="C:/Develop/data/boot_strap.RData")

# calculate variance estimates for SPY
vari_ance <- run_variance(SPY)


###########
### benchmark HighFreq::roll_apply() to functions xts:::rollapply.xts(), 
# xts::period.apply(), and PerformanceAnalytics::apply.rolling()

# extract a single day of SPY data
price_s <- SPY["2012-02-13"]
win_dow <- 11
library(PerformanceAnalytics)  # load package PerformanceAnalytics
# define aggregation function that returns a vector of values
agg_function <- function(x_ts)  c(max(x_ts[, 1]), min(x_ts[, 4]))
agg_function <- function(x_ts)  sapply(x_ts, mean)
agg_function(price_s)
# perform aggregations for vector-valued aggregation function
foo <- xts:::rollapply.xts(price_s, FUN=agg_function, width=win_dow, align="right", by.column=FALSE)
bar <- HighFreq::roll_apply(price_s, agg_fun=agg_function, win_dow=win_dow)
# perform aggregations for scalar-valued aggregation function
foo <- xts:::rollapply.xts(price_s, FUN=sum, width=win_dow, align="right")
bar <- HighFreq::roll_apply(price_s, agg_fun=sum, win_dow=win_dow, by_columns=TRUE)

# define end points at 10-minute intervals (SPY is minutely bars)
end_points <- rutils::end_points(price_s, inter_val=win_dow)
foo <- xts::period.apply(price_s[, 1], FUN=sum, INDEX=end_points)
bar <- HighFreq::roll_apply(price_s[, 1], agg_fun=sum, win_dow=2, end_points=end_points)
foo <- xts::period.apply(price_s, FUN=agg_function, INDEX=end_points)
bar <- HighFreq::roll_apply(price_s, agg_fun=agg_function, win_dow=2, end_points=end_points)

# perform aggregations over length of end_points
foo_bar <- PerformanceAnalytics::apply.rolling(price_s, FUN=agg_function, width=win_dow)
head(cbind(foo, bar, foo_bar), 14)

# benchmark the speed of the functionals
library(microbenchmark)
summary(microbenchmark(
  roll_apply=HighFreq::roll_apply(price_s, agg_fun=sum, win_dow=win_dow, by_columns=TRUE),
  rollapply=xts:::rollapply.xts(price_s, FUN=sum, width=win_dow, align="right"), 
  apply.rolling=apply.rolling(price_s, FUN=sum, width=win_dow), 
  times=10))[, c(1, 4, 5)]

summary(microbenchmark(
  roll_apply=HighFreq::roll_apply(price_s, agg_fun=agg_function, win_dow=win_dow),
  rollapply=xts:::rollapply.xts(price_s, FUN=agg_function, width=win_dow, align="right", by.column=FALSE), 
  times=10))[, c(1, 4, 5)]


