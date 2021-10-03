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
# install.packages("devtools")
# library(devtools)
devtools::install_github(repo="algoquant/HighFreq")

# install HighFreq from local directory using install.packages()
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
### process single day of high frequency data

# load list of symbols from file in cwd
# sym_bols <- read.csv(file="etf_list_hf.csv")
# sym_bols <- sym_bols[[1]]
# the file "hf_data.RData" is part of "HighFreq" package, and contains "sym_bols"
data("hf_data")

# define sym_bol
sym_bol <- "SPY"
inter_val <- "2013-11-11/2013-11-15"

# load a single day of TAQ data
sym_bol <- load(file.path(data_dir, 
            paste0(sym_bol, "/2014.05.02.", sym_bol, ".RData")))

### scrub a single day of TAQ data (don't aggregate)
ta_q <- scrub_taq(ta_q=SPY)

# calculate returns from the secondly TAQ data
returns_running <- run_returns(x_ts=ta_q)

### scrub and aggregate a single day of TAQ data to OHLC
oh_lc <- scrub_agg(ta_q=ta_q)
chart_Series(oh_lc, name=sym_bol)


###########
### process TAQ data using package HighFreq: load TAQ data, aggregate to OHLC, and save to file

# aggregate TAQ data to 15-min OHLC bar data, for a single symbol, and save to file
save_scrub_agg(sym_bol, data_dir=data_dir, output_dir=output_dir, period="15 min")

save_taq(sym_bol, data_dir=data_dir, output_dir=output_dir)

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
chart_Series(SPY, name=sym_bol)


###########
### estimating volatility, skewness, and kurtosis using package TTR

# function volatility() from package TTR
vol_at <- volatility(OHLC=SPY[inter_val], 
                     calc="yang_zhang", n=10)
vol_at <- volatility(OHLC=SPY[inter_val], 
                     calc="rogers_satchell", n=10)
vol_at <- volatility(OHLC=SPY[inter_val], 
                     calc="garman_klass", n=10)


# estimating rolling aggregations and moments using package TTR
vol_at <- runMedian(x=SPY[inter_val], n=100)
vol_at <- runSD(x=SPY[inter_val], n=100)
chart_xts(vol_at)


###########
### estimating rolling moments using package HighFreq
library(HighFreq)

look_back <- 10

### daily open to close variance and skew
volume_daily <- xts::apply.daily(x=Vo(SPY), FUN=sum)
colnames(volume_daily) <- paste0(sym_bol, ".Volume")
var_daily <- (6.5*60*60)*xts::apply.daily(x=SPY, FUN=agg_stats_r, calc_bars="roll_variance", method="rogers_satchell")
colnames(var_daily) <- paste0(sym_bol, ".Var")
skew_daily <- xts::apply.daily(x=SPY, FUN=agg_stats_r, calc_bars="run_skew")
skew_daily <- skew_daily/(var_daily)^(1.5)
colnames(skew_daily) <- paste0(sym_bol, ".Skew")

# daily Sharpe
sharpe_daily <- xts::apply.daily(x=SPY, FUN=agg_stats_r, calc_bars="run_sharpe")
colnames(sharpe_daily) <- paste0(sym_bol, ".Sharpe")
chart_Series(sharpe_daily[inter_val], name=paste(sym_bol, "Sharpe"))

# simple autocorrelation
returns_running <- run_returns(x_ts=SPY)
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
chart_Series(rutils::roll_sum(hurst_daily, look_back=10)[-(1:10)]/10, name=paste(sym_bol, "Hurst"))
abline(h=0.5, col="blue", lwd=2)

# calculate minutely returns of OHLC data
returns_running <- run_returns(x_ts=SPY[inter_val])

# rolling over end_points
n_row <- NROW(returns_running)
num_agg <- n_row %/% 10
end_points <- c(0, n_row-10*num_agg + 10*(0:num_agg))

# rolling average prices
look_back <- 10
prices_rolling <- rutils::roll_sum(SPY[, 1], look_back=look_back)/look_back
colnames(prices_rolling) <- paste0(sym_bol, ".Prices")
chart_Series(SPY["2013-11-12"], name=paste(sym_bol, "Prices"))
add_TA(prices_rolling["2013-11-12"], on=1, col="red", lwd=2)
legend("top", legend=c("SPY prices", "average prices"), 
       bg="white", lty=c(1, 1), lwd=c(2, 2), 
       col=c("black", "red"), bty="n")


### rolling volume-weighted returns
returns_rolling <- roll_vwap(oh_lc=SPY, x_ts=returns_running, look_back=look_back)
returns_rolling <- returns_rolling[end_points, ]

### calculate rolling SPY variance without overnight jumps
var_rolling <- (6.5*60*60)*roll_variance(oh_lc=SPY, method="rogers_satchell")
# calculate rolling SPY variance
# var_rolling <- roll_vwap(oh_lc=SPY, x_ts=var_running, look_back=look_back)
var_rolling <- var_rolling[end_points, ]

### calculate running and rolling SPY skew
skew_running <- (6.5*60*60)*run_skew(oh_lc=SPY)
skew_rolling <- roll_vwap(oh_lc=SPY, x_ts=skew_running, look_back=look_back)

# calculate intraday seasonality of skew
skew_seasonal <- season_ality(skew_running)

inter_val <- "2014-03-01/"
chart_Series(var_rolling[inter_val], name=colnames(var_rolling))


### calculate running and rolling Sharpe ratio
sharpe_running <- run_sharpe(oh_lc=SPY)
sharpe_rolling <- roll_vwap(oh_lc=SPY, x_ts=sharpe_running, look_back=look_back)


### calculate rolling Hurst exponents using different methods

# create synthetic OHLC time series of random prices
oh_lc <- HighFreq::random_ohlc(in_dex=zoo::index(SPY["2009-03"]))

# calculate rolling Hurst for random prices using ratio of price range divided by standard deviation of returns
hurst_rolling <- roll_hurst(oh_lc=oh_lc, look_back=look_back)

# remove overnight price jumps from SPY in March 2009
oh_lc <- remove_jumps(SPY["2009-03"])

# legacy code for identifying overnight price jumps
# find biggest close-to-open price jumps
foo <- log(oh_lc[, 1])-rutils::lag_xts(log(oh_lc[, 4]))
foo <- sort(as.vector(foo))
plot(foo, t="l")
bar <- quantile(foo, probs=c(0.997, 0.003))
bar <- which(foo > bar[1] | foo < bar[2])

# calculate rolling Hurst for SPY in March 2009
hurst_rolling <- roll_hurst(oh_lc=oh_lc, look_back=look_back)
hurst_rolling <- roll_hurst(oh_lc=SPY["2009-03"], look_back=look_back)
chart_Series(hurst_rolling["2009-03-10/2009-03-12"], name=paste("SPY", "hurst_rolling"))

# old version of roll_hurst(): calculate Hurst over end points of SPY
hurst_rolling <- roll_hurst(oh_lc=SPY, look_back=11, off_set=0, roll_end_points=TRUE)
# calculate a series of rolling hurst values using argument off_set
hurst_rolling <- lapply(0:9, roll_hurst, oh_lc=SPY, look_back=11, roll_end_points=TRUE)
hurst_rolling <- rutils::do_call_rbind(hurst_rolling)
# remove daily warmup periods
hurst_rolling <- hurst_rolling["T09:41:00/T16:00:00"]

# calculate rolling Hurst using ratio of range variance estimators
# old version of roll_hurst(): calculate hurst_rolling by applying roll_hurst() over argument off_set
hurst_rolling <- lapply(0:(look_back-1), roll_hurst, oh_lc=SPY, look_back=look_back)
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
roll_autocorr <- rutils::roll_sum(roll_autocorr[, 1], look_back=10)/10
roll_autocorr[1:10, ] <- 0
colnames(roll_autocorr) <- paste0(sym_bol, ".autocorr")
roll_autocorr <- roll_autocorr[end_points, ]


### design matrix called SPY_design containing columns of aggregations

look_back <- 5
returns_running <- run_returns(x_ts=SPY)
returns_rolling <- roll_vwap(oh_lc=SPY, x_ts=returns_running, look_back=look_back)
colnames(returns_running) <- "returns"
colnames(returns_rolling) <- "returns.WA5"

var_running <- (6.5*60*60)*run_variance(oh_lc=SPY)
# var_running <- run_variance(oh_lc=SPY, method="rogers_satchell")
# var_diff <- rutils::diff_xts(x_ts=var_running)
var_rolling <- roll_vwap(oh_lc=SPY, x_ts=var_running, look_back=look_back)
colnames(var_running) <- "variance"
# colnames(var_diff) <- "variance.diff"
colnames(var_rolling) <- "variance.WA5"

skew_running <- (6.5*60*60)*run_skew(oh_lc=SPY)
# skew_diff <- rutils::diff_xts(x_ts=skew_running)
skew_rolling <- roll_vwap(oh_lc=SPY, x_ts=skew_running, look_back=look_back)
colnames(skew_running) <- "skew"
# colnames(skew_diff) <- "skew.diff"
colnames(skew_rolling) <- "skew.WA5"

sharpe_running <- run_sharpe(oh_lc=SPY)
sharpe_rolling <- roll_vwap(oh_lc=SPY, x_ts=sharpe_running, look_back=look_back)
colnames(sharpe_running) <- "sharpe_running"
colnames(sharpe_rolling) <- "sharpe_running.WA5"

sharpe_rolling <- roll_sharpe(oh_lc=SPY, look_back=look_back)
# sharpe_diff <- rutils::diff_xts(x_ts=sharpe_rolling)
colnames(sharpe_rolling) <- "sharpe_rolling"
# colnames(sharpe_diff) <- "sharpe.diff"
# x11()
# chart_Series(sharpe_rolling["2013-11-12", ], name=paste("SPY", "sharpe_rolling"))

hurst_rolling <- roll_hurst(oh_lc=SPY, look_back=look_back)
# hurst_diff <- rutils::diff_xts(x_ts=hurst_rolling)
colnames(hurst_rolling) <- "hurst_rolling"
# colnames(hurst_diff) <- "hurst.diff"
chart_Series(hurst_rolling["2009-11-12", ], name=paste("SPY", "hurst_rolling"))

# vol_diff <- rutils::diff_xts(x_ts=Vo(SPY))
# colnames(vol_diff) <- "volume.diff"
# SPY_design <- cbind(SPY_design, Vo(SPY), vol_diff)

SPY_design <- cbind(returns_running, returns_rolling, var_running, 
                    var_rolling, skew_running, sharpe_running, 
                    sharpe_rolling, hurst_rolling)
colnames(SPY_design)
# SPY_design <- cbind(SPY_design[, -1], SPY_design[, "SPY.Volume"])

# select most significant factors plus interaction terms
SPY_design <- cbind(returns_running, returns_rolling, var_running, skew_running, 
                    hurst_rolling, returns_running*var_running, returns_running*skew_running)
colnames(SPY_design) <- c(colnames(SPY_design)[1:4], "hurst", "rets_var", "rets_skew")

head(SPY_design["2013-11-12", ])
save(SPY_design, file="C:/Develop/data/SPY_design.RData")
load("C:/Develop/data/SPY_design.RData")


### intraday seasonality

# calculate intraday seasonality of volume
volume_seasonal <- season_ality(Vo(SPY[inter_val]))
colnames(volume_seasonal) <- paste0(sym_bol, ".volume_seasonal")
season_data <- volume_seasonal

# calculate intraday seasonality of variance
var_seasonal <- season_ality((6.5*60*60)*run_variance(oh_lc=SPY[inter_val, 1:4], method="rogers_satchell"))
colnames(var_seasonal) <- paste0(sym_bol, ".var_seasonal")
season_data <- var_seasonal

# calculate intraday seasonality of variance
# calculate minutely variance of OHLC data scaled to daily units
x_ts <- (6.5*60*60)*run_variance(SPY)
# remove overnight variance spikes at "09:31"
x_ts <- x_ts["T09:32:00/T16:00:00"]
# calculate intraday seasonality of variance
var_seasonal <- season_ality(x_ts=x_ts)
chart_Series(x=var_seasonal, 
             name=paste(colnames(var_seasonal), "intraday seasonality"))

season_illiquid <- 1e6*sqrt(var_seasonal/volume_seasonal)
foo <- 1e6*ifelse(Vo(SPY)==0, 0, sqrt(run_variance(oh_lc=SPY)/Vo(SPY)))
season_illiquid <- season_ality(foo)
season_illiquid <- season_ality(sqrt(run_variance(oh_lc=SPY)/Vo(SPY)))
colnames(season_illiquid) <- paste0(sym_bol, ".season_illiquid")
season_data <- season_illiquid

# calculate intraday seasonality of skew
skew_seasonal <- season_ality(run_skew(oh_lc=SPY[inter_val, 1:4], method="rogers_satchell"))
# skew_seasonal <- skew_seasonal/(var_seasonal)^(1.5)
colnames(skew_seasonal) <- paste0(sym_bol, ".skew_seasonal")
season_data <- skew_seasonal

# calculate intraday seasonality of Hurst exponent
hurst_seasonal <- season_ality(hurst_rolling)
colnames(hurst_seasonal) <- paste0(colnames(hurst_rolling), ".seasonal")
# plot without daily warmup period
chart_Series(x=hurst_seasonal[-(1:10), ], 
             name=paste(colnames(hurst_seasonal), "intraday seasonality"))
# below is for run_hurst() which isn't really true Hurst
hurst_seasonal <- season_ality(run_hurst(oh_lc=SPY[inter_val, 1:4]))
hurst_seasonal <- hurst_seasonal[-NROW(hurst_seasonal)]
colnames(hurst_seasonal) <- paste0(sym_bol, ".seasonal")
season_data <- hurst_seasonal

season_autocorr <- returns_running*lag(returns_running)
season_autocorr[1, ] <- 0
season_autocorr <- returns_running*lag(rutils::roll_sum(returns_running, look_back=5)/5)
season_autocorr[1:5, ] <- 0
season_autocorr <- season_ality(season_autocorr)
colnames(season_autocorr) <- paste0(sym_bol, ".season_autocorr")
season_data <- rutils::roll_sum(season_autocorr, look_back=5)/5

# daily Hurst exponents
hurst_daily <- xts::apply.daily(x=SPY, FUN=agg_stats_r, calc_bars="run_hurst")
hurst_daily <- xts::apply.daily(x=SPY, FUN=function(x, ...) abs(agg_stats_r(oh_lc=x, ...)), calc_bars="run_hurst")
colnames(hurst_daily) <- paste0(sym_bol, ".Hurst.daily")
chart_Series(roll_sum(hurst_daily, 10)[-(1:10)]/10, name=paste(sym_bol, "Hurst"))
abline(h=0.5, col="blue", lwd=2)


## aggregate SPY to 10-minute bars
SPY_10min <- to.minutes10(SPY)
colnames(SPY_10min) <- colnames(SPY)
head(SPY_10min)

## plot
season_data <- season_data[-(NROW(season_data))]
x11(width=6, height=4)
chart_Series(x=season_data, 
             name=paste(colnames(season_data), "intraday seasonality"))
# or, adjust y-axis range
rutils::chart_xts(season_data,
          name=paste(colnames(season_data), "intraday seasonality"), 
          ylim=c(min(season_data), max(season_data)/2))
# or, adjust y-axis range by hand
plot_theme <- chart_theme()
plot_theme$format.labels <- "%H:%M"
ch_ob <- chart_Series(x=season_data, 
                      name=paste(colnames(season_data), "intraday seasonality"), 
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
skew_rolling <- roll_vwap(oh_lc=SPY[inter_val], x_ts=run_skew(oh_lc=SPY[inter_val]), look_back=10)
skew_rolling <- skew_rolling/(var_rolling)^(1.5)
skew_rolling[1, ] <- 0
skew_rolling <- xts:::na.locf.xts(skew_rolling)
tail(skew_rolling, 11)


## Create minutely synthetic OHLC time series of random prices

set.seed(1121)  # reset random number generator

# use function random_ohlc()
# oh_lc <- HighFreq::random_ohlc()

# create minutely synthetic OHLC time series of random prices by hand
# create time index of one second intervals
in_dex <- seq(from=as.POSIXct("2016-01-01 00:00:00"),
              to=as.POSIXct("2016-01-30 00:00:00"), by="1 sec")
# create xts of random prices
vol_at <- 0.001
x_ts <- xts(exp(cumsum(rnorm(length(in_dex), sd=vol_at) - vol_at^2/2)), order.by=in_dex)
colnames(x_ts) <- "random"
chart_Series(x=x_ts["2016-01-10 09/2016-01-10 10"], name="random prices")
# aggregate to minutes OHLC data
oh_lc <- xts::to.period(x=x_ts, period="minutes", name="random")
chart_Series(x=oh_lc["2016-01-10"], name="random OHLC prices")
# add volume
oh_lc <- cbind(oh_lc, sample(x=10*(2:18), size=NROW(oh_lc), replace=TRUE))
colnames(oh_lc)[ 5] <- "random.volume"
tail(oh_lc)


## calculate variance estimates and their standard errors using bootstrap and function roll_variance()

# calculate variance estimates using function roll_variance()
var_rolling <- (6.5*60*60)*roll_variance(oh_lc)
sum(var_rolling)

# calculate variance using all the different estimators in roll_variance()
meth_ods <- c("close", "garman_klass", "rogers_satchell", "garman_klass_yz", "yang_zhang")
sapply(meth_ods, function(meth_od) {
  sum((6.5*60*60)*roll_variance(oh_lc, method=meth_od))
})

# calculate secondly returns from TAQ data
returns_running <- run_returns(x_ts=SPY_TAQ)
returns_running <- returns_running - sum(returns_running)/NROW(returns_running) - drop(var(returns_running))/2
in_dex <- index(SPY_TAQ)

# calculate standard errors of variance estimators using bootstrap
boot_strap <- sapply(1:100, function(x) {
  boot_sample <- sample(returns_running, replace=TRUE)
  x_ts <- xts::xts(exp(cumsum(boot_sample)), order.by=in_dex)
# aggregate to minutes OHLC data
  oh_lc <- xts::to.period(x=x_ts, period="minutes")
# calculate variance estimates
  sapply(meth_ods, function(meth_od) {
    sum(HighFreq::roll_variance(oh_lc, look_back=look_back, method=meth_od))
  })  # end sapply
})  # end sapply

# calculate standard errors of variance estimators using TTR and bootstrap
meth_ods <- c("close", "garman.klass", "rogers.satchell", "gk.yz", "yang.zhang")
boot_strap <- sapply(1:100, function(x) {
  boot_sample <- sample(returns_running, replace=TRUE)
  x_ts <- xts(exp(cumsum(boot_sample)), order.by=in_dex)
  # aggregate to minutes OHLC data
  oh_lc <- xts::to.period(x=x_ts, period="minutes")
  # calculate variance estimates
  sapply(meth_ods, function(meth_od) {
    sum(na.omit(TTR::volatility(oh_lc, n=look_back, N=1, calc=meth_od))^2)
  })  # end sapply
})  # end sapply

head(t(boot_strap))
boot_errors <- apply(boot_strap, MARGIN=1, sd)


# calculate standard errors of variance estimators using bootstrap - parallel version under Windows
library(parallel)  # load package parallel
num_cores <- detectCores() - 1  # number of cores
clus_ter <- makeCluster(num_cores)  # initialize compute cluster

boot_strap <- parLapply(clus_ter, 1:1000, 
                        function(x, re_turns, in_dex, meth_ods) {
  boot_sample <- sample(re_turns, replace=TRUE)
  x_ts <- xts::xts(exp(cumsum(boot_sample)), order.by=in_dex)
  # aggregate to minutes OHLC data
  oh_lc <- xts::to.period(x=x_ts, period="minutes")
  # calculate variance estimates
  sapply(meth_ods, function(meth_od) {
    sum((6.5*60*60)*HighFreq::roll_variance(oh_lc, method=meth_od))
  })  # end sapply
}, re_turns=returns_running, in_dex=in_dex, meth_ods=meth_ods)  # end parLapply

# analyze bootstrapped variance
boot_strap <- do.call(rbind, boot_strap)
head(boot_strap)
sum(is.na(boot_strap))
# means and standard errors from bootstrap
apply(boot_strap, MARGIN=2, 
      function(x) c(mean=mean(x), std_error=sd(x)))

save(boot_strap, file="C:/Develop/data/boot_strap.RData")

stopCluster(clus_ter)

# calculate variance estimates for SPY
var_rolling <- (6.5*60*60)*roll_variance(SPY)


###########
### benchmark HighFreq::roll_apply() to functions xts:::rollapply.xts(), 
# xts::period.apply(), and PerformanceAnalytics::apply.rolling()

# extract a single day of SPY data
price_s <- SPY["2012-02-13"]
look_back <- 11
library(PerformanceAnalytics)  # load package PerformanceAnalytics
# define aggregation function that returns a vector of values
agg_function <- function(x_ts)  c(max(x_ts[, 1]), min(x_ts[, 4]))
agg_function <- function(x_ts)  sapply(x_ts, mean)
agg_function(price_s)
# perform aggregations for vector-valued aggregation function
foo <- xts:::rollapply.xts(price_s, FUN=agg_function, width=look_back, align="right", by.column=FALSE)
bar <- HighFreq::roll_apply(price_s, agg_fun=agg_function, look_back=look_back)
# perform aggregations for scalar-valued aggregation function
foo <- xts:::rollapply.xts(price_s, FUN=sum, width=look_back, align="right")
bar <- HighFreq::roll_apply(price_s, agg_fun=sum, look_back=look_back, by_columns=TRUE)

# define end points at 10-minute intervals (SPY is minutely bars)
end_points <- rutils::calc_endpoints(price_s, inter_val=look_back)
foo <- xts::period.apply(price_s[, 1], FUN=sum, INDEX=end_points)
bar <- HighFreq::roll_apply(price_s[, 1], agg_fun=sum, look_back=2, end_points=end_points)
foo <- xts::period.apply(price_s, FUN=agg_function, INDEX=end_points)
bar <- HighFreq::roll_apply(price_s, agg_fun=agg_function, look_back=2, end_points=end_points)

# perform aggregations over length of end_points
foo_bar <- PerformanceAnalytics::apply.rolling(price_s, FUN=agg_function, width=look_back)
head(cbind(foo, bar, foo_bar), 14)

# benchmark the speed of the functionals
library(microbenchmark)
summary(microbenchmark(
  roll_apply=HighFreq::roll_apply(price_s, agg_fun=sum, look_back=look_back, by_columns=TRUE),
  rollapply=xts:::rollapply.xts(price_s, FUN=sum, width=look_back, align="right"), 
  apply.rolling=apply.rolling(price_s, FUN=sum, width=look_back), 
  times=10))[, c(1, 4, 5)]

summary(microbenchmark(
  roll_apply=HighFreq::roll_apply(price_s, agg_fun=agg_function, look_back=look_back),
  rollapply=xts:::rollapply.xts(price_s, FUN=agg_function, width=look_back, align="right", by.column=FALSE), 
  times=10))[, c(1, 4, 5)]


