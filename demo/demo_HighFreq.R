################################################
###
### Demos for managing high frequency data using package 'HighFreq'
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
xts_rets <- calc_rets(x_ts=taq_data)

### scrub and aggregate a single day of TAQ data to OHLC
ohlc_data <- scrub_agg(taq_data=get(sym_bol))
chartSeries(ohlc_data, name=sym_bol, theme=chartTheme("white"))


###########
# process TAQ data using package 'HighFreq': load TAQ data, aggregate to OHLC, and save to file

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

## daily open to close variance and skew
daily_volume <- apply.daily(x=get(sym_bol)[, 5], FUN=sum)
colnames(daily_volume) <- paste0(na_me(get(sym_bol)), ".Volume")
daily_var <- apply.daily(x=get(sym_bol), FUN=agg_regate, esti_mator="vari_ance")
colnames(daily_var) <- paste0(na_me(get(sym_bol)), ".Var")
daily_skew <- apply.daily(x=get(sym_bol), FUN=agg_regate, esti_mator="skew_ohlc")
daily_skew <- daily_skew/(daily_var)^(1.5)
colnames(daily_skew) <- paste0(na_me(get(sym_bol)), ".Skew")

# daily autocorrelation or Hurst
foo <- hurst_ohlc(get(sym_bol)["2013-11-12"])
daily_hurst <- apply.daily(x=get(sym_bol), FUN=agg_regate, esti_mator="hurst_ohlc")
colnames(daily_hurst) <- paste0(na_me(get(sym_bol)), ".Hurst")
chart_Series(daily_hurst[inter_val], name=paste(sym_bol, "Hurst"))

re_turns <- calc_rets(x_ts=get(sym_bol))
daily_hurst <- re_turns[, 1]*lag(re_turns[, 1])/re_turns[, 1]^2
daily_hurst[1, ] <- 0
daily_hurst[is.nan(daily_hurst), ] <- 0
daily_hurst <- apply.daily(x=daily_hurst, FUN=sum)

library(PerformanceAnalytics)
daily_hurst <- apply.daily(x=daily_hurst, FUN=HurstIndex)

daily_hurst <- apply.daily(x=daily_hurst, FUN=sum)

colnames(daily_hurst) <- paste0(sym_bol, ".autocorr")

x11(width=6, height=4)
inter_val <- "2013-06-01/"
chart_Series(daily_var[inter_val], name=paste(sym_bol, "variance"))
chart_Series(roll_sum(daily_hurst, win_dow=10)[-(1:10)]/10, name=paste(sym_bol, "Hurst"))
abline(h=0.5, col="blue", lwd=2)

# minutely variance and skew
# minutely returns
re_turns <- calc_rets(x_ts=get(sym_bol)[inter_val])

## rolling over end_points
n_row <- NROW(re_turns)
num_agg <- n_row %/% 10
end_points <- c(0, n_row-10*num_agg + 10*(0:num_agg))

# rolling prices
roll_prices <- roll_sum(get(sym_bol)[, 1], win_dow=10)/10
colnames(roll_prices) <- paste0(sym_bol, ".Rets")
chart_Series(get(sym_bol)["2013-11-12", ], name=paste(sym_bol, "Prices"))
add_TA(roll_prices["2013-11-12"], on=1, col="red", lwd=2)

# rolling returns
roll_returns <- roll_sum(re_turns[, 1], win_dow=10)/10
colnames(roll_returns) <- paste0(sym_bol, ".Rets")
roll_returns[1:10, ] <- 0
roll_returns <- roll_returns[end_points, ]

# rolling volume-weighted variance
roll_var <- roll_agg(ohlc=get(sym_bol), n=10)
roll_var <- roll_var[end_points, ]

# rolling autocorrelation
roll_autocorr <- re_turns[, 1]*lag(re_turns[, 1])
roll_autocorr[1, ] <- 0
roll_autocorr <- roll_sum(roll_autocorr[, 1], win_dow=10)/10
roll_autocorr[1:10, ] <- 0
colnames(roll_autocorr) <- paste0(sym_bol, ".autocorr")
roll_autocorr <- roll_autocorr[end_points, ]

inter_val <- "2014-03-01/"
chart_Series(roll_var[inter_val], name=colnames(roll_var))

# define end points at 10-minute intervals
inter_val <- 10
# rolling Hurst exponent using ratio of range variance estimators
foo <- lapply(0:(inter_val-1), function(off_set) {
  end_points <- rutils::end_points(SPY, inter_val=inter_val, off_set=off_set)
# aggregate over 10-minute end_points:
  SPY10 <- to_period(x_ts=SPY, end_points=end_points)
  vari_ance10 <- vari_ance(SPY10)
# vari_ance <- vari_ance(SPY)
  vari_ance <- diff(cumsum(vari_ance(SPY))[index(vari_ance10)])/inter_val
  roll_hurst <- ifelse((vari_ance==0) | (vari_ance10==0), 
                       NA, 
                       log(vari_ance10/vari_ance)/log(inter_val))
  roll_hurst[1] <- roll_hurst[2]
  roll_hurst <- na.locf(roll_hurst)
#  log((x_ts10[, 2]-x_ts10[, 3])/vol_xts)/log(10))
  colnames(roll_hurst) <- paste0(sym_bol, ".Hurst")
  roll_hurst
})  # end lapply

bar <- do_call_rbind(foo)
which(.index(bar)==(.index(bar)[-1]))
bar[625424] <- mean(bar[625424:625434])
bar <- bar[-(625424:625434)]
bar[625421:625434]
season_hurst <- season_ality(bar)

season_hurst <- season_ality(roll_hurst)
colnames(season_hurst) <- paste0(sym_bol, ".season_hurst")
x11()
chart_Series(x=season_hurst, 
             name=paste(colnames(season_hurst), "daily seasonality"))


agg_regate
season_var <- season_ality(vari_ance(ohlc=get(sym_bol)[inter_val, 1:4], calc_method="rogers.satchell"))
colnames(season_var) <- paste0(na_me(get(sym_bol)), ".season_var")


## daily seasonality

library(HighFreq)
season_volume <- season_ality(Vo(get(sym_bol)[inter_val]))
colnames(season_volume) <- paste0(na_me(get(sym_bol)), ".season_volume")
season_data <- season_volume

season_var <- season_ality(vari_ance(ohlc=get(sym_bol)[inter_val, 1:4], calc_method="rogers.satchell"))
colnames(season_var) <- paste0(na_me(get(sym_bol)), ".season_var")
season_data <- season_var

season_illiquid <- 1e6*sqrt(season_var/season_volume)
foo <- 1e6*ifelse(Vo(get(sym_bol))==0, 0, sqrt(vari_ance(ohlc=get(sym_bol))/Vo(get(sym_bol))))
season_illiquid <- season_ality(foo)
season_illiquid <- season_ality(sqrt(vari_ance(ohlc=get(sym_bol))/Vo(get(sym_bol))))
colnames(season_illiquid) <- paste0(na_me(get(sym_bol)), ".season_illiquid")
season_data <- season_illiquid

season_skew <- season_ality(skew_ohlc(ohlc=get(sym_bol)[inter_val, 1:4], calc_method="rogers.satchell"))
# season_skew <- season_skew/(season_var)^(1.5)
colnames(season_skew) <- paste0(na_me(get(sym_bol)), ".season_skew")
season_data <- season_skew

season_hurst <- season_ality(hurst_ohlc(ohlc=get(sym_bol)[inter_val, 1:4]))
colnames(season_hurst) <- paste0(na_me(get(sym_bol)), ".season_hurst")
season_data <- season_hurst

season_data <- season_ality(get(sym_bol)[inter_val, 4])



season_autocorr <- re_turns[, 1]*lag(re_turns[, 1])
season_autocorr[1, ] <- 0
season_autocorr <- re_turns[, 1]*lag(roll_sum(re_turns[, 1], win_dow=5)/5)
season_autocorr[1:5, ] <- 0
season_autocorr <- season_ality(season_autocorr)
colnames(season_autocorr) <- paste0(na_me(get(sym_bol)), ".season_autocorr")
season_data <- roll_sum(season_autocorr, win_dow=5)/5


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
da_ta <- cbind(roll_autocorr, roll_var)
da_ta <- cbind(foo, roll_var)

da_ta <- cbind(daily_hurst, daily_var)
da_ta <- cbind(season_hurst, season_var)
da_ta <- cbind(daily_hurst, daily_skew)

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
roll_skew <- roll_agg(ohlc=get(sym_bol)[inter_val], esti_mator="skew_ohlc", n=10)
roll_skew <- roll_skew/(roll_var)^(1.5)
roll_skew[1, ] <- 0
roll_skew <- na.locf(roll_skew)
tail(roll_skew, 11)


## Bootstrap of standard errors of all the methods in function vari_ance()

set.seed(1121)  # reset random number generator

# create time index of one second intervals for a single day
in_dex <- seq(from=as.POSIXct("2016-01-01 00:00:00"),
              to=as.POSIXct("2016-01-30 00:00:00"), by="1 sec")
# create xts of random prices
re_turns <- rnorm(length(in_dex), sd=0.001)
x_ts <- xts(exp(cumsum(re_turns)), order.by=in_dex)
x_ts <- xts(exp(cumsum(rnorm(length(in_dex), sd=0.001))), order.by=in_dex)
# aggregate to minutes OHLC data
oh_lc <- xts::to.period(x=x_ts, period="minutes")
# calculate variance estimates
vari_ance <- vari_ance(oh_lc)
sum(vari_ance)

# calculate variance using all the different estimators in vari_ance()
meth_ods <- c("close", "garman.klass", "rogers.satchell", "garman.klass_yz", "yang.zhang")
sapply(meth_ods, function(meth_od) {
  sum(vari_ance(oh_lc, calc_method=meth_od))
})

# calculate standard errors of variance estimators using bootstrap
boot_strap <- sapply(1:100, function(x) {
  boot_sample <- sample(re_turns, replace=TRUE)
  x_ts <- xts(exp(cumsum(boot_sample)), order.by=in_dex)
# aggregate to minutes OHLC data
  oh_lc <- xts::to.period(x=x_ts, period="minutes")
# calculate variance estimates
  sapply(meth_ods, function(meth_od) {
    sum(vari_ance(oh_lc, calc_method=meth_od))
  })  # end sapply
})  # end sapply

# calculate standard errors of variance estimators using bootstrap - parallel version
boot_strap <- sapply(1:100, function(x) {
  boot_sample <- sample(re_turns, replace=TRUE)
  x_ts <- xts(exp(cumsum(boot_sample)), order.by=in_dex)
  # aggregate to minutes OHLC data
  oh_lc <- xts::to.period(x=x_ts, period="minutes")
  # calculate variance estimates
  sapply(meth_ods, function(meth_od) {
    sum(vari_ance(oh_lc, calc_method=meth_od))
  })  # end sapply
})  # end sapply


head(t(boot_strap))
boot_errors <- apply(boot_strap, MARGIN=1, sd)

save(boot_strap, file="C:/Develop/data/boot_strap.RData")

# calculate variance estimates for SPY
vari_ance <- vari_ance(SPY)




