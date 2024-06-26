[![Build
Status](https://travis-ci.org/algoquant/HighFreq.svg?branch=master)](https://travis-ci.org/algoquant/HighFreq)

### Overview

The motivation for the package *HighFreq* is to create a library of
functions designed for managing trade and quote (*TAQ*) and *OHLC* data,
and for efficiently estimating various statistics, like volatility,
skew, Hurst exponent, and Sharpe ratio, from that data.

The are several other packages which offer much of this functionality,
like for example:

-   package
    [xts](https://cran.r-project.org/web/packages/xts/index.html)

-   package
    [TTR](https://cran.r-project.org/web/packages/TTR/index.html)

-   package
    [PerformanceAnalytics](https://cran.r-project.org/web/packages/PerformanceAnalytics/index.html)

-   package
    [highfrequency](https://cran.r-project.org/web/packages/highfrequency/index.html)

Unfortunately many of the functions in these packages are either too
slow, or lack some critical functionality, or produce data in
inconsistent formats (with *NA* values, etc.) The package *HighFreq*
aims to create a unified framework, with consistent data formats and
naming conventions.

The package *HighFreq* relies on *OHLC* price and volume data formatted
as *xts* time series, because the *OHLC* data format provides an
efficient way of compressing *TAQ* data, while preserving information
about price levels, volatility (range), and trading volumes. Most
existing packages don’t rely on *OHLC* data, so their statistical
estimators are much less efficient than those in package *HighFreq*.

### Running and Rolling Statistics Over Time Series Data

Definitions of running and rolling statistics (aggregations):

-   A statistic is some function of *OHLC* data. For example, the
    difference between the *High* minus the *Low* prices is a simple
    statistic. The estimators of volatility, skew, and higher moments
    are also statistics.

-   The functions called *calc\_\** calculate aggregations over columnar
    data, and produce a single number or a vector. For example, the
    function *calc\_var()* calculates the variance of the columns of
    returns data, and produces a row vector.

-   The functions called *agg\_\** calculate various other data
    aggregations over columnar data. For example, the function
    *agg\_ohlc()* aggregates a time series of data into a single bar of
    *OHLC* data.

-   The package *HighFreq* calculates two different types of statistics
    time series: *running* and *rolling* statistics. The *running*
    statistics time series are calculated using individual bars of
    *OHLC* data, while *rolling* statistics time series are calculated
    using multiple bars of *OHLC* data.

-   An example of a *running* statistic is a time series of squared
    differences of the *High* minus *Low* prices. Each point in a
    *running* statistic depends only on a single bar of *OHLC* data (and
    possibly on the neighboring bars.)

-   An example of a *rolling* statistic is a time series of average
    *Close* prices calculated over a rolling look-back interval. Some
    *rolling* statistics can be calculated from *running* statistics by
    calculating weighted averages.

-   The functions called *run\_\** calculate *running* statistics based
    on each bar of *OHLC* data, and produce a single-column *xts* time
    series with the same number of rows as the *OHLC* time series.

-   The functions called *roll\_\** calculate *rolling* statistics based
    on multiple bars of *OHLC* data taken from a rolling look-back
    interval, and often produce a single-column *xts* time series with
    the same number of rows as the *OHLC* time series. The *roll\_\**
    functions perform loops over the rows of the data, and they can
    apply the *calc\_\** functions to subsets of the data over look-back
    intervals.

-   *Running* and *rolling* statistics also represent a form of data
    aggregation, and can include many standard technical indicators,
    like *VWAP*, *Bollinger Bands*, etc.

<br>

### Functions for data scrubbing, formatting, and aggregation

The package *HighFreq* contains several categories of functions designed
for:

-   managing *TAQ* and *OHLC* time series,

-   estimating running and rolling statistics over time series,

The package *HighFreq* contains functions for:

-   chaining and joining time series,

-   scrubbing bad data from time series,

-   managing time zones and alligning time indices,

-   converting *TAQ* data to *OHLC* format,

-   aggregating data to lower frequency (periodicity),

-   estimating running statistics from *OHLC* data, such as volatility,
    skew, and higher moments (functions called *run\_\**),

-   calculating rolling aggregations (VWAP, Hurst exponent, Sharpe
    ratio, etc.),

-   calculating seasonality aggregations,

-   creating random *TAQ* and *OHLC* time series,

### Installation and loading

Install package *HighFreq* from github:

``` r
install.packages("devtools")
devtools::install_github(repo="algoquant/HighFreq")
library(HighFreq)
```

<br>

Install package *HighFreq* from source on local drive:

``` r
install.packages(pkgs="C:/Develop/R/HighFreq", repos=NULL, type="source")
# Install package from source on local drive using R CMD
R CMD INSTALL C:\Develop\R\HighFreq
library(HighFreq)
```

<br>

Build reference manual for package *HighFreq* from *.Rd* files:

``` r
system("R CMD Rd2pdf C:/Develop/R/HighFreq")
R CMD Rd2pdf C:\Develop\R\HighFreq
```

<br>

### Data

Trade and Quote (*TAQ*) data contains intraday trades and quotes on
exchange-traded stocks and futures. *TAQ* data is spaced irregularly in
time, with data recorded each time a new trade or quote arrives. The
rows of *TAQ* data contain the quoted and traded prices, and the
corresponding quote size or trade volume.

*TAQ* data can be aggregated into Open-High-Low-Close (*OHLC*) data.
*OHLC* data is evenly spaced in time, with each row containing the
*Open*, *High*, *Low*, and *Close* prices, and the trade *Volume*,
recorded over the past time interval (called a *bar* of data). The
*Open* and *Close* prices are the first and last trade prices recorded
in the time bar. The *High* and *Low* prices are the highest and lowest
trade prices recorded in the time bar. The *Volume* is the total trading
volume recorded in the time bar.

Aggregating *TAQ* data into *OHLC* data provides data compression, while
preserving information about price levels, volatility (range), and
trading volumes. In addition, evenly spaced data allows analysis of
multiple time series, since all the prices are given at the same moments
of time.

The package *HighFreq* includes three *xts* time series called *SPY*,
*TLT*, and *VXX*, containing intraday 1-minute *OHLC* data for the
*SPY*, *TLT*, and *VXX* ETFs. The package *HighFreq* also includes an
*xts* time series called *SPY\_TAQ* with a single day of *TAQ* data for
the *SPY* ETF. The data is set up for lazy loading, so it doesn’t
require calling `data(hf_data)` to load it before being able to call it.

The data source is the [Wharton Research Data
Service](https://wrds-web.wharton.upenn.edu/wrds/)

List all the data sets included in the *HighFreq* package:

``` r
# list all datasets in package HighFreq
data(package="HighFreq")
```

<br>

### Examples

More examples can be found in the vignettes titled
*managing\_time\_series* and *estimating\_statistics*.

Aggregate *TAQ* data into a 1-minute bar *OHLC* time series:

``` r
# aggregate TAQ data to 1-min OHLC bar data, for a single symbol, and save to file
sym_bol <- "SPY"
save_scrub_agg(sym_bol, 
               data_dir="E:/mktdata/sec/", 
               output_dir="E:/output/data/")
```

<br>

Calculate daily trading volume:

``` r
daily_volume <- apply.daily(x=Vo(SPY), FUN=sum)
colnames(daily_volume) <- "SPY.Volume")
chart_Series(x=daily_volume, name="daily trading volumes for SPY")
```

<br>

Calculate daily average open to close variance from minutely *OHLC*
prices:

``` r
# calculate daily average open to close variance
var_daily <- (6.5*60*60^2)*xts::apply.daily(x=SPY, FUN=agg_stats_r, 
                              calc_bars="run_variance", calc_method="rogers_satchell")
colnames(var_daily) <- "SPY.Var"
chart_Series(100*sqrt(var_daily["/2010"]), name="SPY daily standard deviation")
```

<br>

Calculate daily skew from minutely *OHLC* prices:

``` r
skew_daily <- apply.daily(x=SPY, FUN=agg_stats_r, calc_bars="run_skew")
skew_daily <- skew_daily/(var_daily)^(1.5)
colnames(skew_daily) <- "SPY.Skew")
chart_Series(x=skew_daily, name="daily skew for SPY")
```

<br>

Calculate rolling prices:

``` r
roll_prices <- rutils::roll_sum(Op(SPY), win_dow=10)/10
colnames(roll_prices) <- "SPY.Rets"
# plot candle chart
chart_Series(SPY["2013-11-12", ], name="SPY Prices")
add_TA(roll_prices["2013-11-12"], on=1, col="red", lwd=2)
```

<br>

Calculate rolling volume-weighted variance:

``` r
var_rolling <- roll_stats(oh_lc=SPY["2012"], calc_stats="run_variance", win_dow = 10)
# plot without overnight jump
chart_Series(var_rolling["2012-11-12", ][-(1:11)], name="SPY rolling volume-weighted variance")
```

<br>

Calculate daily seasonality of variance:

``` r
var_seasonal <- season_ality((24*60*60^2)*run_variance(oh_lc=SPY))
colnames(var_seasonal) <- "SPY.var_seasonal"
chart_Series(x=var_seasonal, name="SPY variance daily seasonality")
```
