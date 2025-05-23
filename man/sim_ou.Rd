% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{sim_ou}
\alias{sim_ou}
\title{Simulate an \emph{Ornstein-Uhlenbeck} process using \emph{Rcpp}.}
\usage{
sim_ou(prici, priceq, theta, innov)
}
\arguments{
\item{\code{prici}}{The initial price.}

\item{\code{priceq}}{The equilibrium price.}

\item{\code{theta}}{The strength of mean reversion.}

\item{\code{innov}}{A single-column \emph{matrix} of innovations (random
numbers).}
}
\value{
A single-column \emph{matrix} of simulated prices, with the same
  number of rows as the argument \code{innov}.
}
\description{
Simulate an \emph{Ornstein-Uhlenbeck} process using \emph{Rcpp}.
}
\details{
The function \code{sim_ou()} simulates the following
  \emph{Ornstein-Uhlenbeck} process:
  \deqn{
    r_i = p_i - p_{i-1} = \theta \, (\mu - p_{i-1}) + \xi_i
  }
  \deqn{
    p_i = p_{i-1} + r_i
  }
  Where \eqn{r_i} and \eqn{p_i} are the simulated returns and prices,
  \eqn{\theta}, \eqn{\mu}, and \eqn{\sigma} are the
  \emph{Ornstein-Uhlenbeck} parameters, and \eqn{\xi_i} are the standard
  \emph{innovations}.
  The recursion starts with: \eqn{r_1 = \xi_1} and \eqn{p_1 = init\_price}.

  The function \code{sim_ou()} simulates the percentage returns as equal to
  the difference between the equilibrium price \eqn{\mu} minus the latest
  price \eqn{p_{i-1}}, times the mean reversion parameter \eqn{\theta}, plus
  a random normal innovation. The log prices are calculated as the sum of
  returns (not compounded), so they can become negative.

  The function \code{sim_ou()} simulates the \emph{Ornstein-Uhlenbeck}
  process using fast \emph{Rcpp} \code{C++} code.

  The function \code{sim_ou()} returns a single-column \emph{matrix}
  representing the \emph{time series} of simulated prices.
}
\examples{
\dontrun{
# Define the Ornstein-Uhlenbeck model parameters
prici <- 0.0
priceq <- 1.0
sigmav <- 0.01
thetav <- 0.01
innov <- matrix(rnorm(1e3))
# Simulate Ornstein-Uhlenbeck process using Rcpp
pricev <- HighFreq::sim_ou(prici=prici, priceq=priceq, volat=sigmav, theta=thetav, innov=innov)
plot(pricev, t="l", main="Simulated Ornstein-Uhlenbeck Prices", ylab="prices")
}  # end dontrun

}
