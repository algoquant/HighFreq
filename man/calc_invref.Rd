% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calc_invref}
\alias{calc_invref}
\title{Calculate the inverse of a square \emph{matrix} in place, without copying
the data in memory.}
\usage{
calc_invref(matrixv)
}
\arguments{
\item{\code{matrixv}}{A \emph{matrix} of data to be inverted.  (The argument
is interpreted as a \emph{pointer} to a \emph{matrix}, and it is
overwritten with the inverse matrix.)}
}
\value{
No return value.
}
\description{
Calculate the inverse of a square \emph{matrix} in place, without copying
the data in memory.
}
\details{
The function \code{calc_invref()} calculates the inverse of a square
  \emph{matrix} in place, without copying the data in memory. It accepts a
  \emph{pointer} to the argument \code{matrixv} (which is the \code{matrix}
  to be inverted), and it overwrites the old matrix values with the inverse
  matrix values. It performs the calculation in place, without copying the
  data in memory, which can significantly increase the computation speed for
  large matrices.

  The function \code{calc_invref()} doesn't return a value.
  The function \code{calc_invref()} calls the \code{C++} \code{Armadillo}
  function \code{arma::inv()} to calculate the matrix inverse.
}
\examples{
\dontrun{
# Calculate a random matrix
matrixv <- matrix(rnorm(100), nc=10)
# Copy matrixv to a matrix in a different memory location
invmat <- matrixv + 0
# Calculate the inverse in place using RcppArmadillo
HighFreq::calc_invref(invmat)
# Multiply the matrix times its inverse
multmat <- matrixv \%*\% invmat
round(multmat, 4)
# Calculate the sum of the off-diagonal elements
sum(multmat[upper.tri(multmat)])
# Compare RcppArmadillo with R
all.equal(invmat, solve(matrixv))
# Compare the speed of RcppArmadillo with R code
library(microbenchmark)
summary(microbenchmark(
   rcode=solve(matrixv),
   cppcode=calc_invref(matrixv),
   times=10))[, c(1, 4, 5)]
}  # end dontrun

}
