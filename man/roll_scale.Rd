% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{roll_scale}
\alias{roll_scale}
\title{Perform a rolling scaling (standardization) of the columns of a
\emph{matrix} of data using \code{RcppArmadillo}.}
\usage{
roll_scale(matrix, look_back, use_median = FALSE)
}
\arguments{
\item{use_median}{A \emph{Boolean} argument: if \code{TRUE} then the 
centrality (central tendency) is calculated as the \emph{median} and the 
dispersion is calculated as the \emph{median absolute deviation}
(\emph{MAD}).
If \code{use_median} is \code{FALSE} then the centrality is calculated as 
the \emph{mean} and the dispersion is calculated as the \emph{standard
deviation} (the default is \code{use_median = FALSE})}

\item{\code{matrix}}{A \emph{matrix} of data.}

\item{\code{look_back}}{The length of the look-back interval, equal to the number 
of rows of data used in the scaling.}
}
\value{
A \emph{matrix} with the same dimensions as the input argument
  \code{matrix}.
}
\description{
Perform a rolling scaling (standardization) of the columns of a
\emph{matrix} of data using \code{RcppArmadillo}.
}
\details{
The function \code{roll_scale()} performs a rolling scaling
  (standardization) of the columns of the \code{matrix} argument using
  \code{RcppArmadillo}.
  The function \code{roll_scale()} performs a loop over the rows of 
  \code{matrix}, subsets a number of previous (past) rows equal to 
  \code{look_back}, and scales the subset matrix.  It assigns the last row
  of the scaled subset \emph{matrix} to the return matrix.
  
  If the argument \code{use_median} is \code{FALSE} (the default), then it
  performs the same calculation as the function \code{roll::roll_scale()}.
  If the argument \code{use_median} is \code{TRUE}, then it calculates the
  centrality as the \emph{median} and the dispersion as the \emph{median
  absolute deviation} (\emph{MAD}).
}
\examples{
\dontrun{
mat_rix <- matrix(rnorm(20000), nc=2)
look_back <- 11
rolled_scaled <- roll::roll_scale(data=mat_rix, width = look_back, min_obs=1)
rolled_scaled2 <- roll_scale(matrix=mat_rix, look_back = look_back, use_median=FALSE)
all.equal(rolled_scaled[-1, ], rolled_scaled2[-1, ])
}

}
