% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calc_endpoints}
\alias{calc_endpoints}
\title{Calculate a vector of end points that divides an integer time sequence of
time periods into equal time intervals.}
\usage{
calc_endpoints(length, step = 1L, stub = 0L, stubs = TRUE)
}
\arguments{
\item{\code{length}}{An \emph{integer} equal to the length of the time
sequence to be divided into equal intervals.}

\item{\code{step}}{The number of time periods in each interval between
neighboring end points (the default is \code{step = 1}).}

\item{\code{stub}}{An \emph{integer} equal to the first non-zero end point
(the default is \code{stub = 0}).}

\item{\code{stubs}}{A \emph{Boolean} specifying whether to include stub
intervals (the default is \code{stubs = TRUE}).}
}
\value{
A vector of equally spaced \emph{integers} representing the end
  points.
}
\description{
Calculate a vector of end points that divides an integer time sequence of
time periods into equal time intervals.
}
\details{
The end points are a vector of integers which divide the sequence of time
  periods of length equal to \code{length} into equally spaced time
  intervals.
  The number of time periods between neighboring end points is equal to the
  argument \code{step}.
  If a whole number of intervals doesn't fit over the whole sequence, then
  \code{calc_endpoints()} adds a stub interval at the end.
  A stub interval is one where the number of periods between neighboring end
  points is less than the argument \code{step}.
  
  If \code{stubs = TRUE} (the default) then the first end point is 
  equal to \code{0} (since indexing in \code{C++} code starts at \code{0}).
  The first non-zero end point is equal to \code{step} or \code{stub} (if
  it's not zero).
  If \code{stub = 0} (the default) then the first end point is equal to
  \code{0} (even if \code{stubs = FALSE}).
  If \code{stubs = TRUE} (the default) then the last end point is always
  equal to \code{length-1}.
  The argument \code{stub} should be less than the \code{step}: \code{stub <
  step}.
  
  If \code{step = 1} and \code{stub = 0} (the default), then the vector of
  end points is simply equal to:
  \deqn{
    \{ 0, 1, 2, ..., length - 1 \}
  }

  If \code{stub = 0} (the default) and \code{stubs = TRUE} (the default)
  then the vector of end points is equal to:
  \deqn{
    \{ 0, step, 2*step, ..., length - 1 \}
  }
  
  If \code{stub = 0} (the default) and \code{stubs = FALSE} then the vector
  of end points is equal to:
  \deqn{
    \{ 0, step, 2*step, ..., n*step \}
  }
  
  If \code{stub > 0} and \code{stubs = TRUE} (the default), then the vector
  of end points is equal to:
  \deqn{
    \{ 0, stub, stub + step, ..., length - 1 \}
  }
  
  For example, the end points for \code{length = 20}, divided into intervals
  of \code{step = 5} are equal to: \code{0, 5, 10, 15, 19}.
  
  If \code{stub = 1} then the first non-zero end point is equal to \code{1}
  and the end points are equal to: \code{0, 1, 6, 11, 16, 19}.
  The stub interval at the beginning is equal to \code{2} (including
  \code{0} and \code{1}).
  The stub interval at the end is equal to \code{3 = 19 - 16}.
  
  The end points for \code{length = 21} divided into intervals of length
  \code{step = 5}, with \code{stub = 0}, are equal to: \code{0, 5, 10, 15,
  20}.
  The beginning interval is equal to \code{5}.
  The end interval is equal to \code{5 = 20 - 15}.
  
  If \code{stub = 1} then the first non-zero end point is equal to \code{1}
  and the end points are equal to: \code{0, 1, 6, 11, 16, 20}.
  The beginning stub interval is equal to \code{2}.
  The end stub interval is equal to \code{4 = 20 - 16}.
  
  The function \code{calc_endpoints()} is similar to the function
  \code{rutils::calc_endpoints()} from package
  \href{https://github.com/algoquant/rutils}{rutils}.
  
  But the end points are shifted by \code{-1} compared to \code{R} code
  because indexing starts at \code{0} in \code{C++} code, while it starts at
  \code{1} in \code{R} code. So if \code{calc_endpoints()} is used in
  \code{R} code then \code{1} should be added to it.
}
\examples{
# Calculate the end points without a stub interval
HighFreq::calc_endpoints(length=20, step=5)
# Calculate the end points with a final stub interval
HighFreq::calc_endpoints(length=23, step=5)
# Calculate the end points with initial and final stub intervals
HighFreq::calc_endpoints(length=20, step=5, stub=2)

}
