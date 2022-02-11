\name{wtd.var}
\alias{wtd.var}
\alias{wtd.sd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted Variance and Standard Deviation
}
\description{
Generic function for calculating weighted variance or standard deviation of a vector.
}
\usage{
wtd.var(x, weights, na.rm = FALSE)

wtd.sd(x, weights, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numerical or logical vector.
}
  \item{weights}{
Vector of non-negative weights.  
  }
  \item{na.rm}{
Logical indicating whether \code{NA} values should be removed before calculation.
  }
}
\details{
The weighted variance is defined as

\code{(n / (n - 1)) * sum(weights * (x - xbar)^2) / sum(weights)}

where \code{n} is the number of observations with non-zero weights, and \code{xbar} is the weighted mean computed via the \code{\link{wtd.mean}} function.

The weighted standard deviation is the square root of the weighted variance.
}
\value{
Returns the weighted variance or standard deviation.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
If \code{weights} are missing, the weights are defined to be a vector of ones (which is the same as the unweighted variance or standard deviation).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{wtd.mean}} for weighted mean calculations

\code{\link{wtd.quantile}} for weighted quantile calculations
}
\examples{
# generate data and weights
set.seed(1)
x <- rnorm(10)
w <- rpois(10, lambda = 10)

# weighted mean
xbar <- wtd.mean(x, w)

# weighted variance
wtd.var(x, w)
(10 / 9) * sum(w * (x - xbar)^2) / sum(w)

# weighted standard deviation
wtd.sd(x, w)
sqrt((10 / 9) * sum(w * (x - xbar)^2) / sum(w))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ univar }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
