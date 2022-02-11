\name{wtd.quantile}
\alias{wtd.quantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted Quantiles
}
\description{
Generic function for calculating weighted quantiles.
}
\usage{
wtd.quantile(x, weights, probs = seq(0, 1, 0.25), 
             na.rm = FALSE, names = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numerical or logical vector.
}
  \item{weights}{
Vector of non-negative weights.  
  }
  \item{probs}{
Numeric vector of probabilities with values in [0,1].
  }
  \item{na.rm}{
Logical indicating whether \code{NA} values should be removed before calculation.
  }
  \item{names}{
Logical indicating if the result should have names corresponding to the probabilities.  
  }  
}
\details{
If \code{weights} are missing, the weights are defined to be a vector of ones (which is the same as the unweighted quantiles).

The weighted quantiles are computed by linearly interpolating the empirical cdf via the \code{\link{approx}} function. 
}
\value{
Returns the weighted quantiles corresponding to the input probabilities.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
If the weights are all equal (or missing), the resulting quantiles are equivalent to those produced by the \code{\link{quantile}} function using the 'type = 4' argument.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{wtd.mean}} for weighted mean calculations

\code{\link{wtd.var}} for weighted variance calculations
}
\examples{
# generate data and weights
set.seed(1)
x <- rnorm(10)
w <- rpois(10, lambda = 10)

# unweighted quantiles
quantile(x, probs = c(0.1, 0.9), type = 4)
wtd.quantile(x, probs = c(0.1, 0.9))

# weighted quantiles
sx <- sort(x, index.return = TRUE)
sw <- w[sx$ix]
ecdf <- cumsum(sw) / sum(sw)
approx(x = ecdf, y = sx$x, xout = c(0.1, 0.9), rule = 2)$y
wtd.quantile(x, w, probs = c(0.1, 0.9))

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
