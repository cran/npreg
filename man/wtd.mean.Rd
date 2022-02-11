\name{wtd.mean}
\alias{wtd.mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted Arithmetic Mean
}
\description{
Generic function for calculating the weighted (and possibly trimmed) arithmetic mean.
}
\usage{
wtd.mean(x, weights, trim = 0, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numerical or logical vector.
}
  \item{weights}{
Vector of non-negative weights.  
  }
\item{trim}{
Fraction [0, 0.5) of observations trimmed from each end before calculating mean.
  }
\item{na.rm}{
Logical indicating whether \code{NA} values should be removed before calculation.
  }  
}
\details{
If \code{weights} are missing, the weights are defined to be a vector of ones (which is the same as the unweighted arithmetic mean).

If \code{trim} is non-zero, then \code{trim} observations are deleted from each end before the (weighted) mean is computed. The quantiles used for trimming are defined using the \code{\link{wtd.quantile}} function.
}
\value{
Returns the weighted and/or trimmed arithmetic mean.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The weighted (and possible trimmed) mean is defined as:

\code{sum(weights * x) / sum(weights)}

where \code{x} is the (possibly trimmed version of the) input data.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{wtd.var}} for weighted variance calculations

\code{\link{wtd.quantile}} for weighted quantile calculations
}
\examples{
# generate data and weights
set.seed(1)
x <- rnorm(10)
w <- rpois(10, lambda = 10)

# weighted mean
wtd.mean(x, w)
sum(x * w) / sum(w)

# trimmed mean
q <- quantile(x, probs = c(0.1, 0.9), type = 4)
i <- which(x < q[1] | x > q[2])
mean(x[-i])
wtd.mean(x, trim = 0.1)

# weighted and trimmed mean
q <- wtd.quantile(x, w, probs = c(0.1, 0.9))
i <- which(x < q[1] | x > q[2])
wtd.mean(x[-i], w[-i])
wtd.mean(x, w, trim = 0.1)

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
