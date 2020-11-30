\name{theta.mle}
\alias{theta.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MLE of Theta for Negative Binomial
}
\description{
Computes the maximum likelihood estimate of the size (theta) parameter for the Negative Binomial distribution via a Newton-Raphson algorithm.
}
\usage{
theta.mle(y, mu, theta, wt = 1, 
          maxit = 100, tol = .Machine$double.eps^0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
response vector
}
  \item{mu}{
mean vector
}
  \item{theta}{
initial theta (optional)
}
  \item{wt}{
weight vector
}
  \item{maxit}{
max number of iterations  
}
  \item{tol}{
convergence tolerance  
}
}
\details{
Based on the \code{glm.nb} function in the \bold{MASS} package. If \code{theta} is missing, the initial estimate of theta is given by 

\code{theta <- 1 / mean(wt * (y / mu - 1)^2)} 

which is motivated by the method of moments estimator for the dispersion parameter in a quasi-Poisson model.
}
\value{
Returns estimated theta with attributes
\item{SE}{standard error estimate}
\item{iter}{number of iterations}
}
\references{
Venables, W. N. and Ripley, B. D. (1999) Modern Applied Statistics with S-PLUS. Third Edition. Springer.

https://www.rdocumentation.org/packages/MASS/versions/7.3-51.6/topics/negative.binomial

https://www.rdocumentation.org/packages/MASS/versions/7.3-51.6/topics/glm.nb
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{NegBin}} for details on the Negative Binomial distribution
}
\examples{
# generate data
n <- 1000
x <- seq(0, 1, length.out = n)
fx <- 3 * x + sin(2 * pi * x) - 1.5
mu <- exp(fx)

# simulate negative binomial data
set.seed(1)
y <- rnbinom(n = n, size = 1/2, mu = mu)

# estimate theta
theta.mle(y, mu)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ distribution }
\keyword{ regression }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
