\name{plot.ss}
\alias{plot.ss}
\alias{plot.boot.ss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for Smoothing Spline Fit and Bootstrap
}
\description{
Default plotting methods for \code{\link{ss}} and \code{\link{boot.ss}} objects.
}
\usage{
\method{plot}{ss}(x, n = 201, ci = TRUE, xseq = NULL, ...)

\method{plot}{boot.ss}(x, n = 201, ci = TRUE, xseq = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class 'ss' or 'boot.ss'
}
  \item{n}{
number of points used to plot smoothing spline estimate
}
  \item{ci}{
logical indicating whether to include a confidence interval
}
  \item{xseq}{
ordered sequence of points at which to plot smoothing spline estimate
}
  \item{...}{
optional additional argument for the \code{\link{plotci}} function, e.g., \code{level}, \code{col}, etc.
}
}
\details{
Unless a sequence of points is provided via the \code{xseq} arugment, the plots are created by evaluating the smoothing spline fit at an equidistant sequence of \code{n} values that span the range of the training data.
}
\value{
Plot of the function estimate and confidence interval with the title displaying the effective degrees of freedom.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The \code{plot.ss} and \code{plot.boot.ss} functions produce plots that only differ in terms of their confidence intervals: \code{plot.ss} uses the Bayesian CIs, whereas \code{plot.boot.ss} uses the bootstrap CIs.


}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ss}} and \code{\link{boot.ss}}
}
\examples{
# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
fx <- 2 + 3 * x + sin(2 * pi * x)
y <- fx + rnorm(n, sd = 0.5)

# fit smoothing spline
ssfit <- ss(x, y, nknots = 10)

# plot smoothing spline fit
plot(ssfit)

\dontrun{

# bootstrap smoothing spline
ssfitboot <- boot(ssfit)

# plot smoothing spline bootstrap
plot(ssfitboot)
}

}
\keyword{ aplot }
\keyword{ dplot }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
