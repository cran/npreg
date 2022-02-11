\name{boot.ss}
\alias{boot.ss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap a Fit Smoothing Spline
}
\description{
Bootstraps a fit smoothing spline to form confidence intervals (BCa or percentile) and standard error estimates. Three bootstrapping methods are implemented: nonparametric bootstrapping cases, nonparametric bootsrapping residuals, and parameteric bootstrapping residuals.
}
\usage{
\method{boot}{ss}(object, statistic, ..., 
     R = 9999, level = 0.95, bca = TRUE, 
     method = c("cases", "resid", "param"), 
     fix.lambda = TRUE, cov.mat = FALSE, 
     boot.dist = FALSE, verbose = TRUE,
     parallel = FALSE, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fit from \code{\link{ss}}
}
  \item{statistic}{
a function to compute the statistic (see Details and Note)
}
  \item{...}{
additional arguments to \code{statistic} function (optional)
}
  \item{R}{
number of bootstrap resamples used to form bootstrap distribution
}
  \item{level}{
confidence level for bootstrap confidence intervals
}
  \item{bca}{
logical indicating whether to calculate BCa (default) or percentile intervals
}
  \item{method}{
resampling method used to form bootstrap distribution
}
  \item{fix.lambda}{
logical indicating whether the smoothing parameter should be fixed (default) or re-estimated for each bootstrap sample
}
  \item{cov.mat}{
logical indicating whether the bootstrap estimate of the covariance matrix should be returned
}
  \item{boot.dist}{
logical indicating whether the bootstrap distribution should be returned
}
  \item{verbose}{
logical indicating whether the bootstrap progress bar should be printed
}
  \item{parallel}{
logical indicating if the \code{\link{parallel}} package should be used for parallel computing (of the bootstrap distribution). Defaults to FALSE, which implements sequential computing.
}
  \item{cl}{
cluster for parallel computing, which is used when \code{parallel = TRUE}. Note that if \code{parallel = TRUE} and \code{cl = NULL}, then the cluster is defined as \code{makeCluster(detectCores())}.
}

}
\details{
If the \code{statistic} function is not provided, then the function is internally defined as follows

\verb{
statistic <- 
   function(object, ...){
      xmax <- object$fit$min + object$fit$range
      xseq <- seq(object$fit$min, xmax, length.out = 201)
      predict(object, x = xseq)$y
   }
}

which produces predictions at an equidistant sequence of values spanning the range of the training data.

}
\value{
Produces an object of class 'boot.ss' with the following elements:
\item{t0 }{Observed statistic, computed using \code{statistic(object, ...)}}

\item{se }{Bootstrap estimate of the standard error}

\item{bias }{Bootstrap estimate of the bias}

\item{cov }{Bootstrap estimate of the covariance (if \code{cov.mat = TRUE})}

\item{ci }{Bootstrap estimate of the confidence interval}

\item{boot.dist }{Bootstrap distribution of statistic (if \code{boot.dist = TRUE})}

\item{x0 }{Sequence of 'x' values used to calculate statistic. If \code{statistic} is provided, this is a vector of integers ranging from 1 to \code{length(t0)}. Otherwise it is an equidistant sequence of 201 values spanning the range of the data.}

\item{bias.correct }{Bias correction factor for BCa confidence interval.}

\item{acceleration }{Acceleration parameter for BCa confidence interval.}

The output list also contains the elements \code{object}, \code{R}, \code{level}, \code{bca}, \code{method}, and \code{fix.lambda}, all of which are the same as the corresponding input arguments.

}
\references{
Davison, A. C., & Hinkley, D. V. (1997). \emph{Bootstrap Methods and Their Application}. Cambridge University Press. \doi{10.1017/CBO9780511802843}

Efron, B., & Tibshirani, R. J. (1994). \emph{An Introduction to the Boostrap}. Chapman & Hall/CRC. \doi{10.1201/9780429246593}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The \code{statistic} function must satisfy the following two requirements:
(1) the first input must be a fit smoothing spline object, and
(2) the output must be a scalar or vector calculated from the 'ss' object.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ss}} and \code{\link{plot.boot.ss}}
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

# nonparameteric bootstrap cases (R = 99 for speed)
set.seed(0)
boot.cases <- boot(ssfit, R = 99)

# nonparameteric bootstrap residuals (R = 99 for speed)
set.seed(0)
boot.resid <- boot(ssfit, R = 99, method = "resid")

# parameteric bootstrap residuals (R = 99 for speed)
set.seed(0)
boot.param <- boot(ssfit, R = 99, method = "param")

# plot results
par(mfrow = c(1, 3))
plot(boot.cases, main = "Cases")
plot(boot.resid, main = "Residuals")
plot(boot.param, main = "Parametric")

\dontrun{

# nonparameteric bootstrap cases (recommended default of R = 9999)
boot.cases <- boot(ssfit)

# nonparameteric bootstrap residuals (recommended default of R = 9999)
boot.resid <- boot(ssfit, method = "resid")

# parameteric bootstrap residuals (recommended default of R = 9999)
boot.param <- boot(ssfit, method = "param")

# plot results
par(mfrow = c(1, 3))
plot(boot.cases, main = "Cases")
plot(boot.resid, main = "Residuals")
plot(boot.param, main = "Parametric")

}

}

\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{univar}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
