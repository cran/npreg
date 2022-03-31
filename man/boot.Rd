\name{boot}
\alias{boot}
\alias{boot.ss}
\alias{boot.sm}
\alias{boot.gsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap a Fit Smooth
}
\description{
Bootstraps a fit nonparametric regression model to form confidence intervals (BCa or percentile) and standard error estimates. 
}
\usage{
\method{boot}{ss}(object, statistic, ..., R = 9999, level = 0.95, bca = TRUE, 
     method = c("cases", "resid", "param"), fix.lambda = TRUE, cov.mat = FALSE, 
     boot.dist = FALSE, verbose = TRUE, parallel = FALSE, cl = NULL)

\method{boot}{sm}(object, statistic, ..., R = 9999, level = 0.95, bca = TRUE, 
     method = c("cases", "resid", "param"), fix.lambda = TRUE, 
     fix.thetas = TRUE, cov.mat = FALSE, boot.dist = FALSE, 
     verbose = TRUE, parallel = FALSE, cl = NULL)
     
\method{boot}{gsm}(object, statistic, ..., R = 9999, level = 0.95, bca = TRUE, 
     method = c("cases", "resid", "param"), fix.lambda = TRUE, 
     fix.thetas = TRUE, cov.mat = FALSE, boot.dist = FALSE, 
     verbose = TRUE, parallel = FALSE, cl = NULL)     
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fit from \code{\link{ss}} (smoothing spline), \code{\link{sm}} (smooth model), or \code{\link{gsm}} (generalized smooth model)
}
  \item{statistic}{
a function to compute the statistic (see Details)
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
  \item{fix.thetas}{
logical indicating whether the "extra" smoothing parameters should be fixed (default) or re-estimated for each bootstrap sample. Only applicable to \code{\link{sm}} and \code{\link{gsm}} objects with multiple penalized terms.
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
The \code{statistic} function must satisfy the following two requirements:

(1) the first input must be the \code{object} of class \code{\link{ss}}, \code{\link{sm}}, or \code{\link{gsm}}

(2) the output must be a scalar or vector calculated from the \code{object}

In most applications, the \code{statistic} function will be the model predictions at some user-specified \code{newdata}, which can be passed to \code{statistic} using the \code{...} argument.

If \code{statistic} is not provided, then the function is internally defined to be the model predictions at an equidistance sequence (for \code{\link{ss}} objects) or the training data predictor scores (for \code{\link{sm}} and \code{\link{gsm}} objects).
}
\value{
Produces an object of class 'boot.ss', 'boot.sm', or 'boot.gsm', with the following elements:
\item{t0 }{Observed statistic, computed using \code{statistic(object, ...)}}

\item{se }{Bootstrap estimate of the standard error}

\item{bias }{Bootstrap estimate of the bias}

\item{cov }{Bootstrap estimate of the covariance (if \code{cov.mat = TRUE})}

\item{ci }{Bootstrap estimate of the confidence interval}

\item{boot.dist }{Bootstrap distribution of statistic (if \code{boot.dist = TRUE})}

\item{bias.correct }{Bias correction factor for BCa confidence interval.}

\item{acceleration }{Acceleration parameter for BCa confidence interval.}

The output list also contains the elements \code{object}, \code{R}, \code{level}, \code{bca}, \code{method}, \code{fix.lambda}, and \code{fix.thetas}, all of which are the same as the corresponding input arguments.
}
\references{
Davison, A. C., & Hinkley, D. V. (1997). \emph{Bootstrap Methods and Their Application}. Cambridge University Press. \doi{10.1017/CBO9780511802843}

Efron, B., & Tibshirani, R. J. (1994). \emph{An Introduction to the Boostrap}. Chapman & Hall/CRC. \doi{10.1201/9780429246593}

Moulton, L. H., & Zeger, S. L. (1991). Bootstrapping generalized linear models. \emph{Computational Statistics & Data Analysis, 11}(1), 53-63. \doi{10.1016/0167-9473(91)90052-4}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

\note{
For \code{\link{gsm}} objects, requesting \code{method = "resid"} uses a variant of the one-step technique described in Moulton and Zeger (1991), which forms the bootstrap estimates of the coefficients without refitting the model.

As a result, when bootstrapping \code{\link{gsm}} objects with \code{method = "resid"}:

(1) it is necessary to set \code{fix.lambda = TRUE} and \code{fix.thetas = TRUE}

(2) any logical \code{statistic} must depend on the model \code{coefficients}, e.g., through the model predictions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ss}} for fitting "ss" (smoothing spline) objects

\code{\link{sm}} for fitting "sm" (smooth model) objects

\code{\link{gsm}} for fitting "gsm" (generalized smooth model) objects
}
\examples{
\dontrun{

##########   EXAMPLE 1   ##########
### smoothing spline

# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
fx <- 2 + 3 * x + sin(2 * pi * x)
y <- fx + rnorm(n, sd = 0.5)

# fit smoothing spline
ssfit <- ss(x, y, nknots = 10)

# nonparameteric bootstrap cases
set.seed(0)
boot.cases <- boot(ssfit)

# nonparameteric bootstrap residuals
set.seed(0)
boot.resid <- boot(ssfit, method = "resid")

# parameteric bootstrap residuals
set.seed(0)
boot.param <- boot(ssfit, method = "param")

# plot results
par(mfrow = c(1, 3))
plot(boot.cases, main = "Cases")
plot(boot.resid, main = "Residuals")
plot(boot.param, main = "Parametric")



##########   EXAMPLE 2   ##########
### smooth model

# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
fx <- 2 + 3 * x + sin(2 * pi * x)
y <- fx + rnorm(n, sd = 0.5)

# fit smoothing spline
smfit <- sm(y ~ x, knots = 10)

# define statistic (to be equivalent to boot.ss default)
newdata <- data.frame(x = seq(0, 1, length.out = 201))
statfun <- function(object, newdata) predict(object, newdata)

# nonparameteric bootstrap cases
set.seed(0)
boot.cases <- boot(smfit, statfun, newdata = newdata)

# nonparameteric bootstrap residuals
set.seed(0)
boot.resid <- boot(smfit, statfun, newdata = newdata, method = "resid")

# parameteric bootstrap residuals (R = 99 for speed)
set.seed(0)
boot.param <- boot(smfit, statfun, newdata = newdata, method = "param")
                   
# plot results
par(mfrow = c(1, 3))
plotci(newdata$x, boot.cases$t0, ci = boot.cases$ci, main = "Cases")
plotci(newdata$x, boot.resid$t0, ci = boot.resid$ci, main = "Residuals")
plotci(newdata$x, boot.param$t0, ci = boot.param$ci, main = "Parametric")



##########   EXAMPLE 3   ##########
### generalized smooth model

# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
fx <- 2 + 3 * x + sin(2 * pi * x)
y <- fx + rnorm(n, sd = 0.5)

# fit smoothing spline
gsmfit <- gsm(y ~ x, knots = 10)

# define statistic (to be equivalent to boot.ss default)
newdata <- data.frame(x = seq(0, 1, length.out = 201))
statfun <- function(object, newdata) predict(object, newdata)

# nonparameteric bootstrap cases
set.seed(0)
boot.cases <- boot(gsmfit, statfun, newdata = newdata)

# nonparameteric bootstrap residuals
set.seed(0)
boot.resid <- boot(gsmfit, statfun, newdata = newdata, method = "resid")

# parameteric bootstrap residuals
set.seed(0)
boot.param <- boot(gsmfit, statfun, newdata = newdata,  method = "param")
                   
# plot results
par(mfrow = c(1, 3))
plotci(newdata$x, boot.cases$t0, ci = boot.cases$ci, main = "Cases")
plotci(newdata$x, boot.resid$t0, ci = boot.resid$ci, main = "Residuals")
plotci(newdata$x, boot.param$t0, ci = boot.param$ci, main = "Parametric")
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
