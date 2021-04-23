\name{predict.ss}
\alias{predict.ss}
%\alias{predict.ssnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict method for Smoothing Spline Fits
}
\description{
\code{predict} method for class "ss". %or class "ssnet".
}
\usage{
\method{predict}{ss}(object, x, deriv = 0, se.fit = TRUE, ...)
}
%\method{predict}{ssnet}(object, x, deriv = 0, se.fit = TRUE, ...)
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fit from \code{\link{ss}}. %or \code{\link{ssnet}}.
} 
  \item{x}{
the new values of x.
}
  \item{deriv}{
integer; the order of the derivative required.
}
  \item{se.fit}{
a switch indicating if standard errors are required.
}
  \item{...}{
additional arguments affecting the prediction produced (currently ignored).
}
}
\details{
Inspired by the \code{\link{predict.smooth.spline}} function in R's \bold{stats} package.
}
\value{
A list with components
  \item{x }{The input \code{x}.}
  \item{y }{The fitted values or derivatives at x.}
  \item{se }{The standard errors of the fitted values or derivatives (if requested).}
}
\references{
https://stat.ethz.ch/R-manual/R-devel/library/stats/html/predict.smooth.spline.html

Craven, P. and Wahba, G. (1979). Smoothing noisy data with spline functions: Estimating the correct degree of smoothing by the method of generalized cross-validation. \emph{Numerische Mathematik, 31}, 377-403. \doi{10.1007/BF01404567}

Gu, C. (2013). Smoothing spline ANOVA models, 2nd edition. New York: Springer. \doi{10.1007/978-1-4614-5369-7}

Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), \emph{SAGE Research Methods Foundations.} \doi{10.4135/9781526421036885885}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ss}} %and \code{\link{ssnet}}
}
\examples{
# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
fx <- 2 + 3 * x + sin(2 * pi * x)
y <- fx + rnorm(n, sd = 0.5)

# GCV selection (default)
ss.GCV <- ss(x, y, nknots = 10)

# get predictions and SEs (at design points)
fit <- predict(ss.GCV, x = x)
head(fit)

# compare to original fit
mean((fit$y - ss.GCV$y)^2)

# plot result (with default 95\% CI)
plotci(fit)

# estimate first derivative
d1 <- 3 + 2 * pi * cos(2 * pi * x)
fit <- predict(ss.GCV, x = x, deriv = 1)
head(fit)

# plot result (with default 95\% CI)
plotci(fit)
lines(x, d1, lty = 2)   # truth

}

%## sparse smoothing spline...
%
%# GCV selection (default)
%ss.GCV <- ssnet(x, y, nknots = 10)
%
%# get predictions and SEs (at design points)
%fit <- predict(ss.GCV, x = x)
%head(fit)
%
%# compare to original fit
%mean((fit$y - ss.GCV$y)^2)
%
%# plot result (with default 95\% CI)
%plotci(fit)
%
%# estimate first derivative
%d1 <- 3 + 2 * pi * cos(2 * pi * x)
%fit <- predict(ss.GCV, x = x, deriv = 1)
%head(fit)
%
%# plot result (with default 95\% CI)
%plotci(fit)
%lines(x, d1, lty = 2)   # truth

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
