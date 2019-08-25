\name{plotci}
\alias{plotci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic X-Y Plotting with Confidence Intervals
}
\description{
Modification to the \code{\link{plot}} function that adds confidence intervals. The CIs can be plotted using polygons (default) or error bars.
}
\usage{
plotci(x, y, se, level = 0.95, crit.val = NULL, 
       add = FALSE, col.ci = "gray80", alpha = 0.5, 
       bars = FALSE, bw = 0.2, linkinv = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a vector of 'x' values (\eqn{n} by 1). If \code{y} is missing, the \code{x} input can be a list or matrix containing the \code{x}, \code{y}, and \code{se} arguments.
}
  \item{y}{
a vector of 'y' values (\eqn{n} by 1).  
}
  \item{se}{
a vector of standard error values (\eqn{n} by 1).
}
  \item{level}{
confidence level for the intervals (between 0 and 1).  
}
  \item{crit.val}{
an optional critical value for the intervals. If provided, the \code{level} input is ignored. See Details.  
}
  \item{add}{
a switch controlling whether a new plot should be created (via a call to \code{\link{plot}}) or if the plot should be added to the current plot (via a call to \code{\link{lines}}).
}
  \item{col.ci}{
a character specifying the color for plotting the intervals.  
}
  \item{alpha}{
a scalar between 0 and 1 controlling the transparency of the intervals.
}
  \item{bars}{
a switch controlling whether the intervals should be plotted as bars or polygons.  
}
  \item{bw}{
a positive scalar controlling the bar width. Ignored if \code{bars = FALSE}.  
}
  \item{linkinv}{
an inverse link function for the plotting. If provided, the function plots \code{x} versus \code{linkinv(y)} and the intervals are similary transformed. 
}
  \item{...}{
extra arguments passed to the \code{plot} or \code{lines} function.
}

}
\details{
This function plots \code{x} versus \code{y} with confidence intervals. The CIs have the form \cr \code{lwr = y - crit.val * se}\cr \code{upr = y + crit.val * se}\cr where \code{crit.val} is the critical value.

If \code{crit.val = NULL}, the critival value is determined from the \code{level} input as\cr \code{crit.val <- qnorm(1-(1-level)/2)}\cr where \code{\link{qnorm}} is the quantile function for the standard normal distribution.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

\seealso{
This function is used by \code{\link{plot.ss}} to plot smoothing spline fits.%and \code{\link{plot.ssnet}} to plot smoothing spline fits.
}
\examples{
# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
fx <- 2 + 3 * x + sin(2 * pi * x)
y <- fx + rnorm(n, sd = 0.5)

# fit smooth model
smod <- sm(y ~ x, knots = 10)

# plot fit with 95\% CI polygon
plotci(x, smod$fitted.values, smod$se.fit)

# plot fit with 95\% CI bars
plotci(x, smod$fitted.values, smod$se.fit, bars = TRUE)

# plot fit +/- 1 SE
plotci(x, smod$fitted.values, smod$se.fit, crit.val = 1, bars = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }% use one of  RShowDoc("KEYWORDS")
\keyword{ dplot }% __ONLY ONE__ keyword per line
