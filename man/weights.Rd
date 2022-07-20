\name{weights}
\alias{weights.ss}
\alias{weights.sm}
\alias{weights.gsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract Smooth Model Weights
}
\description{
Extracts prior weights from a fit smoothing spline (fit by \code{\link{ss}}), smooth model (fit by \code{\link{sm}}), or generalized smooth model (fit by \code{\link{gsm}}).
}
\usage{
\method{weights}{ss}(object, ...)

\method{weights}{sm}(object, ...)

\method{weights}{gsm}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class "gsm" output by the \code{\link{gsm}} function, "sm" output by the \code{\link{sm}} function, or "ss" output by the \code{\link{ss}} function
}
  \item{...}{
other arugments (currently ignored)  
}
}
\details{
Returns the "prior weights", which are user-specified via the \code{w} argument (of the \code{\link{ss}} function) or the \code{weights} argument (of the \code{\link{sm}} and \code{\link{gsm}} functions). If no prior weights were supplied, returns the (default) unit weights, i.e., \code{rep(1, nobs)}.
}
\value{
Prior weights extracted from \code{object}
}
\references{
Chambers, J. M. and Hastie, T. J. (1992) \emph{Statistical Models in S}. Wadsworth & Brooks/Cole.

Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), SAGE Research Methods Foundations. \doi{10.4135/9781526421036885885}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ss}}, \code{\link{sm}}, \code{\link{gsm}}
}
\examples{
# generate weighted data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
w <- rep(5:15, length.out = n)
fx <- 2 + 3 * x + sin(2 * pi * x)
y <- fx + rnorm(n, sd = 0.5 / sqrt(w))

# smoothing spline
mod.ss <- ss(x, y, w, nknots = 10)
w.ss <- weights(mod.ss)

# smooth model
mod.sm <- sm(y ~ x, weights = w, knots = 10)
w.sm <- weights(mod.sm)

# generalized smooth model (family = gaussian)
mod.gsm <- gsm(y ~ x, weights = w, knots = 10)
w.gsm <- weights(mod.gsm)

# note: weights are internally rescaled such as 
w0 <- w / mean(w)
max(abs(w0 - w.ss))
max(abs(w0 - w.sm))
max(abs(w0 - w.gsm))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
