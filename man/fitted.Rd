\name{fitted}
\alias{fitted.ss}
\alias{fitted.sm}
\alias{fitted.gsm}
\alias{fitted.values.ss}
\alias{fitted.values.sm}
\alias{fitted.values.gsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract Model Fitted Values
}
\description{
Extracts the fitted values from a fit smoothing spline ("ss"), smooth model ("sm"), or generalized smooth model ("gsm") object. 
}
\usage{
\method{fitted}{ss}(object, ...)
\method{fitted.values}{ss}(object, ...)

\method{fitted}{sm}(object, ...)
\method{fitted.values}{sm}(object, ...)

\method{fitted}{gsm}(object, ...)
\method{fitted.values}{gsm}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class "ss", "sm", or "gsm"
}
  \item{...}{
other arugments (currently ignored)  
}
}
\details{
For objects of class \code{\link{ss}}, fitted values are predicted via \code{predict(object, object$data$x)$y}

For objects of class \code{\link{sm}}, fitted values are extracted via \code{object$fitted.values}

For objects of class \code{\link{gsm}}, fitted values are computed via \code{ginv(object$linear.predictors)} where \code{ginv = object$family$linkinv}
}
\value{
Fitted values extracted (or predicted) from \code{object}
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
# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
fx <- 2 + 3 * x + sin(2 * pi * x)
y <- fx + rnorm(n, sd = 0.5)

# smoothing spline
mod.ss <- ss(x, y, nknots = 10)
fit.ss <- fitted(mod.ss)

# smooth model
mod.sm <- sm(y ~ x, knots = 10)
fit.sm <- fitted(mod.sm)

# generalized smooth model (family = gaussian)
mod.gsm <- gsm(y ~ x, knots = 10)
fit.gsm <- fitted(mod.gsm)

# compare fitted values
mean((fit.ss - fit.sm)^2)
mean((fit.ss - fit.gsm)^2)
mean((fit.sm - fit.gsm)^2)
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
