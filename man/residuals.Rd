\name{residuals}
\alias{residuals.ss}
\alias{residuals.sm}
\alias{residuals.gsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract Model Residuals
}
\description{
Extracts the residuals from a fit smoothing spline ("ss"), smooth model ("sm"), or generalized smooth model ("gsm") object. 
}
\usage{
\method{residuals}{ss}(object, type = c("working", "response", "deviance", 
                           "pearson", "partial"), ...)

\method{residuals}{sm}(object, type = c("working", "response", "deviance", 
                           "pearson", "partial"), ...)
                           
\method{residuals}{gsm}(object, type = c("deviance", "pearson", "working", 
                           "response", "partial"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class "ss", "sm", or "gsm"
}
  \item{type}{
type of residuals
}
  \item{...}{
other arugments (currently ignored)  
}
}
\details{
For objects of class \code{\link{ss}} and \code{\link{sm}} \cr
* the working and response residuals are defined as 'observed - fitted' \cr
* the deviance and Pearson residuals multiply the working residuals by \code{sqrt(weights(object))}

For objects of class \code{\link{gsm}}, the residual types are the same as those produced by the \code{\link{residuals.glm}} function
}
\value{
Residuals from \code{object}
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
res.ss <- residuals(mod.ss)

# smooth model
mod.sm <- sm(y ~ x, knots = 10)
res.sm <- residuals(mod.sm)

# generalized smooth model (family = gaussian)
mod.gsm <- gsm(y ~ x, knots = 10)
res.gsm <- residuals(mod.gsm)

# y = fitted + residuals
mean((y - fitted(mod.ss) - res.ss)^2)
mean((y - fitted(mod.sm) - res.sm)^2)
mean((y - fitted(mod.gsm) - res.gsm)^2)
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
