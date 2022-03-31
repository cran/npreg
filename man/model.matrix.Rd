\name{model.matrix}
\alias{model.matrix.ss}
\alias{model.matrix.sm}
\alias{model.matrix.gsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct Design Matrix for Fit Model
}
\description{
\code{model.matrix} returns the design (or model) matrix used by the input \code{object} to produce the fitted values (for objects of class \code{\link{ss}} or \code{\link{sm}}) or the linear predictors (for objects of class \code{\link{gsm}}).
}
\usage{
\method{model.matrix}{ss}(object, ...)

\method{model.matrix}{sm}(object, ...)

\method{model.matrix}{gsm}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class \code{\link{ss}}, \code{\link{sm}}, or \code{\link{gsm}}
}
  \item{...}{
additional arguments (currently ignored)  
}
}
\details{
For \code{\link{ss}} objects, the \code{\link{basis.poly}} function is used to construct the design matrix. 

For \code{\link{sm}} objects, the \code{\link{predict.sm}} function with option \code{design = TRUE} is used to construct the design matrix. 

For \code{\link{gsm}} objects, the \code{\link{predict.gsm}} function with option \code{design = TRUE} is used to construct the design matrix.
}
\value{
The design matrix that is post-multiplied by the coefficients to produce the fitted values (or linear predictors).
}
\references{
Chambers, J. M. (1992) Data for models. Chapter 3 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), SAGE Research Methods Foundations. \doi{10.4135/9781526421036885885}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{basis.poly}} for the smoothing spline basis

\code{\link{predict.sm}} for predicting from smooth models

\code{\link{predict.gsm}} for predicting from generalized smooth models
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
X.ss <- model.matrix(mod.ss)
mean((mod.ss$y - X.ss \%*\% mod.ss$fit$coef)^2)


# smooth model
mod.sm <- sm(y ~ x, knots = 10)
X.sm <- model.matrix(mod.sm)
mean((mod.sm$fitted.values - X.sm \%*\% mod.sm$coefficients)^2)

# generalized smooth model (family = gaussian)
mod.gsm <- gsm(y ~ x, knots = 10)
X.gsm <- model.matrix(mod.gsm)
mean((mod.gsm$linear.predictors - X.gsm \%*\% mod.gsm$coefficients)^2)

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
