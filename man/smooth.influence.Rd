\name{smooth.influence}
\alias{smooth.influence}
\alias{influence.ss}
\alias{influence.sm}
\alias{influence.gsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Regression Diagnostics
}
\description{
These functions provide the basic quantities that are used to form a variety of diagnostics for checking the quality of a fit smoothing spline (fit by \code{\link{ss}}), smooth model (fit by \code{\link{sm}}), or generalized smooth model (fit by \code{\link{gsm}}).
}
\usage{
\method{influence}{ss}(model, do.coef = TRUE, ...)
\method{influence}{sm}(model, do.coef = TRUE, ...)
\method{influence}{gsm}(model, do.coef = TRUE, ...)

smooth.influence(model, do.coef = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of class "gsm" output by the \code{\link{gsm}} function, "sm" output by the \code{\link{sm}} function, or "ss" output by the \code{\link{ss}} function
}
  \item{do.coef}{
logical indicating if the changed \code{coefficients} are desired (see Details).   
}
  \item{...}{
additional arguments (currently ignored)  
}
}
\details{
Inspired by \code{\link{influence}} and \code{\link{lm.influence}} functions in R's \bold{stats} package.

The functions documented in \code{\link{smooth.influence.measures}} provide a more user-friendly way of computing a variety of regression diagnostics.

For non-Gaussian \code{gsm} objects, these regression diagnostics are based on one-step approximations, which may be inadequate if a case has high influence. 

For all models, the diagostics are computed assuming that the smoothing parameters are fixed at the given values. 
}
\value{
A list with the components
\item{hat }{a vector containing the leverages, i.e., the diagonals of the smoothing matrix}
\item{coefficients }{if \code{do.coef} is true, a matrix whose i-th row contains the change in the estimated coefficients which results when the i-th case is excluded from the fitting. }
\item{deviance }{ a vector whose i-th entry contains the deviance which results when the i-th case is excluded from the fitting.}
\item{df }{ a vector whose i-th entry contains the effective degrees-of-freedom which results when the i-th case is excluded from the fitting.}
\item{sigma }{a vector whose i-th element contains the estimate of the residual standard deviation obtained when the i-th case is excluded from the fitting. }
\item{wt.res }{a vector of \emph{weighted} (or for class \code{gsm} rather \emph{deviance}) residuals.}
}
\references{
See the list in the documentation for \code{\link{influence.measures}}

Chambers, J. M. (1992) \emph{Linear models.} Chapter 4 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The \code{coefficients} returned by \code{\link{smooth.influence}} (and the corresponding functions S3 \code{influence} methods) are the \emph{change} in the coefficients which result from dropping each case, i.e., \eqn{\theta - \theta_i}, where \eqn{\theta} are the original coefficients obtained from the full sample of \eqn{n} observations and \eqn{\theta_i} are the coefficients that result from dropping the i-th case.
}

\section{Warning }{
The approximations used for \code{gsm} objects can result in \code{sigma} estimates being \code{NaN}.
}

\seealso{
\code{\link{ss}}, \code{\link{sm}}, \code{\link{gsm}} for modeling functions

\code{\link{smooth.influence.measures}} for convenient summary

\code{\link{diagnostic.plots}} for regression diagnostic plots
}
\examples{
# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
fx <- 2 + 3 * x + sin(2 * pi * x)
y <- fx + rnorm(n, sd = 0.5)

# fit models
mod.ss <- ss(x, y, nknots = 10)
mod.sm <- sm(y ~ x, knots = 10)
mod.gsm <- gsm(y ~ x, knots = 10)

# calculate influence
infl.ss <- influence(mod.ss)
infl.sm <- influence(mod.sm)
infl.gsm <- influence(mod.gsm)

# compare hat
mean((infl.ss$hat - infl.sm$hat)^2)
mean((infl.ss$hat - infl.gsm$hat)^2)
mean((infl.sm$hat - infl.gsm$hat)^2)

# compare deviance
mean((infl.ss$deviance - infl.sm$deviance)^2)
mean((infl.ss$deviance - infl.gsm$deviance)^2)
mean((infl.sm$deviance - infl.gsm$deviance)^2)

# compare df
mean((infl.ss$df - infl.sm$df)^2)
mean((infl.ss$df - infl.gsm$df)^2)
mean((infl.sm$df - infl.gsm$df)^2)

# compare sigma
mean((infl.ss$sigma - infl.sm$sigma)^2)
mean((infl.ss$sigma - infl.gsm$sigma)^2)
mean((infl.sm$sigma - infl.gsm$sigma)^2)

# compare residuals
mean((infl.ss$wt.res - infl.sm$wt.res)^2)
mean((infl.ss$wt.res - infl.gsm$dev.res)^2)
mean((infl.sm$wt.res - infl.gsm$dev.res)^2)

# NOTE: ss() coef only comparable to sm() and gsm() after rescaling
scale.sm <- rep(c(1, mod.sm$specs$thetas), times = c(2, 10))
scale.gsm <- rep(c(1, mod.gsm$specs$thetas), times = c(2, 10))
mean((coef(mod.ss) / scale.sm - coef(mod.sm))^2)
mean((coef(mod.ss) / scale.gsm - coef(mod.gsm))^2)
mean((coef(mod.sm) - coef(mod.gsm))^2)

# infl.ss$coefficients are *not* comparable to others
mean((infl.ss$coefficients - infl.sm$coefficients)^2)
mean((infl.ss$coefficients - infl.gsm$coefficients)^2)
mean((infl.sm$coefficients - infl.gsm$coefficients)^2)

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
