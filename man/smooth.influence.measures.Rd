\name{smooth.influence.measures}
\alias{smooth.influence.measures}
\alias{rstandard.ss}
\alias{rstandard.sm}
\alias{rstandard.gsm}
\alias{rstudent.ss}
\alias{rstudent.sm}
\alias{rstudent.gsm}
\alias{dfbeta.ss}
\alias{dfbeta.sm}
\alias{dfbeta.gsm}
\alias{dfbetas.ss}
\alias{dfbetas.sm}
\alias{dfbetas.gsm}
\alias{cov.ratio}
\alias{cooks.distance.ss}
\alias{cooks.distance.sm}
\alias{cooks.distance.gsm}
\alias{hatvalues.ss}
\alias{hatvalues.sm}
\alias{hatvalues.gsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Regression Deletion Diagnostics
}
\description{
These functions compute several regression (leave-one-out deletion) diagnostics for a fit smoothing spline (fit by \code{\link{ss}}), smooth model (fit by \code{\link{sm}}), or generalized smooth model (fit by \code{\link{gsm}}).
}
\usage{
smooth.influence.measures(model, infl = smooth.influence(model))

\method{rstandard}{ss}(model, infl = NULL, sd = model$sigma, 
          type = c("sd.1", "predictive"), ...)
\method{rstandard}{sm}(model, infl = NULL, sd = model$sigma, 
          type = c("sd.1", "predictive"), ...)
\method{rstandard}{gsm}(model, infl = NULL, 
          type = c("deviance", "pearson"), ...)

\method{rstudent}{ss}(model, infl = influence(model, do.coef = FALSE), 
         res = infl$wt.res, ...)
\method{rstudent}{sm}(model, infl = influence(model, do.coef = FALSE), 
         res = infl$wt.res, ...)
\method{rstudent}{gsm}(model, infl = influence(model, do.coef = FALSE), ...)

\method{dfbeta}{ss}(model, infl = NULL, ...)
\method{dfbeta}{sm}(model, infl = NULL, ...)
\method{dfbeta}{gsm}(model, infl = NULL, ...)

\method{dfbetas}{ss}(model, infl = smooth.influence(model, do.coef = TRUE), ...)
\method{dfbetas}{sm}(model, infl = smooth.influence(model, do.coef = TRUE), ...)
\method{dfbetas}{gsm}(model, infl = smooth.influence(model, do.coef = TRUE), ...)

cov.ratio(model, infl = smooth.influence(model, do.coef = FALSE),
          res = weighted.residuals(model))

\method{cooks.distance}{ss}(model, infl = NULL, res = weighted.residuals(model), 
               sd = model$sigma, hat = hatvalues(model), ...)
\method{cooks.distance}{sm}(model, infl = NULL, res = weighted.residuals(model), 
               sd = model$sigma, hat = hatvalues(model), ...)
\method{cooks.distance}{gsm}(model, infl = NULL, res = residuals(model, type = "pearson"), 
               dispersion = model$dispersion, hat = hatvalues(model), ...)

\method{hatvalues}{ss}(model, ...)
\method{hatvalues}{sm}(model, ...)
\method{hatvalues}{gsm}(model, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of class "gsm" output by the \code{\link{gsm}} function, "sm" output by the \code{\link{sm}} function, or "ss" output by the \code{\link{ss}} function
}
  \item{infl}{
influence structure as returned by \code{\link{smooth.influence}}  
}
  \item{res}{
(possibly weighted) residuals with proper defaults    
}
  \item{sd}{
standard deviation to use, see defaults
}
  \item{dispersion}{
dispersion (for \code{\link{gsm}} objects) to use, see defaults  
}
  \item{hat}{
hat values \eqn{S_{ii}}, see defaults  
}
  \item{type}{
type of residuals for \code{rstandard}  
}
  \item{...}{
additional arguments (currently ignored)  
}
}
\details{
Inspired by \code{\link{influence.measures}} and related functions in R's \bold{stats} package.

The function \code{smooth.influence.measures} produces a class "infl" object, which displays the DFBETAS for each coefficient, DFFITS, covariance ratios, Cook's distance, and the diagonals of the smoothing matrix. Cases which are influential with respect to any of these measures are marked with an asterisk.

The S3 methods \code{dfbetas}, \code{dffits}, \code{covratio}, and \code{cooks.distance} provide direct access to the corresponding diagnostic quantities. The S3 methods \code{rstandard} and \code{rstudent} give the standardized and Studentized residuals, respectively. (These re-normalize the residuals to have unit variance, using an overall and leave-one-out measure of the error variance, respectively.)

Values for generalized smoothing models are approximations, as described in Williams (1987) (except that Cook's distances are scaled as \eqn{F} rather than chi-square values). THe approximations can be poor when some cases have large influence.

The optional \code{infl}, \code{res}, and \code{sd} arguments are there to encourage the use of these direct access functions in situations where the underlying basic influence measures, e.g., from \code{\link{smooth.influence}}, are already available.

For \code{ss} and \code{sm} objects, the code \code{rstandard(*, type = "predictive")} returns the leave-one-out (ordinary) cross-validation residuals, and the PRESS (PREdictive Sum of Squares) statistic is defined as

\code{PRESS <- sum(rstandard(model, type = "predictive")^2)}

Note that \code{OCV = PRESS / n}, where OCV = ordinary cross-validation criterion

}
\references{
See references listed in \code{\link{influence.measures}}

Williams, D. A. (1987). Generalized linear model diagnostics using the deviance and single case deletions. \emph{Applied Statistics, 36}, 181-191. \doi{10.2307/2347550}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Note: the \code{\link{dffits}} function in R's \bold{stats} package can be used with the following syntax

\code{dffits(model, infl = smooth.influence(model, do.coef = FALSE),
       res = weighted.residuals(model))}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ss}}, \code{\link{sm}}, \code{\link{gsm}} for modeling functions

\code{\link{smooth.influence}} for some basic influence information

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
infl.ss <- smooth.influence.measures(mod.ss)
infl.sm <- smooth.influence.measures(mod.sm)
infl.gsm <- smooth.influence.measures(mod.gsm)

# standardized residuals
rstan.ss <- rstandard(mod.ss)
rstan.sm <- rstandard(mod.sm)
rstan.gsm <- rstandard(mod.gsm)

# studentized residuals
rstud.ss <- rstudent(mod.ss)
rstud.sm <- rstudent(mod.sm)
rstud.gsm <- rstudent(mod.gsm)

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
