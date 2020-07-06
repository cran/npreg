\name{ss}
\alias{ss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a Smoothing Spline
}
\description{
Fits a smoothing spline with the smoothing parameter selected via one of eight methods: GCV, OCV, GACV, ACV, REML, ML, AIC, or BIC.
}
\usage{
ss(x, y = NULL, w = NULL, df, spar = NULL, lambda = NULL,
   method = c("GCV", "OCV", "GACV", "ACV", "REML", "ML", "AIC", "BIC"), 
   m = 2L, periodic = FALSE, all.knots = FALSE, nknots = .nknots.smspl, 
   knots = NULL, keep.data = TRUE, df.offset = 0, penalty = 1, 
   control.spar = list(), tol = 1e-6 * IQR(x))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Predictor vector of length \eqn{n}. Can also input a list or a two-column matrix specifying x and y.
}
  \item{y}{
Response vector of length \eqn{n}. If y is missing or NULL, the responses are assumed to be specified by x, with x the index vector.  
}
  \item{w}{
Weights vector of length \eqn{n}. Defaults to all 1.  
}
  \item{df}{
Equivalent degrees of freedom (trace of the smoother matrix). Must be in \eqn{(1,nx]}, where \eqn{nx} is the number of unique x values, see below.  
}
  \item{spar}{
Smoothing parameter. Typically (but not always) in the range \eqn{(0,1]}. If specified \code{lambda = 256^(3*(spar-1))}.
}
  \item{lambda}{
Computational smoothing parameter. This value is weighted by \eqn{n} to form the penalty coefficient (see Details). Ignored if \code{spar} is provided.
}
  \item{method}{
Method for selecting the smoothing parameter. Ignored if \code{spar} or \code{lambda} is provided.
}
  \item{m}{
Penalty order (integer). The penalty functional is the integrated squared \eqn{m}-th derivative of the function. Defaults to \eqn{m = 2}, which is a cubic smoothing spline. Set \eqn{m = 1} for a linear smoothing spline or \eqn{m = 3} for a quintic smoothing spline.
}
  \item{periodic}{
Logical. If \code{TRUE}, the estimated function \eqn{f(x)} is constrained to be periodic, i.e., \eqn{f(a) = f(b)} where \eqn{a = \min(x)} and \eqn{b = \max(x)}.
}
  \item{all.knots}{
If \code{TRUE}, all distinct points in x are used as knots. If \code{FALSE} (default), a sequence knots is placed at the quantiles of the unique x values; in this case, the input \code{nknots} specifies the number of knots in the sequence. Ignored if the knot values are input using the \code{knots} argument.
}
  \item{nknots}{
Positive integer or function specifying the number of knots. Ignored if either \code{all.knots = TRUE} or the knot values are input using the \code{knots} argument.
}
  \item{knots}{
Vector of knot values for the spline. Should be unique and within the range of the x values (to avoid a warning). 
}
  \item{keep.data}{
Logical. If \code{TRUE}, the original data as a part of the output object.  
}
  \item{df.offset}{
Allows the degrees of freedom to be increased by \code{df.offset} in the GCV criterion.  
}
  \item{penalty}{
The coefficient of the penalty for degrees of freedom in the GCV criterion.  
}
  \item{control.spar}{
Optional list with named components controlling the root finding when the smoothing parameter spar is computed, i.e., missing or NULL, see below.

\bold{Note} that spar is only searched for in the interval \eqn{[lower, upper]}.
  \describe{
  \item{lower:}{
  lower bound for spar; defaults to 0.
  }
  \item{upper:}{
  upper bound for spar; defaults to 1.
  }
  \item{tol:}{
  the absolute precision (\bold{tol}erance) used by \code{\link{optimize}}; defaults to 1e-8.
  }
  }
}
  \item{tol}{
Tolerance for same-ness or uniqueness of the x values. The values are binned into bins of size tol and values which fall into the same bin are regarded as the same. Must be strictly positive (and finite).
}

}
\details{
Inspired by the \code{\link{smooth.spline}} function in R's \bold{stats} package.

Neither \code{x} nor \code{y} are allowed to containing missing or infinite values.

The \code{x} vector should contain at least \eqn{2m} distinct values. 'Distinct' here is controlled by \code{tol}: values which are regarded as the same are replaced by the first of their values and the corresponding \code{y} and \code{w} are pooled accordingly.

Unless \code{lambda} has been specified instead of \code{spar}, the computational \eqn{\lambda} used (as a function of \code{spar}) is \eqn{\lambda = 256^(3*(spar - 1))}. 

If \code{spar} and \code{lambda} are missing or \code{NULL}, the value of \code{df} is used to determine the degree of smoothing. If \code{df} is missing as well, the specified \code{method} is used to determine \eqn{\lambda}. 

Letting \eqn{f_i = f(x_i)}, the function is represented as \deqn{f = X \beta + Z \gamma} where the basis functions in \eqn{X} span the null space (i.e., functions with \eqn{m}-th derivative of zero), and \eqn{Z} contains the reproducing kernel function of the contrast space evaluated at all combinations of observed data points and knots, i.e., \eqn{Z[i,j] = \rho(x_i, k_j)} where \eqn{\rho} is the kernel function and \eqn{k_j} is the \eqn{j}-th knot. The vectors \eqn{\beta} and \eqn{\gamma} contain unknown basis function coefficients. 
Letting \eqn{M =  (X, Z)} and \eqn{\theta = (\beta', \gamma')'}, the penalized least squares problem has the form
\deqn{
(y - M \theta)' W (y - M \theta) + n \lambda \gamma' Q \gamma
}
where \eqn{W} is a diagonal matrix containg the weights, and \eqn{Q} is the penalty matrix. Note that \eqn{Q[i,j] = \rho(k_i, k_j)} contains the reproducing kernel function evaluated at all combinations of knots. The optimal coefficients are the solution to 
\deqn{
(M' W M + n \lambda P) \theta = M' W y
}
where \eqn{P} is the penalty matrix \eqn{Q} augmented with zeros corresponding to the \eqn{\beta} in \eqn{\theta}.

}
\value{
An object of class "ss" with components:
  \item{x}{the distinct \code{x} values in increasing order; see Note.}
  \item{y}{the fitted values corresponding to \code{x}.}
  \item{w}{the weights used at the unique values of \code{x}.}
  \item{yin}{the \code{y} values used at the unique \code{y} values.}
  \item{tol}{the \code{tol} argument (whose default depends on \code{x}).}
  \item{data}{only if keep.data = TRUE: itself a list with components \code{x}, \code{y} and \code{w} (if applicable). These are the original \eqn{(x_i,y_i,w_i), i = 1, \ldots, n}, values where \code{data$x} may have repeated values and hence be longer than the above \code{x} component; see details.}
  \item{lev}{leverages, the diagonal values of the smoother matrix.}
  \item{cv.crit}{cross-validation score.}
  \item{pen.crit}{the penalized criterion, a non-negative number; simply the (weighted) residual sum of squares (RSS).}
  \item{crit}{the criterion value minimized in the underlying \code{df2lambda} function. When \code{df} is provided, the criterion is \eqn{[tr(S_{\lambda}) - df]^2}.}
  \item{df}{equivalent degrees of freedom used.}
  \item{spar}{the value of \code{spar} computed or given, i.e., \eqn{s = 1 + \log_{256}(\lambda)/3}}
  \item{lambda}{the value of \eqn{\lambda} corresponding to \code{spar}, i.e., \eqn{\lambda = 256^{3*(s-1)}}.}
  \item{fit}{list for use by \code{\link{predict.ss}}, with components
  \describe{
    \item{n:}{number of observations.}
    \item{knot:}{the knot sequence.}
    \item{nk:}{number of coefficients (# knots plus \eqn{m}).}
    \item{coef:}{coefficients for the spline basis used.}
    \item{min, range:}{numbers giving the corresponding quantities of \code{x}}
    \item{m:}{spline penalty order (same as input \code{m})}
    \item{periodic:}{is spline periodic?}
    \item{cov.sqrt}{square root of covariance matrix of \code{coef} such that \code{tcrossprod(coef)} reconstructs the covariance matrix.}
    \item{weighted}{were weights \code{w} used in fitting?}
  }
  }
  \item{call}{the matched call.}
  \item{sigma}{estimated error standard deviation.}
  \item{logLik}{log-likelihood (if \code{method} is REML or ML).}
  \item{aic}{Akaike's Information Criterion (if \code{method} is AIC).}
  \item{bic}{Bayesian Information Criterion (if \code{method} is BIC).}
  \item{penalty}{smoothness penalty \eqn{\gamma' Q \gamma}, which is the integrated squared \eqn{m}-th derivative of the estimated function \eqn{f(x)}.}
  \item{method}{smoothing parameter selection method. Will be \code{NULL} if \code{df}, \code{spar}, or \code{lambda} is provided.}
}
\references{
https://stat.ethz.ch/R-manual/R-devel/library/stats/html/smooth.spline.html

Craven, P. and Wahba, G. (1979). Smoothing noisy data with spline functions: Estimating the correct degree of smoothing by the method of generalized cross-validation. \emph{Numerische Mathematik, 31}, 377-403. https://doi.org/10.1007/BF01404567

Gu, C. (2013). Smoothing spline ANOVA models, 2nd edition. New York: Springer. https://doi.org/10.1007/978-1-4614-5369-7

Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), SAGE Research Methods Foundations. https://doi.org/10.4135/9781526421036885885

Wahba, G. (1985). A comparison of {GCV} and {GML} for choosing the smoothing parameters in the generalized spline smoothing problem. \emph{The Annals of Statistics, 4}, 1378-1402. https://doi.org/10.1214/aos/1176349743
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The number of unique x values, nx, are determined by the tol argument, equivalently to

\code{nx <- sum(!duplicated( round((x - mean(x)) / tol) ))}

In this case where not all unique x values are used as knots, the result is not a smoothing spline in the strict sense, but very close unless a small smoothing parameter (or large \code{df}) is used.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\section{Methods }{
The smoothing parameter can be selected using one of eight methods: \cr
Generalized Cross-Validation (GCV) \cr
Ordinary Cross-Validation (OCV) \cr
Generalized Approximate Cross-Validation (GACV) \cr
Approximate Cross-Validation (ACV) \cr
Restricted Maximum Likelihood (REML) \cr
Maximum Likelihood (ML) \cr
Akaike's Information Criterion (AIC) \cr
Bayesian Information Criterion (BIC)

}

\seealso{
\code{\link{summary.ss}} for summarizing \code{ss} objects.

\code{\link{predict.ss}} for predicting from \code{ss} objects.

\code{\link{sm}} for fitting smooth models with multiple predictors of mixed types (Gaussian response).

\code{\link{gsm}} for fitting generalized smooth models with multiple predictors of mixed types (non-Gaussian response).
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
ss.GCV

# OCV selection
ss.OCV <- ss(x, y, method = "OCV", nknots = 10)
ss.OCV

# GACV selection
ss.GACV <- ss(x, y, method = "GACV", nknots = 10)
ss.GACV

# ACV selection
ss.ACV <- ss(x, y, method = "ACV", nknots = 10)
ss.ACV

# ML selection
ss.ML <- ss(x, y, method = "ML", nknots = 10)
ss.ML

# REML selection
ss.REML <- ss(x, y, method = "REML", nknots = 10)
ss.REML

# AIC selection
ss.AIC <- ss(x, y, method = "AIC", nknots = 10)
ss.AIC

# BIC selection
ss.BIC <- ss(x, y, method = "BIC", nknots = 10)
ss.BIC

# compare results
mean( ( fx - ss.GCV$y )^2 )
mean( ( fx - ss.OCV$y )^2 )
mean( ( fx - ss.GACV$y )^2 )
mean( ( fx - ss.ACV$y )^2 )
mean( ( fx - ss.ML$y )^2 )
mean( ( fx - ss.REML$y )^2 )
mean( ( fx - ss.AIC$y )^2 )
mean( ( fx - ss.BIC$y )^2 )

# plot results
plot(x, y)
rlist <- list(ss.GCV, ss.OCV, ss.GACV, ss.ACV,
              ss.REML, ss.ML, ss.AIC, ss.BIC)
for(j in 1:length(rlist)){
   lines(rlist[[j]], lwd = 2, col = j)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
