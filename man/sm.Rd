\name{sm}
\alias{sm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a Smooth Model
}
\description{
Fits a semi- or nonparametric regression model with the smoothing parameter selected via one of eight methods: GCV, OCV, GACV, ACV, REML, ML, AIC, or BIC.
}
\usage{
sm(formula, data, weights, types = NULL, tprk = TRUE, knots = NULL,
   update = TRUE, spar = NULL, lambda = NULL, control = list(),
   method = c("GCV", "OCV", "GACV", "ACV", "REML", "ML", "AIC", "BIC"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. Uses the same syntax as \code{\link{lm}} and \code{\link{glm}}.
}
  \item{data}{
Optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{sm} is called.  
}
  \item{weights}{
Optional vector of weights to be used in the fitting process. If provided, weighted least squares is used. Defaults to all 1.
}
  \item{types}{
Named list giving the type of smooth to use for each predictor. If \code{NULL}, the type is inferred from the data. See "Types of Smooths" section for details.
}
  \item{tprk}{
Logical specifying how to parameterize smooth models with multiple predictors. If \code{TRUE} (default), a \bold{t}ensor \bold{p}roduct \bold{r}eproducing \bold{k}ernel function is used to represent the function. If \code{FALSE}, a tensor product of marginal kernel functions is used to represent the function. See the "Multiple Smooths" section for details.
}
  \item{knots}{
Spline knots for the estimation of the nonparametric effects. For models with multiple predictors, the knot specification will depend on the \code{tprk} input. See the "Choosing Knots" section for details
}
  \item{update}{
If \code{TRUE}, steps 1-2 of Gu and Wahba's (1991) algorithm 3.2 are used to update the "extra" smoothing parameters. If \code{FALSE}, only step 1 of algorithm 3.2 is used, so each effect is given equal influence on the penalty. Only applicable when multiple smooth terms are included. 
}
  \item{spar}{
Smoothing parameter. Typically (but not always) in the range \eqn{(0,1]}. If specified \code{lambda = 256^(3*(spar-1))}.
}
  \item{lambda}{
Computational smoothing parameter. This value is weighted by \eqn{n} to form the penalty coefficient (see Details). Ignored if \code{spar} is provided.
}
  \item{control}{
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
  \item{method}{
Method for selecting the smoothing parameter. Ignored if \code{lambda} is provided.
}
}
\details{
Letting \eqn{f_i = f(x_i)} with \eqn{x_i = (x_{i1}, \ldots, x_{ip})}, the function is represented as \deqn{f = X \beta + Z \gamma} where the basis functions in \eqn{X} span the null space (i.e., parametric effects), and \eqn{Z} contains the kernel function(s) of the contrast space (i.e., nonparametric effects) evaluated at all combinations of observed data points and knots. The vectors \eqn{\beta} and \eqn{\gamma} contain unknown basis function coefficients.

Letting \eqn{M =  (X, Z)} and \eqn{\theta = (\beta', \gamma')'}, the penalized least squares problem has the form
\deqn{
(y - M \theta)' W (y - M \theta) + n \lambda \gamma' Q \gamma
}
where \eqn{W} is a diagonal matrix containg the weights, and \eqn{Q} is the penalty matrix. The optimal coefficients are the solution to 
\deqn{
(M' W M + n \lambda P) \theta = M' W y
}
where \eqn{P} is the penalty matrix \eqn{Q} augmented with zeros corresponding to the \eqn{\beta} in \eqn{\theta}.
}
\value{
An object of class "sm" with components:
  \item{fitted.values}{the fitted values, i.e., predictions.}
  \item{se.fit}{the standard errors of the fitted values.}
  \item{sse}{the sum-of-squared errors.}
  \item{cv.crit}{the cross-validation criterion.}
  \item{df}{the estimated degrees of freedom (Df) for the fit model.}
  \item{nsdf}{the degrees of freedom (Df) for the null space.}
  \item{r.squared}{the observed coefficient of multiple determination.}
  \item{sigma}{the estimate of the error standard deviation.}
  \item{logLik}{the log-likelihood (if \code{method} is REML or ML).}
  \item{aic}{Akaike's Information Criterion (if \code{method} is AIC).}
  \item{bic}{Bayesian Information Criterion (if \code{method} is BIC).}
  \item{spar}{the value of \code{spar} computed or given, i.e., \eqn{s = 1 + \log_{256}(\lambda)/3}}
  \item{lambda}{the value of \eqn{\lambda} corresponding to \code{spar}, i.e., \eqn{\lambda = 256^{3*(s-1)}}.}
  \item{penalty}{the smoothness penalty \eqn{\gamma' Q \gamma}.}
  \item{coefficients}{the basis function coefficients used for the fit model.}
  \item{cov.sqrt}{the square-root of the covariance matrix of \code{coefficients}. Note: \code{tcrossprod(cov.sqrt)} reconstructs the covariance matrix. }
  \item{specs}{a list with information used for prediction purposes:
    \describe{
      \item{knots}{the spline knots used for each predictor.}
      \item{thetas}{the "extra" tuning parameters used to weight the penalties.}
      \item{xrng}{the ranges of the predictor variables.}
      \item{xlev}{the factor levels of the predictor variables (if applicable).}
      \item{tprk}{logical controlling the formation of tensor product smooths.}
    }
  }
  \item{data}{the data used to fit the model.}
  \item{types}{the type of smooth used for each predictor.}
  \item{terms}{the terms included in the fit model.}
  \item{method}{the \code{method} used for smoothing parameter selection. Will be \code{NULL} if \code{lambda} was provided.}
  \item{formula}{the formula specifying the fit model.}
  \item{call}{the matched call.}
}
\references{
Craven, P. and Wahba, G. (1979). Smoothing noisy data with spline functions: Estimating the correct degree of smoothing by the method of generalized cross-validation. \emph{Numerische Mathematik, 31}, 377-403. https://doi.org/10.1007/BF01404567

Gu, C. (2013). Smoothing spline ANOVA models, 2nd edition. New York: Springer. https://doi.org/10.1007/978-1-4614-5369-7

Gu, C. and Wahba, G. (1991). Minimizing GCV/GML scores with multiple smoothing parameters via the Newton method. \emph{SIAM Journal on Scientific and Statistical Computing, 12(2)}, 383-398. https://doi.org/10.1137/0912021

Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), SAGE Research Methods Foundations. https://doi.org/10.4135/9781526421036885885

%Helwig, N. E. (in prep). Spectrally sparse nonparametric regression via Elastic Net penalized smoothing splines.

%Zou, H. and Hastie, T. (2005). Regularization and variable selection via the Elastic Net. \emph{Journal of the Royal Statistical Society, Series B., 67}, 301-320. https://doi.org/10.1111/j.1467-9868.2005.00503.x
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\section{Types of Smooths }{

The following codes specify the spline types:
\tabular{ll}{
par \tab Parametric effect (factor, integer, or numeric). \cr
nom \tab Nominal smoothing spline (unordered factor). \cr
ord \tab Ordinal smoothing spline (ordered factor). \cr
lin \tab Linear smoothing spline (integer or numeric). \cr
cub \tab Cubic smoothing spline (integer or numeric). \cr
qui \tab Quintic smoothing spline (integer or numeric). \cr
per \tab Periodic smoothing spline (integer or numeric).\cr  
sph \tab Spherical spline (matrix with \eqn{d = 3} columns). \cr
tps \tab Thin-plate spline (matrix with \eqn{d \ge 1} columns).
}

For finer control of some specialized spline types:
\tabular{ll}{
per.lin \tab Linear periodic spline (\eqn{m = 1}). \cr
per.cub \tab Cubic periodic spline (\eqn{m = 2}). \cr
per.qui \tab Quintic periodic spline (\eqn{m = 3}). \cr
sph.lin \tab Linear spherical spline (\eqn{m = 1}). \cr
sph.cub \tab Cubic spherical spline (\eqn{m = 2}). \cr
sph.qui \tab Quintic spherical spline (\eqn{m = 3}). \cr
tps.lin \tab Linear thin-plate spline (\eqn{m = 1}). \cr
tps.cub \tab Cubic thin-plate spline (\eqn{m = 2}). \cr
tps.qui \tab Quintic thin-plate spline (\eqn{m = 3}). \cr
}

For details on the spline kernel functions, see \code{\link{basis_nom}} (nominal), \code{\link{basis_ord}} (ordinal), \code{\link{basis_poly}} (polynomial), \code{\link{basis_sph}} (spherical), and \code{\link{basis_tps}} (thin-plate).

}

\section{Choosing Knots }{
If \code{tprk = TRUE}, the four options for the \code{knots} input include: 
\tabular{ll}{
1. \tab a scalar giving the total number of knots to sample \cr
2. \tab a vector of integers indexing which rows of data are the knots \cr
3. \tab a list with named elements giving the marginal knot values for each predictor (to be combined via \code{\link{expand.grid}}) \cr
4. \tab a list with named elements giving the knot values for each predictor (requires the same number of knots for each predictor)
}

If \code{tprk = FALSE}, the three options for the \code{knots} input include:
\tabular{ll}{
1. \tab a scalar giving the common number of knots for each continuous predictor \cr
2. \tab a list with named elements giving the number of marginal knots for each predictor \cr
3. \tab a list with named elements giving the marginal knot values for each predictor
}

}

\section{Multiple Smooths }{
Suppose \code{formula = y ~ x1 + x2} so that the model contains additive effects of two predictor variables. 

The \eqn{k}-th predictor's marginal effect can be denoted as \deqn{f_k = X_k \beta_k + Z_k \gamma_k} where \eqn{X_k} is the \eqn{n} by \eqn{m_k} null space basis function matrix, and \eqn{Z_k} is the \eqn{n} by \eqn{r_k} contrast space basis function matrix. 

If \code{tprk = TRUE}, the null space basis function matrix has the form \eqn{X = [1, X_1, X_2]} and the contrast space basis function matrix has the form \deqn{Z = \theta_1 Z_1 + \theta_2 Z_2} where the \eqn{\theta_k} are the "extra" smoothing parameters. Note that \eqn{Z} is of dimension \eqn{n} by \eqn{r = r_1 = r_2}.

If \code{tprk = FALSE}, the null space basis function matrix has the form \eqn{X = [1, X_1, X_2]}, and the contrast space basis function matrix has the form \deqn{Z = [\theta_1 Z_1, \theta_2 Z_2]} where the \eqn{\theta_k} are the "extra" smoothing parameters. Note that \eqn{Z} is of dimension \eqn{n} by \eqn{r = r_1 + r_2}.
}


\seealso{
\code{\link{summary.sm}} for summarizing \code{sm} objects.

\code{\link{predict.sm}} for predicting from \code{sm} objects.

\code{\link{ss}} for fitting a smoothing spline with a single predictor (Gaussian response).

\code{\link{gsm}} for fitting generalized smooth models with multiple predictors of mixed types (non-Gaussian response).
}
\examples{
##########   EXAMPLE 1   ##########
### 1 continuous predictor

# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
fx <- 2 + 3 * x + sin(2 * pi * x)
y <- fx + rnorm(n, sd = 0.5)

# fit sm with 10 knots (tprk = TRUE)
sm.ssa <- sm(y ~ x, knots = 10)

# fit sm with 10 knots (tprk = FALSE)
sm.gam <- sm(y ~ x, knots = 10, tprk = FALSE)

# print both results (note: they are identical)
sm.ssa
sm.gam

# summarize both results (note: they are identical)
summary(sm.ssa)
summary(sm.gam)

# compare true MSE values (note: they are identical) 
mean( ( fx - sm.ssa$fit )^2 )
mean( ( fx - sm.gam$fit )^2 )



##########   EXAMPLE 2   ##########
### 1 continuous and 1 nominal predictor
### additive model

# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
z <- factor(sample(letters[1:3], size = n, replace = TRUE))
fun <- function(x, z){
  mu <- c(-2, 0, 2)
  zi <- as.integer(z)
  fx <- mu[zi] + 3 * x + sin(2 * pi * x)
}
fx <- fun(x, z)
y <- fx + rnorm(n, sd = 0.5)

# define marginal knots
probs <- seq(0, 0.9, by = 0.1)
knots <- list(x = quantile(x, probs = probs),
              z = letters[1:3])

# fit sm with specified knots (tprk = TRUE)
sm.ssa <- sm(y ~ x + z, knots = knots)

# fit sm with specified knots (tprk = FALSE)
sm.gam <- sm(y ~ x + z, knots = knots, tprk = FALSE)

# print both results (note: they are identical)
sm.ssa
sm.gam

# summarize both results (note: they are almost identical)
summary(sm.ssa)
summary(sm.gam)

# compare true MSE values (note: they are identical) 
mean( ( fx - sm.ssa$fit )^2 )
mean( ( fx - sm.gam$fit )^2 )



##########   EXAMPLE 3   ##########
### 1 continuous and 1 nominal predictor
### interaction model

# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
z <- factor(sample(letters[1:3], size = n, replace = TRUE))
fun <- function(x, z){
  mu <- c(-2, 0, 2)
  zi <- as.integer(z)
  fx <- mu[zi] + 3 * x + sin(2 * pi * x + mu[zi]*pi/4)
}
fx <- fun(x, z)
y <- fx + rnorm(n, sd = 0.5)

# define marginal knots
probs <- seq(0, 0.9, by = 0.1)
knots <- list(x = quantile(x, probs = probs),
              z = letters[1:3])

# fit sm with specified knots (tprk = TRUE)
sm.ssa <- sm(y ~ x * z, knots = knots)

# fit sm with specified knots (tprk = FALSE)
sm.gam <- sm(y ~ x * z, knots = knots, tprk = FALSE)

# print both results (note: they are slightly different)
sm.ssa
sm.gam

# summarize both results (note: they are slightly different)
summary(sm.ssa)
summary(sm.gam)

# compare true MSE values (note: they are slightly different) 
mean( ( fx - sm.ssa$fit )^2 )
mean( ( fx - sm.gam$fit )^2 )



##########   EXAMPLE 4   ##########
### 4 continuous predictors
### additive model

# generate data
set.seed(1)
n <- 100
fun <- function(x){
  sin(pi*x[,1]) + sin(2*pi*x[,2]) + sin(3*pi*x[,3]) + sin(4*pi*x[,4])
}
data <- as.data.frame(replicate(4, runif(n)))
colnames(data) <- c("x1v", "x2v", "x3v", "x4v")
fx <- fun(data)
y <- fx + rnorm(n)

# define marginal knots
knots <- list(x1v = quantile(data$x1v, probs = seq(0, 1, length.out = 10)),
              x2v = quantile(data$x2v, probs = seq(0, 1, length.out = 10)),
              x3v = quantile(data$x3v, probs = seq(0, 1, length.out = 10)),
              x4v = quantile(data$x4v, probs = seq(0, 1, length.out = 10)))
              
# define ssa knot indices
knots.indx <- c(bin_sample(data$x1v, nbin = 10, index.return = TRUE)$ix,
                bin_sample(data$x2v, nbin = 10, index.return = TRUE)$ix,
                bin_sample(data$x3v, nbin = 10, index.return = TRUE)$ix,
                bin_sample(data$x4v, nbin = 10, index.return = TRUE)$ix)

# fit sm with specified knots (tprk = TRUE)
sm.ssa <- sm(y ~ x1v + x2v + x3v + x4v, data = data, knots = knots.indx)

# fit sm with specified knots (tprk = FALSE)
sm.gam <- sm(y ~ x1v + x2v + x3v + x4v, data = data, knots = knots, tprk = FALSE)

# print both results (note: they are slightly different)
sm.ssa
sm.gam

# summarize both results (note: they are slightly different)
summary(sm.ssa)
summary(sm.gam)

# compare true MSE values (note: they are slightly different) 
mean( ( fx - sm.ssa$fit )^2 )
mean( ( fx - sm.gam$fit )^2 )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
