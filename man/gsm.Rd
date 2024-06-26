\name{gsm}
\alias{gsm}
\alias{family.gsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a Generalized Smooth Model
}
\description{
Fits a generalized semi- or nonparametric regression model with the smoothing parameter selected via one of seven methods: GCV, OCV, GACV, ACV, PQL, AIC, or BIC.
}
\usage{
gsm(formula, family = gaussian, data, weights, types = NULL, tprk = TRUE, 
    knots = NULL, skip.iter = TRUE, spar = NULL, lambda = NULL, control = list(),
    method = c("GCV", "OCV", "GACV", "ACV", "PQL", "AIC", "BIC"),
    xrange = NULL, thetas = NULL, mf = NULL)
    
\method{family}{gsm}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
Arguments for \code{gsm}:
  \item{formula}{
Object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. Uses the same syntax as \code{\link{lm}} and \code{\link{glm}}.
}
  \item{family}{
Description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function, or the result of a call to a family function. See the "Family Objects" section for details.
}
  \item{data}{
Optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{sm} is called.  
}
  \item{weights}{
Optional vector of weights to be used in the fitting process. If provided, weighted (penalized) likelihood estimation is used. Defaults to all 1.
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
  \item{skip.iter}{
Set to \code{FALSE} for deep tuning of the hyperparameters. Only applicable when multiple smooth terms are included. See the "Parameter Tuning" section for details.
}
\item{spar}{
Smoothing parameter. Typically (but not always) in the range \eqn{(0,1]}. If specified \code{lambda = 256^(3*(spar-1))}.
}
  \item{lambda}{
Computational smoothing parameter. This value is weighted by \eqn{n} to form the penalty coefficient (see Details). Ignored if \code{spar} is provided.
}
  \item{control}{
Optional list with named components that control the optimization specs for the smoothing parameter selection routine.

\bold{Note} that spar is only searched for in the interval \eqn{[lower, upper]}.
  \describe{
    \item{lower:}{lower bound for spar; defaults to 0.}
    \item{upper:}{upper bound for spar; defaults to 1.}
    \item{tol:}{the absolute precision (\bold{tol}erance) used by \code{\link{optimize}}; defaults to 1e-8.}
    \item{iterlim:}{the iteration limit used by \code{\link{nlm}}; defaults to 5000.}
    \item{print.level:}{the print level used by \code{\link{nlm}}; defaults to 0 (no printing).}
    \item{epsilon:}{relative convergence tolerance for IRPLS algorithm; defaults to 1e-8}
    \item{maxit:}{maximum number of iterations for IRPLS algorithm; defaults to 25}
    \item{epsilon.out:}{relative convergence tolerance for iterative NegBin update; defaults to 1e-6}
    \item{maxit.out:}{maximum number of iterations for iterative NegBin update; defaults to 10}
  }
}
  \item{method}{
Method for selecting the smoothing parameter. Ignored if \code{lambda} is provided.
}
\item{xrange}{
Optional named list containing the range of each predictor. If \code{NULL}, the ranges are calculated from the input \code{data}.
}
  \item{thetas}{
Optional vector of hyperparameters to use for smoothing. If \code{NULL}, these are tuned using the requested \code{method}.
}
\item{mf}{
Optional model frame constructed from \code{formula} and \code{data} (and potentially \code{weights}).
}
Note: the last two arguments are not intended to be called by the typical user of this function. These arguments are included primarily for internal usage by the \code{\link{boot.gsm}} function.

Arguments for \code{family.gsm}:
\item{object}{an object of class "gsm"}
\item{...}{additional arguments (currently ignored)}
}
\details{
Letting \eqn{\eta_i = \eta(x_i)} with \eqn{x_i = (x_{i1}, \ldots, x_{ip})}, the function is represented as \deqn{\eta = X \beta + Z \alpha} where the basis functions in \eqn{X} span the null space (i.e., parametric effects), and \eqn{Z} contains the kernel function(s) of the contrast space (i.e., nonparametric effects) evaluated at all combinations of observed data points and knots. The vectors \eqn{\beta} and \eqn{\alpha} contain unknown basis function coefficients.

Let \eqn{\mu_i = E(y_i)} denote the mean of the \eqn{i}-th response. The unknown function is related to the mean \eqn{\mu_i} such as \deqn{g(\mu_i) = \eta_i} where \eqn{g()} is a known link function. Note that this implies that \eqn{\mu_i = g^{-1}(\eta_i)} given that the link function is assumed to be invertible. 

The penalized likelihood estimation problem has the form
\deqn{
- \sum_{i = 1}^n [y_i \xi_i - b(\xi_i)] + n \lambda \alpha' Q \alpha
}
where \eqn{\xi_i} is the canonical parameter, \eqn{b()} is a known function that depends on the chosen family, and \eqn{Q} is the penalty matrix. Note that \eqn{\xi_i = g_0(\mu_i)} where \eqn{g_0} is the canonical link function. This implies that \eqn{\xi_i = \eta_i} when the chosen link function is canonical, i.e., when \eqn{g = g_0}.
}
\value{
An object of class "gsm" with components:
  %\item{fitted.values}{the fitted values, obtained by transforming the linear predictors by the inverse of the link function.}
  \item{linear.predictors}{the linear fit on link scale. Use \code{\link{fitted.gsm}} to obtain the fitted values on the response scale.}
  \item{se.lp}{the standard errors of the linear predictors.}
  \item{deviance}{up to a constant, minus twice the maximized log-likelihood. Where sensible, the constant is chosen so that a saturated model has deviance zero.}
  \item{cv.crit}{the cross-validation criterion.}
  \item{nsdf}{the degrees of freedom (Df) for the null space.}
  \item{df}{the estimated degrees of freedom (Df) for the fit model.}
  \item{df.residual}{the residual degrees of freedom = \code{nobs - df}}
  \item{r.squared}{the squared correlation between response and fitted values.}
  \item{dispersion}{the estimated dispersion parameter.}
  \item{logLik}{the log-likelihood.}
  \item{aic}{Akaike's Information Criterion.}
  \item{bic}{Bayesian Information Criterion.}
  \item{spar}{the value of \code{spar} computed or given, i.e., \eqn{s = 1 + \log_{256}(\lambda)/3}}
  \item{lambda}{the value of \eqn{\lambda} corresponding to \code{spar}, i.e., \eqn{\lambda = 256^{3*(s-1)}}.}
  \item{penalty}{the smoothness penalty \eqn{\alpha' Q \alpha}.}
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
  \item{weights}{the weights used for fitting (if applicable)}
  \item{call}{the matched call.}
  \item{family}{the input family evaluated as a function using .}
  \item{iter}{the number of iterations of IRPLS used.}
  \item{residuals}{the working (IRPLS) residuals from the fitted model.}
  \item{null.deviance}{the deviance of the null model (i.e., intercept only).}
}
\references{
Berry, L. N., & Helwig, N. E. (2021). Cross-validation, information theory, or maximum likelihood? A comparison of tuning methods for penalized splines. \emph{Stats, 4}(3), 701-724. \doi{10.3390/stats4030042}

Craven, P. and Wahba, G. (1979). Smoothing noisy data with spline functions: Estimating the correct degree of smoothing by the method of generalized cross-validation. \emph{Numerische Mathematik, 31}, 377-403. \doi{10.1007/BF01404567}

Gu, C. (2013). Smoothing spline ANOVA models, 2nd edition. New York: Springer. \doi{10.1007/978-1-4614-5369-7}

Gu, C. and Wahba, G. (1991). Minimizing GCV/GML scores with multiple smoothing parameters via the Newton method. \emph{SIAM Journal on Scientific and Statistical Computing, 12(2)}, 383-398. \doi{10.1137/0912021}

Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), SAGE Research Methods Foundations. \doi{10.4135/9781526421036885885}

Helwig, N. E. (2021). Spectrally sparse nonparametric regression via elastic net regularized smoothers. \emph{Journal of Computational and Graphical Statistics, 30}(1), 182-191. \doi{10.1080/10618600.2020.1806855}

Helwig, N. E. (2024). Precise tensor product smoothing via spectral splines. \emph{Stats, }
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

\section{Family Objects }{
Supported families and links include:
\tabular{lll}{
\code{family} \tab \tab \code{link} \cr
binomial \tab \tab logit, probit, cauchit, log, cloglog \cr
gaussian \tab \tab identity, log, inverse \cr
Gamma \tab \tab inverse, identity, log \cr
inverse.gaussian \tab \tab 1/mu^2, inverse, identity, log \cr
poisson \tab \tab log, identity, sqrt \cr
NegBin \tab \tab log, identity, sqrt \cr
}
See \code{\link{NegBin}} for information about the Negative Binomial family.
}

\section{Methods }{
The smoothing parameter can be selected using one of seven methods: \cr
Generalized Cross-Validation (GCV) \cr
Ordinary Cross-Validation (OCV) \cr
Generalized Approximate Cross-Validation (GACV) \cr
Approximate Cross-Validation (ACV) \cr
Penalized Quasi-Likelihood (PQL) \cr
Akaike's Information Criterion (AIC) \cr
Bayesian Information Criterion (BIC)

}

\section{Types of Smooths }{

The following codes specify the spline types:
\tabular{ll}{
par \tab Parametric effect (factor, integer, or numeric). \cr
ran \tab Random effect/intercept (unordered factor). \cr
nom \tab Nominal smoothing spline (unordered factor). \cr
ord \tab Ordinal smoothing spline (ordered factor). \cr
lin \tab Linear smoothing spline (integer or numeric). \cr
cub \tab Cubic smoothing spline (integer or numeric). \cr
qui \tab Quintic smoothing spline (integer or numeric). \cr
per \tab Periodic smoothing spline (integer or numeric).\cr  
sph \tab Spherical spline (matrix with \eqn{d = 2} columns: lat, long). \cr
tps \tab Thin plate spline (matrix with \eqn{d \ge 1} columns).
}

For finer control of some specialized spline types:
\tabular{ll}{
per.lin \tab Linear periodic spline (\eqn{m = 1}). \cr
per.cub \tab Cubic periodic spline (\eqn{m = 2}). \cr
per.qui \tab Quintic periodic spline (\eqn{m = 3}). \cr
sph.2 \tab 2nd order spherical spline (\eqn{m = 2}). \cr
sph.3 \tab 3rd order spherical spline (\eqn{m = 3}). \cr
sph.4 \tab 4th order spherical spline (\eqn{m = 4}). \cr
tps.lin \tab Linear thin plate spline (\eqn{m = 1}). \cr
tps.cub \tab Cubic thin plate spline (\eqn{m = 2}). \cr
tps.qui \tab Quintic thin plate spline (\eqn{m = 3}). \cr
}

For details on the spline kernel functions, see \code{\link{basis.nom}} (nominal), \code{\link{basis.ord}} (ordinal), \code{\link{basis.poly}} (polynomial), \code{\link{basis.sph}} (spherical), and \code{\link{basis.tps}} (thin plate).

Note: "ran" is default for unordered factors when the number of levels is 20 or more, whereas "nom" is the default for unordered factors otherwise. 

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

The \eqn{k}-th predictor's marginal effect can be denoted as \deqn{f_k = X_k \beta_k + Z_k \alpha_k} where \eqn{X_k} is the \eqn{n} by \eqn{m_k} null space basis function matrix, and \eqn{Z_k} is the \eqn{n} by \eqn{r_k} contrast space basis function matrix. 

If \code{tprk = TRUE}, the null space basis function matrix has the form \eqn{X = [1, X_1, X_2]} and the contrast space basis function matrix has the form \deqn{Z = \theta_1 Z_1 + \theta_2 Z_2} where the \eqn{\theta_k} are the "extra" smoothing parameters. Note that \eqn{Z} is of dimension \eqn{n} by \eqn{r = r_1 = r_2}.

If \code{tprk = FALSE}, the null space basis function matrix has the form \eqn{X = [1, X_1, X_2]}, and the contrast space basis function matrix has the form \deqn{Z = [\theta_1 Z_1, \theta_2 Z_2]} where the \eqn{\theta_k} are the "extra" smoothing parameters. Note that \eqn{Z} is of dimension \eqn{n} by \eqn{r = r_1 + r_2}.
}

\section{Parameter Tuning }{
When multiple smooth terms are included in the model, there are smoothing (hyper)parameters that weight the contribution of each combination of smooth terms. These hyperparameters are distinct from the overall smoothing parameter \code{lambda} that weights the contribution of the penalty. 

\code{skip.iter = TRUE} (default) estimates the smoothing hyperparameters using Algorithm 3.2 of Gu and Wahba (1991), which typically provides adequate results when the model form is correctly specified. The \code{lambda} parameter is tuned via the specified smoothing parameter selection \code{method}.

\code{skip.iter = FALSE} uses Algorithm 3.2 as an initialization, and then the \code{\link{nlm}} function is used to tune the hyperparameters via the specified smoothing parameter selection \code{method}. Setting \code{skip.iter = FALSE} can (substantially) increase the model fitting time, but should produce better results---particularly if the model \code{formula} is misspecified.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\bold{Related Modeling Functions}: 

\code{\link{ss}} for fitting a smoothing spline with a single predictor (Gaussian response).

\code{\link{sm}} for fitting smooth models with multiple predictors of mixed types (Gaussian response). \cr

\bold{S3 Methods and Related Functions for "gsm" Objects}:

\code{\link{boot.gsm}} for bootstrapping \code{gsm} objects.

\code{\link{coef.gsm}} for extracting coefficients from \code{gsm} objects.

\code{\link{cooks.distance.gsm}} for calculating Cook's distances from \code{gsm} objects.

\code{\link{cov.ratio}} for computing covariance ratio from \code{gsm} objects.

\code{\link{deviance.gsm}} for extracting deviance from \code{gsm} objects.

\code{\link{dfbeta.gsm}} for calculating DFBETA from \code{gsm} objects.

\code{\link{dfbetas.gsm}} for calculating DFBETAS from \code{gsm} objects.

\code{\link{diagnostic.plots}} for plotting regression diagnostics from \code{gsm} objects.

\code{\link{family.gsm}} for extracting \code{family} from \code{gsm} objects.

\code{\link{fitted.gsm}} for extracting fitted values from \code{gsm} objects.

\code{\link{hatvalues.gsm}} for extracting leverages from \code{gsm} objects.

\code{\link{model.matrix.gsm}} for constructing model matrix from \code{gsm} objects.

\code{\link{plot.gsm}} for plotting effects from \code{gsm} objects.

\code{\link{predict.gsm}} for predicting from \code{gsm} objects.

\code{\link{residuals.gsm}} for extracting residuals from \code{gsm} objects.

\code{\link{rstandard.gsm}} for computing standardized residuals from \code{gsm} objects.

\code{\link{rstudent.gsm}} for computing studentized residuals from \code{gsm} objects.

\code{\link{smooth.influence}} for calculating basic influence information from \code{gsm} objects.

\code{\link{smooth.influence.measures}} for convenient display of influential observations from \code{gsm} objects.

\code{\link{summary.gsm}} for summarizing \code{gsm} objects.

\code{\link{vcov.gsm}} for extracting coefficient covariance matrix from \code{gsm} objects.

\code{\link{weights.gsm}} for extracting prior weights from \code{gsm} objects.
}
\examples{
##########   EXAMPLE 1   ##########
### 1 continuous predictor

# generate data
n <- 1000
x <- seq(0, 1, length.out = n)
fx <- 3 * x + sin(2 * pi * x) - 1.5

# gaussian (default)
set.seed(1)
y <- fx + rnorm(n, sd = 1/sqrt(2))
mod <- gsm(y ~ x, knots = 10)
plot(mod)
mean((mod$linear.predictors - fx)^2)

# compare to result from sm (they are identical)
mod.sm <- sm(y ~ x, knots = 10)
plot(mod.sm)
mean((mod$linear.predictors - mod.sm$fitted.values)^2)

# binomial (no weights)
set.seed(1)
y <- rbinom(n = n, size = 1, p = 1 / (1 + exp(-fx)))
mod <- gsm(y ~ x, family = binomial, knots = 10)
plot(mod)
mean((mod$linear.predictors - fx)^2)

# binomial (w/ weights)
set.seed(1)
w <- as.integer(rep(c(10,20,30,40,50), length.out = n))
y <- rbinom(n = n, size = w, p = 1 / (1 + exp(-fx))) / w
mod <- gsm(y ~ x, family = binomial, weights = w, knots = 10)
plot(mod)
mean((mod$linear.predictors - fx)^2)

# poisson
set.seed(1)
y <- rpois(n = n, lambda = exp(fx))
mod <- gsm(y ~ x, family = poisson, knots = 10)
plot(mod)
mean((mod$linear.predictors - fx)^2)

# negative binomial (known theta)
set.seed(1)
y <- rnbinom(n = n, size = 1/2, mu = exp(fx))
mod <- gsm(y ~ x, family = NegBin(theta = 1/2), knots = 10)
plot(mod)
mean((mod$linear.predictors - fx)^2)
mod$family$theta    # fixed theta

# negative binomial (unknown theta)
set.seed(1)
y <- rnbinom(n = n, size = 1/2, mu = exp(fx))
mod <- gsm(y ~ x, family = NegBin, knots = 10)
plot(mod)
mean((mod$linear.predictors - fx)^2)
mod$family$theta    # estimated theta

# gamma
set.seed(1)
y <- rgamma(n = n, shape = 2, scale = (1 / (2 + fx)) / 2)
mod <- gsm(y ~ x, family = Gamma, knots = 10)
plot(mod)
mean((mod$linear.predictors - fx - 2)^2)

# inverse.gaussian (not run; requires statmod)
##set.seed(1)
##y <- statmod::rinvgauss(n = n, mean = sqrt(1 / (2 + fx)), shape = 2)
##mod <- gsm(y ~ x, family = inverse.gaussian, knots = 10)
##plot(mod)
##mean((mod$linear.predictors - fx - 2)^2)



##########   EXAMPLE 2   ##########
### 1 continuous and 1 nominal predictor
### additive model

# generate data
n <- 1000
x <- seq(0, 1, length.out = n)
z <- factor(sample(letters[1:3], size = n, replace = TRUE))
fun <- function(x, z){
  mu <- c(-2, 0, 2)
  zi <- as.integer(z)
  fx <- mu[zi] + 3 * x + sin(2 * pi * x) - 1.5
}
fx <- fun(x, z)

# define marginal knots
probs <- seq(0, 0.9, by = 0.1)
knots <- list(x = quantile(x, probs = probs),
              z = letters[1:3])
              
# gaussian (default)
set.seed(1)
y <- fx + rnorm(n, sd = 1/sqrt(2))
mod <- gsm(y ~ x + z, knots = knots)
plot(mod)
mean((mod$linear.predictors - fx)^2)

# compare to result from sm (they are identical)
mod.sm <- sm(y ~ x + z, knots = knots)
plot(mod.sm)
mean((mod$linear.predictors - mod.sm$fitted.values)^2)

# binomial (no weights)
set.seed(1)
y <- rbinom(n = n, size = 1, p = 1 / (1 + exp(-fx)))
mod <- gsm(y ~ x + z, family = binomial, knots = knots)
plot(mod)
mean((mod$linear.predictors - fx)^2)

# binomial (w/ weights)
set.seed(1)
w <- as.integer(rep(c(10,20,30,40,50), length.out = n))
y <- rbinom(n = n, size = w, p = 1 / (1 + exp(-fx))) / w
mod <- gsm(y ~ x + z, family = binomial, weights = w, knots = knots)
plot(mod)
mean((mod$linear.predictors - fx)^2)

# poisson
set.seed(1)
y <- rpois(n = n, lambda = exp(fx))
mod <- gsm(y ~ x + z, family = poisson, knots = knots)
plot(mod)
mean((mod$linear.predictors - fx)^2)

# negative binomial (known theta)
set.seed(1)
y <- rnbinom(n = n, size = 1/2, mu = exp(fx))
mod <- gsm(y ~ x + z, family = NegBin(theta = 1/2), knots = knots)
plot(mod)
mean((mod$linear.predictors - fx)^2)
mod$family$theta    # fixed theta

# negative binomial (unknown theta)
set.seed(1)
y <- rnbinom(n = n, size = 1/2, mu = exp(fx))
mod <- gsm(y ~ x + z, family = NegBin, knots = knots)
plot(mod)
mean((mod$linear.predictors - fx)^2)
mod$family$theta    # estimated theta

# gamma
set.seed(1)
y <- rgamma(n = n, shape = 2, scale = (1 / (4 + fx)) / 2)
mod <- gsm(y ~ x + z, family = Gamma, knots = knots)
plot(mod)
mean((mod$linear.predictors - fx - 4)^2)

# inverse.gaussian (not run; requires statmod)
##set.seed(1)
##y <- statmod::rinvgauss(n = n, mean = sqrt(1 / (4 + fx)), shape = 2)
##mod <- gsm(y ~ x + z, family = inverse.gaussian, knots = knots)
##plot(mod)
##mean((mod$linear.predictors - fx - 4)^2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
