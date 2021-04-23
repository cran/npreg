\name{polynomial}
\alias{polynomial}
\alias{basis.poly}
\alias{basis_poly}
\alias{penalty.poly}
\alias{penalty_poly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Polynomial Smoothing Spline Basis and Penalty
}
\description{
Generate the smoothing spline basis and penalty matrix for a polynomial spline. Derivatives of the smoothing spline basis matrix are supported.
}
\usage{
basis.poly(x, knots, m = 2, d = 0, xmin = min(x), xmax = max(x), 
           periodic = FALSE, rescale = FALSE, intercept = FALSE, 
           bernoulli = TRUE, ridge = FALSE)

penalty.poly(x, m = 2, xmin = min(x), xmax = max(x), 
             periodic = FALSE, rescale = FALSE, bernoulli = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Predictor variable (basis) or spline knots (penalty). Numeric or integer vector of length \eqn{n}.
}
  \item{knots}{
Spline knots. Numeric or integer vector of length \eqn{r}.
}
  \item{m}{
Penalty order. "m=1" for linear smoothing spline, "m=2" for cubic, and "m=3" for quintic.  
}
  \item{d}{
Derivative order. "d=0" for smoothing spline basis, "d=1" for 1st derivative of basis, and "d=2" for 2nd derivative of basis.  
}
  \item{xmin}{
Minimum value of "x".  
}
  \item{xmax}{
Maximum value of "x".  
}
  \item{periodic}{
If \code{TRUE}, the smoothing spline basis is periodic w.r.t. the interval [\code{xmin}, \code{xmax}].
}
  \item{rescale}{
If \code{TRUE}, the nonparametric part of the basis is divided by the average of the reproducing kernel function evaluated at the \code{knots}.
}
  \item{intercept}{
If \code{TRUE}, the first column of the basis will be a column of ones. 
}
  \item{bernoulli}{
If \code{TRUE}, scaled Bernoulli polynomials are used for the basis and penalty functions.
}
  \item{ridge}{
If \code{TRUE}, the basis matrix is post-multiplied by the inverse square root of the penalty matrix. See Note and Examples.
}
}
\details{
Generates a basis function or penalty matrix used to fit linear, cubic, and quintic smoothing splines (or evaluate their derivatives). 

For non-periodic smoothing splines, the basis function matrix has the form 
\deqn{X = [X_0,  X_1]}
where the matrix \code{X_0} is of dimension \eqn{n} by \eqn{m-1} (plus 1 if an intercept is included), and \code{X_1} is a matrix of dimension \eqn{n} by \eqn{r}. 

The \code{X_0} matrix contains the "parametric part" of the basis, which includes polynomial functions of \code{x} up to degree \eqn{m-1}.

The matrix \code{X_1} contains the "nonparametric part" of the basis, which consists of the \emph{reproducing kernel} function
\deqn{\rho(x, y) = \kappa_m(x) \kappa_m(y) + (-1)^{m-1} \kappa_{2m}(|x-y|)}
evaluated at all combinations of \code{x} and \code{knots}. The \eqn{\kappa_v} functions are scaled Bernoulli polynomials.

For periodic smoothing splines, the \eqn{X_0} matrix only contains the intercept column and the modified reproducing kernel function 
\deqn{\rho(x, y) = (-1)^{m-1} \kappa_{2m}(|x-y|)}
is evaluated for all combinations of \code{x} and \code{knots}.

For non-periodic smoothing splines, the penalty matrix consists of the \emph{reproducing kernel} function
\deqn{\rho(x, y) = \kappa_m(x) \kappa_m(y) + (-1)^{m-1} \kappa_{2m}(|x-y|)}
evaluated at all combinations of \code{x}. For periodic smoothing splines, the modified reproducing kernel function
\deqn{\rho(x, y) = (-1)^{m-1} \kappa_{2m}(|x-y|)}
is evaluated for all combinations of \code{x}.

If \code{bernoulli = FALSE}, the reproducing kernel function is defined as
\deqn{\rho(x, y) = (1/(m-1)!)^2 \int_0^1 (x - u)_+^{m-1} (y - u)_+^{m-1} du }
where \eqn{(.)_+ = \max(., 0)}. This produces the "classic" definition of a smoothing spline, where the function estimate is a piecewise polynomial function with pieces of degree \eqn{2m - 1}.
}
\value{
Basis: Matrix of dimension \code{c(length(x), df)} where \code{df >= length(knots)}. If the smoothing spline basis is not periodic (default), then the number of columns is \code{df = length(knots) + m - !intercept}. For periodic smoothing splines, the basis has \code{m} fewer columns.

Penalty: Matrix of dimension \code{c(r, r)} where \code{r = length(x)} is the number of knots.
}
\references{
Gu, C. (2013). Smoothing Spline ANOVA Models. 2nd Ed. New York, NY: Springer-Verlag.  \doi{10.1007/978-1-4614-5369-7}

Helwig, N. E. (2017). Regression with ordered predictors via ordinal smoothing splines. \emph{Frontiers in Applied Mathematics and Statistics, 3}(15), 1-13. \doi{10.3389/fams.2017.00015}

Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), \emph{SAGE Research Methods Foundations.} \doi{10.4135/9781526421036885885}

Helwig, N. E., & Ma, P. (2015). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. \emph{Journal of Computational and Graphical Statistics, 24}(3), 715-732. \doi{10.1080/10618600.2014.926819}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Inputs \code{x} and \code{knots} should be within the interval [\code{xmin}, \code{xmax}].

If \code{ridge = TRUE}, the penalty matrix is the identity matrix.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{thinplate}} for a thin-plate spline basis and penalty.

See \code{\link{ordinal}} for a basis and penalty for ordered factors.
}
\examples{
######***######   standard parameterization   ######***######

# generate data
set.seed(0)
n <- 101
x <- seq(0, 1, length.out = n)
knots <- seq(0, 0.95, by = 0.05)
eta <- 1 + 2 * x + sin(2 * pi * x)
y <- eta + rnorm(n, sd = 0.5)

# cubic smoothing spline basis
X <- basis.poly(x, knots, intercept = TRUE)

# cubic smoothing spline penalty
Q <- penalty.poly(knots, xmin = min(x), xmax = max(x))

# pad Q with zeros (for intercept and linear effect)
Q <- rbind(0, 0, cbind(0, 0, Q))

# define smoothing parameter
lambda <- 1e-5

# estimate coefficients
coefs <- solve(crossprod(X) + n * lambda * Q) \%*\% crossprod(X, y)

# estimate eta
yhat <- X \%*\% coefs

# check rmse
sqrt(mean((eta - yhat)^2))

# plot results
plot(x, y)
lines(x, yhat)



######***######   ridge parameterization   ######***######

# generate data
set.seed(0)
n <- 101
x <- seq(0, 1, length.out = n)
knots <- seq(0, 0.95, by = 0.05)
eta <- 1 + 2 * x + sin(2 * pi * x)
y <- eta + rnorm(n, sd = 0.5)

# cubic smoothing spline basis
X <- basis.poly(x, knots, intercept = TRUE, ridge = TRUE)

# cubic smoothing spline penalty (ridge)
Q <- diag(rep(c(0, 1), times = c(2, ncol(X) - 2)))

# define smoothing parameter
lambda <- 1e-5

# estimate coefficients
coefs <- solve(crossprod(X) + n * lambda * Q) \%*\% crossprod(X, y)

# estimate eta
yhat <- X \%*\% coefs

# check rmse
sqrt(mean((eta - yhat)^2))

# plot results
plot(x, y)
lines(x, yhat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
