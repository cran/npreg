\name{spherical}
\alias{spherical}
\alias{basis_sph}
\alias{penalty_sph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spherical Spline Basis and Penalty
}
\description{
Generate the smoothing spline basis and penalty matrix for a spherical spline. This basis is designed for a 3D predictor where the values are points on a sphere.
}
\usage{
basis_sph(x, knots, m = 2, rescale = TRUE, intercept = FALSE)

penalty_sph(x, m = 2, rescale = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Predictor variables (basis) or spline knots (penalty). Matrix of dimension \eqn{n} by \eqn{3}.
}
  \item{knots}{
Spline knots. Matrix of dimension \eqn{r} by \eqn{3}.
}
  \item{m}{
Penalty order. "m=1" for linear spherical spline, "m=2" for cubic, and "m=3" for quintic.
}
  \item{rescale}{
If \code{TRUE}, the nonparametric part of the basis is divided by the average of the reproducing kernel function evaluated at the \code{knots}.
}
  \item{intercept}{
If \code{TRUE}, the first column of the basis will be a column of ones. 
}
}
\details{
Generates a basis function or penalty matrix used to fit linear, cubic, and spherical splines. 

With an intercept included, the basis function matrix has the form 
\deqn{X = [X_0,  X_1]}
where matrix \code{X_0} is an \eqn{n} by 1 matrix of ones, and \code{X_1} is a matrix of dimension \eqn{n} by \eqn{r}. 

The \code{X_0} matrix contains the "parametric part" of the basis (i.e., the intercept). 

The matrix \code{X_1} contains the "nonparametric part" of the basis, which consists of the \emph{reproducing kernel} function
\deqn{\rho(x, y) = 1 + [s_{2m}(x.y) - \alpha_m] / \beta_m  }
evaluated at all combinations of \code{x} and \code{knots}. Note that \eqn{\alpha_m = 1/(2m + 1)} and \eqn{\beta_m = 2\pi(2m)!} are constants, \eqn{s_{2m}(.)} is the spherical spline semi-kernel function, and \eqn{x.y} denote the inner product between \eqn{x} and \eqn{y} (see References).

The penalty matrix consists of the \emph{reproducing kernel} function
\deqn{\rho(x, y) = 1 + [s_{2m}(x.y) - \alpha_m] / \beta_m  }
evaluated at all combinations of \code{x}.
}
\value{
Basis: Matrix of dimension \code{c(length(x), df)} where \code{df = nrow(knots) + intercept}.

Penalty: Matrix of dimension \code{c(r, r)} where \code{r = nrow(x)} is the number of knots.
}
\references{
Gu, C. (2013). Smoothing Spline ANOVA Models. 2nd Ed. New York, NY: Springer-Verlag. \url{https://doi.org/10.1007/978-1-4614-5369-7}

Wahba, G (1981). Spline interpolation and smoothing on the sphere. \emph{SIAM Journal on Scientific Computing, 2}(1), 5-16. \url{https://doi.org/10.1137/0902002}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The inputs \code{x} and \code{knots} must have the same dimension.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{thinplate}} for a thin-plate spline basis and penalty.
}
\examples{
# function with three spherical predictors
set.seed(0)
n <- 1000
myfun <- function(x){
  sin(pi*x[,1]) + cos(2*pi*x[,2]) + cos(pi*x[,3])
  }
x <- cbind(runif(n), runif(n), runif(n)) - 0.5
x <- t(apply(x, 1, function(x) x / sqrt(sum(x^2))))
eta <- myfun(x)
y <- eta + rnorm(n, sd = 0.5)
knots <- x[1:100,]

# cubic spherical spline basis
X <- basis_sph(x, knots, intercept = TRUE)

# cubic spherical spline penalty
Q <- penalty_sph(knots)

# padd Q with zeros (for intercept)
Q <- rbind(0, cbind(0, Q))

# define smoothing parameter
lambda <- 1e-5

# estimate coefficients
coefs <- psolve(crossprod(X) + n * lambda * Q) \%*\% crossprod(X, y)

# estimate eta
yhat <- X \%*\% coefs

# check rmse
sqrt(mean((eta - yhat)^2))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
