\name{spherical}
\alias{spherical}
\alias{basis.sph}
\alias{basis_sph}
\alias{penalty.sph}
\alias{penalty_sph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spherical Spline Basis and Penalty
}
\description{
Generate the smoothing spline basis and penalty matrix for a spherical spline. This basis is designed for predictors where the values are points on a sphere.
}
\usage{
basis.sph(x, knots, m = 2, intercept = FALSE, ridge = FALSE)

penalty.sph(x, m = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Predictor variables (basis) or spline knots (penalty). Matrix of dimension \eqn{n} by \eqn{2}. Column 1 is latitude (-90 to 90 deg) and column 2 is longitude (-180 to 180 deg).
}
  \item{knots}{
Spline knots. Matrix of dimension \eqn{r} by \eqn{2}. Column 1 is latitude (-90 to 90 deg) and column 2 is longitude (-180 to 180 deg).
}
  \item{m}{
Penalty order. "m=2" for 2nd order spherical spline, "m=3" for 3rd order, and "m=4" for 4th order.
}
  \item{intercept}{
If \code{TRUE}, the first column of the basis will be a column of ones. 
}
  \item{ridge}{
If \code{TRUE}, the basis matrix is post-multiplied by the inverse square root of the penalty matrix. See Note and Examples.
}
}
\details{
Generates a basis function or penalty matrix used to fit spherical splines of order 2, 3, or 4.

With an intercept included, the basis function matrix has the form 
\deqn{X = [X_0,  X_1]}
where matrix \code{X_0} is an \eqn{n} by 1 matrix of ones, and \code{X_1} is a matrix of dimension \eqn{n} by \eqn{r}. 

The \code{X_0} matrix contains the "parametric part" of the basis (i.e., the intercept). 

The matrix \code{X_1} contains the "nonparametric part" of the basis, which consists of the \emph{reproducing kernel} function
\deqn{\rho(x, y) = [q_{2m-2}(x.y) - \alpha] / \beta  }
evaluated at all combinations of \code{x} and \code{knots}. Note that \eqn{\alpha = 1/(2m - 1)} and \eqn{\beta = 2\pi(2m-2)!} are constants, \eqn{q_{2m-2}(.)} is the spherical spline semi-kernel function, and \eqn{x.y} denotes the cosine of the angle between \eqn{x} and \eqn{y} (see References).

The penalty matrix consists of the \emph{reproducing kernel} function
\deqn{\rho(x, y) = [q_{2m-2}(x.y) - \alpha] / \beta  }
evaluated at all combinations of \code{x}.
}
\value{
Basis: Matrix of dimension \code{c(length(x), df)} where \code{df = nrow(knots) + intercept}.

Penalty: Matrix of dimension \code{c(r, r)} where \code{r = nrow(x)} is the number of knots.
}
\references{
Gu, C. (2013). Smoothing Spline ANOVA Models. 2nd Ed. New York, NY: Springer-Verlag. \doi{10.1007/978-1-4614-5369-7}

Wahba, G (1981). Spline interpolation and smoothing on the sphere. \emph{SIAM Journal on Scientific Computing, 2}(1), 5-16. \doi{10.1137/0902002}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The inputs \code{x} and \code{knots} must have the same dimension.

If \code{ridge = TRUE}, the penalty matrix is the identity matrix.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{thinplate}} for a thin plate spline basis and penalty.
}
\examples{
######***######   standard parameterization   ######***######

# function with spherical predictors
set.seed(0)
n <- 1000
myfun <- function(x){
  sin(pi*x[,1]) + cos(2*pi*x[,2]) + cos(pi*x[,3])
  }
x3d <- cbind(runif(n), runif(n), runif(n)) - 0.5
x3d <- t(apply(x3d, 1, function(x) x / sqrt(sum(x^2))))
eta <- myfun(x3d)
y <- eta + rnorm(n, sd = 0.5)

# convert x latitude and longitude
x <- cbind(latitude = acos(x3d[,3]) - pi/2,
           longitude = atan2(x3d[,2], x3d[,1])) * (180 / pi)

# select first 100 points as knots
knots <- x[1:100,]

# cubic spherical spline basis
X <- basis.sph(x, knots, intercept = TRUE)

# cubic spherical spline penalty
Q <- penalty.sph(knots)

# pad Q with zeros (for intercept)
Q <- rbind(0, cbind(0, Q))

# define smoothing parameter
lambda <- 1e-5

# estimate coefficients
coefs <- psolve(crossprod(X) + n * lambda * Q) \%*\% crossprod(X, y)

# estimate eta
yhat <- X \%*\% coefs

# check rmse
sqrt(mean((eta - yhat)^2))



######***######   ridge parameterization   ######***######

# function with spherical predictors
set.seed(0)
n <- 1000
myfun <- function(x){
  sin(pi*x[,1]) + cos(2*pi*x[,2]) + cos(pi*x[,3])
  }
x3d <- cbind(runif(n), runif(n), runif(n)) - 0.5
x3d <- t(apply(x3d, 1, function(x) x / sqrt(sum(x^2))))
eta <- myfun(x3d)
y <- eta + rnorm(n, sd = 0.5)

# convert x latitude and longitude
x <- cbind(latitude = acos(x3d[,3]) - pi/2,
           longitude = atan2(x3d[,2], x3d[,1])) * (180 / pi)

# select first 100 points as knots
knots <- x[1:100,]

# cubic spherical spline basis
X <- basis.sph(x, knots, intercept = TRUE, ridge = TRUE)

# cubic spherical spline penalty (ridge)
Q <- diag(rep(c(0, 1), times = c(1, ncol(X) - 1)))

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
