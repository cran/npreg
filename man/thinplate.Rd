\name{thinplate}
\alias{thinplate}
\alias{basis_tps}
\alias{penalty_tps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Thin-Plate Spline Basis and Penalty
}
\description{
Generate the smoothing spline basis and penalty matrix for a thin-plate spline.
}
\usage{
basis_tps(x, knots, m = 2, rk = TRUE, intercept = FALSE)

penalty_tps(x, m = 2, rk = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Predictor variables (basis) or spline knots (penalty). Numeric or integer vector of length \eqn{n}, or matrix of dimension \eqn{n} by \eqn{p}.
}
  \item{knots}{
Spline knots. Numeric or integer vector of length \eqn{r}, or matrix of dimension \eqn{r} by \eqn{p}.
}
  \item{m}{
Penalty order. "m=1" for linear thin-plate spline, "m=2" for cubic, and "m=3" for quintic. Must satisfy \eqn{2m > p}.
}
  \item{rk}{
If true (default), the reproducing kernel parameterization is used. Otherwise, the standard thin-plate basis is returned.  
}
  \item{intercept}{
If \code{TRUE}, the first column of the basis will be a column of ones. 
}
}
\details{
Generates a basis function or penalty matrix used to fit linear, cubic, and quintic thin-plate splines. 

The basis function matrix has the form 
\deqn{X = [X_0,  X_1]}
where the matrix \code{X_0} is of dimension \eqn{n} by \eqn{M-1} (plus 1 if an intercept is included) where \eqn{M = {p+m-1 \choose p}}, and \code{X_1} is a matrix of dimension \eqn{n} by \eqn{r}. 

The \code{X_0} matrix contains the "parametric part" of the basis, which includes polynomial functions of the columns of \code{x} up to degree \eqn{m-1} (and potentially interactions).

The matrix \code{X_1} contains the "nonparametric part" of the basis.

If \code{rk = TRUE}, the matrix \code{X_1} consists of the \emph{reproducing kernel} function
\deqn{\rho(x, y) = (I - P_x) (I - P_y) E(|x - y|)}
evaluated at all combinations of \code{x} and \code{knots}. Note that \eqn{P_x} and \eqn{P_y} are projection operators, \eqn{|.|} denotes the Euclidean distance, and the TPS semi-kernel is defined as
\deqn{E(z) = \alpha z^{2m-p} \log(z)}
if \eqn{p} is even and
\deqn{E(z) = \beta z^{2m-p}}
otherwise, where \eqn{\alpha} and \eqn{\beta} are positive constants (see References).

If \code{rk = FALSE}, the matrix \code{X_1} contains the TPS semi-kernel \eqn{E(.)} evaluated at all combinations of \code{x} and \code{knots}. Note: the TPS semi-kernel is \emph{not} positive (semi-)definite, but the projection is.

If \code{rk = TRUE}, the penalty matrix consists of the \emph{reproducing kernel} function
\deqn{\rho(x, y) = (I - P_x) (I - P_y) E(|x - y|)}
evaluated at all combinations of \code{x}. If \code{rk = FALSE}, the penalty matrix contains the TPS semi-kernel \eqn{E(.)} evaluated at all combinations of \code{x}.
}
\value{
Basis: Matrix of dimension \code{c(length(x), df)} where \code{df = nrow(as.matrix(knots)) + choose(p + m - 1, p) - !intercept} and \code{p = ncol(as.matrix(x))}.

Penalty: Matrix of dimension \code{c(r, r)} where \code{r = nrow(as.matrix(x))} is the number of knots.
}
\references{
Gu, C. (2013). Smoothing Spline ANOVA Models. 2nd Ed. New York, NY: Springer-Verlag. doi: https://doi.org/10.1007/978-1-4614-5369-7

Helwig, N. E. (2017). Regression with ordered predictors via ordinal smoothing splines. Frontiers in Applied Mathematics and Statistics, 3(15), 1-13. doi: https://doi.org/10.3389/fams.2017.00015

Helwig, N. E., & Ma, P. (2015). Fast and stable multiple smoothing parameter selection in smoothing spline analysis of variance models with large samples. Journal of Computational and Graphical Statistics, 24(3), 715-732. doi: https://doi.org/10.1080/10618600.2014.926819
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The inputs \code{x} and \code{knots} must have the same dimension.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{polynomial}} for a basis and penalty for numeric variables.

See \code{\link{spherical}} for a basis and penalty for spherical variables.
}
\examples{
# generate data
set.seed(0)
n <- 101
x <- seq(0, 1, length.out = n)
knots <- seq(0, 0.95, by = 0.05)
eta <- 1 + 2 * x + sin(2 * pi * x)
y <- eta + rnorm(n, sd = 0.5)

# cubic thin-plate spline basis
X <- basis_tps(x, knots, intercept = TRUE)

# cubic thin-plate spline penalty
Q <- penalty_tps(knots)

# padd Q with zeros (for intercept and linear effect)
Q <- rbind(0, 0, cbind(0, 0, Q))

# define smoothing parameter
lambda <- 1e-5

# estimate coefficients
coefs <- psolve(crossprod(X) + n * lambda * Q) \%*\% crossprod(X, y)

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
