\name{nominal}
\alias{nominal}
\alias{basis.nom}
\alias{basis_nom}
\alias{penalty.nom}
\alias{penalty_nom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nominal Smoothing Spline Basis and Penalty
}
\description{
Generate the smoothing spline basis and penalty matrix for a nominal spline. This basis and penalty are for an unordered factor.
}
\usage{
basis.nom(x, knots, K = NULL, intercept = FALSE, ridge = FALSE)

penalty.nom(x, K = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Predictor variable (basis) or spline knots (penalty). Factor or integer vector of length \eqn{n}.
}
  \item{knots}{
Spline knots. Factor or integer vector of length \eqn{r}.
}
  \item{K}{
Number of levels of \code{x}. If \code{NULL}, this argument is defined as \code{K = length(unique(x))}.
}
  \item{intercept}{
If \code{TRUE}, the first column of the basis will be a column of ones. 
}
  \item{ridge}{
If \code{TRUE}, the basis matrix is post-multiplied by the inverse square root of the penalty matrix. See Note and Examples.
}
}
\details{
Generates a basis function or penalty matrix used to fit nominal smoothing splines.

With an intercept included, the basis function matrix has the form 
\deqn{X = [X_0,  X_1]}
where matrix \code{X_0} is an \eqn{n} by 1 matrix of ones, and \code{X_1} is a matrix of dimension \eqn{n} by \eqn{r}. 

The \code{X_0} matrix contains the "parametric part" of the basis (i.e., the intercept). The matrix \code{X_1} contains the "nonparametric part" of the basis, which consists of the \emph{reproducing kernel} function
\deqn{\rho(x, y) = \delta_{x y} - 1/K }
evaluated at all combinations of \code{x} and \code{knots}. The notation \eqn{\delta_{x y}} denotes Kronecker's delta function.

The penalty matrix consists of the \emph{reproducing kernel} function
\deqn{\rho(x, y) = \delta_{x y} - 1/K }
evaluated at all combinations of \code{x}. 
}
\value{
Basis: Matrix of dimension \code{c(length(x), df)} where \code{df = length(knots) + intercept}.

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
If the inputs \code{x} and \code{knots} are factors, they should have the same levels.

If the inputs \code{x} and \code{knots} are integers, the \code{knots} should be a subset of the input \code{x}.

If \code{ridge = TRUE}, the penalty matrix is the identity matrix.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{ordinal}} for a basis and penalty for ordered factors.
}
\examples{
######***######   standard parameterization   ######***######

# generate data
set.seed(0)
n <- 101
x <- factor(sort(rep(LETTERS[1:4], length.out = n)))
knots <- LETTERS[1:3]
eta <- 1:4
y <- eta[x] + rnorm(n, sd = 0.5)

# nominal smoothing spline basis
X <- basis.nom(x, knots, intercept = TRUE)

# nominal smoothing spline penalty
Q <- penalty.nom(knots, K = 4)

# pad Q with zeros (for intercept)
Q <- rbind(0, cbind(0, Q))

# define smoothing parameter
lambda <- 1e-5

# estimate coefficients
coefs <- solve(crossprod(X) + n * lambda * Q) \%*\% crossprod(X, y)

# estimate eta
yhat <- X \%*\% coefs

# check rmse
sqrt(mean((eta[x] - yhat)^2))



######***######   ridge parameterization   ######***######

# generate data
set.seed(0)
n <- 101
x <- factor(sort(rep(LETTERS[1:4], length.out = n)))
knots <- LETTERS[1:3]
eta <- 1:4
y <- eta[x] + rnorm(n, sd = 0.5)

# nominal smoothing spline basis
X <- basis.nom(x, knots, intercept = TRUE, ridge = TRUE)

# nominal smoothing spline penalty (ridge)
Q <- diag(rep(c(0, 1), times = c(1, ncol(X) - 1)))

# define smoothing parameter
lambda <- 1e-5

# estimate coefficients
coefs <- solve(crossprod(X) + n * lambda * Q) \%*\% crossprod(X, y)

# estimate eta
yhat <- X \%*\% coefs

# check rmse
sqrt(mean((eta[x] - yhat)^2))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
