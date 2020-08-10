\name{ordinal}
\alias{ordinal}
\alias{basis_ord}
\alias{penalty_ord}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ordinal Smoothing Spline Basis and Penalty
}
\description{
Generate the smoothing spline basis and penalty matrix for an ordinal spline. This basis and penalty are for an ordered factor.
}
\usage{
basis_ord(x, knots, K = NULL, intercept = FALSE)

penalty_ord(x, K = NULL, xlev = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Predictor variable (basis) or spline knots (penalty). Ordered factor or integer vector of length \eqn{n}.
}
  \item{knots}{
Spline knots. Ordered factor or integer vector of length \eqn{r}.
}
  \item{K}{
Number of levels of \code{x}. If \code{NULL}, this argument is defined as \code{K = length(unique(x))}.
}
  \item{xlev}{
Factor levels of \code{x} (for penalty). If \code{NULL}, the levels are defined as \code{levels(as.ordered(x))}.
}
  \item{intercept}{
If \code{TRUE}, the first column of the basis will be a column of ones. 
}
}
\details{
Generates a basis function or penalty matrix used to fit ordinal smoothing splines.

With an intercept included, the basis function matrix has the form 
\deqn{X = [X_0,  X_1]}
where matrix \code{X_0} is an \eqn{n} by 1 matrix of ones, and \code{X_1} is a matrix of dimension \eqn{n} by \eqn{r}. The \code{X_0} matrix contains the "parametric part" of the basis (i.e., the intercept). The matrix \code{X_1} contains the "nonparametric part" of the basis, which consists of the \emph{reproducing kernel} function
\deqn{\rho(x, y) = 1 - (x \vee y) + (1 / 2K) * ( x(x-1) + y(y-1) ) + c }
evaluated at all combinations of \code{x} and \code{knots}. The notation \eqn{(x \vee y)} denotes the maximum of \eqn{x} and \eqn{y}, and the constant is \eqn{c = (K-1)(2K-1) / (6K)}.

The penalty matrix consists of the \emph{reproducing kernel} function
\deqn{\rho(x, y) = 1 - (x \vee y) + (1 / 2K) * ( x(x-1) + y(y-1) ) + c }
evaluated at all combinations of \code{x}.
}
\value{
Basis: Matrix of dimension \code{c(length(x), df)} where \code{df = length(knots) + intercept}.

Penalty: Matrix of dimension \code{c(r, r)} where \code{r = length(x)} is the number of knots.
}
\references{
Gu, C. (2013). Smoothing Spline ANOVA Models. 2nd Ed. New York, NY: Springer-Verlag. \url{https://doi.org/10.1007/978-1-4614-5369-7}

Helwig, N. E. (2017). Regression with ordered predictors via ordinal smoothing splines. \emph{Frontiers in Applied Mathematics and Statistics, 3}(15), 1-13. \url{https://doi.org/10.3389/fams.2017.00015}

Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), \emph{SAGE Research Methods Foundations.} \url{https://doi.org/10.4135/9781526421036885885}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
If the inputs \code{x} and \code{knots} are factors, they should have the same levels.

If the inputs \code{x} and \code{knots} are integers, the \code{knots} should be a subset of the input \code{x}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{nominal}} for a basis and penalty for unordered factors.

See \code{\link{polynomial}} for a basis and penalty for numeric variables.
}
\examples{
# generate data
set.seed(0)
n <- 101
x <- factor(sort(rep(LETTERS[1:4], length.out = n)))
knots <- LETTERS[1:3]
eta <- 1:4
y <- eta[x] + rnorm(n, sd = 0.5)

# ordinal smoothing spline basis
X <- basis_ord(x, knots, intercept = TRUE)

# ordinal smoothing spline penalty
Q <- penalty_ord(knots, K = 4)

# padd Q with zeros (for intercept)
Q <- rbind(0, cbind(0, Q))

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
