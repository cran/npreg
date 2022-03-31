\name{psolve}
\alias{psolve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pseudo-Solve a System of Equations
}
\description{
This generic function solves the equation \code{a \%*\% x = b} for \code{x}, where \code{b} can be either a vector or a matrix. This implementation is similar to \code{\link{solve}}, but uses a pseudo-inverse if the system is computationally singular.
}
\usage{
psolve(a, b, tol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
a rectangular numeric matrix containing the coefficients of the linear system. 
}
  \item{b}{
a numeric vector or matrix giving the right-hand side(s) of the linear system. If missing, \code{b} is taken to be an identity matrix and solve will return the (pseudo-)inverse of \code{a}.
}
  \item{tol}{
the tolerance for detecting linear dependencies in the columns of a. The default is \code{.Machine$double.eps}.
}
}
\details{
If \code{a} is a symmetric matrix, \code{\link{eigen}} is used to compute the (pseudo-)inverse. This assumes that \code{a} is a positive semi-definite matrix. Otherwise \code{svd} is used to compute the (pseudo-)inverse for rectangular matrices.
}
\value{
If \code{b} is missing, returns the (pseudo-)inverse of \code{a}. Otherwise returns \code{psolve(a) \%*\% b}.
}
\references{
Moore, E. H. (1920). On the reciprocal of the general algebraic matrix. \emph{Bulletin of the American Mathematical Society, 26}, 394-395. \doi{10.1090/S0002-9904-1920-03322-7}

Penrose, R. (1955). A generalized inverse for matrices. \emph{Mathematical Proceedings of the Cambridge Philosophical Society, 51(3)}, 406-413. \doi{10.1017/S0305004100030401}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The pseudo-inverse is calculated by inverting the eigen/singular values that are greater than the first value multiplied by \code{tol * min(dim(a))}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{msqrt}}
}
\examples{
# generate X
set.seed(0)
X <- matrix(rnorm(100), 20, 5)
X <- cbind(X, rowSums(X))

# pseudo-inverse of X  (dim = 6 by 20)
Xinv <- psolve(X)

# pseudo-inverse of crossprod(X)  (dim = 6 by 6)
XtXinv <- psolve(crossprod(X))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ algebra }% __ONLY ONE__ keyword per line
