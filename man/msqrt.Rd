\name{msqrt}
\alias{msqrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix (Inverse?) Square Root
}
\description{
Stable computation of the square root (or inverse square root) of a positive semi-definite matrix.
}
\usage{
msqrt(x, inverse = FALSE, symmetric = FALSE, 
      tol = .Machine$double.eps, checkx = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
positive semi-definite matrix
}
  \item{inverse}{
compute inverse square root?
}
  \item{symmetric}{
does the square root need to be symmetric? See Details.
}
  \item{tol}{
tolerance for detecting linear dependencies in \code{x}  
}
  \item{checkx}{
should \code{x} be checked for symmetry using \code{\link{isSymmetric}}?
}
}
\details{
If \code{symmetric = FALSE}, this function computes the matrix \code{z} such that \code{x = tcrossprod(z)}

If \code{symmetric = TRUE}, this function computes the matrix \code{z} such that \code{x = crossprod(z) = tcrossprod(z)}

If \code{inverse = TRUE}, the matrix \code{x} is replaced by the pseudo-inverse of \code{x} in these equations (see \code{\link{psolve}})
}
\value{
The matrix \code{z} that gives the (inverse?) square root of \code{x}. See Details.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
The matrix (inverse?) square root is calculated by (inverting and) square rooting the eigenvalues that are greater than the first value multiplied by \code{tol * nrow(x)}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{psolve}}
}
\examples{
# generate x
set.seed(0)
x <- crossprod(matrix(rnorm(100), 20, 5))

# asymmetric square root (default)
xsqrt <- msqrt(x)
mean(( x - crossprod(xsqrt) )^2)
mean(( x - tcrossprod(xsqrt) )^2)

# symmetric square root
xsqrt <- msqrt(x, symmetric = TRUE)
mean(( x - crossprod(xsqrt) )^2)
mean(( x - tcrossprod(xsqrt) )^2)

# asymmetric inverse square root (default)
xsqrt <- msqrt(x, inverse = TRUE)
mean(( solve(x) - crossprod(xsqrt) )^2)
mean(( solve(x) - tcrossprod(xsqrt) )^2)

# symmetric inverse square root
xsqrt <- msqrt(x, inverse = TRUE, symmetric = TRUE)
mean(( solve(x) - crossprod(xsqrt) )^2)
mean(( solve(x) - tcrossprod(xsqrt) )^2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ algebra }% __ONLY ONE__ keyword per line
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
