\name{vcov}
\alias{vcov.ss}
\alias{vcov.sm}
\alias{vcov.gsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Variance-Covariance Matrix for a Fitted Smooth Model
}
\description{
Returns the variance-covariance matrix for the basis function coefficients from a fit smoothing spline (fit by \code{\link{ss}}), smooth model (fit by \code{\link{sm}}), or generalized smooth model (fit by \code{\link{gsm}}).
}
\usage{
\method{vcov}{ss}(object, ...)

\method{vcov}{sm}(object, ...)

\method{vcov}{gsm}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class "gsm" output by the \code{\link{gsm}} function, "sm" output by the \code{\link{sm}} function, or "ss" output by the \code{\link{ss}} function
}
  \item{...}{
other arugments (currently ignored)  
}
}
\details{
The variance-covariance matrix is calculated using the Bayesian interpretation of a smoothing spline. Unlike the classic treatments (e.g., Wahba, 1983; Nychka, 1988), which interpret the smoothing spline as a Bayesian estimate of a Gaussian process, this treatment applies the Bayesian interpretation directly on the coefficient vector. More specifically, the smoothing spline basis function coefficients are interpreted as Bayesian estimates of the basis function coefficients (see Helwig, 2020). 
}
\value{
Returns the (symmetric) matrix such that cell (\eqn{i,j}) contains the covariance between the \eqn{i}-th and \eqn{j}-th elements of the coefficient vector.
}
\references{
Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), SAGE Research Methods Foundations. \doi{10.4135/9781526421036885885}

Nychka, D. (1988). Bayesian confience intervals for smoothing splines. \emph{Journal of the American Statistical Association, 83(404)}, 1134-1143. \doi{10.2307/2290146}

Wahba, G. (1983). Bayesian "confidence intervals" for the cross-validated smoothing spline. \emph{Journal of the Royal Statistical Society. Series B, 45(1)}, 133-150. \doi{10.1111/j.2517-6161.1983.tb01239.x}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ss}}, \code{\link{sm}}, \code{\link{gsm}} for model fitting

\code{\link{boot.ss}}, \code{\link{boot.sm}}, \code{\link{boot.gsm}} for bootstrapping
}
\examples{
## for 'ss' objects this function is defined as
  function(object, ...){
    Sigma <- tcrossprod(object$fit$cov.sqrt)
    rownames(Sigma) <- colnames(Sigma) <- names(object$fit$coef)
    Sigma
  }

## for 'sm' and 'gsm' objects this function is defined as
  function(object, ...){
    Sigma <- tcrossprod(object$cov.sqrt)
    rownames(Sigma) <- colnames(Sigma) <- names(object$coefficients)
    Sigma
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
