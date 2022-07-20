\name{deviance}
\alias{deviance.gsm}
\alias{deviance.sm}
\alias{deviance.ss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smooth Model Deviance
}
\description{
Returns the deviance from a fit smoothing spline (fit by \code{\link{ss}}), smooth model (fit by \code{\link{sm}}), or generalized smooth model (fit by \code{\link{gsm}}).
}
\usage{
\method{deviance}{gsm}(object, ...)

\method{deviance}{sm}(object, ...)

\method{deviance}{ss}(object, ...)
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
For \code{ss} and \code{sm} objects, the deviance is caculated assuming iid Gaussian errors. 

For \code{gsm} objects, the deviance is calculated by summing the squared deviance residuals, which are calculated using \code{family(object)$dev.resid}
}
\value{
Deviance of the model \code{object}.
}
\references{
Chambers, J. M. and Hastie, T. J. (1992) \emph{Statistical Models in S}. Wadsworth & Brooks/Cole.

Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), SAGE Research Methods Foundations. \doi{10.4135/9781526421036885885}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ss}}, \code{\link{sm}}, \code{\link{gsm}}

\code{\link{fitted.values}} and \code{\link{residuals}}
}
\examples{
## for 'ss' and 'sm' objects, this function is defined as
function(object, ...){
    sum(weighted.residuals(object)^2, na.rm = TRUE)
  }
  
## for 'gsm' objects, this function is defined as
function(object, ...){
    object$deviance
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
