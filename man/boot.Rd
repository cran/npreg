\name{boot}
\alias{boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap a Fit Smooth Model
}
\description{
Bootstraps a fit nonparametric regression model to form confidence intervals (BCa or percentile) and standard error estimates. 
}
\usage{
boot(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class "ss", "sm", or "gsm"
}
  \item{...}{
Additional arugments passed to bootstrap functions.
}
}
\details{
Currently implemented for objects of class "ss" only. Implementations for objects of class "sm" and "gsm" will be included in subsequent package updates.
}
\value{
Bootstrap estimates of confidence intervals and standard errors.

See the \code{\link{boot.ss}}, \code{boot.sm}, or \code{boot.gsm} function for specifics.
}
\references{
Davison, A. C., & Hinkley, D. V. (1997). \emph{Bootstrap Methods and Their Application}. Cambridge University Press. \doi{10.1017/CBO9780511802843}

Efron, B., & Tibshirani, R. J. (1994). \emph{An Introduction to the Boostrap}. Chapman & Hall/CRC. \doi{10.1201/9780429246593}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{boot.ss}} for bootstrapping "ss" objects

\code{boot.sm} for bootstrapping "sm" objects (coming soon)

\code{boot.gsm} for bootstrapping "gsm" objects (coming soon)
}
\examples{
?boot.ss
}

\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{univar}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
