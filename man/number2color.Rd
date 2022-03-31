\name{number2color}
\alias{number2color}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Map Numbers to Colors
}
\description{
Each of the \eqn{n} elements of a numeric vector is mapped onto one of the \eqn{m} specified colors.
}
\usage{
number2color(x, colors, ncol = 21, equidistant = TRUE, xmin = min(x), xmax = max(x))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric vector of observations that should be mapped to colors
}
  \item{colors}{
an optional vector of colors (see Note for default colors)  
}
  \item{ncol}{
number of colors \eqn{m} used for mapping
}
  \item{equidistant}{
if \code{TRUE} (default), the breaks used for binning are an equidistant seqeunce of values spanning the range of \code{x}. Otherwise sample quantiles of \code{x} are used to define the bin breaks.  
}
  \item{xmin}{
minimum \code{x} value to use when defining breaks  
}
  \item{xmax}{
maximum \code{x} value to use when defining breaks  
}
}
\details{
Elements of a numeric vector are binned using either an equidistant sequence (default) or sample quantiles. Each bin is associated with a unique color, so binning the observations is equivalent to mapping the numbers to colors. The \code{colors} are input to the \code{\link{colorRampPalette}} function to create a color palette with length specified by the \code{ncol} argument.
}
\value{
Returns of vector of colors the same length as \code{x}
}

\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
If \code{colors} is missing, the default color palette is defined as

\code{colors <- c("darkblue", rainbow(12)[c(9, 8, 7, 5, 3, 2, 1)], "darkred")}

which is a modified version of the \code{\link{rainbow}} color palette.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{.bincode}} is used to bin the data
}
\examples{
x <- 1:100
xcol <- number2color(x)
plot(x, col = xcol)
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
