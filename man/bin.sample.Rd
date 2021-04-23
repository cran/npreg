\name{bin.sample}
\alias{bin.sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bin Sample a Vector, Matrix, or Data Frame
}
\description{
Bin elements of a vector (or rows of a matrix/data frame) and randomly sample a specified number of elements from each bin. Returns sampled data and (optionally) indices of sampled data and/or breaks for defining bins.
}
\usage{
bin.sample(x, nbin = 5, size = 1, equidistant = FALSE, 
           index.return = FALSE, breaks.return = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector, matrix, or data frame to bin sample. Factors are allowed.
}
  \item{nbin}{
Number of bins for each variable (defaults to 5 bins for each dimension of \code{x}). If \code{length(bins) != ncol(x)}, then \code{nbin[1]} is used for each variable. 
}
  \item{size}{
Size of sample to randomly draw from each bin (defaults to 1).
}
  \item{equidistant}{
Should bins be defined equidistantly for each predictor? If \code{FALSE} (default), sample quantiles define bins for each predictor. If \code{length(equidistant) != ncol(x)}, then \code{equidistant[1]} is used for each variable.
}
  \item{index.return}{
If \code{TRUE}, returns the (row) indices of the bin sampled observations.
}
  \item{breaks.return}{
If \code{TRUE}, returns the (lower bounds of the) breaks for the binning.
}
}
\details{
For a single variable, the unidimensional bins are defined using the \code{\link{.bincode}} function. For multiple variables, the multidimensional bins are defined using the algorithm described in the appendix of Helwig et al. (2015), which combines the unidimensional bins (calculated via \code{\link{.bincode}}) into a multidimensional bin code.
}
\value{
If \code{index.return = FALSE} and \code{breaks.return = FALSE}, returns the bin sampled \code{x} observations.

If \code{index.return = TRUE} and/or \code{breaks.return = TRUE}, returns a list with elements:
\item{x }{bin sampled \code{x} observations.}
\item{ix }{row indices of bin sampled observations (if \code{index.return = TRUE}).}
\item{bx }{lower bounds of breaks defining bins (if \code{breaks.return = TRUE}).}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Helwig, N. E., Gao, Y., Wang, S., & Ma, P. (2015). Analyzing spatiotemporal trends in social media data via smoothing spline analysis of variance. \emph{Spatial Statistics, 14}(C), 491-504. \doi{10.1016/j.spasta.2015.09.002}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
For factors, the number of bins is automatically defined to be the number of levels.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{.bincode}} for binning a numeric vector
}
\examples{
##########   EXAMPLE 1   ##########
### unidimensional binning

# generate data
x <- seq(0, 1, length.out = 101)

# bin sample (default)
set.seed(1)
bin.sample(x)

# bin sample (return indices)
set.seed(1)
xs <- bin.sample(x, index.return = TRUE)
xs$x             # sampled data
x[xs$ix]         # indexing sampled data

# bin sample (return indices and breaks)
set.seed(1)
xs <- bin.sample(x, index.return = TRUE, breaks.return = TRUE)
xs$x             # sampled data
x[xs$ix]         # indexing sampled data
xs$bx            # breaks



##########   EXAMPLE 2   ##########
### bidimensional binning

# generate data
x <- expand.grid(x1 = seq(0, 1, length.out = 101),
                 x2 = seq(0, 1, length.out = 101))

# bin sample (default)
set.seed(1)
bin.sample(x)

# bin sample (return indices)
set.seed(1)
xs <- bin.sample(x, index.return = TRUE)
xs$x             # sampled data
x[xs$ix,]        # indexing sampled data

# bin sample (return indices and breaks)
set.seed(1)
xs <- bin.sample(x, index.return = TRUE, breaks.return = TRUE)
xs$x             # sampled data
x[xs$ix,]        # indexing sampled data
xs$bx            # breaks

# plot breaks and 25 bins
plot(xs$bx, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "x1", ylab = "x2", main = "25 bidimensional bins")
grid()
text(xs$bx + 0.1, labels = 1:25)

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
