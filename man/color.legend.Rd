\name{color.legend}
\alias{color.legend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adds Color Legend to Plot Margin
}
\description{
This function can be used to add a color legend to the margin of a plot produced by \code{\link{image}}.
}
\usage{
color.legend(zlim, side = 4, col = NULL, ncol = NULL, zlab = "z", 
             zline = 2.5, box = TRUE, zcex = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zlim}{
numeric vector of the form \code{c(min, max)} giving the range of values for the color legend.
}
  \item{side}{
which side (margin) should the legend be added to? 1 = bottom, 2 = left, 3 = top, 4 = right (default).
}
  \item{col}{
colors to use for the legend. Can input the name of a color palette (see \code{\link{hcl.colors}}) or a vector of colors to create a palette (see \code{\link{colorRampPalette}}).  
}
  \item{ncol}{
number of colors to use for the legend. Defaults to \code{length(col)}.
}
  \item{zlab}{
axis label for the color legend.
}
  \item{zline}{
line number to draw axis label.
}
  \item{box}{
add a box around the legend?
}
  \item{zcex}{
scale factor for axis label.
}
  \item{...}{
additional arguments passed to \code{\link{image}} function.
}
}
\details{
The \code{\link{colorRampPalette}} function is used to create a vector of colors of length \code{ncol} that span the colors included in \code{col}. Then the \code{\link{image}} function is used to draw a color legend with values spanning \code{zlim}. 
}
\value{
Produces a color legend.
}
\references{
Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), SAGE Research Methods Foundations. https://doi.org/10.4135/9781526421036885885
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
You will likely need to use \code{par()$plt} or \code{par()$fig} to make enough room in the appropriate margin (see example).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.gsm}} for effect plots from \code{\link{gsm}} objects

\code{\link{plot.sm}} for effect plots from \code{\link{sm}} objects
}
\examples{
# define function
fun <- function(x){
  exp(-rowSums(x^2)/2)
}

# define xgrid
nx <- 101
x <- y <- seq(-3, 3, length.out = nx)
xy <- expand.grid(x1 = x, x2 = y)

# evaluate function
z <- matrix(fun(xy), nx, nx)

# define colors
colors <- c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
            "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
col <- colorRampPalette(colors)(21)

# setup par
oplt <- par()$plt
par(plt = c(0.15, 0.8, oplt[3:4]))

# plot image
image(x, y, z, col = col)

# add legend
par(plt = c(0.85, 0.9, oplt[3:4]), new = TRUE)
color.legend(range(z), col = col, ncol = length(col))

# restore original par()$plt
par(plt = oplt)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ aplot }
\keyword{ color }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
