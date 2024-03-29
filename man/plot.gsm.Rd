\name{plot.gsm}
\alias{plot.gsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Effects for Generalized Smooth Model Fits
}
\description{
Plots the main and two-way interaction effects for objects of class "gsm".
}
\usage{
\method{plot}{gsm}(x, terms = x$terms, se = TRUE, n = 201, intercept = FALSE,
     ask = prod(par("mfcol")) < length(terms) && dev.interactive(),
     zero.line = TRUE, zero.lty = 3, zero.col = "black", ncolor = 21, 
     colors = NULL, rev = FALSE, zlim = NULL, lty.col = NULL, 
     legend.xy = "top", main = NULL, xlab = NULL, ylab = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a fit from \code{gsm}.
}
  \item{terms}{
which terms to include in the plot. The default plots all terms.  
}
  \item{se}{
a switch indicating if standard errors are required.  
}
  \item{n}{
number of points to use for line plots. Note \code{sqrt(n)} points are used for image plots.
}
  \item{intercept}{
a switch indicating if an intercept should be added to the effect plot(s).
}
  \item{ask}{
a swith indicating if the user should be prompted before switching plots (if \code{length(terms) > 1})
}
  \item{zero.line}{
a switch indicating if the zero line should be added to the effect plot(s).
}
  \item{zero.lty}{
line type for the zero line (if \code{zero.line = TRUE}).
}
  \item{zero.col}{
color for the zero line (if \code{zero.line = TRUE}).
}
  \item{ncolor}{
number of colors to use for image plot(s).
}
  \item{colors}{
colors to use for image plots. Can input the name of a color palette (see \code{\link{hcl.colors}}) or a vector of colors to create a palette (see \code{\link{colorRampPalette}}).
}
  \item{rev}{
if \code{colors} is the name of a palette, should it be reversed? See \code{\link{hcl.colors}}.
}
  \item{zlim}{
limits to use for image plot(s) when mapping numbers to colors.  
}
  \item{lty.col}{
color(s) to use for lines when plotting effects of continuous predictors.
}
  \item{legend.xy}{
location to place the legend for line plots involving interactions.
}
  \item{main}{
title for plot (ignored unless plotting a single term).
}
  \item{xlab}{
x-axis label for plot (ignored unless plotting a single term).
}
  \item{ylab}{
y-axis label for plot (ignored unless plotting a single term).
}
  \item{...}{
additional arguments passed to \code{\link{plotci}} or \code{\link{image}}
}
}
\details{
Plots main and two-way interaction effects for fit smooth models using either line or image plots. The \code{terms} arugment can be used to plot a specific effect term. Main and interaction effects are plotted by creating predictions from the fit model that only include the requested terms (see \code{\link{predict.sm}}), and then using either the \code{\link{plotci}} function (for line plots) or the \code{\link{image}} function (for heatmaps).
}
\value{
Produces a line or image plot for each requested term in the model.
}
\references{
Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), SAGE Research Methods Foundations. https://doi.org/10.4135/9781526421036885885
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Three-way interaction effects are not plotted.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gsm}} for fitting \code{sm} objects.
}
\examples{
# see examples in gsm() help file
?gsm
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }% use one of  RShowDoc("KEYWORDS")
\keyword{ regression }% __ONLY ONE__ keyword per line
\keyword{ smooth }
