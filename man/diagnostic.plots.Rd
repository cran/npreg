\name{diagnostic.plots}
\alias{diagnostic.plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Nonparametric Regression Diagnostics
}
\description{
Six regression diagnostic plots for a fit smoothing spline (fit by \code{\link{ss}}), smooth model (fit by \code{\link{sm}}), or generalized smooth model (fit by \code{\link{gsm}}).
}
\usage{
diagnostic.plots(x, which = c(1, 2, 3, 5), 
           caption = list("Residuals vs Fitted", 
                          "Normal Q-Q", "Scale-Location", 
                          "Cook's distance", "Residuals vs Leverage", 
                          "Cook's dist vs Variance ratio"), 
           panel = if (add.smooth) function(x, y, ...) 
                      panel.smooth(x, y, iter = iter.smooth, ...) 
                   else points, 
           sub.caption = NULL, main = "", 
           ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
           ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75, cex.pt = 1,
           qqline = TRUE, cook.levels = c(0.5, 1), add.smooth = getOption("add.smooth"), 
           iter.smooth = if (isGlm) 0 else 3, label.pos = c(4, 2), cex.caption = 1, 
           cex.oma.main = 1.25, cex.lab = 1, line.lab = 3, xlim = NULL, ylim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class "gsm" output by the \code{\link{gsm}} function, "sm" output by the \code{\link{sm}} function, or "ss" output by the \code{\link{ss}} function
}
  \item{which}{subset of the integers \code{1:6} indicating which plots to produce}

  \item{caption}{captions to appear above the plots}
  
  \item{panel}{panel function (panel.smooth or points?)}
  
  \item{sub.caption}{common title (for use above multiple figures)}
  
  \item{main}{title to each plot (in addition to \code{caption})}
  
  \item{ask}{if \code{TRUE}, the user is asked before each plot}
  
  \item{...}{other parameters to be passed through to plotting functions}
  
  \item{id.n}{number of points to be labeled in each plot, starting with the most extreme}
  
  \item{labels.id}{vector of labels for extreme observations (\code{NULL} uses the observation numbers)}
  
  \item{cex.id}{magnification of point labels}
  
  \item{cex.pt}{magnification of points}
  
  \item{qqline}{logical indicating if a \code{\link{qqline}} should be added to the normal Q-Q plot}
  
  \item{cook.levels}{levels of Cook's distance at which to draw contours}
  
  \item{add.smooth}{logical indicating if a smoother should be added to most plots}
  
  \item{iter.smooth}{the number of robustness iterations, the argument \code{iter} in \code{\link{panel.smooth}}}
  
  \item{label.pos}{positioning of the labels, for the left hald and right half of the graph respectively, for plots 1-3, 5, and 6}
  
  \item{cex.caption}{controls the size of the \code{caption}}
  
  \item{cex.oma.main}{controls the size of the \code{sub.caption} only if that is above the figures (when there is more than one figure)}
  \item{cex.lab}{character expansion factor for axis labels}
  \item{line.lab}{on which margin line should the axis labels be drawn?}
  \item{xlim}{Limits for x-axis. If \code{length(which) == 1}, a vector of the form \code{c(xmin, xmax)}. Otherwise a list the same length as \code{which} such that each list entry gives the x-axis limits for the corresponding plot.}
  \item{ylim}{Limits for y-axis. If \code{length(which) == 1}, a vector of the form \code{c(ymin, ymax)}. Otherwise a list the same length as \code{which} such that each list entry gives the y-axis limits for the corresponding plot.}
}
\details{
This function is modeled after the \code{\link{plot.lm}} function. The structure of the arguments, as well as the internal codes, mimics the \code{\link{plot.lm}} function whenever possible. By default, only plots 1-3 and 5 are provided, but any subset of plots can be requested using the \code{which} argument.

The six plots include: (1) residuals versus fitted values, (2) normal Q-Q plot, (3) scale-location plot of \eqn{\sqrt{|residuals|}} versus fitted values, (4) Cook's distances, (5) residuals versus leverages, and (6) Cook's distance versus variance ratio = leverage/(1-leverage).
}
\references{
Belsley, D. A., Kuh, E. and Welsch, R. E. (1980). Regression Diagnostics. New York: Wiley.

Cook, R. D. and Weisberg, S. (1982). Residuals and Influence in Regression. London: Chapman and Hall.

McCullagh, P. and Nelder, J. A. (1989). Generalized Linear Models. London: Chapman and Hall.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ss}}, \code{\link{sm}}, \code{\link{gsm}}

\code{\link{smooth.influence.measures}} and \code{\link{smooth.influence}}
}
\examples{
# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
fx <- 2 + 3 * x + sin(2 * pi * x)
y <- fx + rnorm(n, sd = 0.5)

# smoothing spline
mod.ss <- ss(x, y, nknots = 10)
diagnostic.plots(mod.ss)

# smooth model
mod.sm <- sm(y ~ x, knots = 10)
diagnostic.plots(mod.sm)

# generalized smooth model (family = gaussian)
mod.gsm <- gsm(y ~ x, knots = 10)
diagnostic.plots(mod.gsm)

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
