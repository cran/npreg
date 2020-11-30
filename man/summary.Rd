\name{summary}
\alias{summary}
\alias{summary.gsm}
\alias{summary.sm}
\alias{summary.ss}
\alias{print.summary.gsm}
\alias{print.summary.sm}
\alias{print.summary.ss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary methods for Fit Models
}
\description{
\code{summary} methods for object classes "gsm", "sm", and "ss".
}
\usage{
\method{summary}{gsm}(object, ...)

\method{summary}{sm}(object, ...)

\method{summary}{ss}(object, ...)

\method{print}{summary.gsm}(x, digits = max(3, getOption("digits") - 3), 
      signif.stars = getOption("show.signif.stars"), ...)

\method{print}{summary.sm}(x, digits = max(3, getOption("digits") - 3), 
      signif.stars = getOption("show.signif.stars"), ...)

\method{print}{summary.ss}(x, digits = max(3, getOption("digits") - 3), 
      signif.stars = getOption("show.signif.stars"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class "gsm" output by the \code{\link{gsm}} function, "sm" output by the \code{\link{sm}} function, or "ss" output by the \code{\link{ss}} function
}
  \item{x}{
an object of class "summary.gsm" output by the \code{\link{summary.gsm}} function, "summary.sm" output by the \code{\link{summary.sm}} function, or "summary.ss" output by the \code{\link{summary.ss}} function.
}
  \item{digits}{
the minimum number of significant digits to be printed in values.  
}
  \item{signif.stars}{
logical. If \code{TRUE}, `significance stars' are printed for each coefficient.  
}
  \item{...}{
additional arguments affecting the summary produced (currently ignored).
}
}
\details{
Summary includes information for assessing the statistical and practical significance of the model terms. 

Statistical inference is conducted via (approximate) frequentist chi-square tests using the Bayesian interpretation of a smoothing spline (Nychka, 1988; Wahba, 1983). 

With multiple smooth terms included in the model, the inferential results may (and likely will) differ slightly depending on the \code{tprk} argument (when using the \code{\link{gsm}} and \code{\link{sm}} functions). 

If significance testing is of interest, the \code{tprk = FALSE} option may be desirable, given that this allows for unique basis function coefficients for each model term.

In all cases, the inferential results are based on a (pseudo) F or chi-square statistic which fails to consider the uncertainty of the smoothing parameter estimation.
}
\value{
  \item{residuals }{the deviance residuals.}
  \item{fstatistic }{the F statistic for testing all effects (parametric and smooth).}
  \item{dev.expl }{the explained deviance.}
  \item{p.table }{the coefficient table for (approximate) inference on the parametric terms.}
  \item{s.table }{the coefficient table for (approximate) inference on the smooth terms.}
  \item{dispersion }{the estimate of the dispersion parameter.}
  \item{r.squared }{the observed coefficient of multiple determination.}
  \item{adj.r.squared }{the adjusted coefficient of multiple determination.}
  \item{kappa }{the collinearity indices. A value of 1 indicates no collinearity, and higher values indicate more collinearity of a given term with other model terms.}
  \item{pi }{the importance indices. Larger values indicate more importance, and the values satisfy \code{sum(pi) = 1}. Note that elements of \code{pi} can be negative.}
  \item{call }{the original function call.}
  \item{family }{the specified family (for gsm objects).}
}
\references{
Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), \emph{SAGE Research Methods Foundations.} \url{https://doi.org/10.4135/9781526421036885885}

Nychka, D. (1988). Bayesian confience intervals for smoothing splines. \emph{Journal of the American Statistical Association, 83(404)}, 1134-1143. \url{https://doi.org/10.2307/2290146}

Wahba, G. (1983). Bayesian "confidence intervals" for the cross-validated smoothing spline. \emph{Journal of the Royal Statistical Society. Series B, 45(1)}, 133-150. \url{https://doi.org/10.1111/j.2517-6161.1983.tb01239.x}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gsm}}, \code{\link{sm}}, and \code{\link{ss}}
}
\examples{
### Example 1: gsm

# generate data
set.seed(1)
n <- 1000
x <- seq(0, 1, length.out = n)
z <- factor(sample(letters[1:3], size = n, replace = TRUE))
fun <- function(x, z){
  mu <- c(-2, 0, 2)
  zi <- as.integer(z)
  fx <- mu[zi] + 3 * x + sin(2 * pi * x + mu[zi]*pi/4)
}
fx <- fun(x, z)
y <- rbinom(n = n, size = 1, p = 1 / (1 + exp(-fx)))

# define marginal knots
probs <- seq(0, 0.9, by = 0.1)
knots <- list(x = quantile(x, probs = probs),
              z = letters[1:3])

# fit sm with specified knots (tprk = TRUE)
gsm.ssa <- gsm(y ~ x * z, family = binomial, knots = knots)
summary(gsm.ssa)

# fit sm with specified knots (tprk = FALSE)
gsm.gam <- gsm(y ~ x * z, family = binomial, knots = knots, tprk = FALSE)
summary(gsm.gam)


### Example 2: sm

# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
z <- factor(sample(letters[1:3], size = n, replace = TRUE))
fun <- function(x, z){
  mu <- c(-2, 0, 2)
  zi <- as.integer(z)
  fx <- mu[zi] + 3 * x + sin(2 * pi * x + mu[zi]*pi/4)
}
fx <- fun(x, z)
y <- fx + rnorm(n, sd = 0.5)

# define marginal knots
probs <- seq(0, 0.9, by = 0.1)
knots <- list(x = quantile(x, probs = probs),
              z = letters[1:3])

# fit sm with specified knots (tprk = TRUE)
sm.ssa <- sm(y ~ x * z, knots = knots)
summary(sm.ssa)

# fit sm with specified knots (tprk = FALSE)
sm.gam <- sm(y ~ x * z, knots = knots, tprk = FALSE)
summary(sm.gam)


### Example 3: ss

# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
fx <- 2 + 3 * x + sin(2 * pi * x)
y <- fx + rnorm(n, sd = 0.5)

# regular smoothing spline
ss.reg <- ss(x, y, nknots = 10)
summary(ss.reg)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
