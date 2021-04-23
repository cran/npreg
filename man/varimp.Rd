\name{varimp}
\alias{varimp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variable Importance
}
\description{
Computes variable importance indices for terms of a smooth model.
}
\usage{
varimp(object, combine = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class "sm" output by the \code{\link{sm}} function or an object of class "gsm" output by the \code{\link{gsm}} function.
}
  \item{combine}{
a switch indicating if the parametric and smooth components of the importance should be combined (default) or returned separately.  
}
}
\details{
Suppose that the function can be written as
\deqn{\eta = \eta_0 + \eta_1 + \eta_2 + ... + \eta_p}
where \eqn{\eta_0} is a constant (intercept) term, and \eqn{\eta_j} denotes the \eqn{j}-th effect function, which is assumed to have mean zero. Note that \eqn{\eta_j} could be a main or interaction effect function for all \eqn{j = 1, ..., p}.

The variable importance index for the \eqn{j}-th effect term is defined as
\deqn{\pi_j = (\eta_j^\top \eta_*) / (\eta_*^\top \eta_*)}
where \eqn{\eta_* = \eta_1 + \eta_2 + ... + \eta_p}. Note that \eqn{\sum_{j = 1}^p \pi_j = 1} but there is no guarantee that \eqn{\pi_j > 0}. 

If all \eqn{\pi_j} are non-negative, then \eqn{\pi_j} gives the proportion of the model's R-squared that can be accounted for by the \eqn{j}-th effect term. Thus, values of \eqn{\pi_j} closer to 1 indicate that \eqn{\eta_j} is more important, whereas values of \eqn{\pi_j} closer to 0 (including negative values) indicate that \eqn{\eta_j} is less important.
}
\value{
If \code{combine = TRUE}, returns a named vector containing the importance indices for each effect function (in \code{object$terms}).

If \code{combine = FALSE}, returns a data frame where the first column gives the importance indices for the \code{p}arametric components and the second column gives the importance indices for the \code{s}mooth (nonparametric) components.
}
\references{
Gu, C. (2013). Smoothing spline ANOVA models, 2nd edition. New York: Springer. \doi{10.1007/978-1-4614-5369-7}

Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), SAGE Research Methods Foundations. \doi{10.4135/9781526421036885885}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
When \code{combine = FALSE}, importance indices will be equal to zero for non-existent components of a model term. For example, a \code{\link{nominal}} effect does not have a parametric component, so the \code{$p} component of the importance index for a nominal effect will be zero.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{summary.sm}} for more thorough summaries of smooth models.

See \code{\link{summary.gsm}} for more thorough summaries of generalized smooth models.
}
\examples{

##########   EXAMPLE 1   ##########
### 1 continuous and 1 nominal predictor

# generate data
set.seed(1)
n <- 100
x <- seq(0, 1, length.out = n)
z <- factor(sample(letters[1:3], size = n, replace = TRUE))
fun <- function(x, z){
  mu <- c(-2, 0, 2)
  zi <- as.integer(z)
  fx <- mu[zi] + 3 * x + sin(2 * pi * x)
}
fx <- fun(x, z)
y <- fx + rnorm(n, sd = 0.5)

# define marginal knots
probs <- seq(0, 0.9, by = 0.1)
knots <- list(x = quantile(x, probs = probs),
              z = letters[1:3])

# fit correct (additive) model
sm.add <- sm(y ~ x + z, knots = knots)

# fit incorrect (interaction) model
sm.int <- sm(y ~ x * z, knots = knots)

# true importance indices
eff <- data.frame(x = 3 * x + sin(2 * pi * x), z = c(-2, 0, 2)[as.integer(z)])
eff <- scale(eff, scale = FALSE)
fstar <- rowSums(eff)
colSums(eff * fstar) / sum(fstar^2)

# estimated importance indices
varimp(sm.add)
varimp(sm.int)



##########   EXAMPLE 2   ##########
### 4 continuous predictors
### additive model

# generate data
set.seed(1)
n <- 100
fun <- function(x){
  sin(pi*x[,1]) + sin(2*pi*x[,2]) + sin(3*pi*x[,3]) + sin(4*pi*x[,4])
}
data <- as.data.frame(replicate(4, runif(n)))
colnames(data) <- c("x1v", "x2v", "x3v", "x4v")
fx <- fun(data)
y <- fx + rnorm(n)

# define ssa knot indices
knots.indx <- c(bin.sample(data$x1v, nbin = 10, index.return = TRUE)$ix,
                bin.sample(data$x2v, nbin = 10, index.return = TRUE)$ix,
                bin.sample(data$x3v, nbin = 10, index.return = TRUE)$ix,
                bin.sample(data$x4v, nbin = 10, index.return = TRUE)$ix)

# fit correct (additive) model
sm.add <- sm(y ~ x1v + x2v + x3v + x4v, data = data, knots = knots.indx)

# fit incorrect (interaction) model
sm.int <- sm(y ~ x1v * x2v + x3v + x4v, data = data, knots = knots.indx)

# true importance indices
eff <- data.frame(x1v = sin(pi*data[,1]), x2v = sin(2*pi*data[,2]),
                  x3v = sin(3*pi*data[,3]), x4v = sin(4*pi*data[,4]))
eff <- scale(eff, scale = FALSE)
fstar <- rowSums(eff)
colSums(eff * fstar) / sum(fstar^2)

# estimated importance indices
varimp(sm.add)
varimp(sm.int)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ regression }
\keyword{ smooth }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{ importance }
% \concept{ ~cpt2 }
% Use only one concept per line.
