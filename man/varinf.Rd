\name{varinf}
\alias{varinf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance Inflation Factors
}
\description{
Computes variance inflation factors for terms of a smooth model.
}
\usage{
varinf(object, newdata = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class "sm" output by the \code{\link{sm}} function or an object of class "gsm" output by the \code{\link{gsm}} function.
}
  \item{newdata}{
the data used for variance inflation calculation (if \code{NULL} training data are used).  
}
}
\details{
Let \eqn{\kappa_j^2} denote the VIF for the \eqn{j}-th model term. 

Values of \eqn{\kappa_j^2} close to 1 indicate no multicollinearity issues for the \eqn{j}-th term. Larger values of \eqn{\kappa_j^2} indicate that \eqn{\eta_j} has more collinearity with other terms. 

Thresholds of \eqn{\kappa_j^2 > 5} or \eqn{\kappa_j^2 > 10} are typically recommended for determining if multicollinearity is too much of an issue.

To understand these thresholds, note that
\deqn{\kappa_j^2 = \frac{1}{1 - R_j^2}}
where \eqn{R_j^2} is the R-squared for the linear model predicting \eqn{\eta_j} from the remaining model terms.
}
\value{
a named vector containing the variance inflation factors for each effect function (in \code{object$terms}).
}
\references{
Gu, C. (2013). Smoothing spline ANOVA models, 2nd edition. New York: Springer. \doi{10.1007/978-1-4614-5369-7}

Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), SAGE Research Methods Foundations. \doi{10.4135/9781526421036885885}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Suppose that the function can be written as
\deqn{\eta = \eta_0 + \eta_1 + \eta_2 + ... + \eta_p}
where \eqn{\eta_0} is a constant (intercept) term, and \eqn{\eta_j} denotes the \eqn{j}-th effect function, which is assumed to have mean zero. Note that \eqn{\eta_j} could be a main or interaction effect function for all \eqn{j = 1, ..., p}.

Defining the \eqn{p \times p} matrix \eqn{C} with entries
\deqn{C_{jk} = \cos(\eta_j, \eta_k)}
where the cosine is defined with respect to the training data, i.e.,
\deqn{\cos(\eta_j, \eta_k) = \frac{\sum_{i=1}^n \eta_j(x_i) \eta_k(x_i)}{\sqrt{\sum_{i=1}^n \eta_j^2(x_i)} \sqrt{\sum_{i=1}^n \eta_k^2(x_i)}}}

The variane inflation factors are the diagonal elements of \eqn{C^{-1}}, i.e.,
\deqn{\kappa_j^2 = C^{jj}}
where \eqn{\kappa_j^2} is the VIF for the \eqn{j}-th term, and \eqn{C^{jj}} denotes the \eqn{j}-th diagonal element of the matrix \eqn{C^{-1}}. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{summary.sm}} for more thorough summaries of smooth models.

See \code{\link{summary.gsm}} for more thorough summaries of generalized smooth models.
}
\examples{
##########   EXAMPLE 1   ##########
### 4 continuous predictors
### no multicollinearity

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

# fit model
mod <- sm(y ~ x1v + x2v + x3v + x4v, data = data, tprk = FALSE)

# check vif
varinf(mod)


##########   EXAMPLE 2   ##########
### 4 continuous predictors
### multicollinearity

# generate data
set.seed(1)
n <- 100
fun <- function(x){
  sin(pi*x[,1]) + sin(2*pi*x[,2]) + sin(3*pi*x[,3]) + sin(3*pi*x[,4])
}
data <- as.data.frame(replicate(3, runif(n)))
data <- cbind(data, c(data[1,2], data[2:n,3]))
colnames(data) <- c("x1v", "x2v", "x3v", "x4v")
fx <- fun(data)
y <- fx + rnorm(n)

# check collinearity
cor(data)
cor(sin(3*pi*data[,3]), sin(3*pi*data[,4]))

# fit model
mod <- sm(y ~ x1v + x2v + x3v + x4v, data = data, tprk = FALSE)

# check vif
varinf(mod)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ regression }
\keyword{ smooth }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{ multicollinearity }
% \concept{ ~cpt2 }
% Use only one concept per line.
