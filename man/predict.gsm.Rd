\name{predict.gsm}
\alias{predict.gsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict method for Generalized Smooth Model Fits
}
\description{
\code{predict} method for class "gsm".
}
\usage{
\method{predict}{gsm}(object, newdata = NULL, se.fit = FALSE, 
        type = c("link", "response", "terms"), 
        terms = NULL, na.action = na.pass,
        intercept = NULL, combine = TRUE, design = FALSE, 
        check.newdata = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fit from \code{gsm}. %or \code{smnet}.
} 
  \item{newdata}{
an optional list or data frame in which to look for variables with which to predict. If omitted, the original data are used.
}
  \item{se.fit}{
a switch indicating if standard errors are required.
}
  \item{type}{
type of prediction (link, response, or model term). Can be abbreviated.
}
  \item{terms}{
which terms to include in the fit. The default of \code{NULL} uses all terms. This input \bold{is} used regardless of the \code{type} of prediction.
}
  \item{na.action}{
function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.  
}
  \item{intercept}{
a switch indicating if the intercept should be included in the prediction. If \code{NULL} (default), the intercept is included in the fit only when \code{type = "r"} and \code{terms} includes all model terms.
}
  \item{combine}{
a switch indicating if the parametric and smooth components of the prediction should be combined (default) or returned separately.  
}
  \item{design}{
a switch indicating if the model (design) matrix for the prediction should be returned.
}
  \item{check.newdata}{
a switch indicating if the \code{newdata} should be checked for consistency (e.g., class and range). Ignored if \code{newdata} is not provided.
}
  \item{...}{
additional arguments affecting the prediction produced (currently ignored).  
}
}
\details{
Inspired by the \code{\link{predict.glm}} function in R's \bold{stats} package.

Produces predicted values, obtained by evaluating the regression function in the frame \code{newdata} (which defaults to \code{model.frame(object)}). If the logical \code{se.fit} is \code{TRUE}, standard errors of the predictions are calculated. 

If \code{newdata} is omitted the predictions are based on the data used for the fit. Regardless of the \code{newdata} argument, how cases with missing values are handled is determined by the \code{na.action} argument. If \code{na.action = na.omit} omitted cases will not appear in the predictions, whereas if \code{na.action = na.exclude} they will appear (in predictions and standard errors), with value \code{NA}.

Similar to the \code{glm} function, setting \code{type = "terms"} returns a matrix giving the predictions for each of the requested model \code{terms}. Unlike the \code{glm} function, this function allows for predictions using any subset of the model terms. Specifically, the predictions (on both the \code{link} and \code{response} scale) will only include the requested \code{terms}, which makes it possible to obtain estimates (and standard errors) for subsets of model terms. In this case, the \code{newdata} only needs to contain data for the subset of variables that are requested in \code{terms}.
}
\value{
Default use returns a vector of predictions. Otherwise the form of the output will depend on the combination of argumments: \code{se.fit}, \code{type}, \code{combine}, and \code{design}.

\code{type = "link"}: \cr
When \code{se.fit = FALSE} and \code{design = FALSE}, the output will be the predictions on the link scale. When \code{se.fit = TRUE} or \code{design = TRUE}, the output is a list with components \code{fit}, \code{se.fit} (if requested), and \code{X} (if requested). 

\code{type = "response"}: \cr
When \code{se.fit = FALSE} and \code{design = FALSE}, the output will be the predictions on the data scale. When \code{se.fit = TRUE} or \code{design = TRUE}, the output is a list with components \code{fit}, \code{se.fit} (if requested), and \code{X} (if requested). 

\code{type = "terms"}: \cr
When \code{se.fit = FALSE} and \code{design = FALSE}, the output will be the predictions for each term on the link scale. When \code{se.fit = TRUE} or \code{design = TRUE}, the output is a list with components \code{fit}, \code{se.fit} (if requested), and \code{X} (if requested). 

Regardless of the \code{type}, setting \code{combine = FALSE} decomposes the requested result(s) into the \bold{p}arametric and \bold{s}mooth contributions.

}
\references{
https://stat.ethz.ch/R-manual/R-devel/library/stats/html/predict.glm.html

Craven, P. and Wahba, G. (1979). Smoothing noisy data with spline functions: Estimating the correct degree of smoothing by the method of generalized cross-validation. \emph{Numerische Mathematik, 31}, 377-403. \doi{10.1007/BF01404567}

Gu, C. (2013). Smoothing spline ANOVA models, 2nd edition. New York: Springer. \doi{10.1007/978-1-4614-5369-7}

Helwig, N. E. (2020). Multiple and Generalized Nonparametric Regression. In P. Atkinson, S. Delamont, A. Cernat, J. W. Sakshaug, & R. A. Williams (Eds.), \emph{SAGE Research Methods Foundations.} \doi{10.4135/9781526421036885885}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gsm}}
}
\examples{
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

# fit gsm with specified knots (tprk = TRUE)
gsm.ssa <- gsm(y ~ x * z, family = binomial, knots = knots)
pred <- predict(gsm.ssa)
term <- predict(gsm.ssa, type = "terms")
mean((gsm.ssa$linear.predictors - pred)^2)
mean((gsm.ssa$linear.predictors - rowSums(term) - attr(term, "constant"))^2)

# fit gsm with specified knots (tprk = FALSE)
gsm.gam <- gsm(y ~ x * z, family = binomial, knots = knots, tprk = FALSE)
pred <- predict(gsm.gam)
term <- predict(gsm.gam, type = "terms")
mean((gsm.gam$linear.predictors - pred)^2)
mean((gsm.gam$linear.predictors - rowSums(term) - attr(term, "constant"))^2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
