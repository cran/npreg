\name{predict.sm}
\alias{predict.sm}
%\alias{predict.smnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict method for Smooth Model Fits
}
\description{
\code{predict} method for class "sm". %or class "smnet".
}
\usage{
\method{predict}{sm}(object, newdata = NULL, se.fit = FALSE, 
        interval = c("none", "confidence", "prediction"),
        level = 0.95, type = c("response", "terms"), 
        terms = NULL, na.action = na.pass,
        intercept = NULL, combine = TRUE, design = FALSE, 
        check.newdata = TRUE, ...)
}
%\method{predict}{smnet}(object, newdata = NULL, se.fit = FALSE, 
%        interval = c("none", "confidence", "prediction"),
%        level = 0.95, type = c("response", "terms"), 
%        terms = NULL, na.action = na.pass,
%        intercept = NULL, combine = TRUE, design = FALSE, 
%        check.newdata = TRUE, ...)
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fit from \code{sm}. %or \code{smnet}.
} 
  \item{newdata}{
an optional list or data frame in which to look for variables with which to predict. If omitted, the original data are used.
}
  \item{se.fit}{
a switch indicating if standard errors are required.
}
  \item{interval}{
type of interval calculation. Can be abbreviated.
}
  \item{level}{
tolerance/confidence level.
}
  \item{type}{
type of prediction (response or model term). Can be abbreviated.  
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
Inspired by the \code{\link{predict.lm}} function in R's \bold{stats} package.

Produces predicted values, obtained by evaluating the regression function in the frame \code{newdata} (which defaults to \code{model.frame(object)}). If the logical \code{se.fit} is \code{TRUE}, standard errors of the predictions are calculated. Setting \code{interval}s specifies computation of confidence or prediction (tolerance) intervals at the specified level, sometimes referred to as narrow vs. wide intervals. 

If \code{newdata} is omitted the predictions are based on the data used for the fit. Regardless of the \code{newdata} argument, how cases with missing values are handled is determined by the \code{na.action} argument. If \code{na.action = na.omit} omitted cases will not appear in the predictions, whereas if \code{na.action = na.exclude} they will appear (in predictions, standard errors or interval limits), with value \code{NA}.

Similar to the \code{lm} function, setting \code{type = "terms"} returns a matrix giving the predictions for each of the requested model \code{terms}. Unlike the \code{lm} function, this function allows for predictions using any subset of the model terms. Specifically, when \code{type = "response"} the predictions will only include the requested \code{terms}, which makes it possible to obtain estimates (and standard errors and intervals) for subsets of model terms. In this case, the \code{newdata} only needs to contain data for the subset of variables that are requested in \code{terms}.

}
\value{
Default use returns a vector of predictions. Otherwise the form of the output will depend on the combination of argumments: \code{se.fit}, \code{interval}, \code{type}, \code{combine}, and \code{design}.

\code{type = "response"}: \cr
When \code{se.fit = FALSE} and \code{design = FALSE}, the output will be the predictions (possibly with \code{lwr} and \code{upr} interval bounds). When \code{se.fit = TRUE} or \code{design = TRUE}, the output is a list with components \code{fit}, \code{se.fit} (if requested), and \code{X} (if requested). 

\code{type = "terms"}: \cr
When \code{se.fit = FALSE} and \code{design = FALSE}, the output will be the predictions for each term (possibly with \code{lwr} and \code{upr} interval bounds). When \code{se.fit = TRUE} or \code{design = TRUE}, the output is a list with components \code{fit}, \code{se.fit} (if requested), and \code{X} (if requested). 

Regardless of the \code{type}, setting \code{combine = FALSE} decomposes the requested result(s) into the \bold{p}arametric and \bold{s}mooth contributions. 
}
\references{
https://stat.ethz.ch/R-manual/R-devel/library/stats/html/predict.lm.html

Craven, P. and Wahba, G. (1979). Smoothing noisy data with spline functions: Estimating the correct degree of smoothing by the method of generalized cross-validation. \emph{Numerische Mathematik, 31}, 377-403. https://doi.org/10.1007/BF01404567

Gu, C. (2013). Smoothing spline ANOVA models, 2nd edition. New York: Springer. https://doi.org/10.1007/978-1-4614-5369-7

Helwig, N. E. (in press). Multiple and generalized nonparametric regression. SAGE Research Methods Foundations: An Encyclopedia.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

\seealso{
\code{\link{sm}} %and \code{\link{smnet}}
}
\examples{
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

# fit sm with specified knots
smod <- sm(y ~ x * z, knots = knots)

# get model "response" predictions
fit <- predict(smod)
mean((smod$fitted.values - fit)^2)

# get model "terms" predictions
trm <- predict(smod, type = "terms")
attr(trm, "constant")
head(trm)
mean((smod$fitted.values - rowSums(trm) - attr(trm, "constant"))^2)

# get predictions with "newdata" (= the original data)
fit <- predict(smod, newdata = data.frame(x = x, z = z))
mean((fit - smod$fitted.values)^2)

# get predictions and standard errors
fit <- predict(smod, se.fit = TRUE)
mean((fit$fit - smod$fitted.values)^2)
mean((fit$se.fit - smod$se.fit)^2)

# get 99\% confidence interval
fit <- predict(smod, interval = "c", level = 0.99)
head(fit)

# get 99\% prediction interval
fit <- predict(smod, interval = "p", level = 0.99)
head(fit)

# get predictions only for x main effect
fit <- predict(smod, newdata = data.frame(x = x), 
               se.fit = TRUE, terms = "x")
plotci(x, fit$fit, fit$se.fit)

# get predictions only for each group
fit.a <- predict(smod, newdata = data.frame(x = x, z = "a"), se.fit = TRUE)
fit.b <- predict(smod, newdata = data.frame(x = x, z = "b"), se.fit = TRUE)
fit.c <- predict(smod, newdata = data.frame(x = x, z = "c"), se.fit = TRUE)

# plot results (truth as dashed line)
plotci(x = x, y = fit.a$fit, se = fit.a$se.fit,
       col = "red", col.ci = "pink", ylim = c(-6, 6))
lines(x, fun(x, rep(1, n)), lty = 2, col = "red")
plotci(x = x, y = fit.b$fit, se = fit.b$se.fit,
       col = "blue", col.ci = "cyan", add = TRUE)
lines(x, fun(x, rep(2, n)), lty = 2, col = "blue")
plotci(x = x, y = fit.c$fit, se = fit.c$se.fit,
       col = "darkgreen", col.ci = "lightgreen", add = TRUE)
lines(x, fun(x, rep(3, n)), lty = 2, col = "darkgreen")

# add legends
legend("bottomleft", legend = c("Truth", "Estimate", "CI"), 
       lty = c(2, 1, NA), lwd = c(1, 2, NA), 
       col = c("black", "black","gray80"),
       pch = c(NA, NA, 15), pt.cex = 2, bty = "n")
legend("bottomright", legend = letters[1:3], 
       lwd = 2, col = c("red", "blue", "darkgreen"), bty = "n")

}

%## sparse smooth model
%
%# fit smnet with specified knots
%smod <- smnet(y ~ x * z, knots = knots)
%
%# get model "response" predictions
%fit <- predict(smod)
%mean((smod$fitted.values - fit)^2)
%
%# get model "terms" predictions
%trm <- predict(smod, type = "terms")
%attr(trm, "constant")
%head(trm)
%mean((smod$fitted.values - rowSums(trm) - attr(trm, "constant"))^2)
%
%# get predictions with "newdata" (= the original data)
%fit <- predict(smod, newdata = data.frame(x = x, z = z))
%mean((fit - smod$fitted.values)^2)
%
%# get predictions and standard errors
%fit <- predict(smod, se.fit = TRUE)
%mean((fit$fit - smod$fitted.values)^2)
%mean((fit$se.fit - smod$se.fit)^2)
%
%# get 99\% confidence interval
%fit <- predict(smod, interval = "c", level = 0.99)
%head(fit)
%
%# get 99\% prediction interval
%fit <- predict(smod, interval = "p", level = 0.99)
%head(fit)

%# get predictions only for x main effect
%fit <- predict(smod, newdata = data.frame(x = x), 
%               se.fit = TRUE, terms = "x")
%plotci(x, fit$fit, fit$se.fit)

%# get predictions only for each group
%fit.a <- predict(smod, newdata = data.frame(x = x, z = "a"), se.fit = TRUE)
%fit.b <- predict(smod, newdata = data.frame(x = x, z = "b"), se.fit = TRUE)
%fit.c <- predict(smod, newdata = data.frame(x = x, z = "c"), se.fit = TRUE)

%# plot results (truth as dashed line)
%plotci(x = x, y = fit.a$fit, se = fit.a$se.fit,
%       col = "red", col.ci = "pink", ylim = c(-6, 6))
%lines(x, fun(x, rep(1, n)), lty = 2, col = "red")
%plotci(x = x, y = fit.b$fit, se = fit.b$se.fit,
%       col = "blue", col.ci = "cyan", add = TRUE)
%lines(x, fun(x, rep(2, n)), lty = 2, col = "blue")
%plotci(x = x, y = fit.c$fit, se = fit.c$se.fit,
%       col = "darkgreen", col.ci = "lightgreen", add = TRUE)
%lines(x, fun(x, rep(3, n)), lty = 2, col = "darkgreen")

%# add legends
%legend("bottomleft", legend = c("Truth", "Estimate", "CI"), 
%       lty = c(2, 1, NA), lwd = c(1, 2, NA), 
%       col = c("black", "black", "gray80"),
%       pch = c(NA, NA, 15), pt.cex = 2, bty = "n")
%legend("bottomright", legend = letters[1:3], 
%       lwd = 2, col = c("red", "blue", "darkgreen"), bty = "n")

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ smooth }% __ONLY ONE__ keyword per line
