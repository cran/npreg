cooks.distance.ss <- cooks.distance.sm <- 
  function(model, infl = NULL, res = weighted.residuals(model), 
           sd = model$sigma, hat = hatvalues(model), ...){
    wsqrt <- sqrt(weights(model))
    M <- wsqrt * model.matrix(model)
    Ginv <- vcov(model) / sd^2
    GMtMG <- crossprod(M %*% Ginv)
    dfun0 <- function(x) as.numeric(t(x) %*% GMtMG %*% x)
    si.si <- apply(M, 1, dfun0)
    res <- (res / (sd * (1 - hat)))^2 * si.si / model$df
    res[is.infinite(res)] <- NaN
    res
  }

cooks.distance.gsm <- 
  function(model, infl = NULL, res = residuals(model, type = "pearson"), 
           dispersion = model$dispersion, hat = hatvalues(model), ...){
    mueta <- model$family$mu.eta(model$linear.predictors)
    varx <- model$family$variance(fitted(model))
    wsqrt <- sqrt(weights(model)) * sqrt(mueta^2 / varx)
    M <- wsqrt * model.matrix(model)
    Ginv <- vcov(model) / dispersion
    GMtMG <- crossprod(M %*% Ginv)
    dfun0 <- function(x) as.numeric(t(x) %*% GMtMG %*% x)
    si.si <- apply(M, 1, dfun0)
    res <- (res / (1 - hat))^2 * si.si / (dispersion * model$df)
    res[is.infinite(res)] <- NaN
    res
  }