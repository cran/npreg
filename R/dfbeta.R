dfbeta.ss <- dfbeta.sm <-
  function(model, infl = NULL, ...){
    r.w <- weighted.residuals(model)
    hat <- hatvalues(model)
    M <- sqrt(weights(model)) * model.matrix(model)
    Ginv <- vcov(model) / model$sigma^2
    MGinv <- M %*% Ginv
    (r.w / (1 - hat)) * MGinv
  }

dfbeta.gsm <-
  function(model, infl = NULL, ...){
    r.w <- residuals(model, type = "pearson")
    hat <- hatvalues(model)
    M <- sqrt(weights(model)) * model.matrix(model)
    Ginv <- vcov(model) / model$dispersion
    MGinv <- M %*% Ginv
    (r.w / (1 - hat)) * MGinv
  }