cov.ratio <-
  function(model, infl = smooth.influence(model, do.coef = FALSE),
           res = weighted.residuals(model)){
    clamod <- class(model)
    sigma <- ifelse(clamod == "gsm", sqrt(model$dispersion), model$sigma)
    (1 / (1 - infl$hat)) * (infl$sigma / sigma)^(2 * model$df)
  }