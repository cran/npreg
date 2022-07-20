dfbetas.ss <- dfbetas.sm <- dfbetas.gsm <-
  function(model, infl = smooth.influence(model, do.coef = TRUE), ...){
    clamod <- class(model)
    coefvar <- if(clamod == "ss") rowSums(model$fit$cov.sqrt^2) else rowSums(model$cov.sqrt^2)
    sigma0 <- ifelse(clamod == "gsm", sqrt(model$dispersion), model$sigma)
    infl$coefficients / outer(infl$sigma, sqrt(coefvar) / sigma0)
  }