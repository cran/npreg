deviance.ss <- deviance.sm <-
  function(object, ...){
    sum(weighted.residuals(object)^2, na.rm = TRUE)
  }

deviance.gsm <-
  function(object, ...){
    object$deviance
  }