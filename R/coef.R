coef.ss <- 
  function(object, ...){
    object$fit$coef
  }

coef.sm <- coef.gsm <- 
  function(object, ...){
    object$coefficients
  }