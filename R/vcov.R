vcov.ss <- 
  function(object, ...){
    Sigma <- tcrossprod(object$fit$cov.sqrt)
    rownames(Sigma) <- colnames(Sigma) <- names(object$fit$coef)
    Sigma
  }

vcov.sm <- vcov.gsm <- 
  function(object, ...){
    Sigma <- tcrossprod(object$cov.sqrt)
    rownames(Sigma) <- colnames(Sigma) <- names(object$coefficients)
    Sigma
  }
