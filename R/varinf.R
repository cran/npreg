varinf <- 
  function(object, newdata = NULL){
    
    fit.terms <- scale(predict(object, newdata = newdata, type = "terms"), scale = FALSE)
    Cmat <- crossprod(fit.terms) / tcrossprod(sqrt(colSums(fit.terms^2)))
    vif <- diag(psolve(Cmat))
    names(vif) <- object$terms
    vif
    
  } # end varinf