varimp <-
  function(object, combine = TRUE){
    # variable importance for class "sm" and "gsm"
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2021-04-13
    
    # calculate variable importance
    if(combine){
      fit.terms <- scale(predict(object, type = "terms"), scale = FALSE)
      fitc <- rowSums(fit.terms)
      pi <- as.numeric(crossprod(fit.terms, fitc) / sum(fitc^2))
      names(pi) <- object$terms
    } else {
      pred <- predict(object, type = "terms", combine = combine)
      fit.terms <- scale(cbind(pred$p, pred$s), scale = FALSE)
      fitc <- rowSums(fit.terms)
      pi <- as.numeric(crossprod(fit.terms, fitc) / sum(fitc^2))
      nterms <- length(object$terms)
      pi <- data.frame(p = pi[1:nterms], s = pi[1:nterms + nterms])
      rownames(pi) <- object$terms
    }
    
    # return pi
    pi
    
  } # end varimp