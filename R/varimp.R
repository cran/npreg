varimp <-
  function(object, newdata = NULL, combine = TRUE){
    # variable importance for class "sm" and "gsm"
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2022-05-24
    
    # calculate variable importance
    if(combine){
      fit.terms <- scale(predict(object, newdata = newdata, type = "terms"), scale = FALSE)
      fitc <- rowSums(fit.terms)
      pi <- as.numeric(crossprod(fit.terms, fitc) / sum(fitc^2))
      names(pi) <- object$terms
    } else {
      pred <- predict(object, newdata = newdata, type = "terms", combine = combine)
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