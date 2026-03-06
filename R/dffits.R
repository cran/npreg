dffits.ss <- dffits.sm <- dffits.gsm <-
  function(model, infl = smooth.influence(model), 
           res = weighted.residuals(model), ...){
    res <- res * sqrt(infl$hat)/(infl$sigma * (1 - infl$hat))
    res[is.infinite(res)] <- NaN
    res
  }