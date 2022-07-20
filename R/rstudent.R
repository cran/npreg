rstudent.ss <- rstudent.sm <-
  function(model, infl = influence(model, do.coef = FALSE), 
           res = infl$wt.res, ...){
    res <- res / (infl$sigma * sqrt(1 - infl$hat))
    res[is.infinite(res)] <- NaN
    res
  }

rstudent.gsm <-
  function(model, infl = influence(model, do.coef = FALSE), ...){
    r <- infl$dev.res
    r <- sign(r) * sqrt( r^2 + (infl$hat * infl$pear.res^2) / (1 - infl$hat) )
    r[is.infinite(r)] <- NaN
    if (any(family(model)$family == c("binomial", "poisson", "NegBin"))) 
      r
    else 
      r / infl$sigma
  }