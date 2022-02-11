wtd.sd <-
  function(x, weights, na.rm = FALSE){
    sqrt(wtd.var(x, weights, na.rm))
  }