wtd.var <-
  function(x, weights, na.rm = FALSE){
    # weighted variance
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2021-10-29
    
    # check data and weights
    nobs <- length(x)
    if(missing(weights)){
      weights <- rep(1, nobs)
    } else {
      if(length(weights) != nobs)
        stop("Inputs x and weights must have the same length")
    }
    
    # check for na
    if(na.rm){
      ix <- which(is.na(x))
      x <- x[-ix]
      weights <- weights[-ix]
      nobs <- length(x)
    }
    
    # calculate weighted mean
    xbar <- wtd.mean(x = x, weights = weights, na.rm = na.rm)
    
    # number of non-zero weights
    nw <- sum(weights != 0)
    
    # weighted variance
    (nw / (nw - 1)) * sum(weights * (x - xbar)^2) / sum(weights)
    
  }