wtd.mean <- 
  function(x, weights, trim = 0, na.rm = FALSE){
    # weighted (and possibly trimmed) mean
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2021-10-28
    
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
    
    # check trim
    trim <- as.numeric(trim[1])
    if(trim < 0 | trim >= 0.5) stop("Input trim must satisfy:  0 <= trim < 0.5")
    
    # trim data
    if(trim > 0){
      bounds <- wtd.quantile(x, weights, c(trim, 1 - trim))
      ix <- which(x < bounds[1] | x > bounds[2])
      x <- x[-ix]
      weights <- weights[-ix]
    }
    
    # weighted mean
    sum(weights * x) / sum(weights)
    
  }