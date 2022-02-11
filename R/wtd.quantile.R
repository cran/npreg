wtd.quantile <- 
  function(x, weights, probs = seq(0, 1, 0.25), 
           na.rm = FALSE, names = TRUE){
    # weighted quantiles
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2021-10-27
    
    # check data and weights
    nobs <- length(x)
    if(missing(weights)){
      weights <- rep(1, nobs)
    } else {
      if(length(weights) != nobs)
        stop("Inputs 'x' and 'weights' must have the same length")
    }
    
    # check probs
    probs <- sort(probs)
    p <- length(probs)
    if(probs[1] < 0 | probs[p] > 1) stop("Input 'probs' must be between 0 and 1")
    if(names) names(probs) <- paste0(100 * round(probs, 5), "%")
    
    # check for na
    if(na.rm){
      ix <- which(is.na(x))
      x <- x[-ix]
      weights <- weights[-ix]
      nobs <- length(x)
    }
    
    # sort x
    xs <- sort(x, index.return = TRUE)
    x <- xs$x
    weights <- weights[xs$ix]
    
    # cumulative sum of weights
    wcs <- cumsum(weights) / sum(weights)
    
    # interpolate x at probs
    q <- approx(wcs, x, xout = probs, rule = 2)$y
    
    # return results
    if(names) names(q) <- names(probs)
    return(q)
    
  }