penalty_sph <- 
  function(x, m = 2, rescale = TRUE){
    # Spherical Smoothing Spline Penalty
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2019-04-04
    
    # initializations
    m <- as.integer(m)
    if(m < 1L | m > 3L) stop("Input 'm' must be 1, 2, or 3.")
    x <- as.matrix(x)
    nx <- nrow(x)
    xdim <- ncol(x)
    if(xdim != 3L) stop("Input 'x' should be a matrix with 3 columns.")
    
    # cos theta
    ssx <- rowSums(x ^ 2)
    if(any(ssx == 0)) ssx[ssx == 0] <- 1
    x <- matrix(1 / sqrt(ssx), nrow = nx, ncol = xdim) * x
    ang <- tcrossprod(x)
    
    # constants
    m <- 2L * m
    alpha <- 1 / (2 * m + 1)
    beta <- 2 * pi * factorial(2 * m)
    
    # return kernel
    if(m == 2L){
      X <- (s2fun(ang) - alpha) / beta
    } else if(m == 4L){
      X <- (s4fun(ang) - alpha) / beta
    } else {
      X <- (s6fun(ang) - alpha) / beta
    }
    
    # rescale?
    if(rescale){
      X <- X / mean(diag(X))
    }
    return(X)
    
  }