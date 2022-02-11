penalty.sph <- 
  function(x, m = 2){
    # Spherical Smoothing Spline Penalty
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2021-07-14
    
    # initializations
    m <- as.integer(m)
    if(m < 2L | m > 4L) stop("Input 'm' must be 2, 3, or 4.")
    x <- as.matrix(x)
    nx <- nrow(x)
    xdim <- ncol(x)
    if(xdim != 2L) stop("Input 'x' should be a matrix with 2 columns (latitude, longitude).")
    
    # check range of x
    if(any(abs(x[,1]) > 90)) stop("First column of 'x' must be latitude (-90 to 90 degrees)")
    if(any(abs(x[,2]) > 180)) stop("Second column of 'x' must be longitude (-180 to 180 degrees)")
    
    # convert to radians
    x <- x * (pi / 180)
    
    # cosine of angle
    z <- tcrossprod(cos(x[,1])) * cos(outer(x[,2], x[,2], FUN = "-")) + tcrossprod(sin(x[,1]))
    
    # constants
    alpha <- 1 / (2 * m - 1)
    beta <- 2 * pi * factorial(2 * m - 2)
    
    # return kernel
    if(m == 2L){
      X <- (q2fun(z) - alpha) / beta
    } else if(m == 3L){
      X <- (q4fun(z) - alpha) / beta
    } else {
      X <- (q6fun(z) - alpha) / beta
    }
    return(X)
    
  }