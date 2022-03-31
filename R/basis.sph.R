basis.sph <- 
  function(x, knots, m = 2, intercept = FALSE, ridge = FALSE){
    # Spherical Smoothing Spline Basis
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2022-03-22
    
    # initializations
    m <- as.integer(m)
    if(m < 2L | m > 4L) stop("Input 'm' must be 2, 3, or 4.")
    x <- as.matrix(x)
    knots <- as.matrix(knots)
    nx <- nrow(x)
    nk <- nrow(knots)
    xdim <- ncol(x)
    if(xdim != 2L) stop("Input 'x' should be matrices with 2 columns (latitude, longitude).")
    if(ncol(knots) != xdim) stop("Need inputs 'x' and 'knots' to satisfy:  ncol(x) == ncol(knots)")
    
    # check range of x
    if(any(abs(x[,1]) > 90)) stop("First column of 'x' must be latitude (-90 to 90 degrees)")
    if(any(abs(x[,2]) > 180)) stop("Second column of 'x' must be longitude (-180 to 180 degrees)")
    if(any(abs(knots[,1]) > 90)) stop("First column of 'knots' must be latitude (-90 to 90 degrees)")
    if(any(abs(knots[,2]) > 180)) stop("Second column of 'knots' must be longitude (-180 to 180 degrees)")
    
    # convert to radians
    x <- x * (pi / 180)
    knots <- knots * (pi / 180)
    
    # column names
    knot.names <- paste("knot", 1:nk, sep = ".")
    
    # cosine of angle
    z <- outer(cos(x[,1]), cos(knots[,1])) * cos(outer(x[,2], knots[,2], FUN = "-")) + outer(sin(x[,1]), sin(knots[,1]))
    
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
    
    # check ridge
    if(ridge){
      Q <- penalty.sph(knots * (180 / pi), m = m)
      X <- X %*% msqrt(Q, inverse = TRUE, checkx = FALSE)
      knot.names <- knot.names[1:ncol(X)]
    }
    
    # intercept?
    if(intercept){
      X <- cbind(1, X)
      knot.names <- c("(Intercept)", knot.names)
    }
    
    # return results
    colnames(X) <- knot.names
    return(X)
        
  }

# spherical spline (m = 2)
q2fun <- function(x){
  w <- (1 - x) / 2
  zix <- which(w <= .Machine$double.eps)
  w[zix] <- 0
  a <- log(1 + 1 / sqrt(w))
  a[zix] <- 0
  c <- 2 * sqrt(w)
  (a * (12 * w^2 - 4 * w) - 6 * c * w + 6 * w + 1) / 2
}

# spherical spline (m = 3)
q4fun <- function(x){
  w <- (1 - x) / 2
  zix <- which(w <= .Machine$double.eps)
  w[zix] <- 0
  a <- log(1 + 1 / sqrt(w))
  a[zix] <- 0
  c <- 2 * sqrt(w)
  p1 <- a * (840 * w^4 - 720 * w^3 + 72 * w^2 )
  p2 <- c * (220 * w^2 - 420 * w^3)
  p3 <- 420 * w^3 - 150 * w^2 - 4 * w + 3
  (p1 + p2 + p3) / 12
}

# spherical spline (m = 4)
q6fun <- function(x){
  w <- (1 - x) / 2
  zix <- which(w <= .Machine$double.eps)
  w[zix] <- 0
  a <- log(1 + 1 / sqrt(w))
  a[zix] <- 0
  c <- 2 * sqrt(w)
  p1 <- a * ( 27720 * w^6 - 37800 * w^5 + 12600 * w^4 - 600 * w^3)
  p2 <- c * ( 14280 * w^4 - 13860 * w^5 - 2772 * w^3 )
  p3 <- 13860 * w^5 - 11970 * w^4 + 1470 * w^3 + 15 * w^2 - 3 * w + 5
  (p1 + p2 + p3) / 30
}