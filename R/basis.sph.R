basis.sph <- 
  function(x, knots, m = 2, rescale = TRUE, intercept = FALSE, ridge = FALSE){
    # Spherical Smoothing Spline Basis
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Update: 2021-04-09
    
    # initializations
    m <- as.integer(m)
    if(m < 1L | m > 3L) stop("Input 'm' must be 1, 2, or 3.")
    x <- as.matrix(x)
    knots <- as.matrix(knots)
    nx <- nrow(x)
    nk <- nrow(knots)
    xdim <- ncol(x)
    if(xdim != 3L) stop("Inputs 'x' and 'knots' should be matrices with 3 columns.")
    if(ncol(knots) != xdim) stop("Need inputs 'x' and 'knots' to satisfy:  ncol(x) == ncol(knots)")
    
    # rescale?
    if(rescale){
      theta <- 0
      for(k in 1:nk){
        theta <- theta + c(penalty.sph(knots[k,,drop=FALSE], m = m, rescale = FALSE))
      }
      theta <- nk / theta
    } else {
      theta <- 1
    }
    
    # column names
    knot.names <- paste("knot", 1:nk, sep = ".")
    
    # cos theta
    ssx <- rowSums(x ^ 2)
    if(any(ssx == 0)) ssx[ssx == 0] <- 1
    ssk <- rowSums(knots ^ 2)
    if(any(ssk == 0)) ssk[ssk == 0] <- 1
    x <- matrix(1 / sqrt(ssx), nrow = nx, ncol = xdim) * x
    knots <- matrix(1 / sqrt(ssk), nrow = nk, ncol = xdim) * knots
    ang <- tcrossprod(x, knots)
    
    # constants
    m <- 2L * m
    alpha <- 1 / (2 * m + 1)
    beta <- 2 * pi * factorial(2 * m)
    
    # return kernel
    if(m == 2L){
      X <- theta * (s2fun(ang) - alpha) / beta
    } else if(m == 4L){
      X <- theta * (s4fun(ang) - alpha) / beta
    } else {
      X <- theta * (s6fun(ang) - alpha) / beta
    }
    
    # check ridge
    if(ridge){
      Q <- penalty.sph(knots, m = m / 2, rescale = rescale)
      Qeig <- eigen(Q, symmetric = TRUE)
      Qrnk <- sum(Qeig$values > Qeig$values[1] * ncol(Q) * .Machine$double.eps)
      Qisqrt <- Qeig$vectors[,1:Qrnk,drop=FALSE] %*% diag(1/sqrt(Qeig$values[1:Qrnk]), nrow = Qrnk, ncol = Qrnk)
      X <- X %*% Qisqrt
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

# linear spherical spline
s2fun <- function(x){
  w <- (1 - x) / 2
  zix <- which(w <= .Machine$double.eps)
  w[zix] <- 0
  a <- log(1 + 1 / sqrt(w))
  a[zix] <- 0
  c <- 2 * sqrt(w)
  (a * (12 * w^2 - 4 * w) - 6 * c * w + 6 * w + 1) / 2
}

s4fun <- function(x){
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

s6fun <- function(x){
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