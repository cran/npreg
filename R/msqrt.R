msqrt <- 
  function(x, inverse = FALSE, symmetric = FALSE, 
           tol = .Machine$double.eps, checkx = TRUE){
    # matrix square root (or inverse square root)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2024-03-28
    
    x <- unname(as.matrix(x))
    if(checkx && !isSymmetric(x)) stop("Input 'x' must be a symmetric matrix")
    m <- nrow(x)
    i <- ifelse(inverse, -1/2, 1/2)
    if(m == 1L) return(x^i)
    if(max(abs(x[lower.tri(x)])) < .Machine$double.eps) return(diag(diag(x)^i))
    xeig <- eigen(x, symmetric = TRUE)
    xrnk <- sum(xeig$values > xeig$values[1] * tol * m)
    dmat <- diag(xeig$values[1:xrnk]^i, nrow = xrnk, ncol = xrnk)
    xsqrt <- xeig$vectors[,1:xrnk,drop=FALSE] %*% dmat
    if(symmetric) xsqrt <- xsqrt %*% t(xeig$vectors[,1:xrnk,drop=FALSE])
    return(xsqrt)
    
  }